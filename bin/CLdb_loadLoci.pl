#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_loadLoci.pl -- adding/updating loci entries in to CLdb

=head1 SYNOPSIS

CLdb_loadLoci.pl [flags] < loci_table.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -forks  <int>

Number of files to process in parallel. [1]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_loadLoci.pl

=head1 DESCRIPTION

Add and/or update CLdb 'loci' table entries.
The 'Locus_ID' column must be unique! 
Existing entries that have the same Locus_ID
will be replaces (entry update). 

Leader start-stop and sequence information can be provided also.
This information will be added to the 'leaders' table.

PAM start-stop and/or sequence information can be provided also.
This information will be added to the 'PAM' table.

Scaffold names will try to be obtained from the genome genbank/fasta
if not provided.

=head2 REQUIRED columns/values in loci file

=over

=item * locus_id

=item * taxon_id

=item * taxon_name

=item * locus_start

=item * locus_end

=item * operon_status

=item * array_status

=item * genbank_file

=item * array_file

=item * author

=back

=head2 WARNING!

The loci table must be in tab-delimited format.

Extra columns (not in CLdb Loci table) can exist in the 
input loci table. They just won't be added to CLdb.

=head1 EXAMPLES

=head2 Usage:

CLdb_loadLoci.pl -d CLdb.sqlite < loci.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use File::Copy;
use Bio::SeqIO;
use DBI;
use List::Util qw/max/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::query qw/
		    table_exists
		    list_columns/;
use CLdb::load qw/
		   load_db_table/;
use CLdb::load::loadLoci qw/
			     unix_line_breaks
			     just_table_columns
			     make_external_file_dirs
			     get_leader_seq
			     get_pam_seq
			   /;
use CLdb::seq qw/
		  read_fasta
		  seq_from_genome_fasta/;
use CLdb::utilities qw/
			file_exists 
			connect2db
			lineBreaks2unix
			get_file_path/;
use CLdb::genbank::genbank2fasta qw/
				     genbank2fasta
				   /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my $forks = 0;
GetOptions(
	   "database=s" => \$database_file,
	   "forks=i" => \$forks,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );


#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# database path 
my $db_path = get_file_path($database_file);

# getting input loci table 
my ($loci_r, $header_r) = get_loci_table();
unix_line_breaks($loci_r, $db_path, $forks) unless $^O =~ /win/i;   # line breaks to unix unless windows OS

# checks
check_locus_id($loci_r);
make_external_file_dirs($loci_r, $header_r, $db_path);	 # copying array files & genbanks if not in ./genbank & ./array #

# connect 2 db 
my $dbh = connect2db($database_file);

# database metadata 
table_exists($dbh, "loci"); 


# inferring from table
## getting genome fasta from genbank
genbank2fasta($loci_r, $db_path, $header_r); 	
## getting scaffold name 
get_scaffold_name($loci_r, $db_path, $header_r);

# updating / loading_db 
## loci
print STDERR "\n### loading entries into CLdb ###\n";
my $loci_header_r = just_table_columns($header_r, 'loci');
load_db_table($dbh, "loci", $loci_header_r, $loci_r);

## loading leader table 
# TODO: update
my $leader_header_r = just_table_columns($header_r, 'leader');
$leader_header_r->{"leader_sequence"} = 1;
my $leader_loci_r = get_leader_seq($loci_r, $db_path);
load_db_table($dbh, "leaders", $leader_header_r, $leader_loci_r);

## loading pam table
# TODO: update
#my $pam_header_r = just_table_columns($header_r, 'pam');
#my $pam_loci_r = get_pam_seq($loci_r, $db_path, $pam_header_r);
#$pam_header_r->{"pam_sequence"} = 1;
#load_db_table($dbh, "pam", $pam_header_r, $pam_loci_r);


# disconnect to db #
$dbh->disconnect();
exit;


#--- Subroutines ---#
sub check_locus_id{
# no '|' in locus_ID #
  my ($loci_r) = @_;
  
  print STDERR "### checking locus_ID values ###\n";
  foreach my $locus_id (keys %$loci_r){
    die "ERROR: locus_id '$locus_id' contains '|', which is not allowed!\n"
			if $locus_id =~ /\|/;
  }
  print STDERR "...locus_ID values are OK\n";
}



sub get_scaffold_name{
# if no scaffold value, try to get from fasta #
  my ($loci_r, $db_path, $header_r) = @_;
  
  print STDERR "### checking for genome fasta for scaffold names ###\n";
  
  # getting all fasta files needed #
  my %fasta_need;
  foreach my $locus_id (keys %$loci_r){
    next unless exists $loci_r->{$locus_id}{'fasta_file'};		# cannot do w/out genome fasta
    $fasta_need{$loci_r->{$locus_id}{'fasta_file'}} = 1
      unless $loci_r->{$locus_id}{'scaffold'};				# scaffold already exists
		}
  
  # if needed fastas, add 'scaffold' to header #
  $header_r->{'scaffold'} = 'X';			# fake column index
  
  # extracting scaffolds #
  my %fasta_scaf;
  foreach my $file (keys %fasta_need){
    my $path_file = "$db_path/fasta/$file";
    die " ERROR: cannot find $path_file!\n" unless -e $path_file;
    
    open IN, $path_file or die $!;
    while(<IN>){
      chomp;
      if(/^>/){
	s/^>//;
	die " ERROR: multiple scaffolds in $file! You must designate the scaffold name yourself!\n"
	  if exists $fasta_scaf{$file};
	$fasta_scaf{$file} = $_; 	# file => scaffold_name
      }
      die " ERROR: count not find a scaffold name in the genome fasta: $file! You must designate the scaffold name yourself!\n"
	unless exists $fasta_scaf{$file};
    }
    close IN;
  }
  
  # adding scaffold values #
  foreach my $locus_id (keys %$loci_r){
    $loci_r->{$locus_id}{'scaffold'} = $fasta_scaf{$loci_r->{$locus_id}{'fasta_file'}}
			unless exists $loci_r->{$locus_id}{'scaffold'};
    die " ERROR: could not add a scaffold to locus: $locus_id! Add it yourself!\n"
      unless exists $loci_r->{$locus_id}{'scaffold'};
		}
  
  #print Dumper %$loci_r; exit;
  #print Dumper @fasta_need; exit;
}



sub get_loci_table{
# loading tab-delimited table from STDIN #
  
  my $tbl_r = lineBreaks2unix(\*STDIN);
  
  # loading into a hash #
  my %loci;
  my %header;
  my %header_rev;
  my $line_cnt = 0;
  foreach (@$tbl_r){
    chomp;
    $line_cnt++;
    next if /^\s*$/;
    
    if($line_cnt == 1){ 					# loading header
      tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
      tr/ /_/;
      my @line = split /\t/;
      for my $i (0..$#line){
	#next unless grep(/^$line[$i]$/, @column_list);		# only including columns in loci DB table
	next if exists $header{$line[$i]}; 					# only using 1st column found with a particular name
	$header{$line[$i]} = $i;		# column_name => index
	$header_rev{$i} = $line[$i];		# column_name => index
      }
      
      # checking for locus_id #
      die " ERROR: 'locus_id' column not found!\n"
	unless exists $header{'locus_id'};
      
    }
    else{
      my @line = split /\t/;	
      
      my $locus_id = $line[$header{'locus_id'}];
      die " ERROR: locus_id: '$locus_id' is not unique in loci table!\n"
	if exists $loci{$locus_id};
      
      foreach my $i (0..$#line){
	next unless $header_rev{$i};
	$loci{$locus_id}{$header_rev{$i}} = $line[$i];
      }
    }
  }
  # sanity check #
  die " ERROR: entries found in loci table!\n" unless %loci;	
  #print Dumper %loci; exit; 
  #print Dumper %header; exit;	
  return (\%loci, \%header);
}


