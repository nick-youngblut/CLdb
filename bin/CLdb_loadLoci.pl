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
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	list_columns/;
use CLdb::load qw/
	load_db_table/;
use CLdb::seq qw/
	read_fasta
	seq_from_genome_fasta/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	lineBreaks2unix
	get_file_path/;
use CLdb::load::loadLoci qw/
			     unix_line_breaks
			   /;
use CLdb::genbank::genbank2fasta qw/
				     genbank2fasta
				   /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $quiet);
GetOptions(
	   "database=s" => \$database_file,
	   "quiet" => \$quiet,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );


#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# database path #
my $db_path = get_file_path($database_file);

# database metadata #
table_exists($dbh, "loci"); 
#my $column_list_r = list_columns($dbh, "loci", 1);

# getting input loci table #
my ($loci_r, $header_r) = get_loci_table();

# checks#
check_locus_id($loci_r);
make_external_file_dirs($loci_r, $header_r, $db_path);	 # copying array files & genbanks if not in ./genbank & ./array #
unix_line_breaks($loci_r, $db_path);						# line breaks of all copied files to unix


# inferring from table #
## getting genome fasta from genbank ##
genbank2fasta($loci_r, $db_path, $header_r); 	
## getting scaffold name 
get_scaffold_name($loci_r, $db_path, $header_r);

# updating / loading_db #
## loci ##
print STDERR "\n### loading entries into CLdb ###\n";
my $loci_header_r = just_table_columns($header_r, 'loci');
load_db_table($dbh, "loci", $loci_header_r, $loci_r);

## loading leader table ##
my $leader_header_r = just_table_columns($header_r, 'leader');
$leader_header_r->{"leader_sequence"} = 1;
my $leader_loci_r = get_leader_seq($loci_r, $db_path);
load_db_table($dbh, "leaders", $leader_header_r, $leader_loci_r);

## loading pam table ##
my $pam_header_r = just_table_columns($header_r, 'pam');
my $pam_loci_r = get_pam_seq($loci_r, $db_path, $pam_header_r);
$pam_header_r->{"pam_sequence"} = 1;
load_db_table($dbh, "pam", $pam_header_r, $pam_loci_r);

# disconnect to db #
$dbh->disconnect();
exit;


### Subroutines
sub check_locus_id{
# no '|' in locus_ID #
	my ($loci_r) = @_;
	
	warn "### checking locus_ID values ###\n";
	
	foreach my $locus_id (keys %$loci_r){
		die "ERROR: locus_id '$locus_id' contains '|', which is not allowed!\n"
			if $locus_id =~ /\|/;
		}
	}

sub get_pam_seq{
# getting pam seq if needed #
	my ($loci_r, $db_path, $pam_header_r) = @_;
	
	if(exists $pam_header_r->{"pam_start"} && exists $pam_header_r->{"pam_end"}
		&& exists $pam_header_r->{"pam_sequence"}){
		print STDERR "### PAM sequence & start-end columns provided. Getting PAM sequence if needed. Loading table ###\n"
		}
	elsif(exists $pam_header_r->{"pam_start"} && exists $pam_header_r->{"pam_end"}){
		print STDERR "### PAM start-end columns provided. Getting PAM sequence if needed. Loading table ###\n"
		}
	elsif(exists $pam_header_r->{"pam_sequence"}){
		print STDERR "### PAM sequence column provided. Loading values into PAM table ###\n"
		}
	else{
		print STDERR "### no PAM info provided. Skipping PAM loading ###\n"
		}
	
	# getting pam sequence if possible #
	my %cp;
	foreach my $locus_id (keys %$loci_r){
		if(exists $loci_r->{$locus_id}{'pam_sequence'}){
			$cp{$locus_id} = $loci_r->{$locus_id};
			}
		elsif( exists $loci_r->{$locus_id}{'pam_start'}
			   && exists $loci_r->{$locus_id}{'pam_start'} 
			   && exists $loci_r->{$locus_id}{'scaffold'}
			   && exists $loci_r->{$locus_id}{'fasta_file'}){		# geting sequence 
			my $fasta_r = read_fasta("$db_path/fasta/$loci_r->{$locus_id}{'fasta_file'}");
			
			$loci_r->{$locus_id}{'pam_sequence'} = 
				seq_from_genome_fasta( $fasta_r, 
						[$loci_r->{$locus_id}{'scaffold'},
						$loci_r->{$locus_id}{'pam_start'}, 
						$loci_r->{$locus_id}{'pam_end'}]
						);
			}
		else{
			next;
			}
		
	
		# checking for existence of genome fasta, if yes, extract leader sequence #
		if(exists $loci_r->{$locus_id}{'fasta_file'}){
			
			my $fasta_r = read_fasta("$db_path/fasta/$loci_r->{$locus_id}{'fasta_file'}");
			
			$loci_r->{$locus_id}{'leader_sequence'} = 
				seq_from_genome_fasta( $fasta_r, 
						[$loci_r->{$locus_id}{'scaffold'},
						$loci_r->{$locus_id}{'leader_start'}, 
						$loci_r->{$locus_id}{'leader_end'}]
						) unless exists $loci_r->{$locus_id}{'leader_sequence'};
			
			if($loci_r->{$locus_id}{'pam_sequence'} eq ""){
				print STDERR "WARNING: '", $loci_r->{$locus_id}{'scaffold'}, 
					"' not found for ", $loci_r->{$locus_id}{'fasta_file'}, 
					". Not loading pam sequence!\n";
				}
			else{		# just entries w/ leader sequence #
				$cp{$locus_id} = $loci_r->{$locus_id};
				}
			}
		}

		#print Dumper %cp; exit;
	return \%cp;
	
	}

sub get_leader_seq{
# loading leader; pulling out sequence from genome if available #
	my ($loci_r, $db_path) = @_;
	
	# status #
	print STDERR "### Leader start-end provided. Loading values into leader table ###\n"
		unless $verbose;
	
	# getting leader sequence if possible #
	my %cp;
	foreach my $locus_id (keys %$loci_r){
		# next unless leader_start && leader_end #
		next unless exists $loci_r->{$locus_id}{'leader_start'} 
					&& exists $loci_r->{$locus_id}{'leader_end'}
					&& exists $loci_r->{$locus_id}{'scaffold'}
					&& exists $loci_r->{$locus_id}{'fasta_file'};
	
		# checking for existence of genome fasta, if yes, extract leader sequence #
		if(exists $loci_r->{$locus_id}{'fasta_file'}){
			
			my $fasta_r = read_fasta("$db_path/fasta/$loci_r->{$locus_id}{'fasta_file'}");
			
			$loci_r->{$locus_id}{'leader_sequence'} = 
				seq_from_genome_fasta( $fasta_r, 
						[$loci_r->{$locus_id}{'scaffold'},
						$loci_r->{$locus_id}{'leader_start'}, 
						$loci_r->{$locus_id}{'leader_end'}]
						) unless exists $loci_r->{$locus_id}{'leader_sequence'};
			
			if($loci_r->{$locus_id}{'leader_sequence'} eq ""){
				print STDERR "WARNING: '", $loci_r->{$locus_id}{'scaffold'}, 
					"' not found for ", $loci_r->{$locus_id}{'fasta_file'}, 
					". Not loading leader sequence!\n";
				}
			else{		# just entries w/ leader sequence #
				$cp{$locus_id} = $loci_r->{$locus_id};
				}
			}
		}

		#print Dumper %cp; exit;
	return \%cp;
	}

sub just_table_columns{
#-- Description --#
# loading just the columns of interest #
  my ($header_r, $table) = @_;
  $table =~ tr/A-Z/a-z/;
  
  my %col_lists;
  @{$col_lists{'loci'}} = qw/Locus_ID Taxon_ID Taxon_Name 
			     Subtype Scaffold 
			     Locus_start Locus_end
			     CAS_Start CAS_End  
			     Array_Start Array_End
			     CAS_Status Array_Status 
			     Genbank_File Fasta_File Array_File 
			     Scaffold_count 
			     File_Creation_Date Author/;
  
  @{$col_lists{'leader'}} = qw/Locus_ID Leader_start Leader_end Leader_end Leader_sequence/;
  @{$col_lists{'pam'}} = qw/Locus_ID PAM_seq PAM_start PAM_end PAM_sequence/;
  
  die "ERROR: '$table' not found in column lists!\n"
    unless exists $col_lists{$table};
  
  map{ $_ =~ tr/A-Z/a-z/ } @{$col_lists{$table}};		# table of interest
	
  my %header_parse;
  map{$header_parse{$_} = $header_r->{$_} if exists $header_r->{$_}} @{$col_lists{$table}};
  
  #print Dumper %header_parse; exit;
  return \%header_parse; 
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

sub make_external_file_dirs{
  # moving files to directory in CLdb_home #
  my ($loci_r, $header_r, $dir) = @_;
  
  foreach my $locus_id (keys %$loci_r){			
    $loci_r->{$locus_id}{"genbank_file"} = 
      make_CLdb_dir($dir, 'genbank', $loci_r->{$locus_id}{"genbank_file"}) 
	if exists $loci_r->{$locus_id}{"genbank_file"};
    $loci_r->{$locus_id}{"array_file"} = 
      make_CLdb_dir($dir, 'array', $loci_r->{$locus_id}{"array_file"}) 
	if exists $loci_r->{$locus_id}{"array_file"};
    $loci_r->{$locus_id}{"fasta_file"} = 
      make_CLdb_dir($dir, 'fasta', $loci_r->{$locus_id}{"fasta_file"})
	if exists $loci_r->{$locus_id}{"fasta_file"};
  }
  
  #print Dumper %$loci_r; exit;	
  return $dir;
}

sub make_CLdb_dir{
  my ($dir, $name, $infile) = @_;
  my @parts = File::Spec->splitpath( $infile );
  
  # making dir; copying files #
  if(File::Spec->rel2abs($parts[1]) ne "$dir/$name"){
    mkdir "$dir/$name" unless -d "$dir/$name";
    unless(-e "$dir/$name/$parts[2]"){		# can't find in needed directory, is it in specified directory?
      die " ERROR: cannot find $infile\n" 
	unless -e $infile;			# specified file can't be found anywhere
      
      print STDERR "'$infile' not in CLdb_HOME/$name/\n";
      copy($infile, "$dir/$name/$parts[2]") or die $!;
      print STDERR "\tCopied $infile to $dir/$name/$parts[2]\n" unless $quiet;
    }
  }
  
  return $parts[2];

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


