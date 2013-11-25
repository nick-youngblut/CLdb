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

Load and/or update CRISPR loci table entries.

Loci entries will be added if no Locus_id is provided;
otherwise, the entrie will be updated.

=head2 WARNING!

The loci table must be in tab-delimited format.

Loci entries must be unique by 'taxon_name', 'taxon_id', 'scaffold', 'locus_start', 'locus_end'.
Only the 1st of duplicated entries will be kept.

A scaffold name will be added to any entries lacking one. 

File names in the loci table should match files in the 
$CLdb_HOME/genbank/ & $CLdb_HOME/array/ directories.

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
use DBI;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	list_columns/;
use CLdb::load qw/
	load_db_table/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	lineBreaks2unix/;



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

# database metadata #
table_exists($dbh, "loci"); 
my $column_list_r = list_columns($dbh, "loci", 1);

# getting input loci table #
my ($loci_r, $header_r) = get_loci_table($column_list_r);

# adding scaffold values if not found #
add_scaffold_ID($loci_r, $header_r);

# copying array files & genbanks if not in ./genbank & ./array #
make_genbank_array_dirs($loci_r, $header_r);

# updating / loading_db #
load_db_table($dbh, "loci", $header_r, $loci_r);


# loading leader table #
load_leader($dbh, $loci_r, $db_path) 
	if exists $header_r->{'leader_start'} && exists $header_r->{'leader_end'};

# loading pam table #
#load_pam($dbh, $loci_r) if exists $header_r->{'pam_seq'};

# loading arrays if found #


# status report #

# disconnect to db #
$dbh->disconnect();
exit;


### Subroutines
sub make_genbank_array_dirs{
	my ($loci_r, $header_r) = @_;
	
	# current directory #
	my $dir = File::Spec->rel2abs(File::Spec->curdir());
	
	my @cp_warning;
	foreach my $locus_id (keys %$loci_r){
		my $row = $loci_r->{$locus_id};
		
		# genbank file #
		if( $$row[$header_r->{"genbank_file"}] ){
			my @parts = File::Spec->splitpath( $$row[$header_r->{"genbank_file"}] );
				# making genbank dir; copying file #
			if(File::Spec->rel2abs($parts[1]) ne "$dir/genbank"){
				mkdir "$dir/genbank" unless -d "$dir/genbank";
				unless(-e "$dir/genbank/$parts[2]"){
					die " ERROR: cannot find ", $$row[$header_r->{"genbank_file"}], "\n"
						unless -e $$row[$header_r->{"genbank_file"}];
						
					copy($$row[$header_r->{"genbank_file"}], "$dir/genbank/$parts[2]") or die $!;
					print STDERR "...Copied ", $$row[$header_r->{"genbank_file"}], 
						" to $dir/genbank/$parts[2]\n" unless $quiet;
					}
				}
			# stripping path from genbank value #
			$$row[$header_r->{"genbank_file"}] = $parts[2];
			# sanity check #
			my $file_chk = join("/", $dir, "genbank", $$row[$header_r->{"genbank_file"}]);
			die " ERROR: cannot find ", $file_chk, "\n"
				unless -e $file_chk;
			} 
		
		# array file #
		if( exists $header_r->{"array_file"} && $$row[$header_r->{"array_file"}] ){
			my @parts = File::Spec->splitpath( $$row[$header_r->{"array_file"}] );

			# making genbank dir; copying file #
			if(File::Spec->rel2abs($parts[1]) ne "$dir/array"){
				mkdir "$dir/array" unless -d "$dir/array";
				unless(-e "$dir/array/$parts[2]"){
					die " ERROR: cannot find ", $$row[$header_r->{"array_file"}], "\n"
						unless -e $$row[$header_r->{"array_file"}];
						
					copy($$row[$header_r->{"array_file"}], "$dir/array/$parts[2]") or die $!;
					print STDERR "...Copied ", $$row[$header_r->{"array_file"}], 
						" to $dir/array/$parts[2]\n" unless $quiet;
					}
				}
			# stripping path from array value #
			$$row[$header_r->{"array_file"}] = $parts[2];
			# sanity check #
			my $file_chk = join("/", $dir, "array", $$row[$header_r->{"array_file"}]);
			die " ERROR: cannot find ", $file_chk, "\n"
				unless -e $file_chk;
			}
		# fasta file #
		if( exists $header_r->{"fasta_file"} && $$row[$header_r->{"fasta_file"}] ){
			my @parts = File::Spec->splitpath( $$row[$header_r->{"fasta_file"}] );

			# making genbank dir; copying file #
			if(File::Spec->rel2abs($parts[1]) ne "$dir/fasta"){
				mkdir "$dir/fasta" unless -d "$dir/fasta";
				unless(-e "$dir/fasta/$parts[2]"){
					die " ERROR: cannot find ", $$row[$header_r->{"fasta_file"}], "\n"
						unless -e $$row[$header_r->{"fasta_file"}];
						
					copy($$row[$header_r->{"fasta_file"}], "$dir/fasta/$parts[2]") or die $!;
					print STDERR "...Copied ", $$row[$header_r->{"fasta_file"}], 
						" to $dir/fasta/$parts[2]\n" unless $quiet;
					}
				}
			# stripping path from array value #
			$$row[$header_r->{"fasta_file"}] = $parts[2];
			# sanity check #
			my $file_chk = join("/", $dir, "fasta", $$row[$header_r->{"fasta_file"}]);
			die " ERROR: cannot find ", $file_chk, "\n"
				unless -e $file_chk;
			}					
		}
	}

sub add_scaffold_ID{
# adding scaffold ID if not found; needed to make sure entries are unique #
	my ($loci_r, $header_r) = @_;
	
	# adding scaffold to header if not found #
	$header_r->{"scaffold"} = scalar keys %$header_r unless exists $header_r->{"scaffold"};
	
	my $scaf_add_bool = 0;
	foreach my $locus_id (keys %$loci_r){
			#foreach my $row (@{$loci_r->{$entry_type}}){
			#print Dumper ${$loci_r->{$locus_id}}[$header_r->{"scaffold"}]; exit;
		unless( defined ${$loci_r->{$locus_id}}[$header_r->{"scaffold"}]){
			${$loci_r->{$locus_id}}[$header_r->{"scaffold"}] = "CLDB__ONE_CHROMOSOME";
			$scaf_add_bool = 1;
			}
		}

	# status #
	print STDERR "...Did not find 'scaffold' values for some entries. Adding scaffold names as: 'CLDB__ONE_CHROMOSOME'\n"
		if $scaf_add_bool and ! $verbose;
	}

sub get_loci_table{
# loading tab-delimited table from STDIN #
	my $column_list_r = shift;
	
	# checking line breaks #
	my $tbl_r = lineBreaks2unix(\*STDIN);
	
	# loading into a hash #
	my %loci;
	my %header;
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
				next unless exists $column_list_r->{$line[$i]};		# only including columns in loci DB table
				next if exists $header{$line[$i]}; 					# only using 1st column found with a particular name
				$header{$line[$i]} = $i;		# column_name => index
				}
			}
		else{
			check_headers(\%header) if $.==2;	
			my @line = split /\t/;
			
			if(exists $loci{$line[$header{"locus_id"}]} ){
				die "ERROR: '", $line[$header{"locus_id"}], "' is duplicated!\n";
				}
			$loci{$line[$header{"locus_id"}]} = \@line;
			}
		}
	# sanity check #
	die " ERROR: entries found in loci table!\n" unless %loci;	
	
		#print Dumper %loci; exit; 
		#print Dumper %header; exit;
	return (\%loci, \%header);
	}

sub check_headers{
	my ($header_r) = @_;
	my @req = qw/locus_id taxon_id taxon_name locus_start locus_end operon_status array_status genbank_file array_file author/;
	
	# checking for required headers not found #
	my (@not_found, @found);
	foreach (@req){
		if (exists $header_r->{$_}){
			push @found, $_;
			}
		else{
			push @not_found, $_;
			}
		}
	
	if(@not_found){
		print STDERR "ERROR: Required columns not found in loci table!\n\n";
		print STDERR "### Required columns not found (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort @not_found), "\n";
		print STDERR "\n### Correct headers found in the loci table (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort @found), "\n";
		exit;
		}
	
	}
	



