#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_loci2GFF.pl -- make GFF3 file from CRISPR loci features

=head1 SYNOPSIS

CLdb_loci2GFF.pl [flags] > loci.gff

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 options

=over

=item -id  <char>

Refine query to specific a subject taxon_id(s) (>1 argument allowed).

=item -name  <char>

Refine query to specific a subject taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine the query.

=item -loci  <char>

Write CRISPR loci features. [FALSE]

=item -array <char>

Write CRISPR array features. [FALSE]

=item -CAS <char>

Write CRISPR CAS features. [FALSE]

=item -spacer  <char>

Write CRISPR spacer features. [FALSE]

=item -dr  <char>

Write CRISPR direct-repeat features. [FALSE]

=item -leader  <char>

Write CRISPR leader features. [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_loci2GFF.pl

=head1 DESCRIPTION

Convert CRISPR loci data in CLdb 
to GFF3 file format for viewing features 
with a genome viewer (eg. jbrowse)

=head2 Output columns

=over

=item seqid = scaffold (subject_genome)

=item source = 'CLdb'

=item feature = 'region'

=item start = 'feature start'

=item end = 'feature end'

=item score = '.'

=item strand = '+' or '-' (determined by loci start-end)

=item phase = 0

=back 

=head3 'Attributes' column

=over

=item ID = taxon_ID

=item name = taxon_name

=item alias = locus_ID/feature_ID

=item note = info specific to the feature

=back

=head2 Important

'-name' or '-id' should be used
to limit hits to just 1 subject genome.

=head1 EXAMPLES

=head2 GFF3 of all loci in E.col

CLdb_loci2GFF.pl -d CLdb.sqlite -name "E.coli" -loci > ecoli_loci.gff

=head2 GFF3 of all CRISPR arrays & CASs in E.col

CLdb_loci2GFF.pl -d CLdb.sqlite -name "E.coli" -o a > ecoli_array-CAS.gff

=head2 GFF3 of all spacers in taxon_ID: 'FIG|2209.27'

CLdb_loci2GFF.pl -d CLdb.sqlite -id 2209.27 -s > 2209.27_spacers.gff

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
use Set::IntervalTree;
use DBI;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my (@taxon_id, @taxon_name); 				# blast subject query refinement
my ($write_loci, $write_spacers, $write_drs, $write_leaders, $write_CAS, $write_array);
my $extra_query = "";
my $len_cutoff = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "id=s{,}" => \@taxon_id,
	   "name=s{,}" => \@taxon_name,
	   "loci" => \$write_loci,
	   "spacers" => \$write_spacers,
	   "drs" => \$write_drs,
	   "leaders" => \$write_leaders,
	   "CAS" => \$write_CAS,
	   "array" => \$write_array,
	   "query=s" => \$extra_query,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# joining query options (for table join) #
my $join_sqls = join_query_opts_or(\@taxon_id, \@taxon_name);

# getting blast hits #
## just loci table ##
loci2gff($dbh, $join_sqls, $extra_query) if $write_loci;
array2gff($dbh, $join_sqls, $extra_query) if $write_array;
CAS2gff($dbh, $join_sqls, $extra_query) if $write_CAS;

## 2 table join ##
spacers2gff($dbh, $join_sqls, $extra_query) if $write_spacers;
drs2gff($dbh, $join_sqls, $extra_query) if $write_drs;
leaders2gff($dbh, $join_sqls, $extra_query) if $write_leaders;


### Subroutines
sub leaders2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.locus_start,
a.locus_end,
a.array_status,
a.CAS_status,
b.locus_id,
b.leader_start,
b.leader_end
FROM loci a, leaders b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" leader";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[10],			# feature start (sstart)
			$$row[11],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}
	
sub drs2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.locus_start,
a.locus_end,
a.array_status,
a.CAS_status,
b.dr_id,
b.dr_start,
b.dr_end
FROM loci a, drs b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" array_status:\"$$row[7]\" CAS_status:\"$$row[8]\"";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[10],			# feature start (sstart)
			$$row[11],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\__$$row[9]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}

sub spacers2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.locus_start,
a.locus_end,
a.array_status,
a.CAS_status,
b.spacer_id,
b.spacer_start,
b.spacer_end
FROM loci a, spacers b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" array_status:\"$$row[7]\" CAS_status:\"$$row[8]\"";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[10],			# feature start (sstart)
			$$row[11],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\__$$row[9]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}

sub CAS2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.CAS_start,
a.CAS_end,
a.CAS_status
FROM loci a, loci b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" CAS_status:\"$$row[7]\"";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[5],			# feature start (sstart)
			$$row[6],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}

sub array2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.array_start,
a.array_end,
a.array_status
FROM loci a, loci b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" array_status:\"$$row[7]\"";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[5],			# feature start (sstart)
			$$row[6],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}

sub loci2gff{
# writing loci as gff feature #
	my ($dbh, $join_sqls, $query) = @_;

	# selecting loci #
	my $q = "
SELECT
a.locus_id,
a.taxon_name,
a.taxon_id,
a.subtype,
a.scaffold,
a.locus_start,
a.locus_end,
a.array_status,
a.CAS_status
FROM loci a, loci b
where a.locus_id = b.locus_id
$join_sqls
$query
";
	$q =~ s/\n/ /g;
	print STDERR "$q\n" if $verbose;
	
	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
	
	# writing gff #
# seqid = locus_id
# source = 'CLdb'
# feature = 'region'
# start = 'locus_start'
# end = 'locus_end'
# score = .
# strand = strand (based on start-end)
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = array_status, CAS_status
	
	foreach my $row (@$ret){
		my $strand = "+";
		$strand = "-" if $$row[5] > $$row[6]; 		# '-' if start > end
		my $note = "subtype:\"$$row[3]\" array_status:\"$$row[7]\" CAS_status:\"$$row[8]\"";
		#$$row[0] =~ s/^/cli./;
		
		print join("\t",
			$$row[4],			# subject scaffold
			"CLdb",				# source
			"region",
			$$row[5],			# feature start (sstart)
			$$row[6],			# feature end	(send)
			".",				# score
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$row[2]\"",			# query taxon_name
				"Name=\"$$row[1]\"",		# query taxon_id
				"Alias=\"$$row[0]\"",			# locus_id
				"Note='$note'"
				)), "\n";
		
		}
	}

sub join_query_opts_or{
	my ($taxon_id_r, $taxon_name_r) = @_;
	
	return "" unless @$taxon_id_r || @$taxon_name_r;
	
	# adding quotes #
	map{ s/"*(.+)"*/'$1'/ } @$taxon_id_r;
	map{ s/"*(.+)"*/'$1'/ } @$taxon_name_r;	
	
	if(@$taxon_id_r && @$taxon_name_r){
		return join("", " AND (a.taxon_id IN (", join(", ", @$taxon_id_r),
						") OR a.taxon_name IN (", join(", ", @$taxon_name_r),
						"))");
		}
	elsif(@$taxon_id_r){
		return join("", " AND a.taxon_id IN (", join(", ", @$taxon_id_r), ")");
		}
	elsif(@$taxon_name_r){
		return join("", " AND a.taxon_name IN (", join(", ", @$taxon_name_r), ")");	
		}
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND a.$cat IN (", join(", ", @$vals_r), ")");
	}


