#!/usr/bin/env perl

=pod

=head1 NAME

getGeneIDs.pl -- getting Gene_ID values (fig|peg) for particular CRISPR loci

=head1 SYNOPSIS

getGeneIDs.pl [flags] > Gene_IDs.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -verbose  <bool>

Verbose output

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc getGeneIDs.pl

=head1 DESCRIPTION

Get Gene_IDs (FIG|PEG IDs) from the CRISPR database
and write them to a table. Output columns are:

=over

=item * Gene_ID

=item * Taxon_name

=item * Taxon_ID

=item * Locus_ID

=item * Subtype

=back

By default, all Gene_IDs will be written. Use the 
flags for query refinement to get Gene_IDs for
particular loci.

=head1 EXAMPLES

=head2 All Gene_IDs in Subtype I-B loci

getGeneIDs.pl -data CRISPR.sqlite -sub "I-B" 

=head2 All Gene_IDs in Subtype I-B & 2 particular Taxon_ID's

getGeneIDs.pl -da CRISPR.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head2 All Gene_IDs in Subtype I-B loci & pipe to ITEP to get clusters

getGeneIDs.pl -data CRISPR.sqlite -sub "I-B" | db_getClustersContainingGenes.py | less

=head2 All Gene_IDs in Subtype I-B & 2 particular Taxon_ID's

getGeneIDs.pl -da CRISPR.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head2 Using '-q' to pick genes in operons (and in subtype 'I-B')

getGeneIDs.pl -da CRISPR.sqlite -sub I-B -q "AND genes.In_Operon='yes'"

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
use DBI;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
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
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# query db #
get_figs_join($dbh, $extra_query, $join_sql);

# disconnect to db #
$dbh->disconnect();
exit;


### Subroutines
sub get_figs_join{
	my ($dbh, $extra_query, $join_sql) = @_;

	
	# make query #
	my $query = "SELECT  genes.gene_id, loci.taxon_name, loci.taxon_id, loci.locus_id, loci.subtype
	from genes, loci where genes.locus_id = loci.locus_id $join_sql";
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	# writing out figs #
	foreach my $row (@$ret){
		next unless $$row[0];
		map{$_ = "" unless $_} @$row;
		print join("\t", @$row), "\n";
		}
	
	}

sub join_query_opts_OLD{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND b.$cat IN (", join(", ", @$vals_r), ")");
	}



