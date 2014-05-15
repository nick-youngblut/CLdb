#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_array2fasta.pl -- write CRISPR array spacers or direct repeats to fasta

=head1 SYNOPSIS

CLdb_array2fasta.pl [flags] > array.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -repeat  <bool>

Get direct repeats instead of spacers. [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -cluster  <bool>

Get just cluster representatives of array elements (from clustering)? [FALSE]

=item -cutoff  <float>

Get array elements at the specified sequence identity cutoff. [1]

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_array2fasta.pl

=head1 DESCRIPTION

Get spacer or direct repeat sequences from the CRISPR database
and write them to a fasta.

By default, all spacers or direct repeats (if '-r') will be written.
The query can be refined with many of the flags. 

The sequences will be oriented base on the array_sense_strand
field in the loci table. To alter this field, see "CLdb_setSenseStrand.pl".

=head2 Output sequence naming

If -cluster: ">NA|spacer/DR|NA|clusterID|cluster_cutoff"

Else: ">locus_id|spacer/DR|elementID|NA|NA"

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite 

=head2 Write all direct repeats to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite -r

=head2 Write all unique spacers

CLdb_array2fasta.pl -d CLdb.sqlite -clust

=head2 Write all spacers clustered at 80% seqence ID

CLdb_array2fasta.pl -d CLdb.sqlite -clust -cut 0.8

=head2 Refine spacer sequence query:

CLdb_array2fasta.pl -d CLdb.sqlite -q "AND loci.Locus_ID=1" 

=head2 Refine spacer query to a specific subtype & 2 taxon_id's

CLdb_array2fasta.pl -d CLdb.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
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
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
		    table_exists
		    n_entries
		    join_query_opts
		    get_array_seq
		    get_array_seq_preCluster
		  /;
use CLdb::utilities qw/
			file_exists 
			connect2db/;
use CLdb::seq qw/
		  revcomp/;



#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my ($spacer_DR_b, $by_cluster);
my (@subtype, @taxon_id, @taxon_name);
my $strand_spec;
my $extra_query = "";
my $cluster_cutoff = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_DR_b,			 # spacers or repeats? [spacers]
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "cutoff=f" => \$cluster_cutoff, 		# clustering cutoff [1.00]
	   "cluster" => \$by_cluster,
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "spacers") unless $spacer_DR_b;
table_exists($dbh, "DRs") unless $spacer_DR_b;


# determining table of interest
my $tbl_oi = "spacers";
$tbl_oi = "DRs" if $spacer_DR_b;
table_exists($dbh, $tbl_oi) or die "ERROR: $tbl_oi not found in CLdb\n";
die "ERROR: no entries in $tbl_oi table!\n" 
	unless n_entries($dbh, $tbl_oi);


# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting arrays of interest from database #
## defining options for query subroutines ##
my %opts = (
	refine_sql => join("\n", $join_sql, $extra_query),
	by_cluster => $by_cluster,
	spacer_DR_b => $spacer_DR_b,
	cutoff => $cluster_cutoff,
	);

## querying CLdb ##
my $arrays_r;
my $cluster_tbl = $spacer_DR_b ? "DR_clusters" : "spacer_clusters";
if( table_exists($dbh, $cluster_tbl) ){
  $arrays_r =  get_array_seq($dbh,\%opts);
}
else{
  $arrays_r = get_array_seq_preCluster($dbh, \%opts);
}

# writing fasta #
map{ print join("\n", ">$_", $arrays_r->{$_}),"\n" } keys %$arrays_r;


#--- disconnect ---#
$dbh->disconnect();
exit;


