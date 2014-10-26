#!/usr/bin/env perl

=pod

=head1 NAME

setSenseStrand.pl -- set the sense stand for CRISPR arrays

=head1 SYNOPSIS

setSenseStrand.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -locus_id  <char>

Refine query to specific >=1 locus_id (>1 argument allowed).

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -leader  <bool>

Use the leader region (if determined) to determine the sense
strand (leader assumed to be upstream)? [TRUE]

=item -sense_file  <char>

File designating sense for >=1 locus. Format: 'locus_id\t[1|-1]'

=item -x  <bool>

Write loci to table STDOUT to view sense strand designations.
Query refinement flags can be used.

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc setSenseStrand.pl

=head1 DESCRIPTION

The CRISPR array could potentially be
transcribe on either strand.

The sense strand is determined by either:

=over

=item * The leader region position (only for loci where leader has been determined)

=item * array_start & array_end in the loci table (start > end = negative strand)

=item * The sense file (if provided)

=back

=head2 sense_file

Optionally, a file designating sense strand
for 1 or more locus can be provided.
Format: tab-delimited, 2 column (locus_id, strand [1|-1])

=head1 EXAMPLES


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
use CLdb::query qw/
		    table_exists
		    n_entries
		    join_query_opts
		    get_array_seq
		    queryLociTable
		  /;

use CLdb::utilities qw/
			file_exists 
			connect2db/;
use CLdb::senseStrand qw/
			  parse_sense_file
			  setSenseByArraySE
			  getLeaderCoords
			  senseByLeaderLoc
			  updateSenseStrand
			/;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $sense_file);
my (@locus_id, @subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my ($use_leaders, $write_loci_table);
GetOptions(
	   "database=s" => \$database_file,
	   "locus_id=s{,}" => \@locus_id,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "leader" => \$use_leaders,
	   "sense_file=s", \$sense_file,
	   "x" => \$write_loci_table,
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
file_exists($sense_file, "sense_file") if
  defined $sense_file;


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "leaders") if
  ! defined $sense_file || ! defined $use_leaders;  # using leaders unless told not to


# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@locus_id, "locus_id");
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


# if $write_loci_table: writing out locus table & calling it a day
if( $write_loci_table ){
  my $loci_tbl_r = queryLociTable(dbh => $dbh, 
				  refine => join("\t", $join_sql, $extra_query)
				  );
  map{ print join("\t", @$_),"\n" } @$loci_tbl_r; 
  exit;
}


# either loading loci (& leader table), or loading sense strand table
my $sense_r;
if(defined $sense_file){  # using sense_file
  $sense_r = parse_sense_file($sense_file);
}
else{
  # getting loci table
  my $loci_tbl_r = queryLociTable(dbh => $dbh,
				  refine => join("\t", $join_sql, $extra_query),
				  columns => ['locus_id', 'array_start', 'array_end'],
				  );
  # determing sense strand from array start-end
  $sense_r = setSenseByArraySE( $loci_tbl_r,
				{locus_id => 0, 
				 array_start => 1, 
				 array_end => 2}
			      );
  
  # unless -leader: define by leader position (if leader provided)
  unless( $use_leaders ){
    my $leaders_r = getLeaderCoords( dbh => $dbh,
		     loci_tbl => $loci_tbl_r, 
		     locus_id_column => 0 );
    # editing sense strand based on location of leader sequences
    $sense_r = senseByLeaderLoc( $sense_r, $leaders_r );
  }
}

# updating loci table with sense_strand information
updateSenseStrand($dbh, $sense_r);


#--- disconnect ---#
$dbh->disconnect();
exit;


