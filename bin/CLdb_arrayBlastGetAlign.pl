#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastGetAlign.pl -- making fasta alignment of crDNA & protospacer

=head1 SYNOPSIS

CLdb_arrayBlastGetAlign.pl [flags] < blast_hits_crDNA_proto_aln.srl > aln.fasta

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -database  <char>

CLdb sqlite file name 
(if getting CLdb metadata on spacer sequences. eg., taxon_name)

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).
Only works with -database.

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).
Only works with -database.

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).
Only works with -database.

=item -query  <char>

Extra sql to refine CLdb query (must start with 'AND').
Only works with -database.

=item -outfmt  <char>

Output columns added to spacer-protospacer alignments.
See DESCRIPTION for more details.

=item -crDNA  <bool>

Orient by crDNA instead of protospacer? [FALSE]

=item -array  <bool>

Write out hits to spacers in CRISPR arrays (instead of
hits to protospacers)? [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastGetAlign.pl

=head1 DESCRIPTION

CLdb_arrayBlastAlignProto.pl must be
run before this script!

Get alignments of crDNA & protospacer.

The alignment is oriented to the protospacer by
default.

=head2 -outfmt

This flag designates the format
of each sequence in the output fasta.
The first 4 columns are always:
'crDNA/protospacer', 'locus_id', 'spacer_id', 'hsp_id'.
These are the unique ID for the spacer blast hit. 

Additional columns can be designated with a
comma-seperated list.

=head3 Example output (no '-outfmt'): 

">crRNA|1|10|l2puealx43xU"

'crRNA' = crRNA (not protospacer)

'1' = locus_id

'10' = spacer_id (for that locus)

'l2puealx43xU' = blast hsp unique ID


=head3 Support columns:

=over

=item * subtype

=item * taxon_name

=item * taxon_id

=item * scaffold   (of spacer/crDNA)

=item * crDNA_start

=item * crDNA_end

=item * array_sense_strand

=item * blastdb

=item * subject_scaffold

=item * subject_def

=item * subject_accession

=item * proto_start

=item * proto_end

=item * protoX_start  (X = proto + extension)

=item * protoX_end

=item * proto_strand

=item * identity

=item * evalue

=item * bitscore

=back

=head1 EXAMPLES

=head2 Basic Usage:


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
use Sereal qw/ encode_sereal /;
use FindBin;
use lib "$FindBin::RealBin/../lib/";

### CLdb
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::GetAlign qw/ get_alignProto 
				   parse_outfmt/;
use CLdb::query qw/ table_exists
		    n_entries
		    join_query_opts
		    getLociSpacerInfo/;
use CLdb::utilities qw/ file_exists
			connect2db/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $database_file;
my $outfmt;
my ($array, $crDNA_ori);
my $query = "";
my (@subtype, @taxon_id, @taxon_name);
my $evalue_cut;
GetOptions(
	   "database=s" => \$database_file,
	   "array" => \$array,
	   "crDNA" => \$crDNA_ori,
	   "outfmt=s" => \$outfmt,
           "subtype=s{,}" => \@subtype,
           "taxon_id=s{,}" => \@taxon_id,
           "taxon_name=s{,}" => \@taxon_name,	   
	   "query=s" => \$query,
	   "evalue=f" => \$evalue_cut,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, 'database') if defined $database_file;
my $outfmt_r = parse_outfmt($outfmt);

#--- MAIN ---#
# decoding spacer and DR srl
print STDERR "Decoding .srl file...\n" unless $verbose;
my $spacer_r = decode_file( fh => \*STDIN );


# if database: connect & query
my $queries_r;   # 
if(defined $database_file){
  # status
  print STDERR "Getting metadata info from CLdb...\n"
    unless $verbose;

  # connecting to CLdb
  my $dbh = connect2db($database_file);
  table_exists($dbh, 'loci');
  table_exists($dbh, 'spacers');

  my $join_sql = "";
  $join_sql .= join_query_opts(\@subtype, 'subtype');
  $join_sql .= join_query_opts(\@taxon_id, 'taxon_id');
  $join_sql .= join_query_opts(\@taxon_name, 'taxon_name');
  
  $queries_r = getLociSpacerInfo(dbh => $dbh, 
		    extra_sql => join(" ", $join_sql, $query), 
		    columns => [keys %{$outfmt_r->{CLdb}}]
		   );

  $dbh->disconnect;
}


# querying blastDBs for proteospacers
print STDERR "Getting alignments...\n" unless $verbose;
get_alignProto( blast => $spacer_r,
		outfmt => $outfmt_r,
		queries => $queries_r,
		array => $array,
		crDNA_ori => $crDNA_ori,
		evalue_cut => $evalue_cut,
		verbose => $verbose );


