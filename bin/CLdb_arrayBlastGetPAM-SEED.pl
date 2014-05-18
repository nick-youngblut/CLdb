#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastGetPAM-SEED.pl -- getting PAM for each protospacer

=head1 SYNOPSIS

CLdb_arrayBlastGetPAM-SEED.pl [flags] < blast_hits_proto.srl > PAMs.txt

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -PAM  <int>

start-stop of the PAM region (2 values required). 
See DESCRIPTION for details. [1 3]

=item -SEED  <int>

start-end of the seed region (1 or 2 values required).
See DESCRIPTION for details. [-8]

=item -database  <char>

CLdb sqlite file name 
(if getting metadata on spacer sequences. eg., taxon_name)

=head3 If -database provided:

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine CLdb query (must start with 'AND').

=head3 Output

=item -outfmt  <char>

Output columns added to spacer-protospacer alignments.
See DESCRIPTION for more details.

=head3 Other

=item -array  <bool>

Write out hits to spacers in CRISPR arrays (instead of
hits to protospacers)? [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastGetPAM-SEED.pl

=head1 DESCRIPTION

Pull out the PAM regions for each hit to a protospacer. 
The output is a tab-delimited table (with header).

=head2 -PAM

The flag designates the region (relative to the
protospacer) containing the PAM. Negative 
values = upstream from the proto, while
positive values = downstream. 
So, the default [1 3] will select the 3bp
immediately downstream of the protospacer.
should not extend beyond the length of extensions
on the protospacer (default: 10bp)! 

=head2 -SEED

The flag designates the region (relative to
the protospacer) containing the SEED sequence.
Negative values bp from the END of the protospacer.
So, the default [-8] will select the last 8bp
of the protospacer.

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

### CLdb
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::GetAlign qw/ get_alignProto 
				   parse_outfmt/;
use CLdb::arrayBlast::PAM_SEED qw/ make_pam_index
				   make_seed_index /;
use CLdb::query qw/ table_exists
		    n_entries
		    join_query_opts
		    getLociSpacerInfo/;
use CLdb::utilities qw/ file_exists
			connect2db/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $verbose;
my $database_file;
my $outfmt;
my $array;
my @PAM = (1,3);
my @SEED = (-8);
my $query = "";
my (@subtype, @taxon_id, @taxon_name);
GetOptions(
	   "database=s" => \$database_file,
	   "array" => \$array,
	   "outfmt=s" => \$outfmt,
	   "PAM=i{2,2}" => \@PAM,
	   "SEED=i{1,2}" => \@SEED,
           "subtype=s{,}" => \@subtype,
           "taxon_id=s{,}" => \@taxon_id,
           "taxon_name=s{,}" => \@taxon_name,	   
	   "query=s" => \$query,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, 'database') if defined $database_file;
my $outfmt_r = parse_outfmt($outfmt);
my $pam_index_r = make_pam_index(\@PAM);
my $seed_index_r = make_seed_index(\@SEED);

#--- MAIN ---#
# decoding spacer and DR srl
my $spacer_r = decode_file( fh => \*STDIN );


# if database: connect & query
my $queries_r;   # 
if(defined $database_file){
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
get_PAM( blast => $spacer_r,
		outfmt => $outfmt_r,
		queries => $queries_r,
		array => $array,
		crDNA_ori => $crDNA_ori,
		verbose => $verbose );

# encoding
#print encode_sereal( $spacer_r );

