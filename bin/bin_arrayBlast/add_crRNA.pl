#!/usr/bin/env perl

=pod

=head1 NAME

add_crRNA -- recreate crRNA extracting spacer & part of its adjacent DRs

=head1 SYNOPSIS

add_crRNA [flags] < blast_hits.srl > blast_hits_crRNA.srl

=head2 Required flags

=over

=item -database  <char>

CLdb database file

=back

=head2 Optional flags

=over

=item -extension  <int>

Number of bp to include on either side of spacer. [10].

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb_perldoc add_crRNA

=head1 DESCRIPTION

For each spacer, extracting spacer & adjacent 
sequence (defined by '-extension') 
from its genome location in order
to make the DNA template of the crRNA.
Note: this script finds the DNA template of the crRNA, 
not the crRNA itself. This is for easier downstream analysis.

The DNA template of the crRNA can then be
aligned to the protospacer region
in order to determine the PAM, protospacer,
and SEED sequence.

The spacer sequence ID is used to query
CLdb and get the genome position (and genome
fasta file location). 
The 'crRNA' sequence is then extracted from the genome
and added to the blast hit srl file.

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
use Sereal qw/ encode_sereal /;
use FindBin;
use lib "$FindBin::RealBin/../../lib/";

### CLdb
use CLdb::utilities qw/
			file_exists
			connect2db
			get_file_path
		      /;
use CLdb::query qw/ table_exists /;
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::AddcrRNA qw/
				  get_query_IDs
				  detectClusteredSpacers
				  queryBySpacer
				  queryBySpacerCluster
				  getSpacerRegion
				  addcrDNAtoBlast
				 /;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my $ext = 10;
GetOptions(
	   "database=s" => \$database_file,
	   "extension=i" => \$ext,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
my $CLdb_HOME = get_file_path($database_file);
die "'-extension' must be an integer >=0"
  unless $ext >= 0;


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking for existence of tables
table_exists($dbh, 'loci');
table_exists($dbh, 'spacers');

# decoding spacer and DR srl
print STDERR "Decoding .srl file...\n" unless $verbose;
my $spacer_r = decode_file( fh => \*STDIN );

# getting spacer IDs (either sequence or rep sequence name)
my $spacerIDs_r = get_query_IDs($spacer_r);

# detecting which spacers are clustered reps & which are not 
$spacerIDs_r = detectClusteredSpacers($spacerIDs_r);

# querying CLdb for info on spacer start-end & genome_file
print STDERR "Getting info from CLdb on each blast query spacer...\n" 
  unless $verbose;
my %info;
## single spacers
if( $spacerIDs_r->{single} ){
  # TODO: test this subroutine
  queryBySpacer($dbh, $spacerIDs_r->{single});
}
## cluster spacers
if( $spacerIDs_r->{cluster} ){
  table_exists($dbh, 'spacer_clusters');
  queryBySpacerCluster($dbh, 
		       $spacerIDs_r->{cluster}, 
		       \%info,
		       -verbose => $verbose);
}


# getting spacer region from each genome
print STDERR "Getting spacer regions from each genome...\n";
my $byQuery_r = getSpacerRegion(
		CLdb => \%info, 
		CLdb_HOME => $CLdb_HOME, 
		extension => $ext);

# adding crRNA info to *srl data structure
addcrDNAtoBlast($spacer_r, $byQuery_r);

# encoding
print encode_sereal($spacer_r);

# disconnecting from CLdb
$dbh->disconnect;

