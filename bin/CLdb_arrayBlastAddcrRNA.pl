#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddcrRNA.pl -- recreate crRNA extracting spacer & part of its adjacent DRs

=head1 SYNOPSIS

CLdb_arrayBlastAddcrRNA.pl [flags] spacer_blast.txt > spacer_blast_filtered.txt

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

perldoc CLdb_arrayBlastAddcrRNA.pl

=head1 DESCRIPTION

For each spacer, extracting spacer & adjacent 
sequence from its genome location in order
to make the crRNA (actually the crDNA sequence).
This can then be aligned to the protospacer region
in order to determine the PAM, protospacer,
and SEED sequence.

The spacer sequence ID is used to query
CLdb and get the genome position (and genome
file location). The sequence is then extracted 
from the genome.

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
use DBI;
use Sereal qw/ encode_sereal /;

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
				  detect_clustered_spacers
				  queryBySpacer
				  queryBySpacerCluster
				  groupByFastaFile
				  getSpacerRegion
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

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking for existence of tables
table_exists($dbh, 'loci');
table_exists($dbh, 'spacers');

# decoding spacer and DR srl
my $spacer_r = decode_file( fh => \*STDIN );

# getting spacer IDs
my $spacerIDs_r = get_query_IDs($spacer_r);

# detecting which spacers are clustered reps & which are not 
$spacerIDs_r = detect_clustered_spacers($spacerIDs_r);

# querying CLdb for info on spacer start-end & genome_file
print STDERR "Getting info on each spacer from CLdb...\n";
my %ret;
## single spacers
if( $spacerIDs_r->{single} ){
  $ret{single} = queryBySpacer($dbh, $spacerIDs_r->{single});
}
## cluster spacers
if( $spacerIDs_r->{cluster} ){
  table_exists($dbh, 'spacer_clusters');
  $ret{cluster} = queryBySpacerCluster($dbh, $spacerIDs_r->{cluster});
}

# adding spacer info back to blast hit srl DS
#addSpacerInfo

# grouping spacers by fasta_file entry
#my $spacers_byGenome = groupByFastaFile( \%ret );

# getting spacer region from each genome
#print STDERR "Getting spacer regions each blast db...\n";
#getSpacerRegion( #blast => $spacer_r, 
#		byGenome => $spacers_byGenome, 
#		CLdb_HOME => $CLdb_HOME, 
#		extension => $ext );

# adding to *srl data structure
## {BlastOutput_iterations}=>{Iteration}=>
##  [ 'Iteration_crRNA'=>{fasta_file}=>
##    {scaffold|strand|spacer_start|spacer_end|spacer_seq|
##     region_start|region_end|region_seq|ext} ]



# encoding
print encode_sereal( $spacer_r );

# disconnecting from CLdb
$dbh->disconnect;


#--- Subroutines ---#
