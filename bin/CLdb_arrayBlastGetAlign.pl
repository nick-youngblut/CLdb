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
(if getting metadata on spacer sequences. eg., taxon_name)

=item -outfmt  <char>

Output columns added to spacer-protospacer alignments.
The first 3 columns are always 'locus_id', 'spacer_id' 

=item array  <bool>

Write out alignments from the array? 

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastGetAlign.pl

=head1 DESCRIPTION

SCRIPTS TO RUN PRIOR TO THIS ONE: 
CLdb_arrayBlastAddcrRNA.pl &
CLdb_arrayBlastAddProto.pl 

Any query-hit without a crDNA 
sequence for the query and a protospacer 
sequence for the hit will be skipped.

The crDNA & protospacer will be aligned
with clustalw, and the extension sequences
flanking the protospacer will be set to lower
case, which will help for determining the PAM.

The alignment will be added to the *.srl
data structure.

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
use CLdb::arrayBlast::Align qw/ alignProto /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
GetOptions(
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#


#--- MAIN ---#
# decoding spacer and DR srl
my $spacer_r = decode_file( fh => \*STDIN );

# querying blastDBs for proteospacers
alignProto( blast => $spacer_r,
	    verbose => $verbose );

# encoding
print encode_sereal( $spacer_r );

