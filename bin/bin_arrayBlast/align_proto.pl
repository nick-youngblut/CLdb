#!/usr/bin/env perl

=pod

=head1 NAME

align_proto -- aligning protospacer & crDNA 

=head1 SYNOPSIS

align_proto [flags] < blast_hits_crDNA_proto.srl > spacer_blast_crDNA_proto_align.srl

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb_perldoc align_proto

=head1 DESCRIPTION

SCRIPTS TO RUN PRIOR TO THIS ONE: 
arrayBlastAddcrRNA.pl &
arrayBlastAddProto.pl 

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
use FindBin;
use lib "$FindBin::RealBin/../../lib/";

### CLdb
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::Align qw/ alignProto /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database);
GetOptions(
	   "database=s" => \$databse, #unused
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

