#!/usr/bin/env perl

=pod

=head1 NAME

get_proto -- getting protospacer from blast srl (assuming no crDNA info)

=head1 SYNOPSIS

get_proto [flags] < blast_hits.srl > spacer_blast_proto.txt

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

perldoc get_proto

=head1 DESCRIPTION

Run arrayBlastAddProto.pl 
prior to this script!

Extracting the protospacer sequence
(extension from proto lower case)
and writing as a tab-delimited table (with header).
Protospacer sequences are written to their orientation
(ie., sequence start = 5' of protospacer; 
sequence end = 3')

This script can be used without a CLdb.
Just blast your sequences, run arrayBlastAddProto.pl,
then this script.

The output is a tab-delimited table (with header).
Convert this to a fasta with 'table2fasta.pl'

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
use CLdb::arrayBlast::Proto qw/ getProto /;


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
getProto( blast => $spacer_r,  
	  verbose => $verbose);


