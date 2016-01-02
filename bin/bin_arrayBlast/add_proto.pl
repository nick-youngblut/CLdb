#!/usr/bin/env perl

=pod

=head1 NAME

add_proto -- adding protospacer to blast srl file

=head1 SYNOPSIS

add_proto [flags] < blast_hits.srl > spacer_blast_proto.srl

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -extension  <int>

Number of bp to include on either side of the protospacer. [10].

=item -workers  <int>

Max cpus to utilize. [1]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb -- arrayBlast --perldoc -- add_proto

=head1 DESCRIPTION

For each spacer blast hit, extracting 
the full length protospacer & adjacent 
sequence (defined by '-extension') 
from the corresponding blast database.

The protospacer can then be aligned to the spacer 
crRNA region in order to determine 
the PAM, protospacer, and SEED sequence.

'full length protospacer' means that
any partial spacer blast hits are 
extended to the full length of the
spacer query sequence. The extension sequences
(defined by '-extension') extend from the
full length protospacer.

=head1 DEPENDENCIES

blastdbcmd (part of the BLAST+ Toolkit)

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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
use FindBin;
use lib "$FindBin::RealBin/../../lib/";
use Sereal qw/ encode_sereal /;

### CLdb
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::Proto qw/ queryBlastDBs /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database);
my $ext = 10;
my $workers = 1;
GetOptions(
	   "extension=i" => \$ext,
	   "workers=i" => \$workers,
	   "database=s" => \$database,  # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#


#--- MAIN ---#
# decoding spacer and DR srl
my $spacer_r = decode_file( fh => \*STDIN );

# querying blastDBs for proteospacers
queryBlastDBs( blast => $spacer_r,  
	       extension => $ext, 
	       verbose => $verbose,
	       workers => $workers );

# encoding
print encode_sereal( $spacer_r );

