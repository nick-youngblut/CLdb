#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddProto.pl -- adding protospacer to blast srl file

=head1 SYNOPSIS

CLdb_arrayBlastAddProto.pl [flags] < blast_hits.srl > spacer_blast_proto.srl

=head2 Required flags

=over

=back

=head2 Optional flags

=over

=item -extension  <int>

Number of bp to include on either side of the protospacer. [10].

=item -fork  <int>

Number of parallel blast db queriesa. [1]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddProto.pl

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
use CLdb::utilities qw/
			file_exists
			get_file_path
		      /;
use CLdb::query qw/ table_exists /;
use CLdb::arrayBlast::sereal qw/ decode_file /;
use CLdb::arrayBlast::Proto qw/ queryBlastDBs /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $ext = 10;
my $fork = 0;
GetOptions(
	   "extension=i" => \$ext,
	   "fork=i" => \$fork,
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
	       fork => $fork);

# encoding
print encode_sereal( $spacer_r );

