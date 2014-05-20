#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastGetPAM.pl -- getting PAM for each protospacer

=head1 SYNOPSIS

CLdb_arrayBlastGetPAM.pl [flags] < proto.fasta > PAMs.fasta

=head2 Required flags

=over


=back

=head2 Optional flags

=over

=item -SEED  <char>

start-stop of the SEED region. (2 values required) 
See DESCRIPTION for details. [-8 -1]

=item -revcomp  <bool>

Reverse complement protospacer sequence
before extracting PAM (if protospacer
is not oriented on the correct strand). [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastGetPAM.pl

=head1 DESCRIPTION

=head2 -SEED

The flag designates the region of the protospacer
that contains the SEED sequence. The values
indicate the start-end of the SEED sequence 
compared to the start of the protospacer
(if values are > 0) or from the end of the
protospacer (if values are < 0). 
Negative values bp from the END of the protospacer.
So, the default [-8 -1] will select the last 8bp
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


### CLdb
use CLdb::arrayBlast::PAM qw/ make_pam_index/;
use CLdb::arrayBlast::SEED qw/ read_proto_aln
			       parseProtoBySEED/;
use CLdb::seq qw/read_fasta/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $revcomp_b);
my @SEED;
my ($fasta_in, $table_in); 
GetOptions(
	   "SEED=i{2,2}" => \@SEED,	   
	   "fasta=s" => \$fasta_in,
	   "table=s" => \$table_in,
	   "revcomp" => \$revcomp_b,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# pam index
@SEED = (-8, -1) unless @SEED;
my $seed_index_r = make_pam_index(\@SEED);  # actually making SEED index


#--- MAIN ---#
# getting fasta if fasta
my $fasta_r = read_proto_aln(fh => \*STDIN);

# parsing proto-crDNA alignment by SEED
my $aln_r = parseProtoBySEED( $fasta_r, $seed_index_r);


