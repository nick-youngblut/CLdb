#!/usr/bin/env perl

=pod

=head1 NAME

get_SEEDstats -- getting mismatch stats on the seed region

=head1 SYNOPSIS

get_SEEDstats [flags] < proto.fasta

=head2 Required flags

=over


=back

=head2 Optional flags

=over

=item -prefix  <char>

Output file prefix. ['']

=item -SEED  <char>

start-stop of the SEED region. (2 values required) 
See DESCRIPTION for details. [-8 -1]

=item -gap  <bool>

Define a gaps in the alignments as mismatches? [TRUE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb -- arrayBlast --perldoc -- get_SEEDstats

=head1 DESCRIPTION

Sum up the mismatches in all
protospacer-crDNA alignments provided
(as fasta via STDIN). 

The mismatches are normalized by total
number of alignment positions 
('mismatch_norm' columns). This accounts
for varying length of alignments.

Mismatches are grouped by region:

=over

=item SEED: just the SEED region

=item nonSEED: positions in the protospacer that are not SEED

=item protospacer: the entire protospacer

=back

=head2 Output

2 output files are produced. Both are
tab-delimited tables (with headers). The prefix for
both file names is defined by -prefix.

=head3 *_sum.txt

Summing the number of mismatches by region for each
alignment.

=head3 *_byPos.txt

Summing the number of mismatches by position across
all alignments. The field 'pos_rel_SEED' is the
alignment position defined relative to the start 
of the SEED region (accounts for alignments with 
different lengths). The field 'pos_rel_aln'
is absolute position in the alignment.
The 'mismatch_norm' field is number of mismatches
divided by the number of characters in that alignment
position (accounts for missing characters in the aligment).

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

### CLdb
use CLdb::seq qw/read_fasta
		 revcomp/;
use CLdb::arrayBlast::PAM qw/ make_pam_index/;
use CLdb::arrayBlast::SEED qw/ read_proto_aln
			       parseProtoBySEED
			       write_sum_table
			       write_byPos_table/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $no_gap, $database);
my @SEED;
my ($fasta_in, $table_in); 
my $prefix = "";
GetOptions(
	   "SEED=i{2,2}" => \@SEED,	   
	   "prefix=s" => \$prefix,
	   "fasta=s" => \$fasta_in,
	   "table=s" => \$table_in,
	   "gap" => \$no_gap,
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# pam index
@SEED = (-8, -1) if ! defined $SEED[0] and ! defined $SEED[1];
my $seed_index_r = make_pam_index(\@SEED);  # actually making SEED index


#--- MAIN ---#
# getting fasta if fasta
my $fasta_r = read_proto_aln(fh => \*STDIN);

# revcomp all alignments to orient by protospacer
foreach my $pair (keys %$fasta_r){
  map{ revcomp($fasta_r->{$pair}{$_}) } keys %{$fasta_r->{$pair}};
}


# parsing proto-crDNA alignment by SEED
my ($mismatchSum_r, $mismatchByPos_r) = parseProtoBySEED( $fasta_r, 
							  $seed_index_r, 
							  $no_gap);

# output
## writing out summary table
write_sum_table($prefix, $mismatchSum_r);
## writing out count by position table
write_byPos_table($prefix, $mismatchByPos_r);
