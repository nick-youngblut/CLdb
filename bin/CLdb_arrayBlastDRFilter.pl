#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastDRFilter.pl -- filter out spacers hitting CRISPR arrays based on adjacency to DR blast hits

=head1 SYNOPSIS

CLdb_arrayBlastDRFilter.pl [options] spacer_blast.txt DR_blast.txt > spacer_blast_filtered.txt

=head2 options

=over

=item -range  <int>

Range allowable between spacer & DR blast hit (bp). [30]

=item -DR  <int>

A spacer blast is considered in an array if number of adjacent DR blasts is >= '-DR'. [1]

=item -length  <float>

DR blast hits must be >= fraction of total DR sequence. [0.66]

=item -evalue  <float>

DR blast hits must be < 'evalue'. [10] 

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_arrayBlastDRFilter.pl

=head1 DESCRIPTION

Filter out spacer blast hits to CRISPR arrays by
removing all spacer blast hits that hit adjacent
to >= '-DR' direct repeat blast hits.

=head1 EXAMPLES

=head2 Basic Usage:

CLdb_arrayBlastDRFilter.pl spacer_blast.txt repeat_blast.txt > spacer_blast_filter.txt

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
use CLdb::arrayBlast::DRfilter qw/
				   decode_file
				   make_DR_itree
				   DR_filter_blast
				 /;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $evalue_cut = 10;
my $len_cut = 0.66;
my $DR_cnt = 1;
my $range = 30;
GetOptions(
	   "range=i" => \$range,
	   "length=f" => \$len_cut,			# length cutoff of DR hits
	   "evalue=f" => \$evalue_cut,			# evalue cutoff of DR hits
	   "DR=i" => \$DR_cnt, 					# number of adjacent DR hits needed to call 'array'
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die "ERROR: provide a spacer blast file\n"
  unless defined $ARGV[0];
die "ERROR: provide a DR blast file\n"
  unless defined $ARGV[1];
map{ die "ERROR: cannot find '$_'\n" unless -e $_ } @ARGV[0..1];

### MAIN
# decoding spacer and DR srl
my $spacer_r = decode_file($ARGV[0]);
my $DR_r = decode_file($ARGV[1]);

# DR blast hit interval tree
my $itrees_r = make_DR_itree( $DR_r,
			      {range => $range, 
			       len_cut => $len_cut,
			       evalue_cut => $evalue_cut
			      });

# adding 'array_hit' to spacer hash
DR_filter_blast( $spacer_r, $itrees_r, 
		{DR_cnt => $DR_cnt}
	       );

# encoding
print encode_sereal( $spacer_r );



#--- Subroutines ---#
