#!/usr/bin/env perl

=pod

=head1 NAME

srl2csv -- converting blast output (sereal format) to blast '-outfmt 6' or '-outfmt 7'

=head1 SYNOPSIS

srl2csv [flags] < blast_output.srl > blast_output.txt

=head2 Required flags

NONE

=head2 Optional flags

=over

=item -outfmt  <str>

blast fields as in '-outfmt' for blast. 
Only '6' or '7' formats are supported. 
['6 qseqid sseqid pident length mismatch 
gapopen qstart qend sstart send evalue bitscore']

=item -path  <bool>

Keep file path for blastdb field? [FALSE] 

=item -v

Verbose output. [TRUE]

=item -h

This help message

=back

=head2 For more information:

perldoc srl2csv

=head1 DESCRIPTION

Simple script for converting blast output in 
binary data serialization format (Sereal)
to blast tabular format: either
with comments '-outfmt 7' or without '-outfmt 6'

Currnently, only some of the standard blastn fields
are supported:

=head2 Supported fields:

=head3 Standard blast fields (see blastn help for more info)

=over

=item qseqid

=item qlen

=item sseqid

=item sacc

=item slen

=item qstart

=item qend

=item sstart

=item send

=item qseq

=item sseq

=item evalue

=item bitscore

=item score

=item length

=item pident

=item mismatch

=item gapopen

=item gaps

=item frames

=item qframe

=item sframe

=back

=head3 CLdb fields


=head4 #-- hsp --#

=over

=item subjectScaffold -- scaffold/chromo of blast the subject

=item subjectStrand -- strand of blast hit

=item protoFullStart -- start of protospacer (full length)

=item protoFullEnd -- end of protospacer (full length)

=item protoFullSeq -- protospacer sequence (full length)

=item protoFullXStart -- start of protospacer + extension

=item protoFullXEnd -- end of protospacer + extension

=item protoFullXSeq -- protospacer + extension sequence

=item protoX -- extension length (bp)

=item arrayHit -- spacer hits a CRISPR array? (1 = True)

=back


=head4 #-- crRNA (DNA) --#

=over

=item query_id -- unique ID for query

=item locus_id -- locus ID

=item spacer_id -- spacer ID (unique to locus ID)

=item cluster_id -- cluster ID (spacer cluster)

=item genome_fasta -- fasta file containing the spacer

=item spacer_scaffold -- scaffold/chromosome containing the spacer

=item spacer_start -- start position of spacer

=item spacer_end -- end position of spacer

=item spacer_seq -- spacer sequence

=item region_start -- start position of spacer + extension (crRNA)

=item region_end -- end position of spacer + extension (crRNA)

=item crDNA_seq -- crRNA sequence (DNA nucleotides)

=item array_sense_strand -- reading strand of crRNA

=item subtype -- CRISPR subtype classification

=item cas_status -- cas_status field from loci table

=item array_status -- array_status field from loci table

=back


=head3 Misc fields

=over

blastdb -- blast db file

=back

=head2 Undefined values

'undef' = value missing in blast srl file.

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use Sereal qw/ decode_sereal /;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../../lib";
use CLdb::arrayBlast::sereal qw/
				 parse_outfmt
				 classify_fields
				 blast_xml2txt
			       /;

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $keep_path);
my $outfmt = '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore'; 
GetOptions(
	   "outfmt=s" => \$outfmt,
	   "path" => \$keep_path,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide '-outfmt'\n"
  unless defined $outfmt;
my $fields_r = parse_outfmt($outfmt); 
classify_fields($fields_r); 

#--- MAIN ---#
# loading serealized blast output
print STDERR "Decoding the srl file...\n" unless $verbose;
my $srl;
$srl .= $_ while <>;
my $decoder = Sereal::Decoder->new();
my $blast_r =  $decoder->decode( $srl );

# making table
foreach my $blast_run (keys %$blast_r){

  blast_xml2txt(blast => $blast_r->{$blast_run}, 
		fields => $fields_r,
		keep_path => $keep_path);
}

