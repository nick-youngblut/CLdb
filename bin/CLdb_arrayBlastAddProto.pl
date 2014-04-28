#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddProto.pl -- getting the protospacers for spacer blast hits (& adjacent regions)

=head1 SYNOPSIS

CLdb_arrayBlastAddProto.pl [flags] < spacerBlast.txt > spacerBlast_proto.txt

=head2 flags

=over

=item -x  <int>

Extension beyond spacer blast to check for PAMs (bp). [10]

=item -length  <float>

Length cutoff for blast hit (>=; fraction of spacer length). [0.66]

=item -pam  <int>

Columns containing just he supposed PAM region. This is designated with 4 values to determine 
5' & 3' PAM region adjacent to the protospacer. Example: -pam -3 -1 1 3 will get the 3 bp 
adjacent up and downstream of the protospacer. [-3 -1 1 3] 

=item -revcomp  <bool>

Protospacer as opposite strand of spacer blast hit? [FALSE] (default: protospacer oriented to subject + strand)

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>	

This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddProto.pl

=head1 DESCRIPTION

Get sequences of protospacers (found via 
blast hits) plus adjacent nucleotides
to look for PAMs.

The input is a blast table with comment lines ('-outfmt 7').
The output is the same blast table with appended fields & 
spacers not meeting the cutoffs removed.

Adjacent regions to the protospacer ('proto_seq_fullx' field)
are in lower case, while the protospacer is upper case.

By default, the protospacer is oriented to subject + strand.

=head2 Extra fields appended to the spacer blast table:

=over

=item * "proto_seq_rel2query" = protospacer sequence (just blast match); oriented to query + strand

=item * "proto_seq_rel2subject" = protospacer sequence (just blast match); oriented to subject + strand 

=item * "q_start_full" = query start (full length spacer)

=item * "q_end_full" = query end (full length spacer)

=item * "s_start_full" = protospacer start (full length spacer hit)

=item * "s_end_full" = protospacer end (full length spacer hit)

=item * "proto_seq_full" = protospacer sequence (full length spacer hit) 

=item * "s_start_fullx" = protospacer sequence start (full length spacer hit) & extension 

=item * "s_end_fullx" = protospacer sequence end (full length spacer hit) & extension 

=item * "proto_seq_fullx" = protospacer sequence (full length spacer hit) & extension 

=item * "x_length" = length of extension beyond full protospacer sequence (each side), ('-x' flag)

=item * "5p_pam_region" = 5' pam sequence (by subject + strand). ('-pam' flag)

=item * "3p_pam_region" = 3' pam sequence (by subject + strand). ('-pam' flag)

=back 

If using '-revcomp', the protospacer sequences will be 
the reverse-complement of spacer blast hits.
Their start-end values should reflect this.

If a spacer blast hits the end of a scaffold/chromosome, 
extending the spacer hit to the full spacer length ("proto_seq_full")
or protospacer+extension length ("proto_seq_fullx") will be 
limited to the end of the genomic sequence.

=head2 WARNING

The spacer blast table must have comment lines (blastn -outfmt 7)!

The blast databases must be in the location specified
in the "# Database:" comment lines! 

=head1 EXAMPLES

=head2 Protospacers for all hits (by spacer group)

CLdb_arrayBlastAddProto.pl < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

=head2 Protospacers for just full length spacer blast hits

CLdb_arrayBlastAddProto.pl -l 1 < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

=head2 Protospacers with 5 bp extensions 

CLdb_arrayBlastAddProto.pl -x 5 < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

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
use Bio::SeqIO;
use List::Util qw/max/;

# CLdb libs #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::blast qw/
		    parse_blast_hits/;
use CLdb::seq qw/
		  revcomp/;
use CLdb::arrayBlast::addProto qw/
				   get_proto_seq
				   add_new_fields
				   write_editted_blast
				 /;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $full_len, $revcomp_b, @pam);
my $extend = 10;							# length to extend off of each side
my $len_cut = 0.66;							# cutoff for length of blast hit relative to query length
GetOptions(
	   "x=i" => \$extend,
	   "length=f" => \$len_cut,				# min length of spacer hit (fraction of total)
	   "revcomp" => \$revcomp_b,			# reverse complement protospacer relative to blast hit? [FALSE]
	   "pam=i{4,4}" => \@pam, 				# pam region to write out
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if(@pam){
	map{die " ERROR: the pam start-end should be negative values!\n" unless $_ < 0} @pam[0..1];
	map{die " ERROR: the pam start-end should be positive values!\n" unless $_ > 0} @pam[2..3];
	die " ERROR: the 5' pam start should be <= pam end (default: -3 -1)\n" unless $pam[0] <= $pam[1];
	die " ERROR: the 3' pam start should be <= pam end (default: 1 3)\n" unless $pam[2] <= $pam[3];
	}
else{ @pam = (-3,-1,1,3); }

### main
# parsing blast table #
my ($blast_r, $fields_r) = parse_blast_hits();

# getting protospacer & extensions #
## making list of protospacer/extension start-stop, & strand ##
my $new_fields_r = get_proto_seq($blast_r, $fields_r, $extend, 
				 $len_cut, \@pam, $verbose);

# writing editted table #
## adding new fields ##
add_new_fields($fields_r, $new_fields_r);

#print Dumper $blast_r; exit;

## writing table ##
write_editted_blast($blast_r, $fields_r);


#--- Subroutines ---#
sub load_fasta{
  # loading fasta file as a hash #
  my $fasta_in = shift;
  open IN, $fasta_in or die $!;
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;	
    if(/>.+/){
      $_ =~ s/^>//;
      $fasta{$_} = "";
      $tmpkey = $_;	# changing key
    }
    else{$fasta{$tmpkey} .= $_; }
  }
  close IN;
  #print Dumper %fasta; exit;
  return \%fasta;
} #end load_fasta


