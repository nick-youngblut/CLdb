#!/usr/bin/env perl

use strict;
use warnings;

=pod

=head1 NAME

genbank_get_region.pl -- Get all CDS features in >=1 specified genome region

=head1 SYNOPSIS

genbank_get_region.pl [options] -genbank > region.txt

=head2 required arguments

=item -genbank

Genbank file name.

=head2 options

=over

=item -region

region file name. Region file: comma-delimited; each line = contig, region_start, region_stop

=item -list_regions  (<char> <int> <int>)...

List of region contig,region_start,region_stop values (use instead of region_file). 
>=1 region can be specified in this manner.
Example: '-r contig_1 1 1000' for bp 1 to bp 1000 of 'contig_1'

=item -tRNA  <bool>

Include rRNA features? [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc genbank_get_region.pl

=head1 DESCRIPTION

Get all CDS features in >=1 specified genome region.
Regions can be provided either by '-region' flag or
by a region file.

=head2 Output

A tab-delimited table.

The tag values of all features falling at least 
partially with a region are written out in a table.
"NA" is written if the tag value is not present in
the feature. Multiple values for the same tag are 
combined with "::"

The 1st 3 columns of the table are the user-define regions
(column names start with 'region_*'). The rest are values
associated with the features.

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/annotation/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Set::IntervalTree;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::genbank::genbank_get_region qw/genbank_get_regions
					 load_regions/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, @regions, $tRNA, $strand_b, $genbank_file, $region_file);
GetOptions(
	   "genbank=s" => \$genbank_file,
	   "region=s" => \$region_file,
	   "list_regions=s{3,}" => \@regions,
	   "tRNA" => \$tRNA,
	   "strand" => \$strand_b,			# all + strand? [false]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults

### MAIN
# loading/editing region(s) #
my $regions_r;
if(@regions){ $regions_r = load_region_args(\@regions); }
elsif($region_file){ $regions_r = load_regions($region_file); }
else{ die " ERROR: provide at least 1 region (scaffold,start,end). Either via '-region' or '-list_regions'\n"; }


# getting region features
my $all_feats_r = genbank_get_regions($regions_r, $genbank_file);


# getting hash of all tags (finding uniques)
my $all_tags_r = get_all_unique_tags($all_feats_r);

# writing
write_feat_tags($all_feats_r, $all_tags_r); 



#-- Subroutines --#
sub write_feat_tags{
  my $all_feats_r = shift or die "Provide all_feats object";
  my $all_tags_r = shift or die "Provide all_tags object";

  # unique tags; header
  my @utags = sort keys %$all_tags_r;
  print join("\t", 'region_contig', 'region_start', 'region_end', @utags), "\n";

  # each feature
  foreach my $region_r (@$all_feats_r){
    foreach my $feat (keys %{$region_r->{features}}){
      my @line = @{$region_r->{region}};

      foreach my $tag (@utags){
	if (exists $region_r->{features}{$feat}{$tag}){
	  if (ref $region_r->{features}{$feat}{$tag} eq 'ARRAY'){
	    push @line, join("::", @{$region_r->{features}{$feat}{$tag}});
	  }
	  else{
	    push @line, $region_r->{features}{$feat}{$tag};
	  }
	}
	else{
	  push @line, 'NA';
	}
      }
      
      print join("\t", @line), "\n";
    }
  }
}

sub get_all_unique_tags{
# getting all tag features from feature-tag object
  my $all_feats_r = shift or die "Provide ref: [{feature=>feat=>tag=>value}]";


  my %tags;
  foreach my $region_r (@$all_feats_r){
    foreach my $feat (keys %{$region_r->{features}}){      
      map{$tags{$_} = 1} keys %{$region_r->{features}{$feat}};
    }
  }
  return \%tags;
}

