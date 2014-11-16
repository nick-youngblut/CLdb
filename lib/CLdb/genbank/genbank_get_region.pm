#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Set::IntervalTree;
use Regexp::Common;
use Carp  qw( carp confess croak );

# export #
use base 'Exporter';
our @EXPORT_OK = qw/
genbank_get_region/;


=head2 genbank_get_regions

Get features from multiple regions of genbank.

=head3 Input

$regions_r --  regions in @@ format
$genbank_file -- genbank file
$tRNA -- include tRNAs? [false]
$strand_b -- all on + strand? [false]

=head3 Return
[{region:[], feature:[]}, ...]

=cut

sub genbank_get_regions{
  my $regions_r = shift or die "Provide a regions object";
  my $genbank_file = shift or die "Provide a genbank file";
  my $tRNA = shift;
  my $strand_b = shift;
  my $verbose = shift;

  
  # loading features
  my $feats_r = load_genbank_features($genbank_file,$tRNA);    
  my $itrees_r = load_itrees($feats_r);

  # getting regions
  my @all_feats;
  for my $region_r (@$regions_r){
    unless (scalar @$region_r == 3){
      confess "ERROR: region must contain 'scaffold, start, stop'\n"
    }

    my ($feat_tags_r, $tags_r) = get_regions_scaffold($itrees_r, $region_r, $strand_b, 0, 1);
    join_tag_values($feat_tags_r, $tags_r);
    push @all_feats, {region => $region_r, features => $feat_tags_r};
  }


  #print Dumper @all_feats; exit;
  return \@all_feats;
}



sub genbank_get_region{
#--- Description ---#
# getting features in one select region of genbank
#--- Input ---#
# @regions = scaffold,start,end
# $genbank_file = genbank file
# $tRNA = include tRNAs? [false]
# $strand_b = all on + strand? [false]
  my @regions = @_[0..2];		# start, end, scaffold
  my ($genbank_file, $tRNA, $strand_b, $verbose) = @_[3..$#_];
  
  # I/O check #
  my $regions_r;
  if(@regions){ $regions_r = \@regions; }
  else{ confess " ERROR: provide a region (scaffold,start,end)\n"; }
  confess " ERROR: you must provide: scaffold ,start, & end for the region!\n"
    unless scalar @$regions_r == 3;
  
  # loading features
  my $feats_r = load_genbank_features($genbank_file,$tRNA);  
  my $itrees_r = load_itrees($feats_r);

  # getting regions
  my ($feat_tags_r, $tags_r) = get_regions_scaffold($itrees_r, $regions_r, $strand_b);
  join_tag_values($feat_tags_r, $tags_r);
  #my $ret_r = feat_tag_table($feat_tags_r, $tags_r);

  #print Dumper $feat_tags_r; exit;
  #return $ret_r;
  return $feat_tags_r
}


#--- accessory subroutines ---#
sub get_regions_scaffold{
# getting CDS tag info that falls into regions #
# Input:
## regions_r = @{scaf, start, end}
## add_tags_b -- if True, feature tags will be 1st row of output list of lists
# Output:
## list of lists
## if tag  not found for feature, using 'NA'

  my $itrees_r = shift or die "Provide hashref of interval trees\n";
  my $regions_r = shift or die "Provide regions object\n";
  my $strand_b = shift;
  
  my %feat_tags;
  my %tags;
  my $res;
  if( exists $itrees_r->{$$regions_r[0]} ){
    if($$regions_r[1] <= $$regions_r[2]){			# if 'start' is <= 'end'
      $res = $itrees_r->{$$regions_r[0]}->fetch($$regions_r[1], $$regions_r[2]);
    }
    else{			# flippng start and end
      $res = $itrees_r->{$$regions_r[0]}->fetch($$regions_r[2], $$regions_r[1]);
    }
  }
  else{ warn "WARNING: '$$regions_r[0]' not found in".
	  "genbank or no features for that scaffold in genbank!\n"; }

  my $feat_id = 0;
  foreach my $feat (@$res){
    $feat_id++;
    foreach my $tag ($feat->get_all_tags){
      $feat_tags{$feat_id}{$tag} = [$feat->get_tag_values($tag)];
      $tags{$tag} = 1;
    }
    $feat_tags{$feat_id}{"start"} = [($feat->location->start)];
    $feat_tags{$feat_id}{"end"} = [($feat->location->end)];
    $feat_tags{$feat_id}{"contig"} = [($feat->location->seq_id)];
    $feat_tags{$feat_id}{"strand"} = [($feat->location->strand)];
  }
  
  # flipping start-end if strand = -1 #
  unless($strand_b){
    foreach my $feat (keys %feat_tags){
      if(${$feat_tags{$feat}{"strand"}}[0] == -1){
	(${$feat_tags{$feat}{"start"}}[0],${$feat_tags{$feat}{"end"}}[0]) =
	  (${$feat_tags{$feat}{"end"}}[0], ${$feat_tags{$feat}{"start"}}[0]);
      }
    }	
  }

  return \%feat_tags, \%tags; 
}


sub feat_tag_table{
  # create a feature-tag table (list of lists)
  my $feat_tags_r = shift or die "Provide feat_tags hashref";
  my $tags_r = shift or die "Provide tags hashref";
  my $add_tags_b = shift;   # adding tags (header) to output

  # creating list of features to return
  ## if tag  not found for feature, using 'NA'
  my @ret;
  my @tags = ("contig", "start", "end");
  push(@tags, sort keys %$tags_r);
  
  if ($add_tags_b){ # tags added as first row in table
    push @ret, ["Feature_num", @tags];
  }
  
  foreach my $feat (keys %$feat_tags_r){
    my @line;
    foreach my $tag (@tags){
      # values for tag joined with '::'
      if(exists $feat_tags_r->{$feat}{$tag}){
	push(@line, join("::", @{$feat_tags_r->{$feat}{$tag}}));
      }
      # if tag not in feature, return 'NA'
      else{
	push(@line, "NA");
      }
    }
    push @ret, [$feat, @line];
  }
  
  return \@ret;
}



sub join_tag_values{
# joining tag values (defaut: '::')
  my $feat_tags_r = shift or die "Provide feat_tags hash";
  my $tags_r = shift or die "Provide tags hashref";
  my $join_val = shift // '::';
  

  # creating list of features to return
  ## if tag  not found for feature, using 'NA'
  my @tags = ("contig", "start", "end");
  push(@tags, sort keys %$tags_r);
  
  
  #my %ret;
  foreach my $feat (keys %$feat_tags_r){
    my @line;
    foreach my $tag (@tags){
      # values for tag joined with '::'
      if(exists $feat_tags_r->{$feat}{$tag}){
	$feat_tags_r->{$feat}{$tag} =  join("::", @{$feat_tags_r->{$feat}{$tag}});
      }
      # if tag not in feature, return 'NA'
      else{
	$feat_tags_r->{$feat}{$tag} = 'NA';
      }
    }
  }
}



=head2 get_regions

Args:
itree -- interval tree object
regions -- regions object
strand -- strand of region

=cut

sub get_regions{
  # getting CDS tag info that falls into regions #
  my ($itree, $regions_r, $strand_b) = @_;
  
  my %feat_tags;
  my %tags;
  for(my $i=0; $i<=$#$regions_r; $i+=2){
    my $res;
    if($$regions_r[$i] <= $$regions_r[$i+1]){			# if 'start' is <= 'end'
      $res = $itree->fetch($$regions_r[$i], $$regions_r[$i+1]);
    }
    else{			# flippng start and end
      $res = $itree->fetch($$regions_r[$i+1], $$regions_r[$i]);
    }
    
    my $feat_id = 0;
    foreach my $feat (@$res){
      $feat_id++;
      foreach my $tag ($feat->get_all_tags){
	$feat_tags{$feat_id}{$tag} = [$feat->get_tag_values($tag)];
	$tags{$tag} = 1;
      }
      $feat_tags{$feat_id}{"start"} = [($feat->location->start)];
      $feat_tags{$feat_id}{"end"} = [($feat->location->end)];
      $feat_tags{$feat_id}{"strand"} = [($feat->location->strand)];
    }
  }
  
  # flipping start-end if strand = -1 #
  unless($strand_b){
    foreach my $feat (keys %feat_tags){
      if(${$feat_tags{$feat}{"strand"}}[0] == -1){
	(${$feat_tags{$feat}{"start"}}[0],${$feat_tags{$feat}{"end"}}[0]) =
	  (${$feat_tags{$feat}{"end"}}[0], ${$feat_tags{$feat}{"start"}}[0]);
      }
    }	
  }
  
  # writing out features #
  my @ret;
  my @tags = ("start", "end");
  push(@tags, sort keys %tags);
  push @ret, ["Feature_num", @tags];
  
  foreach my $feat (keys %feat_tags){
    my @line;
    foreach my $tag (@tags){
      if(exists $feat_tags{$feat}{$tag}){
	push(@line, join("::", @{$feat_tags{$feat}{$tag}}));
      }
      else{
	push(@line, "NA");
      }
    }
    push @ret, [$feat, @line];
  }
  return \@ret;
}


=head2 load_itrees

Loading hash of interval trees with feature start-end info.
{scaffold}=>itree

=cut

sub load_itrees{
  my ($feats_r) = @_;
  
  my %itrees;
  foreach my $feat (@$feats_r){
    my $scaf = $feat->location->seq_id;
    
    # sanity checks
    unless( defined $scaf ){
      print STDERR " WARNING: no scaffold ID found for feature. Skipping feature\n";
      next;
    }	
    if( $feat->location->start !~ /$RE{num}{real}/ or
	$feat->location->start !~ /$RE{num}{real}/ ){
	print STDERR " WARNING: feature start-end not numeric. Skipping feature\n";
	next;
      }	  

    # loading itree
    $itrees{$scaf} = Set::IntervalTree->new() unless exists $itrees{$scaf};	
    $itrees{$scaf} -> insert($feat, $feat->location->start, $feat->location->end);
  }
  
  return \%itrees;
}

	
sub load_itree{
# loading interval tree of genbank features
# Input: list of genbank feature objects
  my ($feats_r) = @_;  
  
  my $itree = Set::IntervalTree->new();
  foreach my $feat (@$feats_r){
    $itree -> insert($feat, $feat->location->start, $feat->location->end);
  }
  
  return $itree;
}


sub load_genbank_features{
  my ($genbank_file, $tRNA) = @_;

  # check that $genbank_file exists
  die "ERROR: cannot find '$genbank_file'\n"
    if ! -e $genbank_file or -d $genbank_file;

  # loading genbank as Bio::SeqIO object
  my $seqio = Bio::SeqIO->new(-file => $genbank_file, -format => "genbank");
  my @feats;
  while ( my $seqo = $seqio->next_seq){
    push(@feats, grep { $_->primary_tag eq 'CDS' } $seqo->get_SeqFeatures);
    push(@feats, grep { $_->primary_tag eq 'tRNA' } $seqo->get_SeqFeatures) if $tRNA;
  }
  
  #print Dumper @feats; exit;
  return \@feats;
}


sub load_regions{
# loading file w/ regions of interest #
  my ($infile) = @_;
  open IN, $infile or confess $!;
  my @regions;
  while(<IN>){
    chomp;
    s/#.+//;
    next if /^\s*$/;

    my @line = split /\t/;
    die "ERROR: row $. of region file does not contain >=3 columns\n"
      unless scalar @line >= 3;
    push(@regions, [@line[0..2]]);
  }
  close IN or confess $!;

  return \@regions;
}

sub load_region_args{
  # loading regions passed by args
  my $regions_r = shift or die "Provide regions as list ref\n";

  # assertions
  die "ERROR: regions must be sets of 'scaffold,start,end'\n"
    unless scalar(@$regions_r) % 3 == 0;

  my @rr;
  for (my $i=0; $i<scalar(@$regions_r); $i+=3){
    push @rr, [@$regions_r[$i..$i+2]];
  }

  return \@rr;
}

