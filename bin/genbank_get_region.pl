#!/usr/bin/env perl

use strict;
use warnings;

=pod

=head1 NAME

genbank_get_region.pl -- Get all CDS features in >=1 specified genome region

=head1 SYNOPSIS

genbank_get_region.pl [options] -genbank file.gbk > region.txt

=head2 arguments

=item -genbank

Genbank file name.

=head2 options

=over

=item region_file 

region file name. Region file: comma-delimited; each line = contig, region_start, region_stop

=item -region  (<char> <int> <int>)...

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

my ($verbose, @regions, $tRNA, $strand_b, $genbank_file);
GetOptions(
	   "genbank=s" => \$genbank_file,
	   "regions=s{3,}" => \@regions,
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
elsif($ARGV[0]){ $regions_r = load_regions($ARGV[0]); }
else{ die " ERROR: provide at least 1 region (scaffold,start,end). Either via '-region' or as file\n"; }


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


  my @utags = sort keys %$all_tags_r;
  print join("\t", @utags), "\n";

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

### Subroutines
# sub load_regions{
# # loading file w/ regions of interest #
#   my ($infile) = @_;
#   open IN, $infile or die $!;
#   my @regions;
#   while(<IN>){
#     chomp;
#     s/#.+//;
#     next if /^\s*$/;
    
#     push(@regions, split /[ ,]+/);
#   }
#   return \@regions;
# }


# sub get_regions_scaffold{
# # getting CDS tag info that falls into regions #
# 	my ($itrees_r, $regions_r) = @_;
	
# 	my %feat_tags;
# 	my %tags;
# 	for(my $i=0; $i<=$#$regions_r; $i+=3){
# 		my $res;
# 		if( exists $itrees_r->{$$regions_r[$i]} ){
# 			if($$regions_r[$i+1] <= $$regions_r[$i+2]){			# if 'start' is <= 'end'
# 				$res = $itrees_r->{$$regions_r[$i]}->fetch($$regions_r[$i+1], $$regions_r[$i + 2]);
# 				}
# 			else{			# flippng start and end
# 				$res = $itrees_r->{$$regions_r[$i]}->fetch($$regions_r[$i+2], $$regions_r[$i+1]);
# 				}
# 			}
# 		else{ print STDERR " WARNING: $$regions_r[$i] not found in genbank!\n"; }
		
# 		my $feat_id = 0;
# 		foreach my $feat (@$res){
# 			$feat_id++;
# 			foreach my $tag ($feat->get_all_tags){
# 				$feat_tags{$feat_id}{$tag} = [$feat->get_tag_values($tag)];
# 				$tags{$tag} = 1;
# 				}
# 			$feat_tags{$feat_id}{"start"} = [($feat->location->start)];
# 			$feat_tags{$feat_id}{"end"} = [($feat->location->end)];
# 			$feat_tags{$feat_id}{"contig"} = [($feat->location->seq_id)];
# 			$feat_tags{$feat_id}{"strand"} = [($feat->location->strand)];
# 			}
# 		}

# 	# flipping start-end if strand = -1 #
# 	foreach my $feat (keys %feat_tags){
# 		if(${$feat_tags{$feat}{"strand"}}[0] == -1){
# 			(${$feat_tags{$feat}{"start"}}[0],${$feat_tags{$feat}{"end"}}[0]) =
# 				(${$feat_tags{$feat}{"end"}}[0], ${$feat_tags{$feat}{"start"}}[0]);
# 			}
# 		}	
	
# 	# writing out features #
# 	my @tags = ("contig", "start", "end");
# 	push(@tags, sort keys %tags);
# 	print join("\t", "Feature_num", @tags), "\n";
	
# 	foreach my $feat (keys %feat_tags){
# 		my @line;
# 		foreach my $tag (@tags){
# 			if(exists $feat_tags{$feat}{$tag}){
# 				push(@line, join("::", @{$feat_tags{$feat}{$tag}}));
# 				}
# 			else{
# 				push(@line, "NA");
# 				}
# 			}
# 		print join("\t", $feat, @line), "\n";
# 		}
# 	}
	
# sub get_regions{
# # getting CDS tag info that falls into regions #
# 	my ($itree, $regions_r) = @_;
	
# 	my %feat_tags;
# 	my %tags;
# 	for(my $i=0; $i<=$#$regions_r; $i+=2){
# 		my $res;
# 		if($$regions_r[$i] <= $$regions_r[$i+1]){			# if 'start' is <= 'end'
# 			$res = $itree->fetch($$regions_r[$i], $$regions_r[$i+1]);
# 			}
# 		else{			# flippng start and end
# 			$res = $itree->fetch($$regions_r[$i+1], $$regions_r[$i]);
# 			}
		
# 		my $feat_id = 0;
# 		foreach my $feat (@$res){
# 			$feat_id++;
# 			foreach my $tag ($feat->get_all_tags){
# 				$feat_tags{$feat_id}{$tag} = [$feat->get_tag_values($tag)];
# 				$tags{$tag} = 1;
# 				}
# 			$feat_tags{$feat_id}{"start"} = [($feat->location->start)];
# 			$feat_tags{$feat_id}{"end"} = [($feat->location->end)];
# 			$feat_tags{$feat_id}{"strand"} = [($feat->location->strand)];
# 			}
# 		}
	
# 	# flipping start-end if strand = -1 #
# 	foreach my $feat (keys %feat_tags){
# 		if(${$feat_tags{$feat}{"strand"}}[0] == -1){
# 			(${$feat_tags{$feat}{"start"}}[0],${$feat_tags{$feat}{"end"}}[0]) =
# 				(${$feat_tags{$feat}{"end"}}[0], ${$feat_tags{$feat}{"start"}}[0]);
# 			}
# 		}	
	
# 	# writing out features #
# 	my @tags = ("start", "end");
# 	push(@tags, sort keys %tags);
# 	print join("\t", "Feature_num", @tags), "\n";
	
# 	foreach my $feat (keys %feat_tags){
# 		my @line;
# 		foreach my $tag (@tags){
# 			if(exists $feat_tags{$feat}{$tag}){
# 				push(@line, join("::", @{$feat_tags{$feat}{$tag}}));
# 				}
# 			else{
# 				push(@line, "NA");
# 				}
# 			}
# 		print join("\t", $feat, @line), "\n";
# 		}
# 	}

# sub load_itrees{
# 	my ($feats_r) = @_;

# 	my %itrees;
# 	foreach my $feat (@$feats_r){
# 		my $scaf = $feat->location->seq_id;
		
# 		unless(exists $itrees{$scaf}){
# 			$itrees{$scaf} = Set::IntervalTree->new();
# 			}
			
# 		$itrees{$scaf} -> insert($feat, $feat->location->start, $feat->location->end);
# 		}
	
# 	return \%itrees;
# 	}
	
# sub load_itree{
# 	my ($feats_r) = @_;

# 	my $itree = Set::IntervalTree->new();
# 	foreach my $feat (@$feats_r){
# 		$itree -> insert($feat, $feat->location->start, $feat->location->end);
# 		}
	
# 	return $itree;
# 	}

# sub load_genbank_features{
# 	my $seqio = Bio::SeqIO->new(-fh => \*STDIN, -format => "genbank");
# 	my @feats;
# 	while ( my $seqo = $seqio->next_seq){
# 		push(@feats, grep { $_->primary_tag eq 'CDS' } $seqo->get_SeqFeatures);
# 		push(@feats, grep { $_->primary_tag eq 'tRNA' } $seqo->get_SeqFeatures) if $tRNA;
# 		}
		
# 		#print Dumper @feats; exit;
# 	return \@feats;
# 	}





