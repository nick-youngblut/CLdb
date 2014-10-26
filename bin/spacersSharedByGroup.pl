#!/usr/bin/env perl

=pod

=head1 NAME

spacersSharedByGroup.pl -- determine spacers specific to groups 

=head1 SYNOPSIS

spacersSharedByGroup.pl [flags] < spacers_shared.txt > summary.txt

=head2 Required flags

=over

=item group_file  <char>

File designating how the columns in the spacers_shared.txt file should be grouped
(2 column, tab-delimited: spacers_shared_column     group_id).

=back

=head2 Optional flags

=over

=item -percent  <float>

Greater than '-percent' of members in group that must contain the spacer. [50]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc spacersSharedByGroup.pl

=head1 DESCRIPTION

How many spacers are specific to a certain group?

'specific' means found in a group but not found in 
any other groups. 

The output is a count of spacers specific to
particular groups.

=head1 EXAMPLES

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

# CLdb #
use FindBin;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $group_file);
my $percent = 50;
GetOptions(	   
	   "group=s" => \$group_file,
	   "percent=f" => \$percent,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: cannot find $group_file\n"
  unless -e $group_file;


#--- MAIN ---#
# loading group file
my $groups_r = load_groups($group_file);

# loading spacers_shared.txt file
my $shared_r = load_spacers_shared($groups_r);

# determining number (& percent) of specific spacers per group
get_specific_spacers($shared_r, $groups_r, $percent);


#--- Subroutines ---#
sub get_specific_spacers{
  my ($shared_r, $groups_r, $percent) = @_;

  # intializing counts in summary
  my %summary;
  foreach my $group_id (keys %{$groups_r->{counts}}){
    $summary{$group_id} = 0;
  }

  # summarizing
  foreach my $spacer_id (keys %{$shared_r->{spacer}}){
    # determining how many groups have >0 counts
    my $pres_cnt = 0;
    map{ $pres_cnt++ if $_ > 0 } @{$shared_r->{spacer}{$spacer_id}}{ keys %{$shared_r->{spacer}{$spacer_id}} };

    # if pres_cnt == 1, check for pres in >= %percent of group
    if($pres_cnt == 1){
      foreach my $group_id (keys %{$shared_r->{spacer}{$spacer_id}}){
	die "LOGIC ERROR: cannot find '$group_id' in group counts\n"
	  unless exists $groups_r->{counts}{$group_id};
	
	if($shared_r->{spacer}{$spacer_id}{$group_id} > 0){ # must be found in > percent of group
	  $summary{$group_id}++ if 
	    $shared_r->{spacer}{$spacer_id}{$group_id} / $groups_r->{counts}{$group_id} * 100 
	      > $percent;
	}
      }
    }
  }

  
  # writing output
  ## header
  print join("\t", qw/group N_specific_spacers total_group_spacers 
		      percent_of_group total_spacers percent_of_total/), "\n";
  ## body
  foreach my $group_id (sort keys %summary){
    print join("\t", $group_id, $summary{$group_id}, 
	       $shared_r->{total}{$group_id},
	       sprintf("%.2f", ($summary{$group_id} / 
				$shared_r->{total}{$group_id}) * 100),
	       scalar keys %{$shared_r->{spacer}},
	       sprintf("%.2f", ($summary{$group_id} /
				scalar keys %{$shared_r->{spacer}}) * 100)
	      ), "\n";
  }

}

sub load_groups{
# loading group file
  my ($group_file) = @_;

  open IN, $group_file or die $!;

  my %groups;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    my @l = split /\t/;
    die "ERROR: line $. of groups file does not have 2 columns (tab-delimited)\n"
      unless scalar @l == 2;

    $groups{counts}{$l[1]}++ unless exists $groups{grouping}{$l[0]};
    $groups{grouping}{$l[0]} = $l[1];    
  }
  
  #print Dumper %groups; exit;
  return \%groups;
}

sub load_spacers_shared{
# loading spacers shared file
## %spacers: spacer_id => group => count
  my ($groups_r) = @_;

  my %spacers;
  #my %header;
  while(<>){
    chomp;
    next if /^\s*$/;    
    my @l = split /\t/;
    
    if( ! exists $groups_r->{index} ){  # header
      for my $i (1..$#l){
	if(exists $groups_r->{grouping}{$l[$i]}){
	  $groups_r->{index}{$i} = $groups_r->{grouping}{$l[$i]};
	}
	else{ warn "$l[$i] not found in group_file\n" }
      }      
    }
    else{   # body
      # initializing group counts
      foreach my $group (keys %{$groups_r->{counts}}){
	$spacers{spacer}{$l[0]}{$group} = 0;
	$spacers{total}{$group} = 0 
	  unless exists $spacers{total}{$group};
      }

      # counting
      my %group_cnt;
      for my $i (1..$#l){
	my $group = $groups_r->{index}{$i};
	if($l[$i] > 0){
	  $spacers{spacer}{$l[0]}{$group}++;
	  $group_cnt{$group} = 1;
	}
      }  
      foreach my $group_id (keys %group_cnt){ # adding to total if spacer present in group
	$spacers{total}{$group_id}++; 
      }
    }
  }

  #print Dumper %spacers; exit;
  return \%spacers;
}
