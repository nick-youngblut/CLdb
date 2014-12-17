#!/usr/bin/env perl

=pod

=head1 NAME

spacerContentDiff.pl -- Create a table of spacer content variation among CRISPRs 

=head1 SYNOPSIS

spacerContentDiff.pl [flags] 

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb_perldoc spacerContentDiff.pl

=head1 DESCRIPTION



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
use File::Spec;

# 3rd party
use DBI;
use Set::IntervalTree;


# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my $nBins = 10;
GetOptions(
	   "bin=i" => \$nBins,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: bin must be >0\n" unless $nBins > 0;

#--- MAIN ---#
# load alignment scores
my $scores_r = parse_scores();
write_scores($scores_r);
# convert scores to numeric
#scores2numeric($scores_r);
# bin scores & calculating seqID for the bin
#binSeqID($scores_r, $nBins);



#--- Subroutines ---#


sub binSeqID{
  my $scores_r = shift or die $!;
  my $nBins = shift or die $!;  

  my $binFrac = 1 / $nBins;

  foreach my $comp (@$scores_r){
    # assertions
    map{ die "KeyError: '$_' not found\n"
      unless exists $comp->{$_} } qw/aln_scores aln_rel_pos/;

    # summing scores by bin
    my $score_sum = 0;
    my $bin_len = 0;
#    for (my $i=0; $i<=1, $i+=$binFrac){
      
#    }
    
  }

}


sub parse_scores{
  # parsing the scores file produced by spacerAcquisition

  my @scores;
  while(<>){
    next if $. == 1;
    chomp;
    s/#.+//;
    next if /^\s*$/;

    # spliting column
    my @l = split /\t/;
    die "ERROR: line $. has < 3 columns!\n"
      unless scalar @l >= 3;

    # spliting alignment scores
    my $aln_scores = [split //, $l[2]];
    my $aln_len = scalar @$aln_scores;

    # scores to numeric
    _scores2numeric($aln_scores);

    # making itree of scores
    $aln_scores = _make_itree($aln_scores);    

    # loading values
    push @scores, { locus_i => $l[0],
		    locus_j => $l[1],
		    aln_scores => $aln_scores,
		    aln_len => $aln_len
		    };
  }

  #print Dumper @scores; exit;
  return \@scores;
}


sub _make_itree{
  my $scores_r = shift or die $!;

  my $aln_len = scalar @$scores_r;

  # making positional array
  my @aln_rel_pos;
  for my $i (1..$aln_len){
    $aln_rel_pos[$i-1] = $i / $aln_len;
  }

  my $itree = Set::IntervalTree->new;
  
  for my $i (0..$aln_len-1){
    $itree->insert($scores_r->[$i], $aln_rel_pos[$i], $aln_rel_pos[$i])
  }
  
  print Dumper $itree; exit;
}


sub _scores2numeric{
  # scoring from character to numeric
  my $scores_r = shift or die $!;

  # index: character => numeric
  my %index = ( 'm' => 1,
		'x' => 0,
		'g' => 'NA' );

  # iterating over alignment scores
  # converting
  @$scores_r = map{  
    exists $index{$_} ? 
      $index{$_} :
	die "ERROR: illegal character in alignment scores -> '$_' \n";      
  } @$scores_r;

  #print Dumper $scores_r; exit;
}
