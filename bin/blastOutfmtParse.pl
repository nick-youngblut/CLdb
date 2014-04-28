#!/usr/bin/env perl

=pod

=head1 NAME

blastOutfmtParse.pl -- parsing a blast file in '-outfmt 7' format while keeping the comments

=head1 SYNOPSIS

=head2 write indexed comments

blastOutfmtParse.pl -i < blast.txt > comment_index.txt

=head2 write indexed blast hits

blastOutfmtParse.pl -h < blast.txt > hit_index.txt

=head2 recombining comments & (modified) blast hits

blastOutfmtParse.pl -r <(comment_index.txt) <(hit_index.txt)

=head2 using tee to do it all in 1 step (selecting just hits to Ecoli)

cat blast.txt | tee >(blastOutfmtParse.pl -i < blast.txt)
>(blastOutfmtParse.pl -h < blast.txt | egrep "E*coli") >/dev/null |
blastOutfmtParse.pl -r 

=head2 Required flags

=over

=item -index 

Write index of comments

=item -blast

Write index of blast hits

=item -recombine

Recombine comments and blast hits

=back

=head2 Optional flags

=over

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc blastOutfmtParse.pl

=head1 DESCRIPTION

Simple method to parse a blast table with comment lines
(eg., using grep), while still retaining the comments.


=head1 EXAMPLES

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
use File::Tail;

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose);
my ($index, $hit, $rec);
GetOptions(
	   "index" => \$index,
	   "blast" => \$hit,
	   "recombine" => \$rec,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide either '-i', '-h', or '-r'\n"
  unless $index || $hit || $rec;


#--- MAIN ---#
if($index){
  parse_blast_comment();
}
elsif($hit){
  parse_blast_hit();
}
elsif($rec){
  recombine_blast();
}


#--- Subroutines ---#
sub recombine_blast{
  
}

sub parse_blast_comment{
  my $index_cnt = 1;
  my $last = '';
  while(<>){
    chomp;
    next if /^\s*$/;

    if(/^#/){
      print join("\t", "#$index_cnt", $_),"\n";
    }
    else{
      $index_cnt++ if $last =~ /^#/;
    }
    $last = $_;
  }
}

sub parse_blast_hit{
  my $index_cnt = 1;
  my $last = '';
  while(<>){
    chomp;
    next if /^\s*$/;

    if(/^#/){
    }
    else{
      print join("\t", "#$index_cnt", $_),"\n";
      $index_cnt++ if $last =~ /^#/;
    }
    $last = $_;
  }
}




sub parse_blast_add_comments{
# parsing blast table and adding back comments from index file
  my ($index_r) = @_;

  my %index_used;
  while(<>){
    chomp;
    next if /^\s*$/;
    my @l = split /\t/, $_, 2;

    # sanity checks
    die "ERROR: index column not found on line $. of blast hits file\n"
      unless scalar @l == 2 and $l[0] =~ /^\d+$/;
    die "ERROR: cannot find index '$l[0]' in index file\n"
      unless exists $index_r->{$l[0]};

    # adding back index
    if(! exists $index_used{$l[0]} ){ # writing comment
      print join("\n", @{$index_r->{$l[0]}} ), "\n";
      $index_used{$l[0]} = 1;
    }
    
    print $l[1], "\n";
  }
}

sub load_index{
# loading index file
  my ($index) = @_;

  my %index;

#  my $file = File::Tail->new($index);
#  while( defined(my $line=$file->read) ){
  
  use IO::Handle;
  
  open IN, $index or die $!;
  for (;;){
    while(<IN>){
      my $line = $_;
      chomp $line;
      
      next if $line =~ /^\s*$/;
      my @l = split /\t/, $line, 2;
      die "ERROR: line $. has < 2 columns\n"
	unless scalar @l == 2;
      
      push @{$index{$l[0]}}, $l[1];
    }
    sleep 1;
    IN->clearerr();
  }
  close IN;
  
  print Dumper %index; exit;
  return \%index;
}

sub parse_blast_write_index{
# removing comments and adding index to blast hits
  my ($output) = @_;

  open OUT, ">$output" or die $!;
  my @cmt;
  my $cmt_cnt = 0;
  while(<>){
    chomp;
    next if /^\s*$/;
    
    if(/^# BLAST.*\d+\.\d+/ || eof){ # new comment
      $cmt_cnt++;      

      # writing out last comment
      map{ print OUT "$cmt_cnt\t$_\n" } @cmt;

      # reseting
      @cmt = ();      
    }
    if(/^#/){ # loading comment
      push @cmt, $_;
    }
    else{ # writing blast hit w/ comment index
      print "$cmt_cnt\t$_\n";
    }

  }
  close OUT;
}

