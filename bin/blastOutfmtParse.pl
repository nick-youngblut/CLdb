#!/usr/bin/env perl

=pod

=head1 NAME

blastOutfmtParse.pl -- parsing a blast file in '-outfmt 7' format while keeping the comments

=head1 SYNOPSIS

blastOutfmtParse.pl < blast.txt -c > blast_edit.txt

=head2 Required flags

=over

=item -command 

Unix command for editing blast hits (eg., 'grep "E*coli"')

=back

=head2 Optional flags

=over

=item -s 

Index hash size (>= number of comments in BLAST file). [1000000]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb --perldoc -- blastOutfmtParse.pl

=head1 DESCRIPTION

Simple method to parse a blast table with comment lines
(eg., using grep), while still retaining the comments.

The blast table is parsed into comments and blast hits,
with a index added to the beginning of the row for each
(eg., 1\tREST_OF_BLAST_HIT).

'-command' is applied to each blast hit. The comments
are then added back by using the index.

=head1 EXAMPLES

=head2 Just keep entries with hits

CLdb -- blastOutfmtParse -c "cat" < blast.txt > blast_edited.txt

=head2 Just hits to scaffold5

CLdb -- blastOutfmtParse -c "grep 'scaffold5'" < blast.txt > blast_scaf5.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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
use IPC::Open2;
use threads;
use threads::shared;

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database);
my $cmd;
my $hash_size = 1000000;
GetOptions(
	   "command=s" => \$cmd,
	   "size=i" => \$hash_size,
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide a unix command\n"
  unless defined $cmd;


#--- MAIN ---#
my $pid = open2(\*READER, \*WRITER, $cmd);

print STDERR "Making shared index hash (size=$hash_size)\n";
my %index;
for my $i (0..$hash_size){
  share( @{$index{$i}} );
}

print STDERR "Reading & writing to command: '$cmd'\n";
my $thread = async{
  close WRITER;

  my %index_used;
  while(<READER>){
    chomp;
    my @l = split /\t/, $_, 2;
    
    # sanity checks
    die "ERROR: index column not found on line $. of blast hits file\n"
    unless scalar @l == 2 and $l[0] =~ /^\d+$/;
    die "ERROR: cannot find index '$l[0]' in index file\n"
      unless exists $index{$l[0]};
    
    #adding back comments
    if(! exists $index_used{$l[0]} ){ # writing comment
      print join("\n", @{$index{$l[0]}} ), "\n";
      $index_used{$l[0]} = 1;
    }
    
    # writing hits
    print $l[1], "\n";
  }
  close READER;
};


{
  my $cmt_cnt = 0;
  #my %index;
  while(<>){
    chomp;
    next if /^\s*$/;
    
    if(/^# BLAST.*\d+\.\d+/ || eof){ # new comment
      $cmt_cnt++;          
    }
    if(/^#/){ # loading comment
      push @{$index{$cmt_cnt}}, $_;
      #$index{$.}{$.} = 1;
    }
    else{ # writing blast hit w/ comment index
      print WRITER join("\t", $cmt_cnt, $_), "\n";    
    }
  }
  close WRITER;
}  

$thread->join();
waitpid($pid, 0);
