#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;


# usage
my $usage = <<HERE;
Simple script for converting feature table (from genbank_get_features) to fasta.
\tUsage: $0 < feature_file > fasta_file
HERE


if ((@ARGV == 0) && (-t STDIN)){
  print $usage;
  exit()
}


# loading table and writing
my %header;
while(<>){
  chomp;
  next if /^\s*$/;
  my @l = split /\t/;

  # header
  if ($.==1){
    for my $i (0..$#l){
      $header{$l[$i]} = $i;
    }
    # checking for header values
    map{ die "ERROR: cannot find $_ in header\n"
	   unless exists $header{$_} } qw/product translation db_xref/;
  }
  else{
    die "ERROR: no header found!\n"
      unless defined keys %header;

    # checking for existing values
    map{ die "ERROR: cannot find $_ value in row\n"
	   unless defined $l[$header{$_}] } qw/product translation db_xref/;

    # sequence name
    my $db_xref = $l[$header{db_xref}];    
    my $product = $l[$header{product}];    
    $db_xref =~ s/.*(ITEP.+peg\.\d+).*/$1/;
    my $seq_name = join("_::_", $db_xref, $product);
    
    # sequence
    my $seq = $l[$header{translation}];


    # printing
    print join("\n", ">$seq_name", $seq), "\n";
  }
}
