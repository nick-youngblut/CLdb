#!/usr/bin/env perl

=pod

=head1 NAME

crt_parse.pl -- parsing output from CRISPR Recognition Tool (CRT); convert to CRISPRFinder format

=head1 SYNOPSIS

crt_parse.pl [flags] < crt.out

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -prefix  <char>

Output file(s) prefix. [array]

=item -coord  <bool>

Array start-end coords added to each array file name? [FALSE]

=item -verbose  <bool>

Verbose output

=item -h	This help message

=back

=head2 For more information:

CLdb --perldoc -- crt_parse.pl

=head1 DESCRIPTION

Parse the putative CRISPR arrays identified by 
the CRISPR Recognition Tool (CRT)
and convert to CRISPRFinder format.

Each array will be written to a separate file.

Default file naming: 'prefix'_'CRISPR#'_'scaffold_name'

start-end is in relation to the positive strand.

=head1 EXAMPLES

=head2 Basic usage (CRT output file = 'a.out'):

CLdb -- crt_parse.pl < a.out

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $coords_b, $database);
my $prefix = "array";
GetOptions(
	   "prefix=s" => \$prefix, 
	   "coords" => \$coords_b,
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#


#--- MAIN ---#
my $res_r = parse_crt();
write_array_files($res_r, $prefix, $coords_b);

#--- Subroutines ---#
sub write_array_files{
  my ($res_r, $prefix, $coords_b) = @_;
  
  foreach my $array (keys %$res_r){
    my $outname = join("_", $prefix, $array, $res_r->{$array}{'organism'});
    $outname = join("_", $outname, $res_r->{$array}{'start'}, 
		    $res_r->{$array}{'end'})
      if $coords_b;
    
    open OUT, ">$outname.txt" or die $!;
    foreach my $line (@{$res_r->{$array}{'array'}}){
      print OUT join("\t", @$line), "\n";
    }
    close OUT;
    
    print STDERR "File written: $outname.txt\n";
  }
}

sub parse_crt{
  # parsing CRT output 
  my %res;
  my $org;
  while(<>){
    chomp;
    
    ($org = $_) =~ s/.+ // if /^ORGANISM/;
    
    if(/^CRISPR \d+/){
      my @head = split / +/;	# 1=num; 3=start; 5=end
      die "ERROR: line $. is not formatted correctly!\n"
	unless scalar @head == 6;
      
      $res{$head[1]}{'organism'} = $org if defined $org;
      $res{$head[1]}{'start'} = $head[3];
      $res{$head[1]}{'end'} = $head[5];
      die "ERROR: start > end at line $.\n"
	unless $head[3] <= $head[5];
			
      while(<>){
	chomp;
	last if /^Repeats:/;
	next if /^(POSITION|^-)/;
	
	my @line = split /\t+/;
	die "ERROR: line $. does not have 2 or 4 tab-delimited columns!\n"
	  unless scalar @line == 2 || scalar @line == 4;
	
	# calculating end #
	my ($DR_len, $spacer_len) = (0,0);
				$DR_len = length $line[1] if defined $line[1];
	$spacer_len = length $line[2] if defined $line[2];				
	$line[3] = $line[0] + $DR_len + $spacer_len - 1;
	
	$line[2] = "\t" unless defined $line[2];
	push @{$res{$head[1]}{'array'}}, \@line;
      }
			}
  }
  
  #print Dumper %res; exit;
  return \%res;
}

