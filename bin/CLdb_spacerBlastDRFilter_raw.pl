#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Set::IntervalTree;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $array_bool);
my $range = 30;
GetOptions(
	   "range=i" => \$range,
	   "array" => \$array_bool, 			# just provide a column stating in_array or not
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a spacer blast table & a direct repeat blast table!\n"
	unless $ARGV[0] && $ARGV[1];		# spacer blast, DR blast
$range++;

### MAIN
my $itrees_r = make_DR_itree($ARGV[1], $range);
filter_spacer_blast($ARGV[0], $itrees_r, $array_bool);

### Subroutines
sub filter_spacer_blast{
	my ($blast_in, $itrees_r, $array_bool) = @_;
	
	open IN, $blast_in or die $!;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		
		if(exists $itrees_r->{$line[1]}){		# hit on same taxon-scaffold of DR & spacer
			my $res;
			if( $line[8] < $line[9] ){
				$res = $itrees_r->{$line[1]}->fetch($line[8], $line[9]);
				}
			else{
				$res = $itrees_r->{$line[1]}->fetch($line[9], $line[8]);
				}
			
			
			if($array_bool){
				if(scalar @$res >=2){			# spacer is adjacent to 2 DR (in array)
					print "$_\tyes\n"; 			# in array
					}
				else{
					print "$_\tno\n";			# not in array
					}
				}
			else{ 
				print "$_\n" unless @$res && scalar @$res >=2; 		# if spacer is adjacent to 2 DR
				}			
			}
		else{			# do not hit on the same scaffold as a DR
			if($array_bool){
				print "$_\tno\n";
				}
			else{
				print "$_\n";
				}
			}
		}
	
	}

sub make_DR_itree{
	my ($blast_in, $range) = @_;
	
	open IN, $blast_in or die $!;
	
	my %itrees;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		
		# making itree #
		unless(exists  $itrees{$line[1]}){  
			$itrees{$line[1]} = Set::IntervalTree->new();
			}
			
		# loading itree #
		if($line[8] < $line[9]){
			$itrees{$line[1]}->insert($., $line[8] - $range,  $line[9] + $range);		
			}
		else{
			$itrees{$line[1]}->insert($., $line[9] - $range,  $line[8] + $range);		
			}
		}
	close IN;
	
	return \%itrees;
	}


__END__

=pod

=head1 NAME

CLdb_spacerBlastDRFilter_raw.pl -- filter spacer blast table using direct-repeat blast table

=head1 SYNOPSIS

CLdb_spacerBlastDRFilter_raw.pl [options] spacer_blast.txt DR_blast.txt > spacer_blast_filtered.txt

=head2 options

=over

=item -range  <int>

Range allowable between spacer & DR blast hit (bp). [30]

=item -array  <bool>

Do not filter hit. Instead, add an In_Array column (yes|no). [FALSE]

=item -verbose  <booL>

Verbose output.

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_spacerBlastDRFilter_raw.pl

=head1 DESCRIPTION

Filter out spacer blast hits to CRISPR arrays by
removing all spacer blast hits that hit adjacent
to 2 direct repeat blast hits.

Use '-a' if you want all hits including spacer
hits to an array (last column will designate whether
spacer hit falls in a CRISPR array).

=head1 EXAMPLES

=head2 Basic Usage:

CLdb_spacerBlastDRFilter_raw.pl spacer_blast.txt repeat_blast.txt > spacer_blast_filter.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

