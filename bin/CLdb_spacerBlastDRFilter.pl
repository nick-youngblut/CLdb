#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_spacerBlastDRFilter_raw.pl -- filter out spacers hitting CRISPR arrays based on adjacency to DR blast hits

=head1 SYNOPSIS

=head2 Separate spacer & DR blast tables 

CLdb_spacerBlastDRFilter_raw.pl [options] spacer_blast.txt DR_blast.txt > spacer_blast_filtered.txt

=head2 Combined spacer & DR blast tables

CLdb_spacerBlastDRFilter_raw.pl [options] spacer_blast.txt DR_blast.txt > spacer_blast_filtered.txt

=head2 options

=over

=item -range  <int>

Range allowable between spacer & DR blast hit (bp). [30]

=item -array  <bool>

Do not filter hit. Instead, add an 'array_hit' column. [FALSE]

=item -DR  <int>

A spacer blast is considered in an array if number of adjacent DR blasts is >= '-DR'. [1]

=item -length  <float>

DR blast hits must be >= fraction of total DR sequence. [0.66]

=item -evalue  <float>

DR blast hits must be < 'evalue'. [ ] 

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_spacerBlastDRFilter_raw.pl

=head1 DESCRIPTION

Filter out spacer blast hits to CRISPR arrays by
removing all spacer blast hits that hit adjacent
to >= '-DR' direct repeat blast hits.

Use '-a' if you want all hits including spacer
hits to an array (last column of output will designate whether
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
my $evalue_cut;
my $len_cut = 0.66;
my $DR_cnt = 1;
my $range = 30;
GetOptions(
	   "range=i" => \$range,
	   "array" => \$array_bool, 			# just provide a column stating in_array or not
	   "length=f" => \$len_cut,			# length cutoff of DR hits
	   "evalue=f" => \$evalue_cut,			# evalue cutoff of DR hits
	   "DR=i" => \$DR_cnt, 					# number of adjacent DR hits needed to call 'array'
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide either a combined spacer/DR blast table or each separately!\n"
	unless $ARGV[0];		# spacer blast, DR blast
$range++;
map{ die " ERROR: cannot find $_!\n" if $_ && ! -e $_} @ARGV[0..1];

### MAIN
my $itrees_r;
if($ARGV[1]){ $itrees_r = make_DR_itree($ARGV[1], $range, $len_cut, $evalue_cut); }
else{ $itrees_r = make_DR_itree($ARGV[0], $range, $len_cut, $evalue_cut); }

my $blast_r = parse_spacer_blast($ARGV[0]);
DR_filter_blast($blast_r, $itrees_r, $array_bool);

### Subroutines
sub DR_filter_blast{
# filtering DR's based on blast hits #
	my ($blast_r, $trees_r, $array_bool) = @_;
	
	my %filter;
	foreach my $query (keys %$blast_r){
		foreach my $db (keys %{$blast_r->{$query}}){
			# parsing header #
			my $header_r;
			foreach my $row (@{$blast_r->{$query}{$db}{"header"}}){
				if($row =~ /Field/){
					$header_r = fields2header($row);
					}
				}
			
			# body #
			foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
				#my @line = split /\t/, $row;
				
				# subject #
				my $db_subj = join("_", $db, $$row[$header_r->{"subject id"}]);

				# subject start-end #
				my $start = $$row[$header_r->{"s. start"}];
				my $end = $$row[$header_r->{"s. end"}];			
				my $strand = 1;
				$strand = 0 if $start > $end;
				($start, $end) = ($end,$start) if $start > $end;

				$filter{"total"}++;

				# itree check #
				my $DR_adjacent;
				if(exists $itrees_r->{$db_subj}{$strand}){
					$DR_adjacent = $itrees_r->{$db_subj}{$strand}->fetch(
							$start, $end);
					
					if(scalar @$DR_adjacent >= $DR_cnt){	# 'array'
						$blast_r->{$query}{$db}{"write"}++ if $array_bool;		# hit being written	
						push @$row, "array";
						$filter{"array"}++;
						}
					else{
						$blast_r->{$query}{$db}{"write"}++;		# hit being written	
						push @$row, "proto";
						$filter{"proto"}++;
						}
					}
				else{		# no DR hits; proto
					$blast_r->{$query}{$db}{"write"}++ ;		# hit being written	
					push @$row, "proto";
					$filter{"proto"}++;
					}				
				}
			
			$blast_r->{$query}{$db}{"write"} = 0 unless exists $blast_r->{$query}{$db}{"write"};
			
			# writing output #
			if($array_bool){				# writing all 
				# writing header #
				foreach my $head ( @{$blast_r->{$query}{$db}{"header"}} ){
					$head .= ", array_hit" if $head =~ /^# Fields:/;
					print $head, "\n";
					}
									
				# body #
				foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
					print join("\t", @$row), "\n";
					}				
				}
			
			elsif(	$blast_r->{$query}{$db}{"write"} > 0){		# >0 protospacer hits
				# writing header #
				print join("\n", @{$blast_r->{$query}{$db}{"header"}}), "\n";
				
				# body #
				foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
					print join("\t", @$row[0..($#$row-1)]), "\n" if $$row[$#$row] eq "proto";
					}
				}

			}
		}
		
	# status #
	map{ $filter{$_} = 0 unless exists $filter{$_}} qw/total array proto/;
	print STDERR "### DR filtering of spacer blast hits ###\n";
	print STDERR "Total spacer blast hits: $filter{'total'}\n";
	print STDERR "Spacer blast hits hitting an array (filtered out): $filter{'array'}\n";
	print STDERR "Spacer blast hits hitting a protospacer: $filter{'proto'}\n";
	}	

sub parse_spacer_blast{
# outfmt = 7;
	my ($blast_in, $itrees_r, $array_bool) = @_;
	
	my %blast;
	my @header;
	my ($db, $query);
	open IN, $blast_in or die $!;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		
		if(/^# Database/){
			($db = $_) =~ s/.+Database: //; 
			}
		if(/# Query/){
			($query = $_) =~ s/.+Query: //; 
			}
		
		if(/^#/){ push @header, $_; }
		else{
			if($query =~ /^DR_/){	# skipping 
				@header = ();
				next;
				}
			elsif(@header){ 
				push @{$blast{$query}{$db}{"header"}}, @header;
				@header = ();
				}
				
			push @{$blast{$query}{$db}{"hits"}}, [split /\t/];
			}

		}
	
		#print Dumper %blast; exit;
	return \%blast;
	}

sub make_DR_itree{
	my ($blast_in, $range, $len_cut, $evalue_cut) = @_;
	
	open IN, $blast_in or die $!;
	
	my %itrees;
	my $header_r;
	my $db;
	my %filter;
	while(<IN>){
		chomp;
		my @line = split /\t/;
		
		# header lines #
		if(/^# Database/){ ($db = $_) =~ s/^# Database: //; }
		elsif(/^# Fields/){ 
			$header_r = fields2header($_);
			check_header($header_r);
			}
		elsif(!/^#/){
			# sanity check #
			die " ERROR: no fields comment line provided in input!\n"
				unless $header_r;
			die " ERROR: no database comment line provided in input!\n"
				unless defined $db;
			
			# skipping any spacers #
			next if $line[$header_r->{"query id"}] =~ /^Spacer/;
			
			# subject #
			my $db_subj = join("_", $db, $line[$header_r->{"subject id"}]);

			# subject start-end #
			my $start = $line[$header_r->{"s. start"}];
			my $end = $line[$header_r->{"s. end"}];			
			my $strand = 1;
			$strand = 0 if $start > $end;
			($end, $start) = ($start, $end) if $start > $end;

			# filtering #
			$filter{"total"}++;
			if( $evalue_cut && $line[$header_r->{"evalue"}] > $evalue_cut){
				$filter{"evalue"}++;
				next;
				}
			my $hit_len = $end - $start + 1;
			if( $len_cut && $hit_len/$line[$header_r->{"query length"}] < $len_cut){
				$filter{"length"}++;
				next;
				} 

			# making itree #
			unless(exists  $itrees{$db_subj}{$strand} ){  
				$itrees{$db_subj}{$strand} = Set::IntervalTree->new();
				}
			
			# loading itree #
			$filter{"used"}++;
			$itrees{$db_subj}{$strand}->insert(			
						1, $start - $range,  $end + $range + 1
						);
			}
		}
	seek IN, 0, 0;
	close IN;
	
	# status #
	map{ $filter{$_} = 0 unless exists $filter{$_}} qw/total evalue length used/;
	print STDERR "### filter of 'bad' DR hits ###\n";
	print STDERR "Total DR hits: $filter{'total'}\n";
	print STDERR "DR hits with evalues < $evalue_cut: $filter{'evalue'}\n" if $evalue_cut;
	print STDERR "DR hits with hit lengths (fraction of total) < $len_cut: $filter{'length'}\n";
	print STDERR "DR hits used for filtering: $filter{'used'}\n";
	
		#print Dumper %itrees; exit;
	return \%itrees;
	}

sub check_header{
# checking to make sure that the correct fields are available #
	my ($header_r) = @_;
	
	my @needed = ("query id", "subject id", "s. start", "s. end", 
				"evalue", "query length", "subject length");
	
	map{die " ERROR: can't find '$_'!" unless exists $header_r->{$_}} @needed;
	}

sub fields2header{
	my ($line) = @_;
	
	$line =~ s/# Fields: //;
	my @l = split /\s*,\s*/, $line;		
	
	my %header;
	for my $i (0..$#l){
		$header{$l[$i]} = $i;
		}
		
		#print Dumper %header; exit;
	return \%header;
	}


