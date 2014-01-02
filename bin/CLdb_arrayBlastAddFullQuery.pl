#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddInfo.pl -- add CLdb info to spacer or DR IDs ('locusID|spacer/DR|spacer/DR_ID|groupID') in Blast file

=head1 SYNOPSIS

CLdb_arrayBlastAddInfo.pl [flags] < blast_results.txt > blast_results_info.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <bool>

Add subtype? [FALSE]

=item -taxon_id  <bool>

Add taxon_id? [FALSE]

=item -taxon_name  <bool>

Add taxon_name? [FALSE]

=item -position  <bool>

Add start-stop position? [FALSE]

=item -order  <bool>

Add spacer-leader order? [FALSE]
# not implemented yet! #

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddInfo.pl

=head1 DESCRIPTION

Add CLdb info to spacer or DR IDs (locusID|spacer/DR|spacer/DR_ID|groupID).
Also, all blast hit entries for a spacer/DR group will be
duplicated for each individual spacer/DR in the spacer/DR group (replication).

=head1 EXAMPLES

=head2 Add spacer subtype and position (BLAST hits will be replicated for each individual element in the group)

CLdb_arrayBlastAddInfo.pl -d CLdb.sqlite -sub -pos < spacer-group_blast.txt > spacers_ALL_blast.txt

=head2 Add DR subtype and position (BLAST hits will be replicated for each individual element in the group)

CLdb_arrayBlastAddInfo.pl -d CLdb.sqlite -sub -pos < DR-group_blast.txt > DR_ALL_blast.txt

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
use DBI;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
use CLdb::seq qw/
	read_fasta/;
use CLdb::blast qw/
	read_blast_file
	write_blast_file/;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $fasta_in);
GetOptions(
	   "fasta=s" => \$fasta_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($fasta_in, "fasta");

#--- MAIN ---#


# loading files #
## fasta of query sequences ##
my $fasta_r = read_fasta($fasta_in);
## blast ##
my ($lines_r) = read_blast_file();

## adding full length query ##
add_full_query_seq($lines_r, $fasta_r);

# writing edited fasta #
write_blast_file($lines_r);




### Subroutines
sub write_blast_file_OLD{
	my ($lines_r) = @_;
	
	foreach my $query ( sort keys %$lines_r ){
		foreach my $db ( sort keys %{$lines_r->{$query}} ){
			foreach my $blast ( keys %{$lines_r->{$query}{$db}} ){
				print $blast, "\n";
				print $query, "\n";
				print $db, "\n";
				print $lines_r->{$query}{$db}{$blast}{'fields'}, "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'fields'};
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
				# printing all hits #
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'hits'}}), "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'hits'};
				}
			}
		}
	}

sub add_full_query_seq{
# adding full length query sequence to blast table #
	my ($lines_r, $fasta_r) = @_;
	
	foreach my $query ( keys %$lines_r ){
		
		# checking for existence of query in fasta #
		(my $name = $query) =~ s/# Query: //i;
		unless (exists $fasta_r->{$name}){
			warn "WARNING: '$name' not found in fasta. Skipping!\n";
			next;
			}
		
		# adding sequence to blast table #
		foreach my $db ( keys %{$lines_r->{$query}} ){
			
			foreach my $blast ( keys %{$lines_r->{$query}{$db}} ){
				next unless exists $lines_r->{$query}{$db}{$blast}{'fields'};
				
				$lines_r->{$query}{$db}{$blast}{'fields'} .= ", query_seq_full";
				$lines_r->{$query}{$db}{$blast}{'fields_sep'}{"query_seq_full"} = 
					scalar keys %{$lines_r->{$query}{$db}{$blast}{'fields_sep'}};
				
				foreach my $hit ( @{$lines_r->{$query}{$db}{$blast}{'hits'}} ){
					$hit .= "\t$fasta_r->{$name}";
					}
				}
			}
		}
		#print Dumper $lines_r; exit
	}
	
sub read_blast_file_OLD{
# reading in each blast entry & extracting names and line numbers #
	my %lines;
	my $blast;
	my $query;
	my $db;
	while(<>){
		chomp;
		
		if(/^# BLAST/i){
			$blast = $_;
			}
		elsif(/^# Query/i){
			$query = $_;
			}
		elsif(/^# Database/i){
			$db = $_;
			}	
		elsif(/^# /){
			push @{$lines{$query}{$db}{$blast}{'comments'}}, $_;
			}
		else{
			push @{$lines{$query}{$db}{$blast}{'hits'}}, $_;
			}
			
		}		
		#print Dumper %lines; exit;	
	return \%lines;
	}








