#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlast.pl -- BLASTn-short of spacers and/or DRs against >=1 genome in CLdb

=head1 SYNOPSIS

CLdb_arrayBlast.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=item -query  <char>

A fasta of either spacer and/or DR group sequences (use: CLdb_array2fasta.pl -g -l)

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -extra  <char>

Extra sql to refine which sequences are returned.

=item -blast  <char>

BLASTn parameters (besides required flags). [-evalue 0.1]

=item -v  <bool>

Verbose output. [TRUE]

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlast.pl

=head1 DESCRIPTION



=head1 EXAMPLES

=head2 Blasting all spacers against all genomes in CLdb

CLdb_array2fasta.pl -g > all_spacer_groups.fna

CLdb_arrayBlast.pl -d CLdb.sqlite -q all_spacer_groups.fna

=head2 Blasting all DRs against all genomes in CLdb

CLdb_array2fasta.pl -r -g > all_DR_groups.fna

CLdb_arrayBlast.pl -d CLdb.sqlite -q all_DR_groups.fna

=head2 Blasting all spacers against 1 genome

CLdb_array2fasta.pl -g > all_spacer_groups.fna

CLdb_arrayBlast.pl -d CLdb.sqlite -q all_spacer_groups.fna -taxon_name "e.coli"

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
use DBI;
use File::Path qw/rmtree/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	get_file_path/;
use CLdb::blast qw/
	read_blast_file/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $out_dir, $by_query, $by_db);
GetOptions(
	   "directory=s" => \$out_dir,
	   "query" => \$by_query, 		# output by query
	   "database" => \$by_db, 		# output by db
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
unless( defined $out_dir ){
	my $curdir = File::Spec->rel2abs(File::Spec->curdir());
	$out_dir = join("/", $curdir, "arrayBlastAlign");
	}
mkdir $out_dir unless -d $out_dir;			# making output directory 

#--- MAIN ---#
my $blast_r = read_blast_file();

blast_2_align($blast_r, $out_dir);

#--- Subroutines ---#
sub blast_2_align{
# writing a fasta for each query #
	my ($blast_r, $out_dir) = @_;
	
	open OUT, ">$out_dir/all.fna" or die $!
		if ! defined $by_query && ! defined $by_db;		# output all together	
	
	foreach my $query ( keys %$blast_r ){
		
		# checking for formating of query #
		(my $name = $query) =~ s/# Query: //i;
		die "ERROR: '$name' not formatted as 'locus_id|spacer/DR|spacerID|spacer_group'\n"
			unless $query =~ /^[^|]+\|/;
		$name =~ s/\|/-/g;
		
		# output file #
		open OUT, ">$out_dir/$name" or die $!
			if $by_query && ! defined $by_db;
		
		foreach my $db ( keys %{$blast_r->{$query}} ){
			
			my @parts = File::Spec->splitpath($db);
			
			# output file #
			open OUT, ">$out_dir/$name\__$parts[2]" or die $!
				if $by_query && $by_db;
			
			foreach my $blast ( keys %{$blast_r->{$query}{$db}} ){
				# checking for 'proto_seq_fullx' & 'query_seq_full_DRx' #
				unless exists $blast_r->{$query}{$db}{$blast}{
				
				foreach my $hit ( @{$blast_r->{$query}{$db}{$blast}}){
				
					
					
					
					}
				}
			
			
			close OUT if $by_query && $by_db;
			}
		close OUT if $by_query && ! defined $by_db;
		}
	close OUT if ! defined $by_query && ! defined $by_db;		# output all together	
	}


