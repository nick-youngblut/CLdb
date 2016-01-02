#!/usr/bin/env perl

=pod

=head1 NAME

Cldb_loadGenes.pl -- adding/updating gene table in CRISPR database

=head1 SYNOPSIS

Cldb_loadGenes.pl [flags] < gene_table.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back 

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc Cldb_loadGenes.pl

=head1 DESCRIPTION

Add/update genes entries in Genes table within the
specified CRISPR database.

A gene table made by getGenesInLoci.pl can be
piped directly into loadGenes.pl, or 
the aliases (or other values) can be currated (manually)
first. 

=head1 EXAMPLES

=head2 Basic usage:

Cldb -- loadGenes -d CLdb.sqlite < genes_table.txt

=head2 Piping from getGenesInLoci.pl 

CLdb -- getGenesInLoci -d CLdb.sqlite | Cldb -- loadGenes -d CLdb.sqlite

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::utilities qw/
	file_exists 
	connect2db
	lineBreaks2unix/;
use CLdb::query qw/
	table_exists/;
use CLdb::load qw/
	load_db_table/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# database metadata #
table_exists($dbh, "genes");
#get_genes_table_info($dbh);

# loading genes table #
my $tbl_r = lineBreaks2unix(\*STDIN);
my ($genes_r, $header_r) = get_gene_table($tbl_r);

# removing 'sequence' feature; not included in genes CLdb table #
delete $header_r->{'sequence'};

# updating / loading_db #
load_db_table($dbh, "genes", $header_r, $genes_r);

# disconnect to db #
$dbh->disconnect();
exit;


### Subroutines
sub get_gene_table{
	my $tbl_r = shift;
	my %genes;
	my %header;
	my %header_rev;
	my $cnt = 0;
	foreach (@$tbl_r){
		$cnt++;

		if(! %header){ # loading header
			tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
			my @line = split /\t/;
			for my $i (0..$#line){
				$header{$line[$i]} = $i;		# column_name => index
				$header_rev{$i} = $line[$i];
				}
			}
		else{
			my @line = split /\t/;
			# loading as locus_id=>category=>value #
			for my $i (0..$#line){
				my $uID = join("_", $line[$header{"locus_id"}], $cnt);
				$genes{$uID}{$header_rev{$i}} = $line[$i];
				}			
			}
		}
		#print Dumper %header; 
		#print Dumper %genes; exit; 
	return (\%genes, \%header);;
	}



