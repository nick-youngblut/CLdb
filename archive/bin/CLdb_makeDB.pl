#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace, $drop_all);
my @tables = ("leaders", "genes");
GetOptions(
	   "replace" => \$replace,
	   "tables=s{,}" => \@tables,
	   "drop" => \$drop_all,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$ARGV[0] = "CLdb.sqlite" unless $ARGV[0];
map{tr/A-Z/a-z/} @tables if @tables;		# lower case table names

### MAIN
my $sql_r = get_sql();
make_db($sql_r, $ARGV[0], \@tables);

### Subroutines
sub make_db{
	my ($sql_r, $db_name, $tables_r) = @_;
	
	# checking if tables specified exists, deleted if yes, dying if no #
	if(-e $db_name){
		foreach my $table (@$tables_r){
			if(exists $sql_r->{$table}){
				unless($drop_all){
					print STDERR "...Not dropping table: \"$table\"\n";
					delete $sql_r->{$table};
					}
				}
			else{
				print STDERR " ERROR: table: \"$table\" not found in sql for making tables\n";
				print STDERR "### tables in sql (ie. the tables that will be created) ###\n";
				print STDERR join(",\n", keys %$sql_r), "\n";
				exit;
				}
			}
		# checking for overwrite of database #
		die " ERROR: $db_name already exists! Use '-r' to replace\n" unless $replace
		}
	
	# adding tables #
	foreach my $table (keys %$sql_r){
		open PIPE, "| sqlite3 $db_name" or die $!;
		print PIPE "BEGIN TRANSACTION;\n";
		print PIPE $sql_r->{$table}; 				# making table
		print PIPE "COMMIT;\n";
		close PIPE;
		}
		
	print STDERR "...sqlite3 database tables created\n";
	}

sub get_sql{
	my %sql; 		# all tables individually 

	$sql{"loci"} = <<HERE;
/* creating tables */
DROP TABLE IF EXISTS Loci;

CREATE TABLE Loci (
Locus_ID	INTEGER	PRIMARY KEY,
Taxon_ID	TEXT	NOT NULL,
Taxon_Name	TEXT	NOT NULL,
Subtype	TEXT,
Scaffold	TEXT	NOT NULL,
Locus_Start	INTEGER	NOT NULL,
Locus_End	INTEGER	NOT NULL,
Operon_Start	INTEGER,
Operon_End	INTEGER,
Array_Start	INTEGER,
Array_End	INTEGER,
Operon_Status	TEXT	NOT NULL,
Array_Status	TEXT	NOT NULL,
Genbank_File	TEXT	NOT NULL,
Fasta_File	TEXT,
Array_File	TEXT,
Scaffold_count	INTEGER,
File_Creation_Date	TEXT,
Author	TEXT	NOT NULL,
UNIQUE (Taxon_ID, Taxon_name, Scaffold, Locus_Start, Locus_End)
ON CONFLICT IGNORE
);

HERE


	$sql{"spacers"} = <<HERE;
DROP TABLE IF EXISTS Spacers;

CREATE TABLE Spacers (
Locus_ID	INTEGER	NOT NULL,
Spacer_ID	TEXT	NOT NULL,
Spacer_Start	INTEGER	NOT NULL,
Spacer_End	INTEGER	NOT NULL,
Spacer_Sequence	TEXT	NOT NULL,
Spacer_Group	TEXT,
UNIQUE (Locus_ID, Spacer_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"drs"} = <<HERE;
DROP TABLE IF EXISTS DRs;

CREATE TABLE DRs (
Locus_ID	INTEGER	NOT NULL,
DR_ID	INTEGER	NOT NULL,
DR_Start	INTEGER	NOT NULL,
DR_End	INTEGER	NOT NULL,
DR_Sequence	TEXT	NOT NULL,
DR_Group	INTEGER,
UNIQUE (LOCUS_ID, DR_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"DR_consensus"} = <<HERE;
DROP TABLE IF EXISTS DR_Consensus;

CREATE TABLE DR_Consensus (
Locus_ID	INTEGER	NOT NULL,
Consensus_Sequence_IUPAC	TEXT	NOT NULL,
Consensus_Sequence_Threshold	TEXT	NOT NULL,
UNIQUE (Locus_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"leaders"} = <<HERE;
DROP TABLE IF EXISTS Leaders;

CREATE TABLE Leaders (
Locus_ID	INTEGER	NOT NULL,
Leader_Start	INTEGER	NOT NULL,
Leader_End	INTEGER	NOT NULL,
Leader_Sequence	TEXT	NOT NULL,
Leader_Group	TEXT,
UNIQUE (Locus_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"genes"} = <<HERE;
DROP TABLE IF EXISTS Genes;

CREATE TABLE Genes (
Locus_ID	INTEGER	NOT NULL,
Gene_ID	TEXT	NOT NULL,
Gene_Start	INTEGER	NOT NULL,
Gene_End	INTEGER	NOT NULL,
Gene_Length__AA	INTEGER	NOT NULL,
In_Operon	TEXT	NOT NULL,
Gene_Alias	TEXT,
UNIQUE (Locus_ID, Gene_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"blast_hits"} = <<HERE;
DROP TABLE IF EXISTS blast_hits;

CREATE TABLE blast_hits (
blast_id	INTEGER	PRIMARY KEY,
spacer_DR    TEXT    NOT NULL,
Group_ID    TEXT    NOT NULL,
S_taxon_ID    TEXT,
S_taxon_name    TEXT,
S_accession    TEXT,
S_GI    TEXT,
sseqid    TEXT    NOT NULL,
pident    REAL    NOT NULL,
mismatch    INTEGER    NOT NULL,
gaps    INTEGER    NOT NULL,
evalue    TEXT    NOT NULL,
bitscore    INTEGER    NOT NULL,
strand	INTEGER	NOT NULL,
len    INTEGER    NOT NULL,
qlen	INTEGER	NOT NULL,
slen    INTEGER    NOT NULL,
qseq	TEXT,
sseq	TEXT,
qstart    INTEGER    NOT NULL,
qend    INTEGER    NOT NULL,
sstart    INTEGER    NOT NULL,
send    INTEGER    NOT NULL,
qseq_full	TEXT,
sseq_full	TEXT,
qseq_full_start	INTEGER,
qseq_full_end	INTEGER,
sseq_full_start INTEGER,
sseq_full_end 	INTEGER,
proto3px	TEXT,
proto3px_start	INTEGER,
proto3px_end	INTEGER,
proto5px	TEXT,
proto5px_start	INTEGER,
proto5px_end	INTEGER,
array_hit    TEXT,
UNIQUE(Group_ID, S_taxon_ID, S_taxon_name, S_accession, sseqid, sstart, send)
ON CONFLICT REPLACE
);

HERE


	$sql{"spacer_blast_subject"} = <<HERE;
DROP TABLE IF EXISTS spacer_blast_subject;

CREATE TABLE spacer_blast_subject (
Blast_subject_ID	TEXT	NOT NULL,
Taxon_ID	TEXT,
Taxon_name	TEXT,
Accession	TEXT,
Scaffold_name	TEXT	NOT NULL,
Scaffold_sequence	TEXT	NOT NULL,
Fragment_start	INTEGER	NOT NULL,
Fragment_end	INTEGER	NOT NULL,
Extension	INTEGER	NOT NULL,
UNIQUE (Blast_subject_ID)
ON CONFLICT REPLACE
);

HERE


	$sql{"spacer_hclust"} = <<HERE;
DROP TABLE IF EXISTS spacer_hclust;

CREATE TABLE spacer_hclust (
Locus_ID	TEXT	NOT NULL,
Spacer_ID	TEXT	NOT NULL,
Cutoff	REAL	NOT NULL,
Cluster_ID	INTEGER	NOT NULL,
UNIQUE (Locus_ID, Spacer_ID, Cutoff)
ON CONFLICT REPLACE
);

HERE


	$sql{"DR_hclust"} = <<HERE;
DROP TABLE IF EXISTS DR_hclust;

CREATE TABLE DR_hclust (
Locus_ID	TEXT	NOT NULL,
DR_ID	TEXT	NOT NULL,
Cutoff	REAL	NOT NULL,
Cluster_ID	INTEGER	NOT NULL,
UNIQUE (Locus_ID, DR_ID, Cutoff)
ON CONFLICT REPLACE
);

HERE


	$sql{"spacer_pairwise_blast"} = <<HERE;
DROP TABLE IF EXISTS spacer_pairwise_blast;

CREATE TABLE spacer_pairwise_blast (
Query_locus_ID	INTEGER	NOT NULL,
Query_spacer_ID	TEXT	NOT NULL,
Subject_locus_ID	INTEGER	NOT NULL,
Subject_spacer_ID	TEXT	NOT NULL,
pident	REAL	NOT NULL,
length	INTEGER	NOT NULL,
mismatch	INTEGER	NOT NULL,
gapopen	INTEGER	NOT NULL,
qstart	INTEGER	NOT NULL,
qend	INTEGER	NOT NULL,
sstart	INTEGER	NOT NULL,
send	INTEGER	NOT NULL,
evalue	TEXT	NOT NULL,
bitscore	INTEGER	NOT NULL,
UNIQUE( Query_locus_ID, Query_spacer_ID, Subject_locus_ID, Subject_spacer_ID )
ON CONFLICT REPLACE
);

HERE

	return \%sql;
	}

__END__

=pod

=head1 NAME

CLdb_makeDB.pl -- Initial DB construction

=head1 SYNOPSIS

CLdb_makeDB.pl [options] [DATABASE_name]

=head2 options

=over

=item -replace  <bool>

Replace existing database.

=item -table  <char>

Table(s) to keep as is (if they exist). ["leaders" "genes"]

=item -drop  <bool>

Drop all tables. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_makeDB.pl

=head1 DESCRIPTION

Make all of the CL_db tables.

The default database name is "CLdb.sqlite"

=head2 '-t' flag

This is needed if you want to keep 'manual' information 
in tables while still being able to remake the rest
of the database. By default, all leader sequences and gene alias info
is saved because manual processing was involved.
Capitalization of table names doesn't matter. 

=head1 EXAMPLES

=head2 Naming database 'CLdb_test'

CLdb_makeDB.pl CLdb_test

=head2 Remaking a database and keeping old "spacers" table

CLdb_makeDB.pl -r -t "spacers"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

