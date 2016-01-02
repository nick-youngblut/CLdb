#!/usr/bin/env perl

=pod

=head1 NAME

run -- BLASTn-short of spacers and/or DRs against >=1 genome in CLdb

=head1 SYNOPSIS

run [flags] > blast_hits.xml

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=item -query  <char>

A fasta of either spacer and/or DR group sequences.
(to make needed file, use: `array2fasta.pl -d CLdb.sqlite -g -l`)

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine blast subject genomes those containing 
CRISPRs with a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine blast subject genomes to a specific taxon_id(s) (>1 argument allowed). 

=item -taxon_name  <char>

Refine blast subject genomes to a specific taxon_name(s) (>1 argument allowed). 

=item -blast  <char>

BLASTn parameters (besides required flags). [-evalue 1e-3]

=item -fork  <int>

Number of parallel blastn calls. [1] 

=item -v  <bool>

Verbose output. [TRUE]

=item -h  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc run.pl

=head1 DESCRIPTION

Run blastn-short against one or more genomes
already in CLdb (word_size=7, reward=1, DUST=off). 
Specific genomes can be
selected by subtype ('-subtype'), taxon_name
('-taxon_name'), taxon_id ('taxon_id'), or
other sql refinements of the query ('-query').

Genomes are obtained from the 'fasta_file'
values in the loci table. Any loci without an
associated fasta file will be skipped.

Genome fasta files should be in the $HOME/fasta/ directory 
(added by loadLoci.pl if genbank files were 
originally used to make the CLdb).

Blast databases for each genome will be created
in $HOME/arrayBlast/blast_db/.

The blast output from the blast against each genome
will be written be in xml format (-outfmt 5) to
STDOUT.

=head2 What if I want to blast against something
else (eg., NCBI's nr)?

Spacers and DRs can be blasted 'manually'
against the genomes in CLdb or any other 
Blast database (see 'IF
PERFORMING BLAST WITHOUT THIS SCRIPT' below).


=head2 IF PERFORMING BLAST WITHOUT THIS SCRIPT:

You can conduct the blastn run without this script
and blast against any nucleotide database (eg., nt).
Just make sure to use the blastn in the blast+ toolkit
and set the output format as '-outfmt 5' (xml output).

=head3 To basically replicate the default blast used 
in this script (blast db is 'nt' in this example):

blastn -task 'blastn-short' -db nt -query spacers.fna
-evalue 1e-3 -outfmt 5 > spacer_blast_hits.xml

=head1 EXAMPLES

=head2 Blasting all spacers against all genomes in CLdb

CLdb -- array2fasta -d CLdb.sqlite -clust > unique_spacers.fna

CLdb -- run -d CLdb.sqlite -q unique_spacers.fna > blast_hits.xml

=head2 Blasting all DRs against all genomes in CLdb

CLdb -- array2fasta -d CLdb.sqlite -r -clust > unique_DRs.fna

CLdb -- run -d CLdb.sqlite -q unique_DRs.fna > blast_hits.xml

=head2 Blasting all spacers against 1 genome (taxon_name = 'e.coli')

CLdb -- array2fasta.pl -d CLdb.sqlite -clust > unique_spacers.fna

CLdb -- run -d CLdb.sqlite -q unique_spacers.fna -taxon_name "e.coli"

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
use lib "$FindBin::RealBin/../../lib";
use CLdb::query qw/
		    table_exists
		    n_entries
		    join_query_opts/;
use CLdb::utilities qw/
			file_exists 
			connect2db
			get_file_path/;
use CLdb::arrayBlast::blast qw/
				make_blast_dir
				make_blast_db
				call_blastn_short/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $blast_params = "-evalue 1e-3";		# 1e-3
my $fork = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "query=s" => \$query_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "extra=s" => \$extra_query,
	   "blast=s" => \$blast_params,
	   "fork=i" => \$fork,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
file_exists($query_file, "query fasta");

$query_file = File::Spec->rel2abs($query_file);
my $db_path = get_file_path($database_file);	


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking fasta query (query_file); determine if repeat or spacer #
check_query_fasta($query_file);	    # query provided


# subject selection #
## joining query options (for table join) ##
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


## getting/making fasta files for all selected subject taxa ##
my $subject_loci_r = get_loci_fasta($dbh, $join_sql, $extra_query);

## making blast directory; blastdb; and blast alias##
print STDERR "Making BLAST databases...\n";
my $blast_dir = make_blast_dir($db_path);
my $blast_dbs_r = make_blast_db($subject_loci_r, $db_path, 
				$blast_dir, $verbose);

## blasting ##
print STDERR "\nBlasting sequences...\n";
call_blastn_short($blast_dbs_r, $query_file, $blast_params, $fork); 


#--- Subroutines ---#
sub get_loci_fasta{
# querying CLdb for fasta files for each select taxon to use for blastdb #
  my ($dbh, $join_sql, $extra_query) = @_;
  
  my $q = "SELECT taxon_name, taxon_id, fasta_file 
FROM loci 
WHERE locus_id=locus_id 
$join_sql 
$extra_query
GROUP BY taxon_name, taxon_id";
  
  my $res = $dbh->selectall_arrayref($q) or die $dbh->err;
  die " ERROR: no matches found to query!\n"
    unless @$res;
  
  #print Dumper @$res; exit;
  return $res;
}


sub check_query_fasta{
# loading fasta file as a hash #
  my $error = "$query_file not formatted correctly! The query input file should be formatted as: >'locusID'|'spacer/DR'|'spacer/DR_ID'|'groupID'";
  
  my ($fasta_in, $query_bool) = @_;
  open IN, $fasta_in or die $!;
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;	
    if(/>.+/){
      # check formating #
      my @l = split /\|/;
      die " ERROR: $error \n" if $query_bool && $l[1] !~ /(spacer|DR)/i;
      die " ERROR: $error \n" if $query_bool && scalar @l < 4; 			
    }
  }
  close IN;
} 

