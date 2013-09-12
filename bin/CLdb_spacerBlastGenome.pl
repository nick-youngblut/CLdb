#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Bio::SeqIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file);
my (@subtype, @taxon_id, @taxon_name);
my @subject_in;
my $blast_params = "-evalue 0.00001";		# 1e-5
my $extra_query = "";
my $range = 30;		# spacer-DR blast hit overlap (bp)
my $Ncpu = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "fasta=s" => \$query_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "blast=s" => \$blast_params,
	   "num_thread=i" => \$Ncpu,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: provdie a fasta of spacers/DRs to use as a BLAST query\n"
	unless $query_file;
die " ERROR: can't find $query_file\n"
	unless -e $query_file;

	
my $db_path = get_database_path($database_file);
$query_file = File::Spec->rel2abs($query_file);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# checking fasta query (query_file); determine if repeat or spacer #
my $spacer_DR = check_query_fasta($query_file);		# spacer|DR

# subject selection #
## making blast directory ##
my $blast_dir = make_blast_dir($db_path);

## joining query options (for table join) ##
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

## getting/making fasta files for all selected subject taxa ##
my $subject_loci_r = get_loci_fasta_genbank($dbh, $join_sql);
my $fasta_dir = make_fasta_dir($db_path);
genbank2fasta($subject_loci_r, $fasta_dir);

## making blast dbs from subject fasta files (foreach subject) ##


# blasting (foreach subject) #
## blastn & loading blastn output ##



### Subroutines
sub genbank2fasta{
# getting fasta from genbank unless -e fasta #
	my ($subject_loci_r, $fasta_dir) = @_;
	
	foreach my $loci (@$subject_loci_r){
		if(! $$loci[3]){		# if no fasta file
			print STDERR " WARNING: no fasta for taxon_name->$$loci[0] taxon_id->$$loci[1]! Trying to extract sequence from genbank...\n";
			$$loci[3] = genbank2fasta_extract($$loci[2], $fasta_dir);
			}
		elsif($$loci[3] && ! -e $$loci[3]){
			print STDERR " WARNING: cannot find $$loci[3]! Trying to extract sequence from genbank...\n";
			$$loci[3] = genbank2fasta_extract($$loci[2], $fasta_dir);
			}
		}
	}
	
sub genbank2fasta_extract{
	my ($genbank_file, $fasta_dir) = @_;
	# sanity check #
	die " ERROR: cannot find $genbank_file!\n"
		unless -e "$db_path/genbank/$genbank_file";
		
	my $seqio = Bio::SeqIO->new(-file => "$db_path/genbank/$genbank_file", -format => "genbank");
	
	# output #
	my @parts = File::Spec->splitpath($genbank_file);
	$parts[2] =~ s/\.[^.]+$|$/.fasta/;
	open OUT, ">$fasta_dir/$parts[2]" or die $!;
	
	# writing #
	while(my $seqo = $seqio->next_seq){
		# seqID #
		my $scafID = $seqo->display_id;
		print OUT ">$scafID\n";
			for my $feato (grep { $_->primary_tag eq 'source' } $seqo->get_SeqFeatures){
			print OUT $feato->seq->seq, "\n";
			}
		}
	close OUT;
	}

sub get_loci_fasta_genbank{
# querying CLdb for fasta & genbank files for each taxon #
	my ($dbh, $join_sql) = @_;
	
	my $q = "SELECT taxon_name, taxon_id, genbank_file, fasta_file FROM loci $join_sql GROUP BY taxon_name, taxon_id";
	
	my $res = $dbh->selectall_arrayref($q);
	die " ERROR: no matches found to query!\n"
		unless @$res;

		# print Dumper @$res;
	return $res;
	}

sub get_database_path{
	my $database_file = shift;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return $parts[1];
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND $cat IN (", join(", ", @$vals_r), ")");
	}

sub make_fasta_dir{
	my $db_path = shift;
	
	my $dir = join("/", $db_path, "fasta");
	mkdir $dir unless -d $dir;

	return $dir;
	}

sub make_blast_dir{
	my $db_path = shift;
	
	my $dir = join("/", $db_path, "spacer_blast");
	mkdir $dir unless -d $dir;

	return $dir;
	}

sub check_query_fasta{
# checking that fasta is of spacers & DRs from CLdb #
# determining if spacer or DR #
	my ($query_file) = @_;

	my $error = "$query_file not formatted correctly! The query input file should contain grouped spacers or DRs. (use '-g' with CLdb_array2fasta)";

	open IN, $query_file or die $!;
	
	my %spacer_DR;
	while(<IN>){
		if(/^>/){
			die " ERROR: $error\n" unless /(spacer|DR)_group\./i;
			
			if(/spacer/i){ $spacer_DR{"spacer"} = 1; }
			elsif(/DR/i){ $spacer_DR{"DR"} = 1;}
			else{ $spacer_DR{"NA"} = 1;}
			}
		}
	close IN;
	
	die " ERROR: $error\n" if exists $spacer_DR{"NA"} || scalar keys %spacer_DR > 1;
	
	foreach my $key (keys %spacer_DR){
		return $key; 
		last;
		}
	}


__END__

=pod

=head1 NAME

CLdb_spacerBlast.pl -- wrapper for spacer blasting

=head1 SYNOPSIS

CLdb_spacerBlast.pl [flags]

=head2 Required flags

=over

=item -database

CLdb database.

=item -subject

Either subject file or 1-3 arguments (see DESCRIPTION).

=back

=head2 Optional flags

=over

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine which sequences are returned.

=item -blast

BLASTn parameters (besides required flags). [-evalue 0.00001]

=item -range

Range allowable between spacer & DR blast hit (bp). [30]

=item -cpu

Use for '-num_threads' parameter in BLASTn. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlast.pl

=head1 DESCRIPTION

A wrapper around other CLdb scripts for blasting all or
a subset of spacer groups.

=head2 The script's procedure is:

=over

=item * Select spacer & direct repeat (DR) groups (can refine to particular taxa & subtypes)

=item * BLASTn-short of spacer & DR groups against provide subjects (e.g. genomes)

=item * Determining which spacer blast hits are hitting CRISPR arrays (if DRs hit adjacent)

=item * Adding the blast hits and subjects (sequences with a hit) to CLdb

=back

=head2 '-subject' flag

Provide either a subject fasta file and a Taxon_ID and/or a Taxon_Name,
or provide a tab-delimited file (3 columns) with subject fasta files and
a Taxon_IDs and/or Taxon_Names (columns: fasta_file, Taxon_id, Taxon_Name).

Example1 "-subject ecoli.fna 666666.452 escherichia_coli"

Example2 "-subject ecoli.fna '' escherichia_coli"

The Taxon_IDs and Taxon_names can be used, for example, to see if CRISPR
spacers are hitting other places in the same genome (and not in other CRISPR arrays).

=head1 EXAMPLES

=head2 Spacer blast of subtype I-B I-C spacers against Ecoli

CLdb_spacerBlast.pl -da CLdb.sqlite -subtype I-B I-C -subject ecoli.fna 666666.452 Escherichia_coli

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

