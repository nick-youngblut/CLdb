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

Run blastn-short against one or more genomes
already in CLdb (word_size=7, reward=1, DUST=off). 
Specific genomes can be
selected by subtype ('-subtype'), taxon_name
('-taxon_name'), taxon_id ('taxon_id'), or
other sql refinements of the query ('-query').

Genomes are obtained from the 'fasta_file'
values in the loci table. Any loci without an
associated fasta file will be skipped.

Fasta files should be in the CLdb_HOME/fasta/ directory 
(copied by CLdb_loadLoci.pl).

Spacer and DR groups can be blasted at the same time (combined fasta)
Athough, it may be best to blast DRs with a higher evalue cutoff.

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
use Forks::Super;

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

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $blast_params = "-evalue 0.1";		# 1e-1
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
check_query_fasta($query_file);		


# subject selection #
## joining query options (for table join) ##
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


## getting/making fasta files for all selected subject taxa ##
my $subject_loci_r = get_loci_fasta($dbh, $join_sql, $extra_query);

## making blast directory; blastdb; and blast alias##
my $blast_dir = make_blast_dir($db_path);
my $blast_dbs_r = make_blast_db($subject_loci_r, $db_path, $blast_dir);

## blasting ##
call_blastn_short($blast_dbs_r, $query_file, $blast_params, $fork); 


#--- Subroutines ---#
sub call_blastn_short{
	my ($blast_dbs_r, $query_file, $blast_params, $fork) = @_;
		
	my $outfmt = join(" ", qw/7 qseqid sseqid pident length mismatch gapopen qstart qend sstart send
evalue bitscore qlen slen btop/);

	# sanity check #
	map{die " ERROR: cannot find $_!\n" unless -e $_} @$blast_dbs_r;

	#print Dumper @$blast_dbs_r; exit;
	my %output;
	foreach my $blast_db (@$blast_dbs_r){
		my $job = fork{
			share => [\%output],
			sub => sub{
				my $cmd = "blastn -task 'blastn-short' -query $query_file -db $blast_db -outfmt '$outfmt' $blast_params";
				open PIPE, "$cmd |" or die $!;
				while(<PIPE>){
					push @{$output{$blast_db}}, $_;
					}
				close PIPE;
				}
			}
		}
	waitall;
	
	# writing output #
	foreach my $out (sort keys %output){
		print join("", @{$output{$out}});
		}
	}

sub make_blast_db{
# foreach fasta, making a blast DB #
	my ($subject_loci_r, $db_path, $blast_dir) = @_;
	
	my $fasta_dir = "$db_path/fasta";

	# sanity check #
	die " ERROR: cannot find $fasta_dir!\n" unless -d $fasta_dir;
	die " ERROR: cannot find $blast_dir!\n" unless -d $blast_dir;
	
	# status #
	print STDERR "...making blast databases in $blast_dir\n";
	
	# making blast dbs #
	my @blastdbs;
	foreach my $row (@$subject_loci_r){
		unless($$row[2]){
			map{$_ = "" unless $_} @$row;
			print STDERR " Skipping: '", join(",", @$row), "'! 'fasta_file' value empty!\n";
			next; 		# next unless fasta file present
			}

		# sanity check #
		die " ERROR: cannnot find $fasta_dir/$$row[2]!"
			unless -e "$fasta_dir/$$row[2]"; 		
		
		# making symlink in blast directory #
		unless(-e "$blast_dir/$$row[2]" || -l "$blast_dir/$$row[2]"){
			symlink("$fasta_dir/$$row[2]", "$blast_dir/$$row[2]") or die $!;
			}
		
		# making blast db #
		my $cmd = "makeblastdb -dbtype nucl -parse_seqids -in $blast_dir/$$row[2]";
		print STDERR "$cmd\n" unless $verbose;
		`$cmd`;
		
		push @blastdbs, "$blast_dir/$$row[2]";
		}
		
	# making blast alias #
#	my $cmd = join("", "blastdb_aliastool -dblist \"",
#			join(" ", @blastdbs),
#			"\" -dbtype nucl -out $blast_dir/CLdb_genomes -title 'CLdb_genomes'");
#	`$cmd`;
	
	return \@blastdbs;
	}

sub get_loci_fasta{
# querying CLdb for fasta files for each select taxon to use for blastdb #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	my $q = "SELECT taxon_name, taxon_id, fasta_file 
FROM loci 
WHERE locus_id=locus_id 
$join_sql 
GROUP BY taxon_name, taxon_id";
	
	$q .= " $extra_query" if $extra_query;
	
	my $res = $dbh->selectall_arrayref($q);
	die " ERROR: no matches found to query!\n"
		unless @$res;

		#print Dumper @$res; exit;
	return $res;
	}

sub get_database_path{
	my $database_file = shift;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return $parts[1];
	}

sub make_blast_dir{
	my $db_path = shift;
	
	my $dir = join("/", $db_path, "spacer_blast");
	mkdir $dir unless -d $dir;

	return $dir;
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

