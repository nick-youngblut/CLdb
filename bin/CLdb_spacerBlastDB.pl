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


my ($verbose, $database_file, $query_file, $blast_db);
my $extra_query = "";
my $blast_params = "-evalue 0.00001";		# 1e-5
my $num_threads = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "fasta=s" => \$query_file,
	   "db=s" => \$blast_db,
	   "query=s" => \$extra_query,
	   "blast=s" => \$blast_params,
	   "num_threads=i" => \$num_threads,
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
check_query_fasta($query_file);		# spacer|DR

# blasting (foreach subject blast DB) #
## blastn & loading blastn output ##
blastn_call_load($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file);


### Subroutines
sub blastn_call_load{
# calling blastn, loading blast results directly into DB #
	my ($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file) = @_;
	
	# status #
	print STDERR "...BLASTing each subject genome database\n";
	
	# preparing sql #
	my @blast_hits_col = qw/blast_id spacer_DR S_taxon_name S_taxon_ID Group_ID sseqid pident len mismatch gapopen qstart qend sstart send evalue bitscore slen/;
		
	# blasting each subject genome #
	my $insert_cnt = 0;				# summing number of entries for all subjects 
	foreach my $row (@$subject_loci_r){
		next unless $$row[3];			# skipping if no fasta & thus no blast db
		my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen";
		my $cmd = "blastn -task 'blastn-short' -outfmt '$outfmt' -db $blast_dir/$$row[3] -query $query_file -num_threads $num_threads $blast_params |";

		# blasting and loading db #
		## blast ##
		open PIPE, $cmd or die $!;
		while(<PIPE>){
			chomp;
			next if /^\s*$?/;
			
			my @line = split /\t/;
			
			# spacer or DR? #
			(my $spacer_DR = $line[0]) =~ s/_.+//;
			$line[0] =~ s/.+\.//;
			
			# adding input values #
			unshift @line, ($spacer_DR, $$row[0], $$row[1]);
			my $blast_id = join("_", @line[(0..4,11..12)] );
			unshift @line, $blast_id;
			
			# quoting #
			$line[0] = $dbh->quote($line[0]);
			$line[1] = $dbh->quote($line[1]);
			$line[2] = $dbh->quote($line[2]);
			$line[5] = $dbh->quote($line[5]);
			
			# making command #
			my $sql = join(" ", "INSERT INTO Blast_hits( ",
					join(",", @blast_hits_col),
					") VALUES( ",
					join(",", @line), 
					")" );
			
			# loading line #
			$dbh->do( $sql );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @line), "\n";
				}
			else{ $insert_cnt++; }
			}
		close PIPE;
		}
	$dbh->commit();
	
	print STDERR "...Number of blast entries added to CLdb: $insert_cnt\n";
	}

sub make_blast_db{
# foreach fasta, making a blast DB #
	my ($subject_loci_r, $fasta_dir, $blast_dir) = @_;
	
	# sanity check #
	die " ERROR: cannot find $fasta_dir!\n" unless -d $fasta_dir;
	die " ERROR: cannot find $blast_dir!\n" unless -d $blast_dir;
	
	# status #
	print STDERR "...making blast databases in $blast_dir\n";
	
	# making blast dbs #
	foreach my $row (@$subject_loci_r){
		next unless $$row[3]; 		# next unless fasta file present

		# sanity check #
		die " ERROR: cannnot find $fasta_dir/$$row[3]!"
			unless -e "$fasta_dir/$$row[3]"; 		
		
		# making symlink in blast directory #
		unless(-e "$blast_dir/$$row[3]" || -l "$blast_dir/$$row[3]"){
			symlink("$fasta_dir/$$row[3]", "$blast_dir/$$row[3]") or die $!;
			}
		
		# making blast db #
		my $cmd = "makeblastdb -dbtype nucl -in $blast_dir/$$row[3]";
		`$cmd`;		
		}
	}

sub update_loci{
# updating loci w/ fasta file info for newly written fasta #
	my ($dbh, $subject_loci_r) = @_;
	
	# status#
	print STDERR "...updating Loci table with newly made fasta files\n";
	
	# sql #
	my $q = "UPDATE loci SET fasta_file=? WHERE genbank_file=?";
	my $sth = $dbh->prepare($q);
	
	foreach my $row (@$subject_loci_r){
		next unless $$row[3]; 						# if no fasta_file; skipping
		my @parts = File::Spec->splitpath($$row[3]);

		$sth->bind_param(1, $parts[2]);				# fasta_file (just file name)
		$sth->bind_param(2, $$row[2]);				# array_file
		$sth->execute( );
		
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$row), "\n";
			}
		}
	$dbh->commit();
	
	print STDERR "...updates committed!\n";
	}

sub genbank2fasta{
# getting fasta from genbank unless -e fasta #
	my ($subject_loci_r, $fasta_dir) = @_;
	
	my $update_cnt = 0;
	foreach my $loci (@$subject_loci_r){
		if(! $$loci[3]){		# if no fasta file
			print STDERR " WARNING: no fasta for taxon_name->$$loci[0] taxon_id->$$loci[1]! Trying to extract sequence from genbank...\n";
			$$loci[3] = genbank2fasta_extract($$loci[2], $fasta_dir);
			$update_cnt++;
			}
		elsif($$loci[3] && ! -e "$fasta_dir/$$loci[3]"){
			print STDERR " WARNING: cannot find $$loci[3]! Trying to extract sequence from genbank...\n";
			$$loci[3] = genbank2fasta_extract($$loci[2], $fasta_dir);
			$update_cnt++;
			}
		}
	return $update_cnt;		# if >0; update loci tbl
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
	my $fasta_out = "$fasta_dir/$parts[2]";
	open OUT, ">$fasta_out" or die $!;
	
	# writing #
	my $seq_cnt = 0;
	while(my $seqo = $seqio->next_seq){
		$seq_cnt++;
		
		# seqID #
		my $scafID = $seqo->display_id;
		print OUT ">$scafID\n";
			for my $feato (grep { $_->primary_tag eq 'source' } $seqo->get_SeqFeatures){
			print OUT $feato->seq->seq, "\n";
			}
		}
	close OUT;
	
	# if genome seq found and fasta written, return fasta #
	if($seq_cnt == 0){
		print STDERR " WARNING: no genome sequnece found in Genbank file: $genbank_file!\nSkipping BLAST!\n";
		unlink $fasta_out;
		return 0;
		}
	else{ 
		print STDERR "...fasta file extracted from $genbank_file written: $fasta_out\n";
		return $fasta_out; 
		}
	}

sub get_loci_fasta_genbank{
# querying CLdb for fasta & genbank files for each taxon #
	my ($dbh, $join_sql) = @_;
	
	my $q = "SELECT taxon_name, taxon_id, genbank_file, fasta_file FROM loci WHERE locus_id=locus_id $join_sql GROUP BY taxon_name, taxon_id";
	
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
		if (/^>/){
			die " ERROR: $error\n" unless /(spacer|DR)_group\./i;
			}
		}
	close IN;
	
	die " ERROR: $error\n" if exists $spacer_DR{"NA"};
	}


__END__

=pod

=head1 NAME

CLdb_spacerBlastGenome.pl -- BLASTn-short of spacers and/or DRs against a genome in CLdb

=head1 SYNOPSIS

CLdb_spacerBlastGenome.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=item -fasta  <char>

A fasta of either spacer or DR group sequences (use: CLdb_array2fasta.pl -g)

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -blast  <char>

BLASTn parameters (besides required flags). [-evalue 0.00001]

=item -num_threads  <int>

Use for '-num_threads' parameter in BLASTn. [1]

=item -v  <bool>

Verbose output

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlastGenome.pl

=head1 DESCRIPTION

Run blastn-short against one or more genomes
already in CLdb. Specific genomes can be
selected by subtype ('-subtype'), taxon_name
('-taxon_name'), taxon_id ('taxon_id'), or
other sql refinements of the query ('-query').

Genomes are obtained from the 'fasta_file'
values in the loci table.
Fasta files should be in the CLdb/fasta/ directory 
(copied by CLdb_loadLoci.pl).
If no fasta file is found, the genome sequence
is extracted from the genbank file ('genbank_file'
in loci table). The genome is skipped if still 
no sequence can be found.

Spacer and DR groups can be blasted at the same time.

It may be best to blast all DRs at once for 
subsequent spacer array screening (spacers that
have adjacent DR hits are considered to be located
in CRISPR arrays and are thus not protospacers).

=head1 EXAMPLES

=head2 Blasting all spacers against all genomes in CLdb

CLdb_array2fasta.pl -g > all_spacer_groups.fna

CLdb_spacerBlastGenome.pl -d CLdb.sqlite -f all_spacer_groups.fna

=head2 Blasting all DRs against all genomes in CLdb

CLdb_array2fasta.pl -r -g > all_DR_groups.fna

CLdb_spacerBlastGenome.pl -d CLdb.sqlite -f all_DR_groups.fna

=head2 Blasting all spacers against 1 genome

CLdb_array2fasta.pl -g > all_spacer_groups.fna

CLdb_spacerBlastGenome.pl -d CLdb.sqlite -f all_spacer_groups.fna -taxon_name "e.coli"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

