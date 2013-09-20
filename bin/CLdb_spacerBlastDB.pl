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
use Bio::DB::GenBank;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file, $blast_db);
my $extra_query = "";
my $blast_params = "-evalue 0.001";		# 1e-3
my $num_threads = 1;		
my $extend = 20;						# number of bp extended beyond hit
GetOptions(
	   "database=s" => \$database_file,
	   "fasta=s" => \$query_file,
	   "db=s" => \$blast_db,
	   "query=s" => \$extra_query,
	   "blast=s" => \$blast_params,
	   "num_threads=i" => \$num_threads,
	   "xtend=i" => \$extend,
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
unless($blast_db){
	print STDERR "...WARNING: no BLAST database provided. Using 'nt' by default.\n";
	$blast_db = "nt";	# genbank's nt db by default
	}


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# checking fasta query (query_file); determine if repeat or spacer #
check_query_fasta($query_file);		# spacer|DR

# blasting (foreach subject blast DB) #
## blastn & loading blastn output ##
blastn_call_load($dbh, $blast_db, $blast_params, $num_threads, $query_file);


### Subroutines
sub blastn_call_load{
# calling blastn, loading blast results directly into DB #
	my ($dbh, $blast_db, $blast_params, $num_threads, $query_file) = @_;
	
	# preparing sql #
	my @blast_hits_col = qw/blast_id spacer_DR S_accession S_GI Group_ID sseqid pident len mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq frag xstart xend/;
		
	my $outfmt = "6 sacc sgi qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq";
	my $cmd = "blastn -task 'blastn-short' -outfmt '$outfmt' -db $blast_db -query $query_file -num_threads $num_threads $blast_params |";
	print STDERR "$cmd\n" unless $verbose;

	# genbank connect object #
	my $gb = new Bio::DB::GenBank;
	
	# blasting and loading db #
	## blast ##
	open PIPE, $cmd or die $!;
	my %insert_cnt;
	while(<PIPE>){
		chomp;
		next if /^\s*$?/;
		
		my @line = split /\t/;
		
		# spacer or DR? #
		(my $spacer_DR = $line[2]) =~ s/_.+//;
		$line[2] =~ s/.+\.//;			# spacer group
		
		# adding input values #
		unshift @line, $spacer_DR;		# adding spacer_DR
		my $blast_id = join("_", @line[(0..4,12..13)] );
		unshift @line, $blast_id;		# adding blast_id	

		# frag & extension #
		## getting frag of sequence from genbank ##
		push @line, get_frag($gb, $line[2], $line[5], $line[12], $line[13], $line[17]);	# gb-object, accession, start, end, slen
		if($spacer_DR eq "DR"){	# DR -> no frag or qseq needed 
			$line[18] = ""; 
			$line[19] = "";
			}

		# quoting #
		$line[0] = $dbh->quote($line[0]);		# blast_ID
		$line[1] = $dbh->quote($line[1]);		# spacer_dr
		$line[2] = $dbh->quote($line[2]);		# S_accession
		$line[3] = $dbh->quote($line[3]);		# S_GI
		$line[5] = $dbh->quote($line[5]);		# sseqid
		$line[18] = $dbh->quote($line[18]);		# qseq
		$line[19] = $dbh->quote($line[19]);		# frag
	
		# blast_hits #
		## making blast_hit sql ##
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
		else{ $insert_cnt{$spacer_DR}++; }
		}
		
	close PIPE;
	
	$dbh->commit(); 

	$insert_cnt{"Spacer"} = 0 unless exists $insert_cnt{"Spacer"};
	print STDERR "...Number of spacer blast entries added to CLdb: $insert_cnt{'Spacer'}\n";
	$insert_cnt{"DR"} = 0 unless exists $insert_cnt{"DR"};
	print STDERR "...Number of DR blast entries added to CLdb: $insert_cnt{'DR'}\n";
	}
	
sub get_frag{
	my ($gb, $acc, $sseqid, $start, $end, $slen) = @_;
	
	my $seqio = $gb->get_Stream_by_acc( $acc );
	
	while( my $seq = $seqio->next_seq ) {
		if($start <= $end){
			$start -= $extend;
			$start = 0 if $start < 0;
			$end += $extend;
			$end = $slen if $end > $slen;
			return (substr($seq->seq(), $start, $end-$start), $start, $end);
			}
		else{
			$end -= $extend;
			$end = 0 if $end < 0;
			$start += $extend;
			$start = $slen if $start > $slen;
			return (revcomp(substr($seq->seq(), $end, $start-$end)), $start, $end);
			}
		}
	}
	
sub revcomp{
	# reverse complements DNA #
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/ACGTNBVDHKMRYSW\.-/TGCANVBHDMKYRSW\.-/;
	return $seq;
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

CLdb_spacerBlastGenome.pl -- BLASTn-short of spacers and/or DRs against 'nt' or other BLAST DB

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

=item -db  <char>

BLAST database. [nt]

=item -blast  <char>

BLASTn parameters (besides required flags). [-evalue 0.001]

=item -num_threads  <int>

Use for '-num_threads' parameter in BLASTn. [1]

=item -x  <int>

Retain blast hit subject sequence fragment + '-x' bp on either side of fragment. [20]

=item -v  <bool>

Verbose output

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlastGenome.pl

=head1 DESCRIPTION

Run blastn-short against a BLAST database
(word_size=7, reward=1, DUST=off).

Accension numbers from BLAST hits are used
to fetch sequences from GenBank.

Spacer and DR groups can be blasted at the same time.

It may be best to blast all DRs at once for 
subsequent spacer array screening (spacers that
have adjacent DR hits are considered to be located
in CRISPR arrays and are thus not protospacers).

=head1 EXAMPLES

=head2 Blasting all spacers against nt

CLdb_spacerBlastGenome.pl -d CLdb.sqlite

=head2 Blasting all DRs against nt

CLdb_spacerBlastGenome.pl -d CLdb.sqlite 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

