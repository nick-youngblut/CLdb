#!/usr/bin/env perl

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

BLASTn parameters (besides required flags). [-evalue 0.001]

=item -num_fork  <int>

Number of parallel blastn calls (one for each subject genome). [1]

=item -overlap  <bool>

Filter out hits falling within a large blast hit (same query-subject? [TRUE]

=item -x  <int>

Retain blast hit subject sequence fragment + '-x' bp on either side of fragment. [20]

=item -v  <bool>

Verbose output. [TRUE]

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlastGenome.pl

=head1 DESCRIPTION

Run blastn-short against one or more genomes
already in CLdb (word_size=7, reward=1, DUST=off). 
Specific genomes can be
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

Any gaps in the spacer-protospacer alignment are include.

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


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Bio::SearchIO;
use Bio::AlignIO;
use Set::IntervalTree;
use Parallel::ForkManager; 

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts
	get_array_seq
	get_arrays_seq_byLeader/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	get_file_path/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file, $entire, $filter_overlap_hits);
my (@subtype, @taxon_id, @taxon_name);
my @subject_in;
my $extra_query = "";
my $blast_params = "-evalue 0.001";		# 1e-3
my $num_threads = 1;
my $extend = 10;
my $fork = 0;
GetOptions(
	   "database=s" => \$database_file,
	   "fasta=s" => \$query_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,		
	   "blast=s" => \$blast_params,
	   "num_fork=i" => \$fork,
	   "xtend=i" => \$extend,
	   "overlap" => \$filter_overlap_hits, 		# filter out overlapping hits? [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
file_exists($database_file, "database");
file_exists($query_file, "fasta");

my $db_path = get_file_path($database_file);
$query_file = File::Spec->rel2abs($query_file);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# checking fasta query (query_file); determine if repeat or spacer #
#check_query_fasta($query_file);		# spacer|DR
my $query_r = load_fasta($query_file, 1);

# subject selection #
## making blast directory ##
my $blast_dir = make_blast_dir($db_path);

## joining query options (for table join) ##
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

print Dumper $join_sql; exit;


## getting/making fasta files for all selected subject taxa ##
my $subject_loci_r = get_loci_fasta_genbank($dbh, $join_sql);

my $fasta_dir = make_fasta_dir($db_path);
my $loci_update = genbank2fasta($subject_loci_r, $fasta_dir);
update_loci($dbh, $subject_loci_r) if $loci_update;

## making blast DBs from subject fasta files (foreach subject) ##
make_blast_db($subject_loci_r, $fasta_dir, $blast_dir);

# blasting (foreach subject blast DB) #
## blastn & loading blastn output ##
my $blast_out_r = blastn_xml_call($dbh, $subject_loci_r, $blast_dir, 
		$blast_params, $num_threads, $query_file, $query_r);
		
## loading blast output into db ##
load_blast_output($dbh, $blast_out_r);


# commit & exit #
$dbh->commit;
#print Dumper "NO commit!\n";
exit;


### Subroutines
sub load_blast_output{
# loading blast output into CLdb #
	my ($dbh, $blast_out_r) = @_;

	# status #
	print STDERR "...loading blast results into CLdb\n" unless $verbose;

	# columns to load into blast_hits DB #
	my @blast_hits_col = qw/
spacer_DR 
Group_ID 
S_taxon_ID 
S_taxon_name
S_accession
sseqid 
pident 
mismatch 
gaps 
evalue 
bitscore
strand 
len 
qlen 
slen 
qseq 
qstart 
qend 
sstart 
send 
qseq_full 
sseq_full 
qseq_full_start 
qseq_full_end 
sseq_full_start 
sseq_full_end 
proto3px 
proto3px_start 
proto3px_end 
proto5px 
proto5px_start 
proto5px_end 
sseq/;

	my %insert_cnt; 		# summing number of entries added/updated in CLdb
	foreach my $infile (@$blast_out_r){
		print STDERR "...loading $infile\n";
		die " ERROR: cannot find $infile!\n" unless -e $infile;
		open IN, $infile or die $!;
		while(<IN>){
			chomp;
			next if /^\s*$/;
			
			my @line = split /\t/;
			(my $spacer_DR = $line[0]) =~ s/'//g;
			map{$_ = "" if $_ eq "NULL"} @line;
			
			# making sql for entry #
			my $sql;
			if($spacer_DR eq "DR"){		# DR: no need for 
				$sql = join(" ", "INSERT INTO Blast_hits( ",
						join(",", @blast_hits_col[0..21]),
						") VALUES( ",
							join(",", @line[0..21]), ")" );
						}
			elsif($spacer_DR eq "Spacer"){
				$sql = join(" ", "INSERT INTO Blast_hits( ",
						join(",", @blast_hits_col),
						") VALUES( ",
							join(",", @line), ")" );
						}
			else{ die " LOGIC ERROR: neither 'DR' nor 'Spacer' found in spacer_DR column in $infile!\n", join(",", @line), "\n"; }
				
			# loading blast hit into CLdb #
			$dbh->do( $sql );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @line), "\n";
				}
			else{ $insert_cnt{$spacer_DR}++; }	
			}
		close IN;
		unlink $infile or die $!;			# delete blast file
		}
	
	# summary of CLdb upload #	
	$insert_cnt{"Spacer"} = 0 unless exists $insert_cnt{"Spacer"};
	print STDERR "...Number of spacer blast entries added/updated to CLdb: $insert_cnt{'Spacer'}\n";
	$insert_cnt{"DR"} = 0 unless exists $insert_cnt{"DR"};
	print STDERR "...Number of DR blast entries added/updated to CLdb: $insert_cnt{'DR'}\n";
	}

sub blastn_xml_call{
# calling blastn (xml output format); parsing output #
## needed to get blast alignment (find potential gaps) ##
## forking blast runs ##
	my ($dbh, $subject_loci_r, $blast_dir, $blast_params, 
		$num_threads, $query_file, $query_r) = @_;
	
	# status #
	print STDERR "...blasting spacers & direct repeats\n" unless $verbose;
	
	# forking #
	my $pm = new Parallel::ForkManager($fork);
	
	
	# making blast output file list # needed to prevent forking errors #
	my @blast_out_files;
	foreach my $row (@$subject_loci_r){
		my $outfile = "$blast_dir/$$row[3].blast.txt";
		push @blast_out_files, $outfile; 
		}
	
	# blasting each subject genome #
	foreach my $row (@$subject_loci_r){		# taxon_name, taxon_id, genbank, fasta
		# forking #
		my $pid = $pm->start and next;
		
		# skipping if no fasta & thus no blast db #
		unless ($$row[3]){
			$pm->finish;
			next;			
			}
		
		# output file #
		my $outfile = "$blast_dir/$$row[3].blast.txt";
		open OUT, ">$outfile" or die $!;
		print STDERR "...writing blast results to: $outfile\n";
		
		# loading fasta #
		my $fasta_r = load_fasta("$blast_dir/$$row[3]");

		# setting blast cmd #
		my $cmd = "blastn -task 'blastn-short' -outfmt 5 -db $blast_dir/$$row[3] -query $query_file -num_threads $num_threads $blast_params |";
		print STDERR "$cmd\n" unless $verbose;
		
		# parsing xml output w/ bioperl #
		my $seqin = Bio::SearchIO->new( -file => $cmd, -format => "blastxml" );
		
		
		# extracting info to load into table & loading #
		while( my $result = $seqin->next_result ) {		# result object 
			
			# getting all start-end for each hit to find overlapping hits #
			#my $hit_se_r = get_hit_se($result) unless $filter_overlap_hits;		# hit_se_r=>spacer/DR_group=>scaffold=>[sstart, send]
			my $hit_itree_r = get_hit_itrees($result) unless $filter_overlap_hits;
			$result->rewind;		# reset iterator
			
			while( my $hit = $result->next_hit ) {		# hit object
				$hit->rewind; 							# reset iterator
				
				while( my $hsp = $hit->next_hsp ) {		# hsp object					
					# query sequence ID #
					my $qseqid = $result->query_description;
					
					# spacer or DR? #
					(my $spacer_DR = $qseqid) =~ s/_.+//;
					$qseqid =~ s/.+\.//;
					
					
					# blastID # spacer_DR qseqid S_taxon_name S_taxon_ID sseqid sstart send #
					#my $blast_id = join("_", $spacer_DR,	# query ID (spacer/DR group)(
					#				$$row[0],				# S_taxon_name
					#				$$row[1],				# S_taxon_ID
					#				$hit->name, 			# sseqid
					#				$hsp->start('hit'),		# sstart
					#				$hsp->end('hit')		# send
					#				);
					
					# filtering overlapping hits #
					next if filter_overlapping_hits($result, $hit, $hsp, $hit_itree_r);
						#print Dumper filter_overlapping_hits($hit, $hsp, $qseqid, $hit_itree_r);
					
					# getting alignment #
					my $alnIO = $hsp->get_aln; 
					my @seqs = $alnIO->each_seq;
					
					# getting full length & proto extension sequence info if Spacer #
					my ($qseq_full, $sseq_full) = ("", "");
					my ($qseq_full_start, $qseq_full_end);
					my ($sseq_full_start, $sseq_full_end);
					my ($proto3px, $proto3px_start, $proto3px_end);
					my ($proto5px, $proto5px_start, $proto5px_end);
					if($spacer_DR eq "Spacer"){
						# getting full length spacer (qseq); keeping any gaps in blast alignment #
						($qseq_full, $qseq_full_start, $qseq_full_end) = 
							get_full_qseq($result, $hit, $hsp, \@seqs);
							#print Dumper ($qseq_full, $qseq_full_start, $qseq_full_end);
						
						# getting full length protospacer sequence; keeping any gaps in blast alignment #
						($sseq_full, $sseq_full_start, $sseq_full_end) = 
							get_full_sseq($result, $hit, $hsp, \@seqs, $fasta_r->{$hit->name}, 
								$qseq_full_start, $qseq_full_end);
							#print Dumper ($sseq_full, $sseq_full_start, $sseq_full_end);
						
						# getting protospacer 3' & 5' extension #
						## if spacer hits + strand; proto on - strand; proto 3' is on left on + strand
						## here, calling proto3x & proto5x by + strand hit 
						## need flip for - strand hits (proto 5' actually on left because proto on + strand)
						($proto3px, $proto3px_start, $proto3px_end) = 
							get_proto_threep($hsp, \@seqs, $fasta_r->{$hit->name}, 
									 		$sseq_full_start, $extend);
							
						($proto5px, $proto5px_start, $proto5px_end) = 
							get_proto_fivep($hit, $hsp, \@seqs, $fasta_r->{$hit->name}, 
									 		$sseq_full_end, $extend);
							
						
							#print Dumper join(" ", "strand: ", $hsp->strand('hit'));
						
						## flipping 5' & 3' proto extension if blast hit to - strand (proto then on + strand)
						($proto5px, $proto5px_start, $proto5px_end, $proto3px, $proto3px_start, $proto3px_end) =
						 	(revcomp($proto3px), $proto3px_start, $proto3px_end,
						 	 revcomp($proto5px), $proto5px_start, $proto5px_end) if $hsp->strand('hit') == -1;
						}
					
						#print Dumper ($proto3px, $proto3px_start, $proto3px_end);
						#print Dumper ($proto5px, $proto5px_start, $proto5px_end);
						
					# loading CLdb #
					## taxon_name & taxon_id = "" if undef #
					map{$_ = "" unless $_ } (@$row[0..1]);
					
					## making blast_hit sql ##
					### start - end encode bioperl-style!!
					my @hit_matches = $hsp->matches('hit');
					my @vals = (
						#$dbh->quote($blast_id),		# blastID
						$dbh->quote($spacer_DR),	# Spacer|DR
						$qseqid,					# groupID
						$dbh->quote($$row[1]),		# taxon_id
						$dbh->quote($$row[0]),		# taxon_name
						$dbh->quote(""),			# S_accession (must be included as '' so not NULL)
						$dbh->quote($hit->name), 	# sseqid
						$hsp->percent_identity,
						$hsp->length('total') - $hit_matches[0] - $hsp->gaps, 		# mismatches (length - matches - gaps)
						$hsp->gaps,					# number of gaps
						$hsp->evalue,
						$hsp->bits,
						$hsp->strand('hit'), 				# strand (encoded like bioperl)
						$hsp->length('total'),		# length (alignment length)
						$hit->query_length,			# qlen
						$hit->hit_length,			# slen
						$dbh->quote($seqs[0]->seq),				# qseq
						$hsp->start('query'),		# qstart
						$hsp->end('query'),			# qend
						$hsp->start('hit'),			# sstart		(>send if strand is '-')
						$hsp->end('hit'),			# send 
						$dbh->quote($qseq_full),	# query sequence (full length)   #$dbh->quote($seqs[0]->seq)
						$dbh->quote($sseq_full),		
						$qseq_full_start,
						$qseq_full_end,
						$sseq_full_start,
						$sseq_full_end,
						$dbh->quote($proto3px),
						$proto3px_start,
						$proto3px_end,
						$dbh->quote($proto5px),
						$proto5px_start,
						$proto5px_end,
						$dbh->quote($seqs[1]->seq),		# sseq
						);
						
					# writing out file #
					map{$_ = "" if ! defined $_ } @vals;
					print OUT join("\t", @vals), "\n";
					
					}
				}
			}
		# closing output for blast #
		close OUT;

		$pm->finish;
		}
	$pm->wait_all_children;		# all forks must complete
	
	return \@blast_out_files;
	}
	
sub get_proto_fivep{
# 5' of protospacer orientation (if spacer hit to + strand) #
	my ($hit, $hsp, $seqs_r, $scaf_seq, $sseq_full_end, $extend) = @_;
	
	my $slen = $hit->hit_length;
	
	# extension beyond protospacer (adding to protospacer sequence)  #
	my $xend  = $sseq_full_end + $extend;			# 3' extension beyond proto
	$xend = $slen if $xend > $slen; 				# ceiling of scaffold length

	# getting sequence #
	my $fivep = substr($scaf_seq, $sseq_full_end,  $xend - $sseq_full_end); 	# 5' extension
	
	#unless($fivep){
	#	print Dumper
	#		$sseq_full_end, 
	#		$xend,
	#		$slen;
	#	exit(1);
	#	}
	
	return $fivep, $sseq_full_end + 1, $xend;
	}
	
sub get_proto_threep{
# 3' of protospacer orientation (if spacer hit to + strand) #
	my ($hsp, $seqs_r, $scaf_seq, $sseq_full_start, $extend) = @_;
	
	# extension beyond protospacer (adding to protospacer sequence)  #
	my $xstart = $sseq_full_start - $extend;	# 5' extension beyond proto
	$xstart = 1 if $xstart < 1;					# floor of 1
	
	# getting sequence 
	my $threep = substr($scaf_seq, $xstart -1, $sseq_full_start - $xstart);	 # 3' extension; amost forgot 0-indexing!
	
	return $threep, $xstart, $sseq_full_start - 1;
	}

sub get_full_sseq{
# 'completing spacer' if only partial hit #
# getting full length of protospacer (sseq ne sseq_full only if partial blast hit) #
## spacer-protospacer alignment retained (ie. any gaps) ##
	my ($result, $hit, $hsp, $seqs_r, $scaf_seq, 
		$qseq_full_start, $qseq_full_end) = @_;

	# max bp that start or end can actually be extended #	
	my ($min_sstart, $max_send) = max_extend($hit, $hsp);
	
	
	# actual start - end (indexing by 1) #
	my ($threep, $fivep) = ("", "");
	my ($act_start, $act_end);
	if( $hsp->strand('hit') == 1 ){		# + strand hit 
		## 3' end of proto ##
		$act_start = $hsp->start('hit') - ($hsp->start('query') - $qseq_full_start);
		$threep = substr_se($scaf_seq, $act_start, $hsp->start('hit') -1);
	
		## 5' end of proto ##
		$act_end = $hsp->end('hit') + $qseq_full_end - $hsp->end('query');
		$fivep = substr_se($scaf_seq, $hsp->end('hit') + 1, $act_end);
		}
	elsif( $hsp->strand('hit') == -1 ){	 # - strand hit: extension reversed #
		## 5' end of proto ##
		$act_start = $hsp->start('hit') - ($qseq_full_end - $hsp->end('query'));
		$fivep = substr_se($scaf_seq, $act_start, $hsp->start('hit') -1);
		$fivep = revcomp($fivep);
	
		## 5' end of proto ##
		$act_end = $hsp->end('hit') + ($hsp->start('query') - $qseq_full_start);
		$threep = substr_se($scaf_seq, $hsp->end('hit') + 1, $act_end);
		$threep = revcomp($threep);
		}
	else{ die " ERROR: logic $!\n"; }
	
	# sequence from scaffold #
	my $sseq_aln_full = join("", $threep, $$seqs_r[1]->seq, $fivep);

	# check #
	#print Dumper "XXX";
#	print Dumper $hsp->start('query');
#	print Dumper $hsp->end('query');
#	print Dumper $hsp->start('hit');
#	print Dumper $hsp->end('hit');	
#	print Dumper $act_start;
#	print Dumper $act_end;
#	print Dumper $seqs_r;
#	print Dumper $$seqs_r[1]->seq;
#	print Dumper $sseq_aln_full;
#	print Dumper join(" ", "strand:", $hsp->strand('hit'));
#	print Dumper "XXX";	
#	exit;

	return $sseq_aln_full, $act_start, $act_end;	# sseq_full, $sseq_full_start, sseq_full_end
	}

sub get_full_qseq{
# 'completing spacer' if only partial hit #
## spacer-protospacer alignment retained (ie. any gaps) ##
	my ($result, $hit, $hsp, $seqs_r) = @_;

	# setting variables #
	#$hsp->start('query');	# qstart
	#$hsp->end('query');	# qend
	#$hit->query_length;	# qlen
	#$hsp->start('hit');	# sstart
	#$hsp->end('hit');		# send
	#$hit->hit_length;		# slen
	
	my $qseq_aln =  $$seqs_r[0]->seq;
	my $qseq_full = $query_r->{$result->query_description};
	die " ERROR: no query sequence found for ", $result->query_description, "!\n"
		unless $qseq_full;
	
	# max bp that start or end can actually be extended #	
	my ($min_sstart, $max_send) = max_extend($hit, $hsp);
	
	## actual start - end (indexing by 1) ##
	### start ###
	my $act_start = 1;
	$act_start = $hsp->start('query') - $min_sstart 
		if $hsp->start('query') - $min_sstart > $act_start; 
	### end ###
	my $act_end = $hit->query_length;
	$act_end = $hsp->end('query') + $max_send 
		if $hsp->end('query') + $max_send < $act_end;
		
		#print Dumper $qseq_full;
		#print Dumper $hit->query_length;
		#print Dumper $hsp->start('query');
		#print Dumper $hsp->end('query');
		#print Dumper $max_send;
		#print Dumper $$seqs_r[0]->seq;
		#print Dumper $act_end;
	
	## extending spacer sequence ##
	my ($sfivep, $sthreep) = ("", "");
	$sfivep = substr_se($qseq_full, $act_start, $hsp->start('query') -1 );
		#print Dumper join("__", $act_start, $hsp->start('query') -1, $sfivep); 
	$sthreep = substr_se($qseq_full, $hsp->end('query') + 1, $act_end);
		#print Dumper join("__", $hsp->end('query') + 1, $act_end, $sthreep); 
	
	my $qseq_aln_full = join("", $sfivep, $$seqs_r[0]->seq, $sthreep);
		#print Dumper "$qseq_aln_full\nXXXX\n";

	#if($sfivep || $sthreep){
		#print Dumper "5p: $sfivep";
		#print Dumper "3p: $sthreep";
		#print Dumper $qseq_aln_full;
		#print Dumper $$seqs_r[0]->seq;
		#print Dumper join(" ", "strand:", $hsp->strand('hit'));
		#}
		
	# sanity check #
	die " ERROR: full length qseq length not matching start-end: $qseq_aln_full == $act_end - $act_start + 1!\n" 
		unless length $qseq_aln_full == $act_end - $act_start + 1;
	
	return $qseq_aln_full, $act_start, $act_end; 		# full length spacer (+ gaps), qseq_full_start, qseq_full_end;
	}
	
sub substr_se{
# substr using start-end indexed by 1 (blast positioning) #
## args: seq, start (index1), end (index1)
## getting from start to end 
	my ($seq, $start, $end) = @_;
	my $sub = substr($seq, $start - 1, $end - $start + 1);
	if($sub){ return $sub;}
	else{ return ""; }
	}

sub max_extend{
	my ($hit, $hsp) = @_;
	
	# determining max extend range; can't extend spacer-protospacer past scaffold #
	my $min_sstart = $hsp->start('hit') - 1;			# max bp to extend backward
	die " ERROR: max_extend_start is negative!\n" unless $min_sstart >= 0;
	
	my $slen = $hit->hit_length;
	my $send = $hsp->end('hit');	
	my $max_xend = 0;
	$max_xend = $slen - $send - 1 if $slen - $send -1 > 0;					# max bp to extend forward
	#die " ERROR: max_extend_end is negative: slen=$slen, send=$send!\n" unless $max_send >= 0;
	
		#print Dumper $hsp->start('hit'), $send, $slen, $min_sstart, $max_send;
	
	# strand of hit dictates start - end #
	if($hsp->strand('hit') == 1){
		return $min_sstart, $max_xend;
		}
	elsif($hsp->strand('hit') == -1){
		return $max_xend, $min_sstart;
		}
	else{ die " LOGIC ERROR: $!\n"; }
	}

sub filter_overlapping_hits{
# skipping hit if another hit from same group hits longer region of hit #
	my ($result, $hit, $hsp, $hit_itree_r) = @_;

	my $qseqid = $result->query_description;

		#print Dumper $hit_itree_r, $hit->name, $qseqid unless exists $hit_itree_r->{$hit->name}{$qseqid};
	die " ERROR: cannot find ", $hit->name, " in interval tree!\n"
		unless exists $hit_itree_r->{$hit->name}{$qseqid};
	
	my $res = $hit_itree_r->{$hit->name}{$qseqid}->fetch($hsp->start('hit'), $hsp->end('hit'));
		die " ERROR: no hits for found in interval tree!\n" unless @$res;
						
	my $next_bool = 0;
	foreach (@$res){		# skipping if another hit from same spacer/DR group hits larger range
		next if scalar @$res == 1;		# only 1 hit to location, don't have to check 
						
		$next_bool = 1 if 
			($$_[1] <= $hsp->start('hit') && $$_[2] >= $hsp->end('hit')) &&	# same range or greater 
			($$_[1] != $hsp->start('hit') && $$_[2] != $hsp->end('hit')) &&	# not skipping if same range
			($$_[1] < $hsp->start('hit') || $$_[2] > $hsp->end('hit'));		# skipping if another range spans larger region 
		}
	print STDERR "...skipping blast hit falling with large blast hit for same query-subject!\n" 
		if $verbose && $next_bool;

	return $next_bool;
	}

sub get_hit_itrees{
# putting hits in interval trees parsed by scaffold
# used to find hits inclusive of others
	my ($result) = @_;
	
	my %itrees;
	while( my $hit = $result->next_hit ) {		# hit object
		while( my $hsp = $hit->next_hsp ) {
			#$hit->name;					# scaffold ID
			#$hsp->start('hit');			# sstart
  			#$send = $hsp->end('hit');	# send
  			my $qseqid = $result->query_description;
  			
  			$itrees{$hit->name}{$qseqid} = Set::IntervalTree->new 
  				unless exists $itrees{$hit->name}{$qseqid};
  			
  			$itrees{$hit->name}{$qseqid}->insert(
  								[1, $hsp->start('hit'), $hsp->end('hit')],
  								$hsp->start('hit') - 1,
  								$hsp->end('hit') + 1
  								);
  			}
  		}
	  	#print Dumper %itrees; exit;
  	return \%itrees;
	}

sub get_seq_frag{
	my ($start, $end, $slen, $extend, $seq) = @_;
	
	#print Dumper $start, $end; exit;
	if($start > $end){			# need rev-comp of string 
		$start += $extend;
		$start = $slen if $start > $slen;
		$end -= $extend;
		$end = 0 if $end < 0;
		return revcomp(substr($seq, $end, $start - $end)), $start, $end;
		}
	else{
		$start -= $extend;
		$start = 0 if $start < 0;
		$end += $extend;
		$end = $slen if $end > $slen;
		return substr($seq, $start, $end - $start), $start, $end;	
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
		print STDERR "$cmd\n" unless $verbose;
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
		$sth->bind_param(2, $$row[2]);				# genbank_file
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
	
	# writing fasta #
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
		return $parts[2]; 			# just fasta name
		}
	}

sub get_loci_fasta_genbank{
# querying CLdb for fasta & genbank files for each taxon #
	my ($dbh, $join_sql) = @_;
	
	my $q = "SELECT taxon_name, taxon_id, genbank_file, fasta_file 
FROM loci 
WHERE locus_id=locus_id 
$join_sql 
GROUP BY taxon_name, taxon_id";
	
	my $res = $dbh->selectall_arrayref($q);
	die " ERROR: no matches found to query!\n"
		unless @$res;

		#print Dumper @$res; exit;
	return $res;
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

sub load_fasta{
# loading fasta file as a hash #
	my $error = "$query_file not formatted correctly! The query input file should be formatted as: >'locusID'|'spacer/DR'|'spacer/DR_ID'|'groupID'";
	
	my ($fasta_in, $query_bool) = @_;
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			# check formating #
 			my @l = split /\|/;
 			die " ERROR: $error \n" if $query_bool && $l[1] !~ /(spacer|DR)/i;
 			die " ERROR: $error \n" if $query_bool && scalar @l < 4;
 			
 			# loading # 
 			s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta

