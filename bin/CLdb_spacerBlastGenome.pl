#!/usr/bin/env perl

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

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $query_file, $entire, $filter_overlap_hits);
my (@subtype, @taxon_id, @taxon_name);
my @subject_in;
my $extra_query = "";
my $blast_params = "-evalue 0.001";		# 1e-3
my $num_threads = 1;
my $extend = 10;
GetOptions(
	   "database=s" => \$database_file,
	   "fasta=s" => \$query_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "blast=s" => \$blast_params,
	   "num_threads=i" => \$num_threads,
	   "xtend=i" => \$extend,
	   "entire" => \$entire,					# protospacer entire length of spacer regardless of blast hit range? [TRUE]
	   "overlap" => \$filter_overlap_hits, 		# filter out overlapping hits? [TRUE]
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
my $loci_update = genbank2fasta($subject_loci_r, $fasta_dir);
update_loci($dbh, $subject_loci_r) if $loci_update;

## making blast DBs from subject fasta files (foreach subject) ##
make_blast_db($subject_loci_r, $fasta_dir, $blast_dir);

# blasting (foreach subject blast DB) #
## blastn & loading blastn output ##
#blastn_call_load($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file);
blastn_xml_call_load($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file);

print Dumper "NO commit!\n";


### Subroutines
sub blastn_xml_call_load{
# calling blastn (xml output format); parsing output #
## needed to get blast alignment (find potential gaps) 
	my ($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file) = @_;
	
	# columns to load into blast_hits DB
	my @blast_hits_col = qw/blast_id spacer_DR S_taxon_name
S_taxon_ID Group_ID sseqid pident 
len mismatch gapopen qstart qend 
sstart send evalue bitscore qlen 
slen qseq frag xstart xend/;
	
	# blasting each subject genome #
	my %insert_cnt; 		# summing number of entries added/updated in CLdb
	foreach my $row (@$subject_loci_r){		# taxon_name, taxon_id, genbank, fasta
		next unless $$row[3];			# skipping if no fasta & thus no blast db
		
		# loading fasta #
		my $fasta_r = load_fasta("$blast_dir/$$row[3]");

		# setting blast cmd #
		#my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq";
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
					
					# filtering overlapping hits #
					## hit must be in $res & have > range (for totally overlapping hits)
					unless($filter_overlap_hits){
						die " ERROR: cannot find ", $hit->name, " in interval tree!\n"
							unless exists $hit_itree_r->{$hit->name}{$qseqid};
						my $res = $hit_itree_r->{$hit->name}{$qseqid}->fetch($hsp->start('hit'), $hsp->end('hit'));
						die " ERROR: no hits for found in interval tree!\n" unless @$res;
						
						my $next_bool;
						foreach (@$res){		# skipping if another hit from same spacer/DR group hits larger range
							next if scalar @$res == 1;
						
							$next_bool = 1 if 
										($$_[1] <= $hsp->start('hit') && $$_[2] >= $hsp->end('hit')) &&	# same range or greater 
										($$_[1] != $hsp->start('hit') && $$_[2] != $hsp->end('hit')) &&	# not skipping if same range
										($$_[1] < $hsp->start('hit') || $$_[2] > $hsp->end('hit'));		# skipping if another range spans larger region 
							}
						next if $next_bool;
						}
	
					# spacer or DR? #
					(my $spacer_DR = $qseqid) =~ s/_.+//;
					$qseqid =~ s/.+\.//;

					# blastID # spacer_DR qseqid S_taxon_name S_taxon_ID sseqid sstart send #
					my $blast_id = join("_", $spacer_DR,	# query ID (spacer/DR group)(
									$$row[0],				# S_taxon_name
									$$row[1],				# S_taxon_ID
									$hit->name, 			# sseqid
									$hsp->start('hit'),		# sstart
									$hsp->end('hit')		# send
									);

					# getting blast info #
					#print Dumper $result->query_description; 	# qseqid
					#print Dumper $hit->name;					# sseqid
					#print Dumper $$row[0]; 					# S_taxon_name
					#print Dumper $$row[1]; 					# S_taxon_id
					#print Dumper $hit->query_length;			# qlen
					#print Dumper $hit->hit_length;				# slen
    				#print Dumper $hsp->percent_identity;		# pident
    				#print Dumper $hsp->bits;					# bitscore
    				#print Dumper $hsp->evalue;					# evalue
    				#print Dumper $hsp->start('query');			# qstart
    				#print Dumper $hsp->end('query');			# qend
    				#print Dumper $hsp->start('hit');			# sstart
    				#print Dumper $hsp->end('hit');				# send
    				#print Dumper $hsp->strand; 				# strand
    				
					# getting alignment #
					my $alnIO = $hsp->get_aln; 
					my @seqs = $alnIO->each_seq;

					# getting protospacer sequence (& start-end) #
					die " ERROR: cannot find ", $hit->name, " in $$row[0], $$row[1]\n"
						unless exists  $fasta_r->{$hit->name};
					my ($proto, $xstart, $xend);
					if($spacer_DR eq "DR"){
						($proto, $xstart, $xend) = ("", $hsp->start('hit'), $hsp->end('hit'));
						}
					elsif($spacer_DR eq "Spacer"){
						($proto, $xstart, $xend) = get_protospacer($hsp, $hit, \@seqs, $fasta_r->{$hit->name});
						}
					else{ die " LOGIC ERROR: $!\n"; }
						#print Dumper $seqs[0]->seq;
						#print Dumper $proto;
						#print Dumper $xstart;
						#print Dumper $xend;
						#print Dumper $hsp->strand('hit');
						#print Dumper $hit->name;

					my @hit_matches = $hsp->matches('hit');

					## making blast_hit sql ##
					my @vals = (
						$dbh->quote($blast_id),		# blastID
						$dbh->quote($spacer_DR),	# Spacer|DR
						$dbh->quote($$row[0]),		# taxon_name
						$dbh->quote($$row[1]),		# taxon_id
						$qseqid,					# groupID
						$dbh->quote($hit->name), 	# sseqid
						$hsp->percent_identity,
						$hit->hit_length,
						$hsp->length('total') - $hit_matches[0] - $hsp->gaps, 		# mismatches (length - matches - gaps)
						$hsp->gaps,
						$hsp->start('query'),
						$hsp->end('query'),
						$hsp->start('hit'),
						$hsp->end('hit'),
						$hsp->evalue,
						$hsp->bits,
						$hit->query_length,
						$hit->hit_length,
						$dbh->quote($seqs[0]->seq), 				# query sequence
						$dbh->quote($proto),
						$xstart,
						$xend
						);
					
					#print join("\t", @vals), "\n";
					# making sql #
					my $sql = join(" ", "INSERT INTO Blast_hits( ",
						join(",", @blast_hits_col),
						") VALUES( ",
						join(",", @vals), 
						")" );
						
					# loading blast hit into CLdb #
					$dbh->do( $sql );
					if($DBI::err){
						print STDERR "ERROR: $DBI::errstr in: ", join("\t", @vals), "\n";
						}
					else{ $insert_cnt{$spacer_DR}++; }
					}
				}
			}
		}
	
	$insert_cnt{"Spacer"} = 0 unless exists $insert_cnt{"Spacer"};
	print STDERR "...Number of spacer blast entries added to CLdb: $insert_cnt{'Spacer'}\n";
	$insert_cnt{"DR"} = 0 unless exists $insert_cnt{"DR"};
	print STDERR "...Number of DR blast entries added to CLdb: $insert_cnt{'DR'}\n";
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

sub get_hit_se{
# getting all start-end for each hit to find overlapping hits #
	my $result = shift;
	my %hit_se;

	while( my $hit = $result->next_hit ) {		# hit object
		while( my $hsp = $hit->next_hsp ) {
			# spacer/DR group => subject_scaffold => (start,end)
			my $qseqid = $result->query_description;
			my $sstart = $hsp->start('hit');			# sstart
  				my $send = $hsp->end('hit');				# send
				if( exists $hit_se{$qseqid}{$hit->name} ){		# hit exists, check for inclusion
					$hit_se{$qseqid}{$hit->name} = [$sstart, $send]
						if $sstart <= ${$hit_se{$qseqid}{$hit->name}}[0] &&		# range must span previous range and be = or larger
						$send >= ${$hit_se{$qseqid}{$hit->name}}[1];
				}
			else{			# no hit, loading into hash
				$hit_se{$qseqid}{$hit->name} = [$sstart, $send];		
				}
			}
		}
		
		#print Dumper %hit_se; 
	return \%hit_se;
	}

sub get_protospacer{
# getting the protospacer sequence from genome genome fasta 
# also extending protospacer to full lenght of query (unless -entire)
# also adding adjacent sequence for finding PAMs
	my ($hsp, $hit, $seqs_r, $scaf_seq) = @_;

	# setting variables #
	my $qstart = $hsp->start('query');	# qstart
	my $qend = $hsp->end('query');		# qend
	my $qlen = $hit->query_length;		# qlen
	my $sstart = $hsp->start('hit');	# sstart
	my $send = $hsp->end('hit');		# send
	my $slen = $hit->hit_length;		# slen
	
	# full-length protospacer #
	unless($entire){	 				# extending protospacer to full length of spacer
			#print Dumper "entire", $qstart, $qend, $qlen if $qstart != 1 || $qlen != $qend;
			$sstart -= $qstart - 1;		# missing 5' end hit
			$send += $qlen - $qend;		# missing 3' end hit
			}
		
	# extend #
	$sstart -= $extend;					# 5' extension beyond proto
	$sstart = 0 if $sstart < 0;			# floor of 0
	$send += $extend;					# 3' extension beyond proto
	$send = $slen if $send > $slen; 	# ceiling of scaffold length

	
	# making protospacer seq 
	## substr from scaff seq for 5' & 3' extensions 
	## revcomp if - strand; proto revcomp already  by bioperl
	my $fivep = substr($scaf_seq, $sstart -1, $hsp->start('hit') - $sstart);	# 5' extension; amost forgot 0-indexing!
	$fivep = revcomp($fivep) if $hsp->strand('hit') ==  -1;
	my $threep = substr($scaf_seq, $hsp->end('hit') ,  $send - $hsp->end('hit')); # 3' extension
	$threep = revcomp($threep) if $hsp->strand('hit') == -1;	

	# returning sequence; 
	## 5' & 3' switch if - strand ##
	if($hsp->strand('hit') == 1){
		return join("", $fivep, $$seqs_r[1]->seq, $threep), 
			$sstart - 1, $send;		# -1 for 0-indexing
		}
	elsif($hsp->strand('hit') == -1){			 # flip 5' <-> 3'
		return join("", $threep, $$seqs_r[1]->seq, $fivep), 
			$sstart - 1, $send;		# -1 for 0-indexing
		}
	else{die " ERROR: cannot recognize strand: ", $hsp->strand('hit'), "\n"; }
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

sub blastn_call_load{
# calling blastn, loading blast results directly into DB #
# also loading blast subject hit #
	my ($dbh, $subject_loci_r, $blast_dir, $blast_params, $num_threads, $query_file) = @_;
	
	# status #
	print STDERR "...BLASTing each subject genome database\n";
	
	# preparing blast_hits insert #
	my @blast_hits_col = qw/blast_id spacer_DR S_taxon_name S_taxon_ID Group_ID sseqid pident len mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq frag xstart xend/;
		
	# blasting each subject genome #
	my %insert_cnt;			# summing number of entries for all subjects 
	foreach my $row (@$subject_loci_r){
		next unless $$row[3];			# skipping if no fasta & thus no blast db
		
		# loading fasta #
		my $fasta_r = load_fasta("$blast_dir/$$row[3]");
		
		# setting blast cmd #
		my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qseq";
		my $cmd = "blastn -task 'blastn-short' -outfmt '$outfmt' -db $blast_dir/$$row[3] -query $query_file -num_threads $num_threads $blast_params |";
		print STDERR "$cmd\n" unless $verbose;

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
			$line[0] = $dbh->quote($line[0]);		# blast_id
			$line[1] = $dbh->quote($line[1]);		# spacer_DR (qseqid)
			$line[2] = $dbh->quote($line[2]);		# S_taxon_name
			$line[3] = $dbh->quote($line[3]);		# S_taxon_ID
			my $scaf = $line[5];
			$line[5] = $dbh->quote($line[5]);		# sseqid (scaffold)
			
			# qseq; not need if direct repeat #
			if($spacer_DR eq "DR"){ $line[18] = $dbh->quote(""); }		
			else{ $line[18] = $dbh->quote($line[18]); }	
			
			# frag & extension #
			## substr of fasta ##
			my ($frag, $xstart, $xend);
			if($line[1] =~ /Spacer/){				# if spacer blast hit
				die " ERROR: '$scaf' not found in fasta!" unless exists $fasta_r->{"$scaf"};
				my $start = $line[12];
				my $end = $line[13];
				my $slen = $line[17];
				($frag, $xstart, $xend) = get_seq_frag($start, $end, $slen, $extend, $fasta_r->{"$scaf"});
				$frag = $dbh->quote($frag);			
				}
			else{ 
				$frag = $dbh->quote(""); 
				$xstart = $line[12];
				$xend = $line[13];
				}
			push @line, $frag, $xstart, $xend;
			
			# blast_hits #
			## making blast_hit sql ##
			my $sql = join(" ", "INSERT INTO Blast_hits( ",
					join(",", @blast_hits_col),
					") VALUES( ",
					join(",", @line), 
					")" );
			
			## loading line ##
			$dbh->do( $sql );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @line), "\n";
				}
			else{ $insert_cnt{$spacer_DR}++; }
			}
		close PIPE;
		}
	$dbh->commit();
	
	$insert_cnt{"Spacer"} = 0 unless exists $insert_cnt{"Spacer"};
	print STDERR "...Number of spacer blast entries added to CLdb: $insert_cnt{'Spacer'}\n";
	$insert_cnt{"DR"} = 0 unless exists $insert_cnt{"DR"};
	print STDERR "...Number of DR blast entries added to CLdb: $insert_cnt{'DR'}\n";
	}

sub revcomp{
	# reverse complements DNA #
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/ACGTNBVDHKMRYSW\.-/TGCANVBHDMKYRSW\.-/;
	return $seq;
	}

sub load_fasta{
# loading fasta file as a hash #
	my $fasta_in = shift;
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if  /^\s*$/;	
 		
 		if(/^>.+/){
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

BLASTn parameters (besides required flags). [-evalue 0.001]

=item -num_threads  <int>

Use for '-num_threads' parameter in BLASTn. [1]

=item -entire  <bool>

Protospacer sequence as the entire length of the spacer, regardless of blast hit length? [TRUE]

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

