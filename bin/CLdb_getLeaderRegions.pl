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
use Set::IntervalTree;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $overlap_check, $degeneracy_bool, $degen_check);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $length = 1000;				# max length of leader region
my $gene_len_cutoff = 150;		# bp minimum to be considered a real gene
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,	   
	   "query=s" => \$extra_query, 
	   "length=i" => \$length,
	   "overlap" => \$overlap_check,
	   "cutoff=i" => \$gene_len_cutoff,	
	   "both" => \$degen_check,			# use degeneracies to determine leader side? [FALSE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);
my $genbank_path = path_by_database($database_file);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting array start-end from loci table #
my $array_se_r = get_array_se($dbh, $join_sql, $extra_query);

# determinig leader end of array based on DR degeneracy #
## pulling out direct repeat sequences for each loci ##
my $leader_loc_r = get_DR_seq($dbh, $array_se_r);

## getting sequences from genbanks ##
get_leader_seq($array_se_r, $leader_loc_r, $genbank_path, $length);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub path_by_database{
	my ($database_file) = @_;
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "genbank");
	}

sub get_leader_seq{
# getting leader seqs from genbanks #
	my ($array_se_r, $leader_loc_r, $genbank_path, $length) = @_;

	# getting all leader regions for each genbank #
	my %genbank_loci_index;
	foreach my $locus (keys %$array_se_r){
		push(@{$genbank_loci_index{ ${$array_se_r->{$locus}}[3] }}, $locus);
		}

	# output of genbank regions #
	my $out = Bio::SeqIO->new(-format => "fasta", -fh => \*STDOUT);
	
	# opening genbanks & getting sequences #
	foreach my $genbank_file (keys %genbank_loci_index){
		
		# loadign genbank #
		my $seqio = Bio::SeqIO->new(-format => "genbank", 
								-file => "$genbank_path/$genbank_file");
		while( my $seqo = $seqio->next_seq ){
			# checking for sequence in genbank file; if not provided, die #
			unless($seqo->seq){ 
				print STDERR " ERROR: no genomic sequence found in $genbank_file! Cannot extract the leader sequence from the genome!\n";
				next;
				}
			foreach my $locus (@{$genbank_loci_index{ $genbank_file }} ){
				foreach my $loc (@{$leader_loc_r->{$locus}}){		# if both ends, must get both
					
					# only comparing on the same scaffold (or just 1 scaffold)
					my $leader_scaf = ${$array_se_r->{$locus}}[4];
					next unless $leader_scaf eq "CLDB__ONE_CHROMOSOME" ||
						$leader_scaf eq $seqo->display_id;
					
					# determining leader region start-end
					my $leader_start;
					my $leader_end;
					my $array_start = ${$array_se_r->{$locus}}[1];
					my $array_end = ${$array_se_r->{$locus}}[2];
					
					if($loc eq "start" ){
						$leader_start = $array_start - $length;
						$leader_end = $array_start; 
						}
					elsif($loc eq "end"){		# if "end"
						$leader_start = $array_end;
						$leader_end = $array_end + $length;
						}
					else{ die " LOGIC ERROR: $!\n"; }
					
					# no negative (or zero) values (due to $length extending beyond scaffold) #
					$leader_start = 1 if $leader_start < 1;
					$leader_end = 1 if $leader_end < 1;
					
					# leader start & end must be < total scaffold length #
					my $scaf_len = length $seqo->seq;
					$leader_start = $scaf_len if $leader_start > $scaf_len;
					$leader_end = $scaf_len if $leader_end > $scaf_len;
					
					# checking for gene overlap #
					($seqo, $leader_start, $leader_end) = check_gene_overlap($seqo, $leader_start, $leader_end, $loc, $locus)
						unless $overlap_check;
					
					# parsing seq from genbank #
					my $leader_seqo = $seqo->trunc($leader_start, $leader_end);					
					$leader_seqo->display_id("cli.$locus\___$loc\___$leader_scaf\___$leader_start\___$leader_end");
					$leader_seqo->desc(" $loc\___$leader_scaf\___$leader_start\___$leader_end");
					$leader_seqo->desc("");
									
					$out->write_seq($leader_seqo);
					}
				}
			}
		}
	}

sub check_gene_overlap{
# checking to see if and genes overlap w/ leader region #
# if yes, truncating leader region #
	my ($seq, $region_start, $region_end, $leader_loc, $locus) = @_;
	
	# making an interval tree #
	my $itree = Set::IntervalTree->new();
	for my $feat ($seq->get_SeqFeatures){
		next if $feat->primary_tag eq "source";

		my $start = $feat->location->start;
		my $end = $feat->location->end;
	
		# filtering out short genes #
		my @tags = grep(/translation/, $feat->get_all_tags);
		my $gene_length = 0;
		if(@tags){
			my @tmp = $feat->get_tag_values("translation"); 	# amino acid 
			$gene_length = (length $tmp[0]) * 3;
			}
		else{ 
			$gene_length = length $feat->entire_seq->seq;		# nucleotide
			}
		next if $gene_length < $gene_len_cutoff;			# excluding short genes (probably not real) 		
		
		# loading tree #
		if($start <= $end){			# all on + strand
			$itree->insert($feat, $start - 1,  $end + 1);
			}
		else{
			$itree->insert($feat, $end - 1,  $start + 1);		
			}
		}
		
	# debuggin itree #	
		#print_itree($itree, $region_start, $region_end);
		
	my $overlap_pos = 0;
	if($leader_loc eq "start"){		# checking for 1st overlap from 3' - 5' end # 
		for (my $i=$region_end; $i>=$region_start; $i--){
			$overlap_pos = check_pos($itree, $i);
			last if $overlap_pos;
			}
		}
	elsif($leader_loc eq "end"){	# checking for 1st overlap from 5' - 3' end;
		for (my $i=$region_start; $i<=$region_end; $i++){
			$overlap_pos = check_pos($itree, $i);
			last if $overlap_pos;
			}
		}
	
	sub check_pos{
		my ($itree, $i) = @_;
		my $res = $itree->fetch($i, $i);
		return $i if $$res[0];
		}
		
	# altering start-end #
	if($overlap_pos){
		print STDERR "WARNING: Gene(s) overlap the leader in cli.$locus array. Truncating\n" unless $verbose;
		
		if($leader_loc eq "start"){
			$region_start = $overlap_pos;
			}
		elsif($leader_loc eq "end"){
			$region_end = $overlap_pos;
			}
		die " LOGIC ERROR: region start is > region end! $!\n"
			if $region_start > $region_end;
		}

	# sanity checks #
	die " ERROR: For locus$locus, start and/or end is negative!: start=$region_start; end=$region_end\n"
		unless $region_start >= 0 && $region_end >= 0;

	return $seq, $region_start, $region_end;
	}

sub check_gene_overlap_old{
# checking to see if and genes overlap w/ leader region #
# if yes, truncating leader region #
	my ($seq, $region_start, $region_end, $leader_loc, $locus) = @_;
	
	# making an interval tree #
	my $itree = Set::IntervalTree->new();
	for my $feat ($seq->get_SeqFeatures){
		next if $feat->primary_tag eq "source";

		my $start = $feat->location->start;
		my $end = $feat->location->end;
		if($start <= $end){
			$itree->insert($feat, $start - 1,  $end + 1);
			}
		else{
			$itree->insert($feat, $end - 1,  $start + 1);		
			}
		}
		
	# debuggin itree #	
	#	print_itree($itree, 500);
	
	# checking for overlap #
	my $trans = 0; my $last;
	for my $i ($region_start..$region_end){
		my $res = $itree->fetch($i, $i);
		if(! $last){
			if($$res[0]){ $last = 2; }			# gene-leader overlap at start
			else{ $last = 1; } 					# gene-leader overlap at end 
			}
		elsif($$res[0] && $last == 1){		# hitting overlap after no previous overlap (leader at end)
			$trans = $i;
			die " LOGIC ERROR: gene overlap at the wrong end of the leader region for cli.$locus!\n"
				unless $leader_loc eq "end"; 
			last;
			}
		elsif( ! $$res[0] && $last ==2 ){	# no more overlap after previous overlap (leader at start)
			$trans = $i - 1;
			die " LOGIC ERROR: gene overlap at the wrong end of the leader region for cli.$locus!\n"
				unless $leader_loc eq "start"; 
			last;
			}
		elsif($$res[0]){
			$last = 2;
			}
		elsif(! $$res[0]){
			$last = 1;
			}
		else{ die $!; }
		}	
		
	# altering start / stop #
	$trans -= $region_start if $trans;			# making index local to region (instead of genome)
	if($trans){
		print STDERR "WARNING: Gene(s) overlap the leader in cli.$locus array. Truncating\n" unless $verbose;
		
		if($leader_loc eq "start"){				# gene overlaps at 5' end
			$region_start += $trans;			
			}
		elsif($leader_loc eq "end"){			# gene overlaps at the 3' end
			$region_end -= $trans;
			}
		}
	return $seq, $region_start, $region_end;
	}
	
sub print_itree{
	my ($itree, $start, $end) = @_;
	
	for my $i ($start..$end){
		my $res = $itree->fetch($i, $i);
		
		print join("\t", "pos:$i", join("::", @$res)), "\n";
		}
	#exit;
	}	

sub get_DR_seq{
	my ($dbh, $array_se_r) = @_;

	my $cmd = "SELECT DR_id, DR_start, DR_end, DR_group from DRs where Locus_ID = ?";
	my $sql = $dbh->prepare($cmd);
	
	my %leader_loc;
	foreach my $locus (keys %$array_se_r){		# each CRISPR locus
		$sql->execute($locus);
		my $ret = $sql->fetchall_arrayref();
		die " ERROR: no direct repeats found for cli.$locus!\n"
			unless $$ret[0];
		die " ERROR: no 'repeat_group' entries found!\n Run CLdb_groupArrayElements.pl before this script!\n\n"
			unless $$ret[0][3];
		
		# determining leader based on degeneracies #
		if($degen_check){			# skipping if opted 
			$leader_loc{$locus} = ["start", "end"];
			}
		else{	
			$leader_loc{$locus} = determine_leader($ret, $array_se_r->{$locus}, $locus);
			}
		}
	
		#print Dumper %leader_loc; exit;
	return \%leader_loc;
	}

sub determine_leader{
# determing leader region #
	my ($ret, $array_r, $locus) = @_;
	
	# getting halfway point in array #
	my $array_half = ($$array_r[1] +  $$array_r[2]) / 2;
	
	# counting groups for each half of the array #
	my %group_cnt;
	foreach my $DR (@$ret){		# each direct repeat
		# $DR = (repeat, start, end, group_ID)
		if($$DR[1] < $array_half && $$DR[2] < $array_half){
			$group_cnt{1}{$$DR[3]} = 1;		
			}
		elsif($$DR[1] > $array_half && $$DR[2] > $array_half){
			$group_cnt{2}{$$DR[3]} = 1;
			}
		}
	my $first_cnt = scalar keys %{$group_cnt{1}};		# number of groups in 1st half
	my $second_cnt = scalar keys %{$group_cnt{2}};		# number of groups in 2nd half

	# writing degeneracy report to STDERR #
	print STDERR join("\t", "degeneracies", "cli.$locus", $first_cnt, $second_cnt), "\n"
			unless $verbose;
	
	# getting leader end based on number of degeneracies #
	## leader end should have less than trailer ##
	## therefore, the side w/ the most groups is the trailer end ##
 	## based on + strand positioning ##
	if($first_cnt < $second_cnt){		# > degeneracy at end, leader at start
		return ["start"];		
		}
	elsif($first_cnt > $second_cnt){	# > degeneracy at start, leader at end
		return ["end"];		
		}
	else{
		print STDERR " WARNING: for cli.$locus, could not use direct repeat degeneracy determine leader region end! Writing both.\n";
		return ["start", "end"];
		}
	}

sub get_array_se{
# getting the array start-end from loci table #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	my $cmd = "SELECT Locus_ID, array_start, array_end, genbank_file, scaffold FROM loci WHERE (array_start is not null or array_end is not null) $join_sql $extra_query";
	my $ret = $dbh->selectall_arrayref($cmd);

	my %array_se;
	foreach my $row (@$ret){
		# moving array start-end to + strand #
		($$row[1], $$row[2]) = pos_strand_se($$row[1], $$row[2]);
		$array_se{$$row[0]} = $row;
		}
		
		#print Dumper %array_se; exit;
	return \%array_se; 
	}
	
sub pos_strand_se{
	my ($start, $end) = @_;
	if($start > $end){ return $end, $start; }
	else{ return $start, $end; }
	}

sub make_leader_dir{
# making a directory for leader region output #
	my ($leader_path, $database_file) = @_;
	
	unless ($leader_path){
		my @parts = File::Spec->splitpath($database_file);
		$leader_path = join("/", $parts[1], "leader");
		}
	mkdir $leader_path unless -d $leader_path;
	return $leader_path;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND $cat IN (", join(", ", @$vals_r), ")");
	}

__END__

=pod

=head1 NAME

CLdb_getLeaderRegion.pl -- getting leader region sequences

=head1 SYNOPSIS

CLdb_getLeaderRegion.pl [flags] > possible_leaders.fna

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 options

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -length  <int>

Max length of the potential leader region (bp). [1000]

=item -overlap  <bool>

Check for overlapping genes (& trim region if overlap)? [TRUE]

=item -cutoff  <int>

Gene length cutoff for counting gene as real in overlap assessment (bp). [150]

=item -both  <bool>

Get sequence on both sides of the array regardless of degeneracies? [FALSE]

=item -verbose  <bool>

Verbose output

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_getLeaderRegions.pl

=head1 DESCRIPTION

Get the suspected leader regions for CRISPR loci 
in the CRISPR database.

Leader regions are determined by direct repeat (DR) degeneracy,
which should occur at the trailer end.
For example: leader half has 1 cluster due to DR conservation,
while the trailer half has 3 clusters due to DR degeneracy). 
This is determined by the 'repeat_group' column in the 
Directrepeats table in the CRISPR database. Therefore, 
CLdb_groupArrayElements.pl has to be run on the database 
before running this script!

Both end regions of a CRISPR array are written if a CRISPR
array has equal numbers of DR clusters on each half.

The default leader region length is 1000bp from the CRISPR array.

The leader region length will be truncated if any genes overlap in 
that region. 

Leader regions in the output fasta are named as: 
"cli.ID" "start|end"__"region_start"__"region_end"

Leader region start & end are according to the + strand. 

The number of repeat_groups on each side (5' & 3') of the
CRISPR array will be printed to STDERR (unless '-v'). 
The output values are: 'degeracies', locus_id', 
'5-prime_number_repeat_groups', '3-prime_number_repeat_groups'

=head1 EXAMPLES

=head2 Leader regions for all loci

CLdb_getLeaderRegion.pl -d CLdb.sqlite > leaders.fna

=head2 Leader regions for just subtype I-B

CLdb_getLeaderRegions.pl -d CLdb.sqlite -sub I-B > leader_I-B.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut
