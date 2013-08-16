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

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# listing tables #
my $tables_r = list_tables($dbh);

# loading genbanks #
my $gen_list_r = get_genbank_list();
my $gen_io_r = load_genbank_io($gen_list_r);

# making merged => unmerged location index #
foreach my $unmerged (keys %$gen_list_r){				
	# getting scaffold lengths & sequence #
	my $unmerged_r = parse_unmerged_genbank($gen_list_r, $gen_io_r, $unmerged);
	my $merged_r = parse_merged_genbank($gen_list_r, $gen_io_r, $gen_list_r->{$unmerged});

	# finding merged-unmerged contig 1-to-1 assocations in contig length #
	my ($scaf_index_r, $needs_blasting_r) = find_same_length($gen_list_r, $unmerged_r, $merged_r);
	
	# determine if sequences are the same (either identical or rev-comp)
	find_same_sequence($scaf_index_r, $needs_blasting_r, $unmerged_r, $merged_r);

	# checking to make sure all scaffolds accounted for #
	if(%{$needs_blasting_r->{"merged"}}){
		print STDERR " ERROR: not all 1-to-1 connects made. Ambiguous contigs:\n\t";
		print STDERR join(",\n\t", keys %{$needs_blasting_r->{"unmerged"}}), "\n\n";
		exit;
		}
	
	# making an interval tree relating absolute location to scaffold-based location #
	my $itree = make_loc_itree($scaf_index_r, $merged_r, $unmerged_r);
	
	# getting values from database #
	my $loci_r = get_loci_oi($dbh, $gen_list_r->{$unmerged});
	
	# getting associated tables #
	my %tables_oi;
	get_other_tables_oi($dbh, $loci_r, \%tables_oi, "spacers", "spacer")
		if exists $tables_r->{"spacers"};
	get_other_tables_oi($dbh, $loci_r, \%tables_oi, "directrepeats", "repeat")
		if exists $tables_r->{"directrepeats"};
	get_other_tables_oi($dbh, $loci_r, \%tables_oi, "leaderseqs", "leader", "locus")
		if exists $tables_r->{"leaderseqs"};
	get_other_tables_oi($dbh, $loci_r, \%tables_oi, "genes", "gene")
		if exists $tables_r->{"genes"};
	
	# converting values #
	merged_to_unmerged_pos($itree, $loci_r, \%tables_oi);
	
	# updating tables #
	update_loci($dbh, $loci_r);
	update_other_table($dbh, $loci_r, \%tables_oi, "spacers", "spacer")
		if exists $tables_r->{"spacers"};
	update_other_table($dbh, $loci_r, \%tables_oi, "directrepeats", "repeat")
		if exists $tables_r->{"directrepeats"};
	update_other_table($dbh, $loci_r, \%tables_oi, "leaderseqs", "leader", "locus")
		if exists $tables_r->{"leaderseqs"};
	update_other_table($dbh, $loci_r, \%tables_oi, "genes", "gene")
		if exists $tables_r->{"genes"};
	}

# updating loci table values: genbank file names #
update_loci_genbank($dbh, $gen_list_r);

# committing updates #
$dbh->commit;
print STDERR "...All updates committed!\n";

# disconnect #
$dbh->disconnect();
exit;



### Subroutines
sub update_loci_genbank{
# updating loci genbank values to unmerged genbank files #
	my ($dbh, $gen_list_r) = @_;
	print STDERR "...updating genbank file names in loci table\n";
	
	my $cmd = "UPDATE loci SET genbank=? where genbank=?";		# merged to unmerged
	my $sql = $dbh->prepare($cmd);

	my $cnt = 0;
	foreach my $unmerged (keys %$gen_list_r){
		my @parts = File::Spec->splitpath($unmerged);
		$sql->bind_param(1, $parts[2]);
		@parts = File::Spec->splitpath($gen_list_r->{$unmerged});
		$sql->bind_param(2, $parts[2]);		
		$sql->execute();
		
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: $unmerged\n";
			}
		else{ $cnt++; }
		}

	print STDERR "...Number of genbank values updated: $cnt\n";
	print STDERR "...BE SURE to move unmerged genbanks into \$CLdb_HOME/genbank!!!\n";
	}

sub find_same_sequence{
# checking contigs by same sequence; should match all contigs #
## sequences should match exactly & 1-to-1 ##
	my ($scaf_index_r, $needs_blasting_r, $unmerged_r, $merged_r) = @_;

	foreach my $merged_contig (keys %{$needs_blasting_r->{"merged"}}){
		foreach my $unmerged_contig (keys %{$needs_blasting_r->{"unmerged"}}){
			if ($merged_r->{$merged_contig}{"seq"} eq $unmerged_r->{$unmerged_contig}{"seq"}){
				push( @{$scaf_index_r->{$merged_contig}}, $unmerged_contig);
				#print Dumper $unmerged_contig;
				#print Dumper $merged_contig;
				#print Dumper $unmerged_contig;
				#print Dumper $merged_r->{$merged_contig}{"seq"};
				#print Dumper $unmerged_r->{$unmerged_contig}{"seq"};
				
				}
		#	else{
				#print Dumper $merged_r->{$merged_contig}{"seq"};
		#		#print Dumper $unmerged_r->{$unmerged_contig}{"seq"}
		#		}
			}
		}
	
	# checking problem contigs for 1 to 1 associations #
	foreach my $merged_contig (keys %{$needs_blasting_r->{"merged"}}){
		if( exists $scaf_index_r->{$merged_contig} ){
			if(scalar @{$scaf_index_r->{$merged_contig}} == 1){
				# deleting merged #
				delete $needs_blasting_r->{"merged"}{$merged_contig};
			
				# deleting unmerged #
				foreach my $unmerged_contig (@{$scaf_index_r->{$merged_contig}}){
					delete $needs_blasting_r->{"unmerged"}{$unmerged_contig};
					}
				}
			else{
				die " ERROR: Scaffold->\"$merged_contig\" has same sequence as >=2 scaffolds in unmerged genbank!\n";
				}
			}
		else{		# no sequenc hit
			die " ERROR: Scaffold->\"$merged_contig\" has no identical sequence in unmerged genbank!\n";
			}
		}
	}

sub list_tables{
	my $dbh = shift;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	map{delete $all->{$_}; tr/A-Z/a-z/; $all->{$_} = 1} keys %$all;
		#print Dumper %$all; exit;
	return $all;
	}
	
sub update_other_table{
# updating position info in other tables #
	my ($dbh, $loci_r, $tables_oi_r, $table, $prefix, $prefix2) = @_;
	$prefix2 = $prefix unless $prefix2;
	
	my $cmd = "UPDATE $table SET $prefix\_start=?, $prefix\_end=? 
WHERE $prefix2\_id=?
AND locus_id=?";

	$cmd =~ s/[\n\r]/ /g;
	
	my $sql = $dbh->prepare($cmd);
	my $cnt = 0;
	foreach my $locus (@$loci_r){
		foreach my $row ( @{$tables_oi_r->{$table}{$$locus[0]}} ) {
			$sql->execute(@$row[0..2], $$locus[0]);
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr for locus: $$locus[0]\n";
				}
			else{ $cnt++; }
			}
		}
	
	print STDERR "...Number of entries to be updated in $table: $cnt\n";	
	}

sub update_loci{
# updating position info in loci table #
	my ($dbh, $loci_r) = @_;
	
	my $cmd = "UPDATE Loci SET scaffold=?, locus_start=?, locus_end=?,
operon_start=?, operon_end=?, CRISPR_array_start=?, CRISPR_array_end=?
WHERE locus_id=?";
	$cmd =~ s/[\n\r]/ /g;
	
	my $sql = $dbh->prepare($cmd);
	my $cnt = 0;
	foreach my $locus (@$loci_r){
		$sql->execute(@$locus[1..$#$locus], $$locus[0]);
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for locus: $$locus[0]\n";
			}
		else{ $cnt++; }
		}
	
	print STDERR "...Number of entries to be updated in Loci: $cnt\n";
	}

sub merged_to_unmerged_pos{
# getting scaffold that positions fall into #
	my ($itree, $loci_r, $tables_oi_r) = @_;
	
	# updating per-loci #
	foreach my $locus (@$loci_r){			# $$locus[0] = locus_id
		my $scaffold_name;
		# locus table #
		for (my $i=2; $i<=6; $i+=2){
			# skipping if no values #
			next unless $$locus[$i] && $$locus[$i+1];
			
			# checking start-end orientation #
			my $flip_bool = 0;
			$flip_bool = 1 if $$locus[$i] > $$locus[$i+1];
			($$locus[$i], $$locus[$i+1]) = flip_se($$locus[$i], $$locus[$i+1]) if $flip_bool;
			
			# getting position info #
			my $ret = $itree->fetch($$locus[$i], $$locus[$i+1]);		# ret = scaffold that positions span
			check_scaffold_span($ret, $locus, $$locus[$i], $$locus[$i+1], "loci");						# positions should span 1 scaffold

			# updating position info #
			## new start = scaf_start + (current_pos - scaf_gw_pos)
			## new end = scaf_end - (scaf_gw_pos - current_pos)
			$$locus[$i] = $$ret[0]->{"start"} + ($$locus[$i] - $$ret[0]->{"gw_start"});
			$$locus[$i+1] = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $$locus[$i+1]);
			
			# flipping start-end back if necessary #
			($$locus[$i], $$locus[$i+1]) = flip_se($$locus[$i], $$locus[$i+1]) if $flip_bool;
			
			# check start stop #
			die " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->loci; start = $$locus[$i], end = ", $$locus[$i+1], "\n"
				#print STDERR " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->loci; start = $$locus[$i], end = ", $$locus[$i+1], "\n"
				if $$locus[$i] < 0 || $$locus[$i+1] < 0;
			
			# scaffold #
			$scaffold_name = $$ret[0]->{"id"} unless $scaffold_name;
			}
		$$locus[1] = $scaffold_name;
		
		# other tables #
		foreach my $table (keys %$tables_oi_r){
			foreach my $row ( @{$tables_oi_r->{$table}{$$locus[0]}} ){		 # 2 value
				# skipping if not present #
				next unless $$row[0] && $$row[1];

					#print Dumper @$row;
				
				# checking start-end orientation #
				my $flip_bool = 0;
				$flip_bool = 1 if $$row[0] > $$row[1];
				($$row[0], $$row[1]) = flip_se($$row[0], $$row[1]) if $flip_bool;
				
				# getting position info #
				my $ret = $itree->fetch($$row[0], $$row[1]);		# ret = scaffold that positions span
				check_scaffold_span($ret, $locus, $$row[0], $$row[1], $table);						# positions should span 1 scaffold

				# updating position info #
				## new start = scaf_start + (current_pos - scaf_gw_pos)
				## new end = scaf_end - (scaf_gw_pos - current_pos)
				$$row[0] = $$ret[0]->{"start"} + ($$row[0] - $$ret[0]->{"gw_start"});
				$$row[1] = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $$row[1]);
			
				# flipping start-end back if necessary #
				($$row[0], $$row[1]) = flip_se($$row[0], $$row[1]) if $flip_bool;

				# check start stop #
				die " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->$table; start = $$row[0], end = $$row[1]\n"
					#print STDERR " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->$table; start = $$row[0], end = $$row[1]\n"
					if $$row[0] < 0 || $$row[1] < 0;
		
					#print Dumper @$row;
				}
			}
		}
	
		#print Dumper @$loci_r;
		#print Dumper %$tables_oi_r; 
		#exit;
	}

sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}
	
sub check_scaffold_span{
	my ($ret, $locus, $start, $end, $table) = @_;
	die " ERROR: no scaffold spans positions start->$start, end->$end for locus->$$locus[0], table->$table\n"
		unless @$ret;

	if(scalar @$ret > 1){		# error & die
		print STDERR " ERROR: start=$start, end=$end spans >1 scaffold in table->$table\n";
		print STDERR " LOCUS: ", join(", ", $$locus[0], @$locus[2..$#$locus]), "\n";	
			#print Dumper @$ret;
		exit;
		}	
	}

sub get_other_tables_oi{
# other tables that need position updating:
## spacers
## directrepeats
## leaderseqs
## genes

	my ($dbh, $loci_r, $tables_oi_r, $table, $prefix, $prefix2) = @_;
	$prefix2 = $prefix unless $prefix2;
	
	my $cmd = "SELECT $prefix\_start, $prefix\_end, $prefix2\_ID FROM $table where locus_id = ?";
	$cmd =~ s/[\n\r]/ /g;
	my $sql = $dbh->prepare($cmd);

	foreach my $loci (@$loci_r){
		$sql->bind_param(1, $$loci[0]);
		$sql->execute();
		my $ret = $sql->fetchall_arrayref();
		$tables_oi_r->{$table}{$$loci[0]} = $ret;	
		}
	}	

sub get_loci_oi{
# getting loci with genbank of interest #
	my ($dbh, $merged_genbank) = @_;

	my @parts = File::Spec->splitpath($merged_genbank);
	
	my $cmd = "SELECT locus_id, scaffold, locus_start, locus_end, 
operon_start, operon_end, CRISPR_array_start, CRISPR_array_end
FROM Loci where genbank = \'$parts[2]\'";
	$cmd =~ s/[\n\r]/ /g;
	my $ret = $dbh->selectall_arrayref($cmd);
	
	return $ret;
	}

sub make_loc_itree{
# making a location itree #
## input: genome-wide loc; output: unmerged scaf loc ##
	my ($scaf_index_r, $merged_r, $unmerged_r) = @_;
	
	my $itree = Set::IntervalTree->new();

	foreach my $merged_contig (keys %$scaf_index_r){
		foreach my $unmerged_contig (@{$scaf_index_r->{$merged_contig}}){
			# addig genome-wide position info to unmerged contigs #
			$unmerged_r->{$unmerged_contig}{"gw_start"} = $merged_r->{$merged_contig}{"start"};
			$unmerged_r->{$unmerged_contig}{"gw_end"} = $merged_r->{$merged_contig}{"end"};
			
			# deleting sequence info #
			delete $unmerged_r->{$unmerged_contig}{"seq"};
			
			# adding to interval tree #
			$itree->insert($unmerged_r->{$unmerged_contig},
						$merged_r->{$merged_contig}{"start"},
						$merged_r->{$merged_contig}{"end"});
			}
		}

	#print_itree($itree, 1, 10000);
	#my $res_r = $itree->fetch(200, 300);
	#print Dumper $res_r; 
	
	return $itree;
	}

sub print_itree{
	my ($itree, $start, $end) = @_;
	
	for my $i ($start..$end){
		my $res = $itree->fetch($i, $i);
		
		print join("\t", "pos:$i", join("::", @$res)), "\n";
		}
	}

sub find_same_length{
# foreach merged-unmerged, find scaffolds that match in length #
## finding 1-to-1 associations in length ##
	my ($gen_list_r, $unmerged_r, $merged_r) = @_;
	
	my %scaf_index;
	foreach my $mcontig (keys %$merged_r){				# merged => unmerged
		# comparing merged & unmerged contig lengths #
		foreach my $ucontig (keys %$unmerged_r){
			push( @{$scaf_index{$mcontig}}, $ucontig) if
				$merged_r->{$mcontig}{"length"} == $unmerged_r->{$ucontig}{"length"};
			}
		# checking that contig had at least on hit #
		unless (exists $scaf_index{$mcontig} ){
			die " ERROR: Contig$mcontig did not have a match in the merged genbank!\n";
			}
		}		

	# making list of taxa that need blasting #
	my %needs_blasting;
	foreach my $contig (keys %scaf_index){
		if(scalar @{$scaf_index{$contig}} > 1){			# if to many associations
			$needs_blasting{"merged"}{$contig} = 1;
			foreach ( @{$scaf_index{$contig}} ){
				$needs_blasting{"unmerged"}{$_} = 1;
				}
			delete $scaf_index{$contig};				# deleting from index; non 1-to-1
			}
		elsif( scalar @{$scaf_index{$contig}} < 1){
			die " LOGIC ERROR: $!\n";
			}
		}
	
		#print Dumper %needs_blasting; exit;
		#print Dumper %scaf_index; exit;	
	return \%scaf_index, \%needs_blasting;
	}

sub parse_unmerged_genbank{
# pulling out start-stop & sequence for all 'source' features #
	my ($gen_list_r, $gen_io_r, $unmerged) = @_;
	
	my %unmerged;			
	print STDERR "...loading: $unmerged\n" unless $verbose;
	while(my $seqo = $gen_io_r->{$unmerged}->next_seq){
		my @sources = grep{ $_->primary_tag eq "source" } $seqo->get_SeqFeatures;
		foreach my $source (@sources){
			my $contig_id = $source->seq->display_id; 
			$unmerged{$contig_id}{"start"} = $source->location->start;
			$unmerged{$contig_id}{"end"} = $source->location->end;
			$unmerged{$contig_id}{"length"} = abs($unmerged{$contig_id}{"end"} -
															$unmerged{$contig_id}{"start"});
			$unmerged{$contig_id}{"seq"} = $source->seq->seq;
			$unmerged{$contig_id}{"seq"} =~ tr/A-Z/a-z/;
			$unmerged{$contig_id}{"id"} = $contig_id;
			}
		}

	
		#print Dumper %unmerged; exit;
	return \%unmerged;
	}

sub parse_merged_genbank{
# pulling out start-stop & sequence for all 'source' features #
## start-end is absolute positioning ##
	my ($gen_list_r, $gen_io_r, $merged) = @_;
	
	my %merged;
	my $source_cnt = 0;
	print STDERR "...loading: $merged\n" unless $verbose;
	while(my $seqo = $gen_io_r->{$merged}->next_seq){
		my @sources = grep{ $_->primary_tag eq "source" } $seqo->get_SeqFeatures;
		foreach my $source (@sources){
			$merged{$source_cnt}{"start"} = $source->location->start;
			$merged{$source_cnt}{"end"} = $source->location->end;
			$merged{$source_cnt}{"length"} = abs($merged{$source_cnt}{"end"} -
														$merged{$source_cnt}{"start"});
			$merged{$source_cnt}{"seq"} = $source->seq->seq;
			$merged{$source_cnt}{"seq"} =~ tr/A-Z/a-z/;
			$source_cnt++;
			}		
		}
		#print Dumper %merged; exit;
	return \%merged;
	}

sub load_genbank_io{
# parsing the merged genbank into scaffold lengths & sequence #
	my ($gen_list_r) = @_;
	
	my %gen_io;
	foreach my $file (%$gen_list_r){
		$gen_io{$file} = Bio::SeqIO->new(-format => "genbank", -file => $file);
		}
	
	return \%gen_io;
	}

sub get_genbank_list{
# getting list of genbanks #
## 2 column: unmerged, merged ##
	
	my %list;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @tmp = split /\t/;
		die " ERROR: the genbank list must have 2 columns!\n"
			unless scalar @tmp == 2;
		
		$list{$tmp[0]} = $tmp[1];	# unmerged => merged
		}
	
		#print Dumper %list; exit;
	return \%list;
	}


__END__

=pod

=head1 NAME

CLdb_positionByScaffold.pl -- Convert positional values from merged (whole-genome) to by-scaffold

=head1 SYNOPSIS

CLdb_positionByScaffold.pl [Flags] < genbank_unmerged_merged.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_positionByScaffold.pl

=head1 DESCRIPTION

Convert merged genbank position information (merged with EMBOSS 'union')
to scaffold-based position info. 

=head2 Input

The input must be a 2 column tab-delimited table:
'unmerged_genbank_file_name', 'merged_genbank_file_name'

Example of unmerged genbank file format: RAST genbank output for a draft genome.

Updates not committed until all of the values have been converted
(in case something goes wrong, you won't have a hybrid-position database)

=head2 Matching unmerged and merged scaffolds

All scaffold start-end info is obtained form the 'source' tags
in both the unmerged and merged genbank files. 1-to-1
connections between scaffolds are determined in a number of steps:

=over 

=item 1) Comparing scaffold lengths.

=item 2) Comparing scaffold sequences. No duplicate scaffolds allowed!

=back

=head1 EXAMPLES

=head2 Normal usage:

CLdb_positionByScaffold.pl -d CLdb.sqlite < unmerged_merged.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

