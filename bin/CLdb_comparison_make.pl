#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_comparison_make.pl -- making comparisons table for plotting

=head1 SYNOPSIS

CLdb_comparison_make.pl [flags] < dna_segs.txt > compare.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -ITEP  <char>

ITEP database. Used to get BLASTp info.

=item -spacer  <char>

Use spacer sequence clustering ('cluster') or pairwise blast ('blast')? ['cluster']

=item -cutoff  <int>

If -spacer 'cluster': Clustering cutoff. [1]

If -spacer 'blast': blast pident cutoff. [100]

=item -blastp  <int>

BLASTp percent identity cutoff. [30]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_comparison_make.pl

=head1 DESCRIPTION

Make a table of comparisons between genes & 
spacers in the plot. A dna_segs table is used
for determining which comparisons to add to
the comparisons table.

Spacer comparisons can be based on clustering 
of spacer sequences (CLdb_hclusterArrays.pl must
be run prior), or by pairwise spacer blasts 
(CLdb_spacerPairwiseBlast.pl must be run prior).

Gene comparisons are obtained from ITEP.

=head1 EXAMPLES

=head2 Just comparisons among spacers (using BLASTn info)

CLdb_comparison_make.pl -da CLdb.sqlite  < dna_segs.txt > compare.txt

=head2 Just comparisons among spacers (using clustering info)

CLdb_comparison_make.pl -da CLdb.sqlite -s cluster < dna_segs.txt > compare.txt

=head2 Comparisons among spacers and genes

CLdb_comparison_make.pl -da CLdb.sqlite -I ITEP_database.sqlite < dna_segs.txt > compare.txt

=head2 

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
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $CLdb_sqlite, $ITEP_sqlite);
my $spacer_cutoff = 1;
my $blastp_cutoff = 30;
my $spacer_opt = "cluster";
GetOptions(
	   "database=s" => \$CLdb_sqlite,
	   "ITEP=s" => \$ITEP_sqlite,
	   "cutoff=f" =>  \$spacer_cutoff,
	   "spacer=s" => \$spacer_opt, 		# 'blast' or 'cluster'
	   "blastp=f" => \$blastp_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# checking CLdb #
file_exists($CLdb_sqlite, "database");
$spacer_cutoff = 100 if $spacer_opt =~ /blast/i;
print STDERR "...WARNING: no ITEP db provided ('-i')! No blastp information will be included!\n"
	unless $ITEP_sqlite;
	
# checking ITEP input #
die " ERROR: cannot find ITEP database file!\n"
	if $ITEP_sqlite && ! -e $ITEP_sqlite;


#--- MAIN ---#
# connect 2 CLdb #
my $dbh = connect2db($CLdb_sqlite);

# connect to ITEP #
my $dbh_ITEP;
if($ITEP_sqlite){
	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	$dbh_ITEP = DBI->connect("dbi:SQLite:dbname=$ITEP_sqlite", '','', \%attr) 
		or die " Can't connect to $ITEP_sqlite!\n";	
	}

# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r) = load_dna_segs();
my %compare;

# check for spacer_cluster entries #
die "ERROR: no entries in spacer_cluster table!\n"
	if exists $dna_segs_r->{"spacer"} && ! n_entries($dbh, "spacer_clusters");

# querying CLdb for spacer comparisons #
if( exists $dna_segs_r->{"spacer"} ){
	if($spacer_opt =~ /cluster/i){
		get_spacer_cluster($dbh, $dna_segs_r, $dna_segs_order_r, 
						\%compare, $spacer_cutoff);
		}
	elsif($spacer_opt =~ /blast/i){
		get_spacer_pairwise_blast($dbh, $dna_segs_r, $dna_segs_order_r, 
						\%compare, $spacer_cutoff);
		}
	}

# querying ITEP #
if( exists $dna_segs_r->{"gene"} && $dbh_ITEP ){
	my $blastp_res_r = get_ITEP_blastp($dbh_ITEP, $dna_segs_r, $dna_segs_order_r, $blastp_cutoff);
	get_CLdb_gene_start_end($dbh, $dna_segs_r, $dna_segs_order_r, $blastp_res_r, \%compare);
	}

# writing out compare table #
write_compare(\%compare, $dna_segs_order_r);

# disconnect #
$dbh->disconnect();
$dbh_ITEP->disconnect() if $dbh_ITEP;
exit;


### Subroutines
sub write_compare{
# writing out compare table #
	my ($compare_r, $dna_segs_order_r) = @_;
	
	# header #
	print join("\t", qw/start1 end1 start2 end2 col dna_segs_id1 dna_segs_id2 feat feat_id sfeat_id/), "\n";
	
	# body #
	foreach my $dna_seg1 (@$dna_segs_order_r){
		foreach my $feat (keys %$compare_r){
			foreach my $dna_seg2 (keys %{$compare_r->{$feat}{$dna_seg1}}){
				foreach my $feat_id (keys %{$compare_r->{$feat}{$dna_seg1}{$dna_seg2}}){
					foreach my $sfeat_id (keys %{$compare_r->{$feat}{$dna_seg1}{$dna_seg2}{$feat_id}}){
						my $row = $compare_r->{$feat}{$dna_seg1}{$dna_seg2}{$feat_id}{$sfeat_id};
						print join("\t", @$row, 
						$dna_seg1, $dna_seg2, $feat, $feat_id, $sfeat_id), "\n";
						}
					}
				}
			}
		}
	}

sub get_CLdb_gene_start_end{
# querying CLdb for geneID start-end #
	my ($dbh, $dna_segs_r, $dna_segs_order_r, $blastp_res_r, $compare_r) = @_;
	
	# preparing query for query gene #
	## 'AND locus_id = ?' newly added; not sure if needed 
	my $query = "
SELECT gene_start, gene_end
FROM genes
WHERE gene_id = ?
AND locus_id = ?				
";
	print STDERR "$query\n" if $verbose;
	$query =~ s/\r|\n/ /g;
	my $sth1 = $dbh->prepare($query);

	# preparing query for subject gene #
	$query = "
SELECT gene_start, gene_end
FROM genes
WHERE gene_id = ?
AND locus_id = ?
";
	print STDERR "$query\n" if $verbose;
	$query =~ s/\r|\n/ /g;
	my $sth2 = $dbh->prepare($query);	
	
	# querying each geneID #
	my %blastp_res;
	for my $i (0..($#$dna_segs_order_r-1)){
		my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];
		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
			
		next unless defined $locus_id1 && exists $dna_segs_r->{gene}{$dna_seg_id1}{$locus_id1};
		my %keep_cnt;
		foreach my $feat_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id1}{$locus_id1}}){
			# seeing if feat_d is in blast res #
			next unless exists $blastp_res_r->{$feat_id};
			
			# getting start-end of query gene #
			$sth1->bind_param(1, $feat_id);
			$sth1->bind_param(2, $locus_id1);			# newly added; not sure if needed 
			$sth1->execute();
			my $res_q = $sth1->fetchall_arrayref();
			
			if (! @$res_q){
				#print STDERR " WARNING: no gene start-end found for $feat_id. Skipping\n";
				next;
				}
			else{ 
				$keep_cnt{'query'}++; 
				}
			
			# getting start-end of subject gene #
				#print Dumper "q $feat_id";
			foreach my $subject_gene_id (keys %{$blastp_res_r->{$feat_id}}){
					#print Dumper "s $subject_gene_id";
				$sth2->bind_param(1, $subject_gene_id);
				$sth2->bind_param(2, $locus_id2);
				$sth2->execute();
				my $res_s = $sth2->fetchall_arrayref();

				if(! @$res_s){
					next;
					}
				else{
					$keep_cnt{'subject'}++;
					}
			
				# all + strand #
				## needed to inversions in blast connections for genoPlotR ##
				#($$res_q[0][0], $$res_q[0][1]) = to_pos_strand($$res_q[0][0], $$res_q[0][1]);
				#($$res_s[0][0], $$res_s[0][1]) = to_pos_strand($$res_s[0][0], $$res_s[0][1]);
			
				# loading hash #
				$compare_r->{"gene"}{$dna_seg_id1}{$dna_seg_id2}{$feat_id}{$subject_gene_id} =
					[$$res_q[0][0], $$res_q[0][1],
					 $$res_s[0][0], $$res_s[0][1], 	
					 $blastp_res_r->{$feat_id}{$subject_gene_id}
					 ];
				}
			}
		
		# status #
		print STDERR "Number of features for locus '$locus_id1':\t\t", 
			scalar keys %{$dna_segs_r->{"gene"}{$dna_seg_id1}{$locus_id1}}, "\n";
		$keep_cnt{'query'}  = 0 unless exists $keep_cnt{'query'};
		print STDERR "Number of query features with CLdb entry:\t",
					$keep_cnt{'query'}, "\n";
		$keep_cnt{'subject'}  = 0 unless exists $keep_cnt{'subject'};
		print STDERR "Number of subject features with CLdb entry:\t",
					$keep_cnt{'subject'}, "\n\n";
		}
		#print Dumper "here", %{$compare_r->{"gene"}}; exit;
	}
	
sub to_pos_strand{
# flipping to + strand #
	my ($start, $end) = @_;
	if($start > $end){ return $end, $start; }
	else{ return $start, $end; }
	}

sub get_ITEP_blastp{
# getting blastp info from ITEP #
	my ($dbh_ITEP, $dna_segs_r, $dna_segs_order_r, $blastp_cutoff) = @_;
	
	# status #
	print STDERR "### Querying ITEP ###\n";
	
	# preparing query #
	my $query = "
SELECT querygene, targetgene, pctid
FROM blastres_selfbit
WHERE querygene=?
AND pctid >= $blastp_cutoff
";
	print STDERR "$query\n" if $verbose;
	
	$query =~ s/\r|\n/ /g;

	my $sth = $dbh_ITEP->prepare($query);
	
	# querying each geneID #
	my %blastp_res;
	for my $i (0..($#$dna_segs_order_r-1)){
	  	my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];

		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
		
		# sanity check 
		unless( defined $locus_id1 && exists $dna_segs_r->{gene}{$dna_seg_id1}{$locus_id1} ){
		  print STDERR "WARNING: no genes found for $dna_seg_id1\n";
		  next;
		}		  

		# querying ITEP
		foreach my $feat_id (keys %{$dna_segs_r->{"gene"}{$dna_seg_id1}{$locus_id1}}){
			$sth->bind_param(1, $feat_id);
			$sth->execute();
			my $res = $sth->fetchall_arrayref();
			
			print STDERR " WARNING: $feat_id not found in database!\n"
				unless @$res;
			next unless @$res;			# should warn if gene not found in db
			
			# loading hash #
			foreach my $row (@$res){
				if(exists $blastp_res{$$row[0]}{$$row[1]}){
					$blastp_res{$$row[0]}{$$row[1]} = $$row[2]
						if $$row[2] > $blastp_res{$$row[0]}{$$row[1]};
					}
				else{ $blastp_res{$$row[0]}{$$row[1]} = $$row[2] }
				}
			}
		}
		#print Dumper %blastp_res; exit;
	return \%blastp_res;
	}

sub get_spacer_pairwise_blast{
# querying CLdb for spacer_hclust info #
	my ($dbh, $dna_segs_r, $dna_segs_order_r, $compare_r, $spacer_cutoff) = @_;
	
	# status #
	print STDERR "### Getting spacer pairwise blast results ###\n";
	
	# preparing query #
	my $query = "
SELECT  
b.spacer_start, b.spacer_end, 
c.spacer_start, c.spacer_end, 
a.pident, 
a.Query_locus_ID, a.Subject_locus_ID,
a.Query_spacer_ID, a.Subject_spacer_ID
FROM spacer_pairwise_blast a, spacers b, spacers c
WHERE a.query_locus_id = b.locus_id 
AND a.subject_locus_id = c.locus_id
AND a.subject_spacer_id = c.spacer_id
AND a.query_spacer_id = ? 
AND b.spacer_id = ? 
AND a.query_locus_id = ?
AND a.subject_locus_id = ?
AND a.pident >= ?
";
	print STDERR "$query\n" if $verbose;
	
	$query =~ s/\r|\n/ /g;

	my $sth = $dbh->prepare($query);
	
	# querying each spacer ID for matches against adjacent locus #
	for my $i (0..($#$dna_segs_order_r-1)){
		my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];
		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
			
		foreach my $feat_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}{$locus_id1}}){
			#print Dumper $feat_id, $spacer_cutoff, $locus_id1, $locus_id2; exit;
			$sth->bind_param(1, $feat_id);
			$sth->bind_param(2, $feat_id);			
			$sth->bind_param(3, $locus_id1);			# query locus
			$sth->bind_param(4, $locus_id2);			# subject locus
			$sth->bind_param(5, $spacer_cutoff);
			$sth->execute();
			my $res = $sth->fetchall_arrayref();
			
			# skipping any spacers not in both loci #
			next unless @$res;
			
			# loading hash #
			foreach my $row (@$res){
				die " ERROR: feature IDs don't match"
					unless $feat_id eq $$row[7];
				my $sfeat_id = $$row[8];
				$compare_r->{"spacer"}{$dna_seg_id1}{$dna_seg_id2}{$feat_id}{$sfeat_id} = [@$row[0..4]];
				}
			}
		}
		#print Dumper %$compare_r; exit;
	}

sub get_spacer_cluster{
# querying CLdb for spacer_cluster info #
	my ($dbh, $dna_segs_r, $dna_segs_order_r, $compare_r, $spacer_cutoff) = @_;
	
	# status #
	print STDERR "### Getting spacer clusters ###\n";
	
	# preparing query #
	my $query = "
SELECT 
c.spacer_start, c.spacer_end,
d.spacer_start, d.spacer_end,
c.spacer_id, d.spacer_id
FROM
(SELECT a.locus_id, a.spacer_id, b.spacer_start, b.spacer_end, a.cluster_id
FROM spacer_clusters a, spacers b
WHERE a.locus_id = b.locus_id
AND a.spacer_id = b.spacer_id
AND a.spacer_id= ?
AND a.cutoff = ?
AND a.locus_id = ?) c,
(SELECT a.locus_id, a.spacer_id, b.spacer_start, b.spacer_end, a.cluster_id
FROM spacer_clusters a, spacers b
WHERE a.locus_id = b.locus_id
AND a.spacer_id = b.spacer_id
AND a.cutoff = ?
AND a.locus_id = ?) d
WHERE c.cluster_id = d.cluster_id
";
	print STDERR "$query\n" if $verbose;
	
	$query =~ s/\r|\n/ /g;

	my $sth = $dbh->prepare($query);
	
	# querying each spacer ID for matches against adjacent locus #
	for my $i (0..($#$dna_segs_order_r-1)){
		my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];
		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
		foreach my $feat_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}{$locus_id1}}){
				#print Dumper $feat_id, $spacer_cutoff, $locus_id1, $locus_id2; exit;
			$sth->bind_param(1, $feat_id);		
			$sth->bind_param(2, $spacer_cutoff);
			$sth->bind_param(3, $locus_id1);
			$sth->bind_param(4, $spacer_cutoff);
			$sth->bind_param(5, $locus_id2);
			$sth->execute();
			my $res = $sth->fetchall_arrayref();
			
			# skipping any spacers not in both loci #
			next unless @$res;
			
			# adding seqID #
			map{ push(@$_, $spacer_cutoff * 100) } @$res;
			
			# loading hash #
			foreach my $row (@$res){
				my $sfeat_id = $$row[5];
				$compare_r->{"spacer"}{$dna_seg_id1}{$dna_seg_id2}{$feat_id}{$sfeat_id} = [@$row[0..3,$#$row]];
				}
			}
		}
		#print Dumper %$compare_r; exit;
		#print "exit\n"; exit;
	}

sub load_dna_segs{
# loading dna_segs from stdin #

	my @dna_segs_order;
	my %dna_segs;
	my %header;
	my %name_index;
	while(<>){
		chomp;
		next if /^\s*$/;

		# header #
		if($.==1){
			# indexing header #
			my @tmp = split /\t/;
			for my $i (0..$#tmp){
				$header{$tmp[$i]} = $i;
				}
			
			# checking for existence of columns #
			my @check = qw/locus_id dna_segs_id feat feat_id/;
			map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
				unless exists $header{$_} } @check;
			}
		else{
			my @line = split /\t/;
			my $locus_id = $line[$header{"locus_id"}];
			my $dna_seg_id = $line[$header{"dna_segs_id"}];
			my $feat = $line[$header{"feat"}];
			my $feat_id = $line[$header{"feat_id"}];						
			
			# order of loci #
			push( @dna_segs_order, $dna_seg_id)
				unless exists $name_index{$dna_seg_id};
			$name_index{$dna_seg_id} = 1;
			
			# dna_segs #
			die " ERROR: $feat -> $dna_seg_id -> $locus_id -> $feat_id is not unique!\n"
				if exists $dna_segs{$feat}{$dna_seg_id}{$locus_id}{$feat_id};
			$dna_segs{$feat}{$dna_seg_id}{$locus_id}{$feat_id} = \@line;
			
			}
		}
		
		#print Dumper @dna_segs_order; exit;
		#print Dumper %dna_segs; exit;
	return \%dna_segs, \@dna_segs_order, \%header;
	}

sub get_spacer_group_OLD{
# querying CLdb for spacer_group info #
	my ($dbh, $dna_segs_r, $dna_segs_order_r, $compare_r) = @_;
	
	# status #
	print STDERR "### Getting spacer groups ###\n";
	
	# preparing query #
	my $query = "
SELECT 
c.spacer_start, c.spacer_end,
d.spacer_start, d.spacer_end,
c.spacer_id, d.spacer_id
FROM
(SELECT a.locus_id, b.spacer_id, b.spacer_start, b.spacer_end, b.spacer_group
FROM loci a, spacers b
WHERE a.locus_id = b.locus_id
AND b.spacer_id= ?
AND a.locus_id = ?) c,
(SELECT a.locus_id, b.spacer_id, b.spacer_start, b.spacer_end, b.spacer_group
FROM loci a, spacers b
WHERE a.locus_id = b.locus_id
AND a.locus_id = ?) d
WHERE c.spacer_group = d.spacer_group
";
	print STDERR "$query\n" if $verbose;
	
	$query =~ s/\r|\n/ /g;

	my $sth = $dbh->prepare($query);
	
	# querying each spacer ID for matches against adjacent locus #
	for my $i (0..($#$dna_segs_order_r-1)){
		my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];
		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
		foreach my $feat_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}{$locus_id1}}){
				#print Dumper $feat_id, $spacer_cutoff, $locus_id1, $locus_id2
				#		unless exists 
			
			$sth->bind_param(1, $feat_id);		
				#$sth->bind_param(2, $spacer_cutoff);
			$sth->bind_param(2, $locus_id1);
				#$sth->bind_param(3, $spacer_cutoff);
			$sth->bind_param(3, $locus_id2);
			$sth->execute();
			my $res = $sth->fetchall_arrayref();
			
			# skipping any spacers not in both loci #
			next unless @$res;
			
			# adding seqID #
			map{ push(@$_, $spacer_cutoff * 100) } @$res;
			
			# loading hash #
			foreach my $row (@$res){
				my $sfeat_id = $$row[5];
				$compare_r->{"spacer"}{$dna_seg_id1}{$dna_seg_id2}{$feat_id}{$sfeat_id} = [@$row[0..3,$#$row]];
				}
			}
		}
		#print Dumper %$compare_r; exit;
		#print "exit\n"; exit;
	}
