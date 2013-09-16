#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Algorithm::NaiveBayes;
use List::Util qw/shuffle/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $quiet);
my $kmer_len = 8;
my $score_cutoff = 0.99;
my $boot_cutoff = 95;
my $iter = 100;
my $frac_arrays = 0.1;
my $nboot = 100;
GetOptions(
	   "database=s" => \$database_file,
	   "kmer=i" => \$kmer_len,
	   "score=f" => \$score_cutoff,			# bayes score cutoff
	   "iterations=i" => \$iter, 			# number of test iterations
	   "fraction=f" => \$frac_arrays, 		# fraction of arrays to use as query
	   "bootstrap=i" => \$nboot,					# number of bootstrapped classifications to perform
	   "cutoff=i" => \$boot_cutoff, 			# bootstrap cutoff
	   "verbose" => \$verbose,
	   "quiet" => \$quiet, 
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);

### MAIN
# iterating #
my %tpfn;			# true-positives, false-negatives
for my $i (1..$iter){	
	# connect 2 db #
	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
		or die " Can't connect to $database_file!\n";

	# making a temporary loci table #
	make_temp_loci_table($dbh);

	# randomly unclassifying a certain fraction of arrays #
	my $array_count = get_array_count($dbh);
	my $true_class_r = unclassify_frac($dbh, $array_count, $frac_arrays);

	# getting training dataset (classified, intact arrays) #
	my $arrays_class_r = get_classified_arrays($dbh, $kmer_len);

	# getting unclassified arrays #
	my ($arrays_unclass_r, $kmers_r) = get_unclassified_arrays($dbh, $kmer_len);

	# training classifier #
	my $nb = train_classifier($arrays_class_r);

	# classifying #
	my $res_r = classify_arrays($nb, $arrays_unclass_r);

	# filter results #
	$res_r = filter_bayes_scores($res_r, $score_cutoff);

	# bootstrapping #
	my %boot_score;
	for my $ii (1..$nboot){
		# getting subset of kmers #
		my $kmer_sub_r = subset_kmer_pool($kmers_r, $kmer_len);
		
		# classifying #
		my $boot_res_r = classify_arrays($nb, $kmer_sub_r);

		# filter results #
		$boot_res_r = filter_bayes_scores($boot_res_r, $score_cutoff, 1);		
		
		# does bootstrap classification match original classification? #
		score_bootstrap($res_r, $boot_res_r, \%boot_score);
		}
		
	# filter by bootstrap score #
	filter_by_boot_scores($res_r, \%boot_score, $boot_cutoff, $nboot) if $nboot;

	# compare to true classifications #
	compare_classifications($true_class_r, $res_r, $i, \%tpfn);

	# disconnect #
	$dbh->disconnect();
	}

write_summary(\%tpfn, $iter);


### Subroutines 
sub write_summary{
	my ($tpfn_r) = @_;

	# making array of values #
	my @tp;
	foreach my $iter (keys %{$tpfn_r->{"true_pos"}}){
		push @tp, $tpfn_r->{"true_pos"}{$iter};
		}
	my @fn;
	foreach my $iter (keys %{$tpfn_r->{"false_neg"}}){
		push @fn, $tpfn_r->{"false_neg"}{$iter};
		}
	my @na;
	foreach my $iter (keys %{$tpfn_r->{"NA"}}){
		push @na, $tpfn_r->{"NA"}{$iter};
		}	
	my @tpfn;
	foreach my $iter (keys %{$tpfn_r->{"true_pos"}}){
		push @tpfn, $tpfn_r->{"true_pos"}{$iter} / $tpfn_r->{"false_neg"}{$iter};
		}
	
	# summary #
	print join("\t", "TP:FN",
		sprintf("%.3f", average( \@tpfn )), 
		sprintf("%.3f", stdev(\@tpfn))
		), "\n";
	print join("\t", "True positives:",
		sprintf("%.3f", average( \@tp )), 
		sprintf("%.3f", stdev(\@tp))
		), "\n";
	print join("\t", "False negatives:",
		sprintf("%.3f", average( \@fn )), 
		sprintf("%.3f", stdev(\@fn))
		), "\n";
	print join("\t", "NA:",
		sprintf("%.3f", average( \@na )), 
		sprintf("%.3f", stdev(\@na))
		), "\n";

	}

sub compare_classifications{
	my ($true_class_r, $res_r, $iter_num, $tpfn_r) = @_;
	
	foreach my $locus_id (keys %$true_class_r){
		die " ERROR: locus \"cli$locus_id\" not found in classifed values\n"
			unless exists $res_r->{$locus_id};
		
		if($res_r->{$locus_id} eq ""){
			$tpfn_r->{"NA"}{$iter_num}++;
			}
		elsif($true_class_r->{$locus_id} eq $res_r->{$locus_id}){
			$tpfn_r->{"true_pos"}{$iter_num}++;
			}
		else{ 
			print STDERR "FALSE negative! Locus_id: $locus_id; True classification '", 
					$true_class_r->{$locus_id},
					"'; Classified as: '", $res_r->{$locus_id}, "'\n"
					if $verbose;
			$tpfn_r->{"false_neg"}{$iter_num}++;
			}
		}
	}

sub filter_by_boot_scores{
# filtering classifications by bootstrap scores #
## scores must be >= cutoff; or no classification given ##
	my ($res_r, $boot_score_r, $boot_cutoff, $nboot) = @_;
	
	foreach my $locus_id (keys %$res_r){
		$boot_score_r->{$locus_id} = 0 unless exists $boot_score_r->{$locus_id};
		if($boot_score_r->{$locus_id} < $boot_cutoff){			# no classification unless meets cutoff
			$res_r->{$locus_id} = "";
			}
		}	
	}

sub score_bootstrap{
	my ($res_r, $boot_res_r, $boot_score_r) = @_;
	#print Dumper $res_r, $boot_res_r; exit;
	
	foreach my $locus_id (keys %$res_r){	
		die " ERROR: Locus$locus_id not found in bootstrapped classifications!\n"
			unless exists $boot_res_r->{$locus_id};
		
		if( $res_r->{$locus_id} eq $boot_res_r->{$locus_id} ){	# same classification
			$boot_score_r->{$locus_id}++
			}
		}
	}

sub	subset_kmer_pool{
# getting subset of kmers (fraction = 1/$kmer_len) #
	my ($kmers_r, $kmer_len) = @_;
	
	my %subset;
	foreach my $locus_id (keys %$kmers_r){
		my $fract = sprintf("%.0f", scalar @{$kmers_r->{$locus_id}} / $kmer_len);
		my @shuffled_kmers = shuffle(@{$kmers_r->{$locus_id}});
		for my $i (0..($fract - 1)){
			$subset{$locus_id}{$shuffled_kmers[$i]}++;
			}
		}

	return \%subset;
	}

sub filter_bayes_scores{
# filtering bayes scores from NB classification #
	my ($res_r, $score_cutoff, $quiet) = @_;
	
	my %passed;
	foreach my $locus_id (keys %$res_r){
		foreach my $label ( sort{$res_r->{$locus_id}{$b}<=>$res_r->{$locus_id}{$a}} 
							keys %{$res_r->{$locus_id}}){
			
			my $top_score = $res_r->{$locus_id}{$label};
			if($top_score < $score_cutoff){	# less than cutoff, cannot be classified
				print STDERR " WARNING! Cannot classify Locus: '$locus_id'!\tTop score only $top_score\n"
					unless $quiet;
				$passed{$locus_id} = "";			# label = no subtype
				}
			else{
				$passed{$locus_id} = $label;		# label = subtype
				}
			last;
			}
		}
	
	return \%passed;
	}

sub classify_arrays{
# classifying subtype of unclassified arrays #
	my ($nb, $arrays_r) = @_;
	
	my %res;
	foreach my $locus_id (keys %$arrays_r){
		$res{$locus_id} = $nb->predict(
			attributes => $arrays_r->{$locus_id}
     		);
     	}

		#print Dumper %res; 
	return \%res;
	}

sub train_classifier{
# training classifier with classified arrays #
	my ($arrays_r) = @_;
	
	my $nb = Algorithm::NaiveBayes->new;

	foreach my $subtype (keys %$arrays_r){
		$nb->add_instance(
			attributes => $arrays_r->{$subtype},
     		label => $subtype
     		);
     	}

	$nb->train;
	
		#print Dumper $nb; exit;
	return $nb;
	}

sub get_unclassified_arrays{
	my ($dbh, $kmer_len) = @_;
	
	# make query #
	my $query = "SELECT 
DRs.DR_sequence,
loci_temp.locus_id
FROM DRs, loci_temp 
WHERE DRs.locus_id = loci_temp.locus_id 
AND (loci_temp.subtype IS NULL
OR loci_temp.subtype = '')
";
	$query =~ s/\r|\n/ /g;
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no unclassifed arrays that need classifying!\n"
		unless $$ret[0];
	
	my (%arrays, %kmers);
	foreach my $row (@$ret){
		for (my $i = 0; $i<=(length $$row[0]) - $kmer_len; $i++){
			my $kmer = substr($$row[0], $i, $kmer_len);
			push( @{$kmers{$$row[1]}}, $kmer );
			$arrays{$$row[1]}{$kmer}++;		# locus_id => kmer => count
			}
		}
	
	# status #
	print STDERR "...Number of arrays to classify:\t", scalar keys %arrays, "\n"
		unless $quiet;
	
		#print Dumper %arrays; exit;
	return \%arrays, \%kmers;
	}

sub get_classified_arrays{
# using arrays that were classified and have 'intact' operons #
	my ($dbh, $kmer_len) = @_;
	
	# make query #
	my $query = "SELECT 
DRs.DR_sequence,
loci_temp.subtype,
loci_temp.locus_id
FROM DRs, loci_temp 
WHERE DRs.locus_id = loci_temp.locus_id 
AND loci_temp.operon_status = 'intact'
AND loci_temp.subtype IS NOT NULL
";
	$query =~ s/\r|\n/ /g;
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no classified arrays to train the classifier!\n"
		unless $$ret[0];
	
	my %arrays;
	foreach my $row (@$ret){
		next if $$row[1] eq "";
		for (my $i = 0; $i<=(length $$row[0]) - $kmer_len; $i++){
			my $kmer = substr($$row[0], $i, $kmer_len);
			$arrays{$$row[1]}{$kmer}++;
			}
		}
		
	# counting arrays used #
	my %array_cnt;
	map{$array_cnt{$$_[2]}++} @$ret;

	# status #
	print STDERR "...Classifier trained with ", 
		scalar keys %array_cnt, " arrays and ",
		scalar @$ret, " direct_repeats\n"
			unless $quiet;
	
		#print Dumper %arrays; exit;
	return \%arrays;
	}
	
sub unclassify_frac{
# unclassifying a fraction of loci #
	my ($dbh, $array_count, $frac_arrays) = @_;
	
	# number of loci to unclassify #
	my $N_unclass = sprintf("%.0f", $array_count * $frac_arrays);
	print STDERR "...Number of loci unclassified: $N_unclass\n"
		unless $quiet;
	
	# getting arrays to unclassify #
my $q = "
SELECT locus_id, subtype
FROM Loci_temp
ORDER BY Random()
LIMIT $N_unclass;
";
	$q =~ s/\n|\r/ /g;
	
	my $ret = $dbh->selectall_arrayref($q);
	die " ERROR: no loci in temporary loci table!\n"
		unless @$ret;	

	## converting to hash ##
	my %true_class;
	map{$true_class{$$_[0]} = $$_[1]} @$ret;
	
	
	# unclassifying #
	$q = "
UPDATE loci_temp 
SET subtype = ''
WHERE locus_id = ?
";
	$q =~ s/\n|\r/ /g;
	
	foreach my $locus_id (keys %true_class){
		my $sql = $dbh->prepare($q);
		$sql->execute($locus_id);
		if($DBI::err){
			print STDERR "ERROR: could unclassify loci in temporary loci table! $DBI::errstr\n";
			}
		}
	
	return \%true_class;
	}

sub get_array_count{
# getting total number of loci in temp loci table #
	my ($dbh) = @_;
	
	my $q = "SELECT count(*) from loci_temp";
	my $ret = $dbh->selectall_arrayref($q);
	die " ERROR: no loci found in temporary loci table!\n"
		unless @$ret;
		
	return $$ret[0][0];
	}

sub make_temp_loci_table{
# making a temporary loci table #
	my ($dbh) = @_;
	
	# making temporary table #
	my $q = "
CREATE TEMPORARY TABLE Loci_temp (
Locus_ID	INTEGER	PRIMARY KEY,
Taxon_ID	TEXT	NOT NULL,
Taxon_Name	TEXT	NOT NULL,
Subtype	TEXT,
Scaffold	TEXT	NOT NULL,
Locus_Start	INTEGER	NOT NULL,
Locus_End	INTEGER	NOT NULL,
Operon_Start	INTEGER,
Operon_End	INTEGER,
Array_Start	INTEGER,
Array_End	INTEGER,
Operon_Status	TEXT	NOT NULL,
Array_Status	TEXT	NOT NULL,
Genbank_File	TEXT	NOT NULL,
Array_File	TEXT,
Fasta_File	TEXT,
Scaffold_count	INTEGER,
File_Creation_Date	TEXT,
Author	TEXT	NOT NULL
)";
	$q =~ s/\n|\r/ /g;
	
	my $sql = $dbh->prepare($q);
	$sql->execute();
	if($DBI::err){
		print STDERR "ERROR: could not make temporary loci table! $DBI::errstr\n";
		}

	# loading temp table with classifed, 'intact' arrays #
	$q = "
INSERT INTO Loci_temp 
SELECT * FROM Loci
WHERE subtype IS NOT NULL
AND operon_status = 'intact'
AND (array_start IS NOT NULL
OR array_start = '')
";
	$q =~ s/\n|\r/ /g;
	
	$sql = $dbh->prepare($q);
	$sql->execute();
	if($DBI::err){
		print STDERR "ERROR: could not load temporary loci table! $DBI::errstr\n";
		}		
	}

sub copy_file{
	my ($file, $dir) = @_;
	
	my @parts = File::Spec->splitpath($file);
	my $tmpfile = "$dir/$parts[2]";
	
	open IN, $file or die $!;
	open OUT, ">$tmpfile" or die $!;
	while (<IN>){ print OUT; }
	close IN; close OUT;
	
	return $tmpfile;
	}	

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}


__END__

=pod

=head1 NAME

CLdb_classifyArraysByDR_validate.pl -- validation of Bayes classifier for arrays

=head1 SYNOPSIS

CLdb_classifyArraysByDR_validate.pl [flags]

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -iterations

Number of test iterations to perform. [100]

=item -fraction

Fraction of dataset to test. [0.1]

=item -kmer

Kmer length using for naive bayes classifier. [8]

=item -score

Bayes score cutoff (0-1) for accepting classification (>=). [0.99]

=item -bootstrap

Number of bootstraps to perform the classification of each array

=item -cutoff 

Bootstrap score cutoff (>=). [95]

=item -v 	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_classifyArraysByDR_validate.pl

=head1 DESCRIPTION

Test accuracy of using direct repeats to 
classify rogue arrays using a naive bayesian classifier.

=head2 Output

2 columns: average, standard_deviation

=head1 EXAMPLES

=head2 Basic usage:

CLdb_classifyArraysByDR_validate.pl -da CLdb.sqlite

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

