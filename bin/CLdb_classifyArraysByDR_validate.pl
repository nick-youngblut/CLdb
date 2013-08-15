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

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $quiet);
my $kmer_len = 8;
my $score_cutoff = 0.99;
my $iter = 100;
my $frac_arrays = 0.1;
GetOptions(
	   "database=s" => \$database_file,
	   "kmer=i" => \$kmer_len,
	   "cutoff=f" => \$score_cutoff,		# bayes score cutoff
	   "iterations=i" => \$iter, 			# number of test iterations
	   "fraction=f" => \$frac_arrays, 		# fraction of arrays to use as query
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
	my $arrays_unclass_r = get_unclassified_arrays($dbh, $kmer_len);

	# training classifier #
	my $nb = train_classifier($arrays_class_r);

	# classifying #
	my $res_r = classify_arrays($nb, $arrays_unclass_r);

	

	# filter results #
	$res_r = filter_bayes_scores($res_r, $score_cutoff);

	# compare to true classifications #
	compare_classifications($true_class_r, $res_r, $i, \%tpfn);

	# disconnect #
	$dbh->disconnect();
	}

write_summary(\%tpfn);


### Subroutines 
sub write_summary{
	my ($tpfn_r) = @_;

	# initializing if zero count #
	$tpfn_r->{"total"}{"true_pos"} = 0 unless exists $tpfn_r->{"total"}{"true_pos"};
	$tpfn_r->{"total"}{"false_neg"} = 0 unless exists $tpfn_r->{"total"}{"false_neg"};
	
	# total #
	print STDERR "Total true positives: ", $tpfn_r->{"total"}{"true_pos"}, "\n";
	print STDERR "Total false negatives: ", $tpfn_r->{"total"}{"false_neg"}, "\n";
	
	# by iteration #
	foreach my $iter_num (keys %$tpfn_r){
		next if $iter_num eq "total";
		
		# intializing if zero count #
		$tpfn_r->{$iter_num}{"true_pos"} = 0 unless exists $tpfn_r->{$iter_num}{"true_pos"};
		$tpfn_r->{$iter_num}{"false_neg"} = 0 unless exists $tpfn_r->{$iter_num}{"false_neg"};
		
		# writing #
		print join("\t", $iter_num, 
				$tpfn_r->{$iter_num}{"true_pos"},
				$tpfn_r->{$iter_num}{"false_neg"}
				), "\n";
		}
	}

sub compare_classifications{
	my ($true_class_r, $res_r, $iter_num, $tpfn_r) = @_;
	
	foreach my $locus_id (keys %$true_class_r){
		die " ERROR: locus \"cli$locus_id\" not found in classifed values\n"
			unless exists $res_r->{$locus_id};
			
		if($true_class_r->{$locus_id} eq $res_r->{$locus_id}){
			$tpfn_r->{$iter_num}{"true_pos"}++;
			$tpfn_r->{"total"}{"true_pos"}++;
			}
		else{ 
			print STDERR "FALSE negative! Locus_id: $locus_id; True classification '", 
					$true_class_r->{$locus_id},
					"'; Classified as: '", $res_r->{$locus_id}, "'\n"
					if $verbose;
			$tpfn_r->{$iter_num}{"false_neg"}++;
			$tpfn_r->{"total"}{"false_neg"}++;
			}
		}
	
	}

sub filter_bayes_scores{
# filtering bayes scores from NB classification #
	my ($res_r, $score_cutoff) = @_;
	
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
directrepeats.Repeat_sequence,
loci_temp.locus_id
FROM directrepeats, loci_temp 
WHERE directrepeats.locus_id = loci_temp.locus_id 
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
	
	my %arrays;
	foreach my $row (@$ret){
		for (my $i = 0; $i<=(length $$row[0]) - $kmer_len; $i++){
			my $kmer = substr($$row[0], $i, $kmer_len);
			$arrays{$$row[1]}{$kmer}++;		# locus_id => kmer => count
			}
		}
	
	# status #
	print STDERR "...Number of arrays to classify:\t", scalar keys %arrays, "\n"
		unless $quiet;
	
		#print Dumper %arrays; exit;
	return \%arrays;
	}

sub get_classified_arrays{
# using arrays that were classified and have 'intact' operons #
	my ($dbh, $kmer_len) = @_;
	
	# make query #
	my $query = "SELECT 
directrepeats.Repeat_sequence,
loci_temp.subtype,
loci_temp.locus_id
FROM directrepeats, loci_temp 
WHERE directrepeats.locus_id = loci_temp.locus_id 
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
CRISPR_Array_Start	INTEGER,
CRISPR_Array_End	INTEGER,
Operon_Status	TEXT	NOT NULL,
CRISPR_Array_Status	TEXT	NOT NULL,
Genbank	TEXT	NOT NULL,
Array_File	TEXT,
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
AND (CRISPR_array_start IS NOT NULL
OR CRISPR_array_start = '')
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

=item -kmer

Kmer length using for naive bayes classifier. [8]

=item -cutoff

Bayes score cutoff (0-1) for accepting classification (>=). [0.99]

=item -iterations

Number of test iterations to perform. [100]

=item -fraction

Fraction of dataset to test. [0.1]

=item -v 	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_classifyArraysByDR_validate.pl

=head1 DESCRIPTION

Test accuracy of naive bayesian classifier for arrays using direct repeats.

=head1 EXAMPLES

=head2 


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

