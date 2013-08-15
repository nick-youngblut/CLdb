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


my ($verbose, $database_file, $rogue_bool);
my (@taxon_id, @taxon_name, @locus_id);
my $extra_query = "";
my $kmer_len = 8;
my $score_cutoff = 0.95;	
GetOptions(
	   "database=s" => \$database_file,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "locus_id=s{,}" => \@locus_id,
	   "query=s" => \$extra_query, 
	   "kmer=i" => \$kmer_len,
	   "cutoff=f" => \$score_cutoff,		# bayes score cutoff
	   "rogue" => \$rogue_bool, 			# just rounge arrays? [TRUE]
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

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");
map{ s/[A-Za-z]//g } @locus_id;
$join_sql .= join_query_opts(\@locus_id, "locus_id");

# getting  #
my $arrays_class_r = get_classified_arrays($dbh, $kmer_len);

# getting unclassified arrays #
my $arrays_unclass_r = get_unclassified_arrays($dbh, $join_sql, $extra_query, $kmer_len, $rogue_bool);

# training classifier #
my $nb = train_classifier($arrays_class_r);

# classify #
my $res_r = classify_arrays($nb, $arrays_unclass_r);

# filter results #
$res_r = filter_bayes_scores($res_r, $score_cutoff);

# updating loci table #
update_db($dbh, $res_r);


# disconnect #
$dbh->disconnect();
exit;


### Subroutines 
sub update_db{
# updating scaffold info #
	my ($dbh, $res_r) = @_;
	
	my $cmd = "Update Loci SET subtype = ? where locus_id = ?";
	my $sql = $dbh->prepare($cmd);
	
	my $update_cnt = 0;
	foreach my $locus_id (keys %$res_r){
		
		$sql->execute( ($res_r->{$locus_id}, $locus_id) );
					
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", $locus_id, $res_r->{$locus_id}), "\n";
			}
		else{ $update_cnt++; }
		}
	$dbh->commit;	
	
	print STDERR "...Number of newly classified arrays:\t$update_cnt\n" unless $verbose;
	print STDERR "...Loci table updated!\n" unless $verbose;
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
				print STDERR " WARNING! Cannot classify Locus: '$locus_id'!\tTop score only $top_score\n";
				}
			else{
				$passed{$locus_id} = $label;		# label = subtype
				}
			last;
			}
		}
	
		#print Dumper %passed; exit;
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

		#print Dumper %res; exit;
	return \%res;
	}

sub train_classifier{
# training classifier with classified arrays #
	my ($arrays_r) = @_;
	
	my $nb = Algorithm::NaiveBayes->new;

	foreach my $subtype (keys %$arrays_r){
		#print Dumper $subtype;
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
	my ($dbh, $join_sql, $extra_query, $kmer_len, $rogue_bool) = @_;
	
	# just rogues?
	my $rogue_sql = "";
	unless ($rogue_bool){
		$rogue_sql = "AND (loci.operon_start IS NULL OR loci.operon_start = \"\")";
		print STDERR "...Just trying to classify rogue arrays (change with '-rogue')\n"
			unless $verbose;
		} 
	
	# make query #
	my $query = "SELECT 
directrepeats.Repeat_sequence,
loci.locus_id
FROM directrepeats, loci 
WHERE directrepeats.locus_id = loci.locus_id 
AND (loci.subtype IS NULL
OR loci.subtype = \"\")
$rogue_sql
$join_sql
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
		unless $verbose;
	
		#print Dumper %arrays; exit;
	return \%arrays;
	}

sub get_classified_arrays{
# using arrays that were classified and have 'intact' operons #
	my ($dbh, $kmer_len) = @_;
	
	# make query #
	my $query = "SELECT 
directrepeats.Repeat_sequence,
loci.subtype,
loci.locus_id
FROM directrepeats, loci 
WHERE directrepeats.locus_id = loci.locus_id 
AND loci.operon_status = 'intact'
AND loci.subtype IS NOT NULL
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
			unless $verbose;
	
		#print Dumper %arrays; exit;
	return \%arrays;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}



__END__

=pod

=head1 NAME

CLdb_classifyArraysByDR.pl -- classify array subtype using direct repeats

=head1 SYNOPSIS

CLdb_classifyArraysByDR.pl [flags]

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -kmer

Kmer length using for naive bayes classifier. [8]

=item -cutoff

Bayes score cutoff (0-1) for accepting classification (>=). [0.95]

=item -rogue

Just classify 'rogue' arrays (no associated operon)? [TRUE]

=item -locus_id

Refine query arrays to specific a locus_id(s) (>1 argument allowed).

=item -taxon_id

Refine query arrays to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query arrays to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine the query array selection

=item -v 	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_classifyArraysByDR.pl

=head1 DESCRIPTION

Classify 'rogue' arrays that lack operons needed for classification. 

=head2 Naive Bayesian Classifier Method

=over

=item *	All classified arrays with 'intact' operons used for the training dataset.

=item *	Direct repeats sequences (training & query datasets) are parsed into all possible kmers
(similar to the RDPClassifier).

=item * The trained bayes classifier gives a score for each possible subtype
to each query (unclassified) array bases on the kmer composition of all 
of its direct repeats. Score range: 0-1 (0 = not likely; 1 = very likely)

=item *	The subtype label with the top score is used if it is >= '-cutoff' (Defaut: 0.95)

=item *	If the top score is < '-cutoff', the array will remain unclassified

=back

=head1 EXAMPLES

=head2 Classifying all rogue arrays

CLdb_classifyArraysByDR.pl -d CLdb.sqlite 

=head2 Classifying 1 rogue array: 'cli103'

CLdb_classifyArraysByDR.pl -d CLdb.sqlite -l cli103

=head2 Classifying all unclassified arrays

CLdb_classifyArraysByDR.pl -d CLdb.sqlite -r 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

