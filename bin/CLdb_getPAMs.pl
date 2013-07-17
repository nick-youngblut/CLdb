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
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name); 				# blast subject query refinement
my $extra_query = "";						# query refinement
my $extend = 20;						# length to extend off of each side
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "x=i" => \$extend,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);


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
$join_sql .= join_query_opts_or(\@staxon_id, \@staxon_name);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sql, $extra_query);

# getting PAMs & spacer sequences #
get_PAMs($dbh, $blast_hits_r, $extend);


# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub get_PAMs{
# getting spacer blast hit sequence & adjacent regions #
	my ($dbh, $blast_hits_r, $extend) = @_;
	
	my $query = "SELECT substr(scaffold_sequence,?,?) FROM blast_subject
WHERE scaffold_name = ?
AND (taxon_id == ?
 OR taxon_name == ?)";
  	$query =~ s/[\n\r]+/ /g; 

	print STDERR "$query\n" if $verbose;

	my $sql = $dbh->prepare($query);
	
	# getting sequence for each blast hit queried #
	foreach my $hit (@$blast_hits_r){
		# start - end #
		my $hit_start = $$hit[4];
		my $hit_end = $$hit[5];
		
		## flipping to positive strand ##
		if($hit_start > $hit_end){
			($hit_start, $hit_end) = flip_se($hit_start, $hit_end);
			}
		my $hit_len = $hit_end - $hit_start;
		
		## adding -extend ##
		$hit_start = $hit_start - $extend;
		
		## accounting for loss of extension on 5' side #
		my $fiveP_loss = 0;
		if($hit_start < 0){		 		#need to account for loss of extension
			$fiveP_loss = abs $hit_start;	# amount negative
			$hit_start = 0;
			}
		
		# undef of taxon_id & taxon_name #
		undef $$hit[1] if $$hit[1] eq "";
		undef $$hit[2] if $$hit[2] eq "";
		
		# querying; getting substr of scaffold #
		$sql->execute( ($hit_start, $hit_len + $extend * 2, $$hit[3], $$hit[1], $$hit[2]) );	
		my $ret = $sql->fetchall_arrayref();
		print STDERR " WARNING: a blast subject sequence could not be found for hit: ",
				join(", ", @$hit), "\n" unless @$ret;
		
		# converting extended region to lower case #
		foreach my $row (@$ret){
			$$row[0] = spacer2lc($$row[0], $extend, $hit_len, $fiveP_loss, $$hit[4], $$hit[5]);
			}
		
		foreach my $row (@$ret){
			$$hit[1] = "NULL" unless $$hit[1];
			$$hit[2] = "NULL" unless $$hit[2];		
			print join("___", ">Group$$hit[0]", @$hit[1..5]), "\n";
			print join("\t", @$row), "\n";
			}
		}  	
	}

sub spacer2lc{
# converting extension beyond spacer hit to lower case #
	my ($seq, $extend, $hit_len, $fiveP_loss, $hit_s, $hit_e) = @_;

	# all to upper case sequence #
	$seq =~ tr/a-z/A-Z/;

	# spliting sequence #
	## 5' end extension = extension - loss_cause_by_negative #
	my $fiveP_seq = substr $seq, 0, $extend - $fiveP_loss;
	my $hit_seq = substr $seq, $extend - $fiveP_loss, $hit_len;
	my $threeP_seq = substr $seq, $extend - $fiveP_loss + $hit_len;
	
		#print Dumper $seq, $fiveP_seq, $hit_seq, $threeP_seq; exit;
	
	# extensions to lower case #
	$fiveP_seq =~ tr/A-Z/a-z/;
	$threeP_seq =~ tr/A-Z/a-z/;
	
	# rev-comp if needed #
	if($hit_s > $hit_e){
		$fiveP_seq = revcomp($fiveP_seq);
		$hit_seq = revcomp($hit_seq);
		$threeP_seq = revcomp($threeP_seq);
		
		return join("", $threeP_seq, $hit_seq, $fiveP_seq);
		}
	else{
		return join("", $fiveP_seq, $hit_seq, $threeP_seq);
		}
	}

sub revcomp{
	# reverse complements DNA #
	my $seq = shift;
	
	$seq = reverse $seq;
	$seq =~ tr/ACGTNBVDHKMRYSWacgtnbvdhkmrysw\.-/TGCANVBHDMKYRSWtgcanvbhdmkyrsw\.-/;
	return $seq;
	}

sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT blast_hits.Spacer_group, blast_hits.Taxon_ID, 
blast_hits.Taxon_name, blast_hits.Subject, 
blast_hits.sstart, blast_hits.send 
from Loci, Spacers, Blast_hits
WHERE Loci.locus_id == Spacers.locus_id
AND blast_hits.spacer_group == spacers.spacer_group
AND blast_hits.CRISPR_array == 'yes'
$join_sql";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	return $ret;
	}

sub join_query_opts_or{
	my ($staxon_id_r, $staxon_name_r) = @_;
	
	return "" unless @$staxon_id_r || @$staxon_name_r;
	
	# adding quotes #
	map{ s/"*(.+)"*/"$1"/ } @$staxon_id_r;
	map{ s/"*(.+)"*/"$1"/ } @$staxon_name_r;	
	
	if(@$staxon_id_r && @$staxon_name_r){
		return join("", " AND (blast_hits.taxon_id IN (", join(", ", @$staxon_id_r),
						") OR blast_hits.taxon_name IN (", join(", ", @$staxon_name_r),
						"))");
		}
	elsif(@$staxon_id_r){
		return join("", " AND blast_hits.staxon_id IN (", join(", ", @$staxon_id_r), ")");
		}
	elsif(@$staxon_name_r){
		return join("", " AND blast_hits.staxon_name IN (", join(", ", @$staxon_name_r), ")");	
		}
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}

sub path_by_database{
	my ($database_file) = @_;
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "genbank");
	}




__END__

=pod

=head1 NAME

CLdb_getPAMs.pl -- getting PAMs from spacer blast hits

=head1 SYNOPSIS

CLdb_getPAMs.pl [flags] > PAM-spacer.fan

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 options

=over

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine the query.

=item -staxon_id

BLAST subject taxon_id (>1 argument allowed)

=item -staxon_name

BLAST subject taxon_name (>1 argument allowed)

=item -x

Extension beyond spacer blast to check for PAMs (bp). [20]

=item -v	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_getPAMs.pl

=head1 DESCRIPTION

Get sequences of spacer blast hits plus adjacent nucleotides
to look for PAMs. Only spacer blast hits that were not found
to fall into a CRISPR array are queried.

=head2 Output

The spacer blast hit sequence is capitalized, while the 
extended region around the hit is lower case. 

Spacer-PAM sequences are named as following:
"Spacer_group"___"Subject_Taxon_ID"___"Subject_Taxon_name"___"Subject_scaffold_name"___"blast_hit_start"___"blast_hit_end"

"Subject" = subject in the spacer blast.

"NULL" will be used if "Subject_Taxon_ID" or "Subject_Taxon_name"
is not present in CLdb.

=head2 Requirements for analysis

=over

=item * Spacer blasting must be done prior

=item * Subject sequences of spacer blast hits must be in CLdb

=back

=head1 EXAMPLES

=head2 Spacer-PAMs for all hits

CLdb_getPAMs.pl -d CLdb.sqlite

=head2 Spacer-PAMs for just subtype I-B 

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" > spacer_PAM_I-B.fna

=head2 Spacer-PAMs for just subtype I-B, smaller extendion

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" -x 10 > spacer_PAM_I-B.fna

=head2 Spacer-PAMs to just 2 subject taxon_names

CLdb_getPAMs.pl -d CLdb.sqlite -staxon_name ecoli salmonela > spacer_PAM_E-S.fna

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

