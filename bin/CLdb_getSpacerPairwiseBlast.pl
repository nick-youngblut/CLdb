#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use List::Util qw/min/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name); 				# blast subject query refinement
my @overlap_frac;
my $extra_query = "";						# query refinement
my $percentID_cutoff = 90; 					# spacers must have >= percent ID
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "overlap=f{,}" => \@overlap_frac,
	   "percentID=f" => \$percentID_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);
@overlap_frac = (0,1) unless @overlap_frac;

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

# getting spacer info from spacers table #
get_spacer_info($dbh, $blast_hits_r, \@overlap_frac, $percentID_cutoff);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub get_spacer_info{
# getting info on spacers (genome-wide location) #
	my ($dbh, $blast_hits_r, $overlap_frac_r, $percentID_cutoff) = @_;
	
	# query #
	my $query = "SELECT * from spacers where Locus_id = ? AND Spacer_id = ?";
	my $sql1 = $dbh->prepare($query);
	my $sql2 = $dbh->prepare($query);
	
	foreach my $hit (@$blast_hits_r){
		# query #
		$sql1->bind_param(1,$$hit[0]);
		$sql1->bind_param(2,$$hit[1]);
		$sql1->execute();
		my $q_ret = $sql1->fetchall_arrayref();

		# subject #
		$sql2->bind_param(1,$$hit[2]);
		$sql2->bind_param(2,$$hit[3]);
		$sql2->execute();
		my $s_ret = $sql2->fetchall_arrayref();	

		# total spacer length (not just blast hit) #
		my $qlen = abs($$q_ret[0][3] - $$q_ret[0][2]);
		my $slen = abs($$s_ret[0][3] - $$s_ret[0][2]);
		
		# applying genome-wide start-stop to blast hits #
		@$hit[8..9] = change_position([@$hit[8..9]], $$q_ret[0], $qlen);
		@$hit[10..11] = change_position([@$hit[10..11]], $$s_ret[0], $slen);
		
		# filtering by overlap & 
		my $pass = 1;
		if(@$overlap_frac_r){			# filtering by min blast hit overlap (% of length of query or subject)
			my $q_overlap = abs($$hit[8] - $$hit[9]) / $qlen;
			my $s_overlap = abs($$hit[10] - $$hit[11]) / $slen;			
			my $min_overlap = min($q_overlap, $s_overlap);
			$pass = 0 unless $min_overlap >= $$overlap_frac_r[0] && 
				$min_overlap <= $$overlap_frac_r[1];

			push(@$hit, sprintf("%.2f", $q_overlap));
			push(@$hit, sprintf("%.2f", $s_overlap));
			}
		if($percentID_cutoff){
			$pass = 0 if $$hit[4] < $percentID_cutoff;
			}
		
		# writing #
		if($pass){
			$$hit[0] =~ s/^/cli./;
			$$hit[2] =~ s/^/cli./;
			print join("\t", 
					join("__", @$hit[0..1]),
					join("__", @$hit[2..3]),
					@$hit[4..$#$hit]), "\n";
			}
		}
	
	}
	
sub change_position{
# local to genome-wide blast hit #
	my ($pos_r, $ret, $slen) = @_;		# pos_r = hit, ret = spacer location in genome
	
	# flipping spacer start-end if needed #
	@$ret[2..3] = flip_se( 	@$ret[2..3] );
	
	# moving blast to genome wide #
	if($$pos_r[0] > $$pos_r[1]){		# if blast to other strand 
		$$pos_r[0] = $$ret[3] - ($slen - $$pos_r[0]); 			# gen_wide_end - spacer_len - blast_end
		$$pos_r[1] = $$ret[2] + $$pos_r[1] - 1;
		}
	else{
		$$pos_r[1] = $$ret[3] - ($slen - $$pos_r[1]); 			# gen_wide_end - spacer_len - blast_end
		$$pos_r[0] = $$ret[2] + $$pos_r[0] - 1;
		}

	return $$pos_r[0], $$pos_r[1];
	}
	
sub flip_se{
	my ($start, $end) = @_;
	if($start > $end){
		return $end, $start
		}
	else{
		return $start, $end;
		}
	}

sub get_blast_hits{
# getting pairwise spacer blast hits #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	## the 'OR' will allow for getting blasts across selected subtypes or taxa ##
	## the grouping prevents multiple copies of a blast hit entry ##
	my $query = "SELECT spacer_pairwise_blast.*
from spacer_pairwise_blast, loci
WHERE (Loci.locus_id == spacer_pairwise_blast.query_locus_id
OR Loci.locus_id == spacer_pairwise_blast.subject_locus_id)
$join_sql
GROUP BY spacer_pairwise_blast.query_locus_id, spacer_pairwise_blast.query_spacer_id, 
spacer_pairwise_blast.subject_locus_id, spacer_pairwise_blast.subject_spacer_id";
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

CLdb_getSpacerPairwiseBlast.pl -- getting pairwise spacer blast hits in CLdb

=head1 SYNOPSIS

CLdb_getSpacerPairwiseBlast.pl [flags] > spacer-spacer_blasts.txt

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

Extra sql to refine the query.

=item -overlap  <float>

Check for spacer overlap. See DESCRIPTION. [0 1]

=item -percentID  <fload>

PercentID cutoff for blast hits (>=). [90]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_getSpacerPairwiseBlast.pl

=head1 DESCRIPTION

Get spacer-spacer BLASTn-short hits among spacers
in CLdb. 

=head2 Overlap

Overlap = blast_hit_length / spacer_length

'-overlap' defines the cutoff; the min overlap for the query
or spacer is used. Two values are required: min_overlap &
max_overlap. For example: '-overlap 0 0.99' will select
all spacer blasts that only partially overlap.

=head2 Output

The output is in standard blast+ -outfmt 6 format, 
except 2 columns are appended: 
query_spacer_overlap, subject_spacer_overlap

=head1 EXAMPLES

=head2 All spacer-spacer blasts

CLdb_getSpacerPairwiseBlast.pl -d CLdb.sqlite

=head2 All spacer-spacer blasts for just subtype I-B

CLdb_getSpacerPairwiseBlast.pl -d CLdb.sqlite -subtype I-B

=head2 All spacer-spacer blasts that only partially overlap

CLdb_getSpacerPairwiseBlast.pl -d CLdb.sqlite -o 0 0.99

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

