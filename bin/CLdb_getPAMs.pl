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
use List::Util qw/max/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $by_group, $align_pams);
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name, @sacc); 		# blast subject query refinement
my $extra_query = "";						# query refinement
my $extend = 10;							# length to extend off of each side
my $len_cutoff = 1;							# cutoff for length of blast hit relative to query length
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "saccession=s{,}" => \@sacc,
	   "group" => \$by_group,					# by spacer group, not individual spacer
	   "x=i" => \$extend,
	   "align" => \$align_pams, 
	   "length=f" => \$len_cutoff,
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
$join_sql .= join_query_opts_or(\@staxon_id, \@staxon_name, \@sacc);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sql, $extra_query);

# writing out PAMs #
get_PAMs($dbh, $blast_hits_r, $extend, $len_cutoff);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines	
sub get_PAMs{
# getting spacer blast hit sequence & adjacent regions #
	my ($dbh, $blast_hits_r, $extend, $len_cutoff) = @_;
	
	# finding max proto length if aligning 5' & 3' #
	my $max_len;
	unless($align_pams){
		my %proto_lens;
		foreach my $entry (@$blast_hits_r){
			$proto_lens{ abs($$entry[8] - $$entry[7]) } = 1;
			}
		$max_len = max keys %proto_lens;
		}
		
	# writing PAMs #
	foreach my $entry (@$blast_hits_r){
		# skipping unless full length blast hit #
		next unless (($$entry[5] - $$entry[4] + 1) / $$entry[6]) >= $len_cutoff;
		
		# making blast hit all caps, extension = lower case #
		my $frag = pam_caps($entry, $max_len);
			#print $frag, "\n";
		
		# writing fasta of PAMs #
		map{$_ = "NULL" unless $_} @$entry[0..3];
		if($by_group){
			print join("\n",
				  join("__", ">Group." . $$entry[0],
					   @$entry[1..3], 
					   $$entry[16], $$entry[7], $$entry[8]),
				  $frag), "\n";
			}
		else{
			print join("\n",
				  join("__", ">Group." . $$entry[0],
				  		@$entry[12..14], @$entry[1..3],
					   $$entry[16],$$entry[7], $$entry[8]),
				  $frag), "\n";			
			}
		}
	}
	
sub pam_caps{
	my ($entry, $max_len) = @_;
	
	# getting substrings of 5', proto, & 3' #
	my $sstart = $$entry[7];
	my $send = $$entry[8];
	my $frag = $$entry[9];
	my $xstart = $$entry[10];
	my $xend = $$entry[11];
	
	if($sstart > $send){		# flippig start & end
		# sanity check #
		die " ERROR: sstart > send but xstart < xend!\n"
			if $xstart < $xend;
		($sstart,$send) = flip_se($sstart,$send);
		($xstart,$xend) = flip_se($xstart,$xend);
		}
	
	## applying -x ##
	my $substr_start =  $sstart - $xstart - $extend;
	$substr_start = 0 if $substr_start < 0;
	$xstart = $sstart - $extend unless $sstart - $extend < $xstart;
	$xend = $send + $extend unless $send + $extend > $xend;

	## substr of frag ##	
	my $fiveP = substr($frag, $substr_start, $sstart - $xstart);
	my $proto = substr($frag, $substr_start + $sstart - $xstart, $send - $sstart);
	my $threeP = substr($frag, $substr_start + $sstart - $xstart + $send - $sstart, $xend - $send);

	# extensions to lower case; proto to upper #
	$fiveP =~ tr/A-Z/a-z/;
	$threeP =~ tr/A-Z/a-z/;
	$proto =~ tr/a-z/A-Z/;
	
	# adding gaps in middle of proto if aligning 5' & 3' #
	my $proto_len = length $proto;
	if($max_len && $proto_len < $max_len){
		my @gaps = "-" x ($max_len - $proto_len + 1);
		my $mid = sprintf("%.0f", $proto_len / 2);
		my $proto1 = substr($proto, 0, $mid);
		my $proto2 = substr($proto, $mid + 1, $proto_len - $mid + 1);
		$proto = join("", $proto1, @gaps, $proto2);
		}
	
	return join("", $fiveP, $proto, $threeP);
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
		
		return [$threeP_seq, $hit_seq, $fiveP_seq];
		}
	else{
		return [$fiveP_seq, $hit_seq, $threeP_seq];
		}
	}

sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}

sub get_blast_hits{
# just selecting by group
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
c.Group_id, c.S_taxon_name, c.S_taxon_id, c.S_accession,
c.qstart, c.qend, c.qlen,
c.sstart, c.send, 
c.frag, c.xstart, c.xend,
b.spacer_id,
a.taxon_name, a.taxon_id, a.subtype,
c.sseqid
from Loci a, Spacers b, Blast_hits c
WHERE a.locus_id == b.locus_id
AND b.spacer_group == c.group_id
AND c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sql";
	
	$query .= " GROUP BY c.Group_ID" if $by_group;

	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
		
		#print Dumper @$ret; exit;
	return $ret;
	}

sub get_blast_hits_by_group{
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
Group_ID,S_taxon_name,S_taxon_ID,S_accession,
qstart, qend, qlen,
sstart, send,
frag, xstart, xend
FROM loci a, Blast_hits c
WHERE c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sql";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;
		#print Dumper $query; exit;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
		#print Dumper @$ret; exit;
	return $ret;
	}

sub join_query_opts_or{
	my ($staxon_id_r, $staxon_name_r, $sacc_r) = @_;
	
	return "" unless @$staxon_id_r || @$staxon_name_r || @$sacc_r;
	
	# adding quotes #
	map{ s/"*(.+)"*/"$1"/ } @$staxon_id_r;
	map{ s/"*(.+)"*/"$1"/ } @$staxon_name_r;	
	map{ s/"*(.+)"*/"$1"/ } @$sacc_r;	
	
	if(@$staxon_id_r && @$staxon_name_r){
		return join("", " AND (c.S_taxon_id IN (", join(", ", @$staxon_id_r),
						") OR c.S_taxon_name IN (", join(", ", @$staxon_name_r),
						"))");
		}
	elsif(@$staxon_id_r){
		return join("", " AND c.S_taxon_id IN (", join(", ", @$staxon_id_r), ")");
		}
	elsif(@$staxon_name_r){
		return join("", " AND c.S_taxon_name IN (", join(", ", @$staxon_name_r), ")");	
		}
	elsif(@$sacc_r){
		return join("", " AND c.S_accession IN (", join(", ", @$sacc_r), ")");	
		}
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND a.$cat IN (", join(", ", @$vals_r), ")");
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

CLdb_getPAMs.pl [flags] > PAM-spacer.fna

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

=item -staxon_id  <char>

BLAST subject taxon_id (>1 argument allowed)

=item -staxon_name  <char>

BLAST subject taxon_name (>1 argument allowed)

=item -saccession  <char>

BLAST accession number (>1 argument allowed)

=item -align  <bool>

Align the 5' & 3' end of protospacers by adding gaps in
the middle of the protospacer. [TRUE]

=item -group  <bool>

Get protospacers by spacer group (i.e. de-replicated spacers). [FALSE]

=item -length  <int>

Length cutoff for blast hit (>=; fraction of spacer length). [1]

=item -x  <int>

Extension beyond spacer blast to check for PAMs (bp). [10]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>	

This help message

=back

=head2 For more information:

perldoc CLdb_getPAMs.pl

=head1 DESCRIPTION

Get sequences of spacer blast hits plus adjacent nucleotides
to look for PAMs.

Only spacer blast hits that were not found
to fall into a CRISPR array are queried.

By default, only full-length spacer blast hits
are used (change with '-l').

=head2 Output

The spacer blast hit sequence is capitalized, while the 
extended region around the hit is lower case. 


=head3 Default sequence naming 

"Spacer_group"__"spacer_ID"__"Query_taxon_name"__"Query_taxon_ID"__
"Subject_Taxon_Name"__"Subject_Taxon_ID"__"Subject_scaffold_name"__"blast_hit_start"__"blast_hit_end"

=head3 Naming for '-g' flag (by spacer group)

"Spacer_group"__"spacer_ID"__"Subject_Taxon_Name"__"Subject_Taxon_ID"
__"Subject_scaffold_name"__"blast_hit_start"__"blast_hit_end"

"Subject" = subject in the spacer blast.

"NULL" will be used if any values are NULL in CLdb.


=head2 Requirements for analysis

=over

=item * Spacer blasting with CLdb_spacerBlastGenome.pl or CLdb_spacerBlastDB.pl must be done prior


=back

=head1 EXAMPLES

=head2 Protospacers for all hits

CLdb_getPAMs.pl -d CLdb.sqlite > all_proto.fasta

=head2 Protospacers for all hits (hits grouped by spacer group)

CLdb_getPAMs.pl -d CLdb.sqlite -g > all_proto_grouped.fasta

=head2 Protospacers for just spacers in subtype I-B 

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" > I-B_proto.fasta

=head2 Protospacers for just subtype I-Bl smaller extention

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" -x 5 > I-B_proto.fasta

=head2 Protospacers from just 2 subject taxon_names

CLdb_getPAMs.pl -d CLdb.sqlite -staxon_name ecoli salmonela > proto_coli_sal.fasta

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

