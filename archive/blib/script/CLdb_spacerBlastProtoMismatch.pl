#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Set::IntervalTree;
use DBI;
use List::Util qw/sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $by_group);
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name, @sacc); 		# blast subject query refinement
my $extra_query = "";
my $len_cutoff = 1;
my $proto_region = 3;	
my ($pam_regex, @seed_region);
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "saccession=s{,}" => \@sacc,
	   "sacc=s{,}" => \@sacc,
	   "length=f" => \$len_cutoff,			# blast hit must be full length of query
	   "query=s" => \$extra_query,
	   "group" => \$by_group,				# just per spacer group
	   "pam=s" => \$pam_regex,				# pam regex for screening
	   "region=s" => \$proto_region,		# proto region to look for pam; 3' by default
	   "seed=i{2,2}" => \@seed_region, 		# seed region (eg. 1,8 or 8-1)
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
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");
$join_sql .= join_query_opts_or(\@staxon_id, \@staxon_name, \@sacc);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sql, $extra_query);

# flipping proto & extension #
## protospacer stored as spacer match; need to revcomp ##
flip_proto($blast_hits_r);

# screen by user provided pam #
#screen_by_pam($blast_hits_r, $pam_regex, $proto_region);

# make mismatch table #
make_mismatch_table($blast_hits_r);


### Subroutines
sub make_mismatch_table{
# making a mismatch table; determining number of mismatches between spacer & protospacer #
## mismatches aligned by 5' end of protospacer ##
## writing out 
	my ($blast_hits_r) = @_;
	
	# getting max spacer-proto length #
	my $max_len = 0;
	foreach my $entry (@$blast_hits_r){
		$max_len = length $$entry[14] if length $$entry[14] > $max_len;
		}
	
	# counting number of mismatches #
	my @mismatch_sum;
	my $total_len = 0;
	foreach my $entry (@$blast_hits_r){
		calc_mismatches($$entry[14], $$entry[15], $$entry[0], $max_len, \@mismatch_sum);
		$total_len += length $$entry[14];		# length of spacer-protospacer (full)
		}
	
	# adding zeros if needed #
	foreach (@mismatch_sum){ $_ = 0 unless $_; } 		# no mismatches = 0;
	
	# parsing by seed & non-seed #
	#my %mm_region; 		# mismatch by region (seed, non-seed)
	print STDERR "Total number of mismatches: ", sum(@mismatch_sum), "\n";
	print STDERR "Total seed region mismatches: ", sum(@mismatch_sum[0..7]), "\n";
	print STDERR "Total non-seed region mismatches: ", sum(@mismatch_sum[8..$#mismatch_sum]), "\n";	

	print STDERR "PercentID total: ", 
			(1 - sum(@mismatch_sum) / $total_len) * 100, "\n";
	print STDERR "PercentID seed: ", 
			(1 - sum(@mismatch_sum[0..7]) / (8 * scalar @$blast_hits_r)) * 100, "\n";
	print STDERR "PercentID non-seed: ", 
			(1 - sum(@mismatch_sum[8..$#mismatch_sum]) / ($total_len - (8 * scalar @$blast_hits_r))) * 100, "\n";	
	
	
	#print Dumper @mismatch_sum; 
	}

sub calc_mismatches{
# calculating mismatches for spacer-protospacers #
	my ($spacer_seq, $proto_seq, $blast_id, $max_len, $mismatch_sum_r) = @_;
	
	# upper case #
	$spacer_seq =~ tr/a-z/A-Z/;
    $proto_seq =~ tr/a-z/A-Z/;

    # sanity check #
    die " ERROR: query & protospacer seq are not the same length!\n"
        unless length $spacer_seq == length $proto_seq;

	my @spacer_seq = split //, $spacer_seq;
	my @proto_seq = split //, $proto_seq;

	# mismatch count; summing by position #
	my @mismatch;
	for my $i (0..$#spacer_seq){
		if ($spacer_seq[$i] ne $proto_seq[$i]){
			$mismatch[$i]++;
			$$mismatch_sum_r[$i]++;
			}
		}
	# adding zeros if neede #
	for my $i (0..($max_len -1)){ $mismatch[$i] = 0 unless $mismatch[$i]; }
		
	# writing out basic mismatch table #
	print join("\t", $blast_id, @mismatch), "\n";
	
	}

sub screen_by_pam{
# screening by pam regex #
	my ($blast_hits_r, $pam_regex, $proto_region);

	# status #
	print STDERR "...Number of protospacers selected: ", scalar @$blast_hits_r, "\n";
	
	# pam regex #
	$pam_regex = make_pam_regex($pam_regex);

	my @passed;
	foreach my $entry (@$blast_hits_r){
		# screening regions for pam #
		my $hit;
		if($proto_region != 5){			# screening 3' region
			$hit++ if $$entry[16] =~ /$pam_regex/;
			}
		elsif($proto_region != 3){		# screen 5' region
			my $proto5px_revcomp = revcomp($$entry[16]);
			$hit++ if $proto5px_revcomp =~ /$pam_regex/;
			}
			
		push(@passed, $entry) if $hit;
		}
	
	print STDERR "...Number of entries with specified PAM: ", scalar @passed, "\n";
	return \@passed;
	}

sub flip_proto{
# flipping proto & extension #
## protospacer stored as spacer match; need to revcomp ##
## also revcomp qseq_full for finding mismatches ##
	my ($blast_hits_r) = @_;

	foreach my $entry (@$blast_hits_r){
		# qseq_full, sseq_full, proto3px, proto5px
		map{ $_ = revcomp($_) } @$entry[14..17];
		}
	}

sub revcomp{
	# reverse complements DNA #
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/ACGTNBVDHKMRYSW\.-/TGCANVBHDMKYRSW\.-/;
	return $seq;
	}

sub make_pam_regex{
	my ($pam_regex) = @_;

	my %IUPAC = (
		"A" => "A",
		"C" => "C",
		"G" => "G",
		"T" => "T",
		"U" => "T",
		"R" => "AT",
		"Y" => "CT",
		"S" => "GC",
		"W" => "AT",
		"K" => "GT",
		"M" => "AC",
		"B" => "CGT",
		"D" => "AGT",
		"H" => "ACT",
		"V" => "ACG",
		"N" => "ATCG"
		);
	
	$pam_regex =~ tr/[a-z]/[A-Z/; 					# all caps
	my @pam_regex = split //, $pam_regex;
	map{$_ = "\[$IUPAC{$_}\]" if exists $IUPAC{$_}} @pam_regex;
	
	$pam_regex = join("", @pam_regex);
	$pam_regex =~ qr/$pam_regex/;
	
	return $pam_regex;
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sqls, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
c.blast_id,
c.group_id, 
a.locus_id,
a.taxon_name,
a.taxon_id,
b.spacer_id,
a.subtype,
c.S_Taxon_ID, 
c.S_Taxon_name,
c.S_accession,
c.sseqid,
c.qstart,
c.qend,
c.qlen,
c.qseq_full,
c.sseq_full,
c.proto3px,
c.proto5px
FROM Loci a, Spacers b, Blast_hits c
WHERE a.locus_id == b.locus_id
AND c.group_id == b.spacer_group
AND c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sqls";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	$query .= " GROUP BY c.blast_id" unless $by_group;
	
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




__END__

=pod

=head1 NAME

CLdb_spacerBlastProtoMismatch.pl -- Make table of spacer-protospacer mismatches 

=head1 SYNOPSIS

CLdb_spacerBlastProtoMismatch.pl [flags] > mismatch.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 options

=over

=item -bin  <int>

Number of bins for grouping mismatches along the spacer-protospacer alignment. [20]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -staxon_id  <char>

Refine query to specific a subject taxon_id(s) (>1 argument allowed).

=item -staxon_name  <char>

Refine query to specific a subject taxon_name(s) (>1 argument allowed).

=item -saccession  <char>

Refine query to specific a accession(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine the query.

=item -length  <float>

Length cutoff for blast hit (>=; fraction of spacer length). [1]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlastProtoMismatch.pl

=head1 DESCRIPTION

Make a table of protospacer mismatches.

Spacers are not necessarily the same length,
so mismatches are quantified by the number 
of mismatches in a certain percentage of the alignment
(e.g. 0%-5%, with 0% & 100% being the 5' & 3' of the alignment,
respectively).
The spacer-protospacer alignment is
binned into percentages of the alignment length.
Bins are rounded down, so a bin of 0%-5% may 
technically be alignment positions 1-3.5, 
but the program will use positions 1-3.
Bins are not inclusive by percentage (e.g. 
0%-4.99%, 5%-9.99%, 10%-15.99%, ...).

=head2 WARNING!

Spacer blasting and DR filtering must be done prior!

=head1 EXAMPLES

=head2 Mismatch table of all protospacers

CLdb_spacerBlastProtoMismatch.pl -d CLdb.sqlite  > all_proto_mismatch.txt

=head2 Mismatch table of all protospacers in subtype I-A

CLdb_spacerBlastProtoMismatch.pl -d CLdb.sqlite -sub I-A > I-A_proto_mismatch.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

