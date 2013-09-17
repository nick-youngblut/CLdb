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

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $by_group);
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name, @sacc); 		# blast subject query refinement
my $extra_query = "";
my $len_cutoff = 1;
my $bin = 20; 								# 20 bins by default
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
	   "bin=i" => \$bin,					# number of bins to parse sequences
	   "group" => \$by_group,				# just per spacer group
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

# binning mismatches #
bin_mismatches($blast_hits_r, $bin);


### Subroutines
sub bin_mismatches{
	my ($blast_hits_r, $bin) = @_;
	
	# header #
	if($by_group){
		print join("\t", qw/bin mismatch group_ID subtype s_taxon_id s_taxon_name s_accession/), "\n";
		}
	else{
		print join("\t", qw/bin mismatch group_ID q_taxon_name q_taxon_id subtype locus_id spacer_id s_taxon_id s_taxon_name s_accession/), "\n";
		}
	
	# body #
	foreach my $entry (@$blast_hits_r){
		# applying len cutoff #
		next if abs($$entry[15] - $$entry[14] + 1) / $$entry[16] < $len_cutoff;
		
		# checking for undefined values #
		map{$_ = "" unless $_} @$entry[(2..3,5..7,21)];
		
		# trimming frag to just proto #
		# frag = $$entry[17];
		my $xstart = $$entry[19];
		my $xend = $$entry[20];
		my $sstart = $$entry[9];
		my $send = $$entry[10];
		
		## flipping start-end if necessary ##
		($xstart, $xend) = flip_se($xstart, $xend) if $xstart > $xend;
		($sstart, $send) = flip_se($sstart, $send) if $sstart > $send;
		
		## getting substr (protospacer) ##
		my $proto;
		if($$entry[9] <= $$entry[10]){
			$proto = substr($$entry[17], $sstart - $xstart - 1, $send - $sstart + 1);
			}
		else{
			$proto = substr($$entry[17], $sstart - $xstart, $send - $sstart + 1);
			}
		
			#print Dumper$$entry[18], $proto;
		
		# binning mismatches #
		my $mismatch_r = bin_align($$entry[18], $proto, $bin);
		
		# writing mismatch table #
		foreach my $mis_bin (sort{$a<=>$b} keys %$mismatch_r){
			if($by_group){
				print join("\t",  $mis_bin, $mismatch_r->{$mis_bin}, @$entry[(0,21,5..7)]), "\n";
				}		
			else{
				print join("\t",  $mis_bin, $mismatch_r->{$mis_bin}, @$entry[(0,2,3,21,1,4..7)]), "\n";
				}
			}
		}
	}

sub bin_align{
# binning mismatches in align #
	my ($query, $proto, $bin) = @_;

	# upper case #
	$query =~ tr/a-z/A-Z/;
	$proto =~ tr/a-z/A-Z/;

	# sanity check #
	die " ERROR: query & protospacer seq are not the same length!\n"
		unless length $query == length $proto;

	
	my @query = split //, $query;
	my @proto = split //, $proto;
	my $bin_size = length($query) / $bin;
	
	my @bins;
	for (my $i=0; $i<=$#query; $i+=$bin_size){
		push @bins, $i; #sprintf("%.0f", $i);
		}
	
	my %mismatch;
	#for (my $i=0;$i<= $#bins -1; $i+=2){
	for my $i (0..($#bins - 1)){			# foreach bin #
		my $bin_start = sprintf("%.0f", $bins[$i]);
		my $bin_end = sprintf("%.0f", $bins[$i+1] - 0.001);
		my $mismatch = 0;
		for my $ii ($bin_start..$bin_end){
			$mismatch++ if $query[$ii] ne $proto[$ii] &&
				($query[$ii] ne "-" && $proto[$ii] ne "-");
				#print "$query[$ii] <-> $proto[$ii] ; $mismatch\n";
			}
		$mismatch{sprintf("%.3f", $i * $bin_size / length($proto) * 100)} = $mismatch;
		}
	
		#print Dumper %mismatch; exit;
	return \%mismatch;
	}

sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sqls, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
c.group_id, 
a.locus_id,
a.taxon_name,
a.taxon_id,
b.spacer_id,
c.S_Taxon_ID, 
c.S_Taxon_name,
c.S_accession,
c.sseqid,
c.sstart, 
c.send, 
c.evalue,
c.mismatch,
c.pident,
c.qstart,
c.qend,
c.qlen,
c.frag,
c.qseq,
c.xstart,
c.xend,
a.subtype
FROM Loci a, Spacers b, Blast_hits c
WHERE a.locus_id == b.locus_id
AND c.group_id == b.spacer_group
AND c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sqls";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	$query .= " GROUP BY c.group_id, c.S_taxon_id, c.S_taxon_name, c.S_accession, c.sseqid, c.sstart, c.send, a.subtype" if $by_group;
	
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

