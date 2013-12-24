#!/usr/bin/env perl

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

Get protospacers by spacer group (ie. de-replicated spacers). [TRUE]

=item -region  <int>

Write just the 3' (-r 3) or 5' (-r 5) region adjacent to protospacer. 

=item -qlength  <int>

Length cutoff for blast hit (>=; fraction of spacer length). [0.66]

=item -slength  <int>

Spacer max total length (<=; bp). [ ]

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

Get sequences of protospacers (found via 
blast hits) plus adjacent nucleotides
to look for PAMs.

Only spacer blast hits that were not found
to fall into a CRISPR array are queried
(must run CLdb_spacerBlastDRFilter.pl 1st!).


=head2 Output

The spacer blast hit sequence is capitalized, while the 
extended region around the hit is lower case. 

=head3 Default sequence naming (by spacer group)

"Spacer_group"__"cli#"__"spacer_ID"__"Subject_Taxon_Name"__
"Subject_Taxon_ID"__"Subject_scaffold_name"__"protospacer_start"__"protospacer_end"__"strand"

=head3 Naming if '-g' flag (not by spacer group; by individual spacer)

"Spacer_group"__"Subject_Taxon_Name"__"Subject_Taxon_ID"
__"Subject_scaffold_name"__"protospacer_start"__"protospacer_end"__"strand"

"Subject" = subject in the spacer blast.

"NULL" will be used if any values are NULL in CLdb.

=head2 Requirements for analysis

=head3 Orientation of the written sequences 

The protospacer (& adjacent) regions as oriented so that 
the 3' end is on the left. The sequence is actually the complement of
the protospacer. 

=over

=item * Spacer blasting with CLdb_spacerBlastGenome.pl or CLdb_spacerBlastDB.pl must be done prior

=back

=head1 EXAMPLES

=head2 Protospacers for all hits (by spacer group)

CLdb_getPAMs.pl -d CLdb.sqlite > all_proto.fasta

=head2 Protospacers for all hits (each identical spacer in each spacer group)

CLdb_getPAMs.pl -d CLdb.sqlite -g > all_proto_each.fasta

=head2 Protospacers for just spacers in subtype I-B 

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" > I-B_proto.fasta

=head2 Protospacers for just subtype I-B smaller extention

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" -x 5 > I-B_proto.fasta

=head2 Protospacers from just 2 subject taxon_names

CLdb_getPAMs.pl -d CLdb.sqlite -staxon_name ecoli salmonela > proto_coli_sal.fasta

=head2 Just the 3' region adjacent to the protospacer

CLdb_getPAMs.pl -d CLdb.sqlite -subtype I-B" -r 3 > I-B_proto_3prime.fasta

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
use Bio::SeqIO;
use Set::IntervalTree;
use List::Util qw/max/;

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

my ($verbose, $database_file, $align_pams, $by_group);
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name, @sacc); 		# blast subject query refinement
my $extra_query = "";						# query refinement
my $extend = 10;							# length to extend off of each side
my $len_cutoff = 0.66;						# cutoff for length of blast hit relative to query length
my $proto_region = 0;
my $slen_cutoff = 50;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "saccession=s{,}" => \@sacc,
	   "group" => \$by_group,					# TRUE
	   "x=i" => \$extend,
	   "align" => \$align_pams, 				# TRUE
	   "qlength=f" => \$len_cutoff,				# full-length blast hits?
	   "slength=i" => \$slen_cutoff,			# spacer length must be <= 
	   "region=i" => \$proto_region,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
$database_file = File::Spec->rel2abs($database_file);


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);


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

print Dumper $blast_hits_r; exit;

# writing out PAMs #
get_PAMs($dbh, $blast_hits_r, $extend, $len_cutoff);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines	
sub get_PAMs{
# writing out protospacers & extensions #
	my ($dbh, $blast_hits_r, $extend, $len_cutoff) = @_;

    # finding max proto length for aligning 5' & 3' #
	my $max_len;
	my $max_3px_len = 0;
	my $max_5px_len = 0;
	unless($align_pams){
		my %proto_lens;
    	foreach my $entry (@$blast_hits_r){
        	$proto_lens{"proto"}{ abs($$entry[9] - $$entry[8] + 1) } = 1;		# qseq_full_end - qseq_full_start
            if($$entry[12]){ $proto_lens{"3p"}{ length $$entry[12] } = 1; }
            else{ $proto_lens{"3p"}{ 0 } = 1; }
            if($$entry[15]){ $proto_lens{"5p"}{ length $$entry[15] } = 1; }
            else{ $proto_lens{"5p"}{ 0 } = 1; }
            }
        $max_len = max keys %{$proto_lens{"proto"}};
        $max_3px_len = max keys %{$proto_lens{"3p"}};
       	$max_5px_len = max keys %{$proto_lens{"5p"}}; 
        }
    
    
    # writing PAMs #
    foreach my $entry (@$blast_hits_r){
    	# filling in partial proto3-5x #
    	my ($proto3px, $proto5px) = fill_partial_protox(
    					$$entry[12], $$entry[15], 
    					$max_3px_len, $max_5px_len);
    	
    	# uncaps for 3' & 5' protospacer extensions #
    	$proto3px =~ tr/[A-Z]/[a-z]/;
    	$proto5px =~ tr/[A-Z]/[a-z]/;    	
    	
    	# trimming adjacent regions if specified #
    	$proto3px = reverse(substr(reverse($proto3px), 0, $extend))
    		if length $proto3px > $extend;
    	$proto5px = substr($proto3px, 0, $extend)
    		if length $proto5px > $extend;


    	# adding gaps to sseq_full & capitalizating #
    	my $proto_seq = $$entry[7];
    	$proto_seq = add_gaps($proto_seq, $max_len) unless $align_pams; 		# adding gaps, unless not wanted
    	$proto_seq =~ tr/[a-z]/[A-Z]/;  
    	
    	# full proto seq (& extensions); oriented by proto (5'-3') #
    	my $full_seq;
    	if($proto_region == 3){ $full_seq = $proto3px; }
	    elsif($proto_region == 5){ $full_seq = $proto5px; }
		else{ $full_seq = join("", $proto3px, $proto_seq, $proto5px); }
    	
    	# writing fasta of PAMs #
		map{$_ = "NULL" unless $_} @$entry[0..3];

		if($by_group){		# all loci hitting the pam
			print join("\n",
				  join("__", ">Group" . $$entry[0],
				  		"cli.$$entry[22]", $$entry[18],
					   @$entry[1..4], 
					   $$entry[10], $$entry[11], $$entry[5]),
				  $full_seq), "\n";		
			}
		else{				# just spacer groups
			print join("\n",
				  join("__", ">Group." . $$entry[0],
					   @$entry[1..4], 
					   $$entry[10], $$entry[11], $$entry[5]),
				  $full_seq), "\n";
			}
    	}
	}
	
sub fill_partial_protox{
	my ($threep, $fivep, $max_3p, $max_5p) = @_;

	$threep = "" unless $threep;
	$fivep = "" unless $fivep;	
	my $ngap =  "-" x ($max_3p - length($threep));
	$threep = $ngap . $threep if $ngap;
	$ngap =  "-" x ($max_5p - length($fivep));
	$fivep = $fivep . $ngap if $ngap;
	
	return ($threep, $fivep);
	}

sub add_gaps{
	my ($proto_seq, $max_len) = @_;
	
	my $proto_len = length $proto_seq;
	if($max_len && $proto_len < $max_len){
		my @gaps = "-" x ($max_len - $proto_len + 1);
		my $mid = sprintf("%.0f", $proto_len / 2);
		my $proto1 = substr($proto_seq, 0, $mid);
		my $proto2 = substr($proto_seq, $mid + 1, $proto_len - $mid + 1);
		$proto_seq = join("", $proto1, @gaps, $proto2);
		}
		
	return $proto_seq;
	}

sub flip_proto{
# flipping proto & extension #
## protospacer stored as spacer match; need to revcomp ##
	my ($blast_hits_r) = @_;

	foreach my $entry (@$blast_hits_r){
		$$entry[7] = revcomp($$entry[7]);		# sseq_full (full proto)
		my $prime3 = revcomp($$entry[12]);
		my $prime5 = revcomp($$entry[15]);
		($$entry[12], $$entry[15]) = ($prime5, $prime3);
		#$$entry[15] = revcomp($$entry[]);		# 5' end of proto
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

sub get_blast_hits{
# just selecting by group
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
c.Group_id, 
c.S_taxon_name, 
c.S_taxon_id, 
c.S_accession,
c.sseqid,
c.strand,
c.qseq_full,
c.sseq_full,
c.qseq_full_start, 
c.qseq_full_end, 
c.sseq_full_start, 
c.sseq_full_end,
c.proto3px,
c.proto3px_start,
c.proto3px_end,
c.proto5px,
c.proto5px_start,
c.proto5px_end,
b.spacer_id,
loci.subtype,
loci.taxon_name, 
loci.taxon_id,
loci.locus_id,
c.blast_id
from Loci, Spacers b, Blast_hits c
WHERE loci.locus_id == b.locus_id
AND b.spacer_group == c.group_id
AND c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sql";
	
	$query =~ s/[\n\r]+/ /g;
	$query .= " AND c.qlen <= $slen_cutoff" if $slen_cutoff;
	$query = join(" ", $query, $extra_query);
	$query .= " GROUP BY c.blast_id" unless $by_group;
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries. Did you filter the spacer hits using CLdb_spacerBlastDRFilter.pl?\n"
		unless $$ret[0];
	
	# checking for subject sequence #
	foreach my $entry (@$ret){
		die " ERROR: no qseq_full found for $$entry[23]\n"
			unless $$entry[6];
		die " ERROR: no sseq_full found for $$entry[23]\n"
			unless $$entry[7];
		}

	
		#print Dumper  @$ret; exit;
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

sub join_query_opts_OLD{
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

sub get_blast_hits_by_group_OLD{
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
Group_ID,
S_taxon_name,
S_taxon_ID,
S_accession,
qstart, 
qend, 
qlen,
sstart, 
send,
frag, 
xstart, 
xend,
strand
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
	
sub get_PAMs_OLD{
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
		# skipping unless protospacer is full length of spacer + extension #
		# next unless length(proto_seq) == qlen + 5' ext + 3' ext
		print Dumper $entry;
		print Dumper length($$entry[16]), $$entry[6], abs($$entry[7] - $$entry[10]), 
											+ abs($$entry[11] - $$entry[8]); exit;
		next unless length($$entry[16]) == $$entry[6]
											+ abs($$entry[7] - $$entry[10]) 
											+ abs($$entry[11] - $$entry[8]);
		
		print Dumper $entry;
		
		# making blast hit all caps, extension = lower case #
		my $frag = pam_caps($entry, $max_len);
			#print $frag, "\n";
		
		# writing fasta of PAMs #
		map{$_ = "NULL" unless $_} @$entry[0..3];
		}
	}
	
sub pam_caps_OLD{
	my ($entry, $max_len) = @_;
	
	# getting substrings of 5', proto, & 3' #
	my $sstart = $$entry[7];
	my $send = $$entry[8];
	my $frag = $$entry[9];
	my $xstart = $$entry[10];
	my $xend = $$entry[11];
	
	#if($sstart > $send){		# flippig start & end
	#	# sanity check #
	#	die " ERROR: sstart > send but xstart < xend!\n"
	#		if $xstart < $xend;
	#	($sstart,$send) = flip_se($sstart,$send);
	#	($xstart,$xend) = flip_se($xstart,$xend);
	#	}
	
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

sub spacer2lc_OLD{
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

