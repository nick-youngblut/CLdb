#!/usr/bin/env perl

=pod

=head1 NAME

spacerAquisition.pl -- Identify instances of recent spacer acquisition

=head1 SYNOPSIS

spacerAquisition.pl [flags] 

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -cutoff  <float>

Compare spacers at the specified sequence identity cutoff. [1]

=item -truncation  <int>

All loci ending within  '-truncation' bp away from scaffold end 
will be considered truncated, and only compare up to end of truncation. 
Negative values to turn option off. [500]

=item -all  <bool>

Compare all CRISPRs regardless of whether they have differing
spacer content at the leader end of the array.
(useful for make spacer content score tables of all pairwise
comparisons of CRISPRs).

=item -out  <char>

Prefix for output files. [psblNewSpacers]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb_perldoc spacerAquisition.pl

=head1 DESCRIPTION

Identify instances of potential recent spacer acquisition.

For all CRISPR arrays selected by query, each will be
compared pairwise with all other selected loci.

The arrays will be aligned based on the cluster_ID of
each spacer. 

CRISPR arrays are oriented by the leader (and thus must
have an identified leader region). 

Mismatches and gaps in the start of the alignment
prior to any matches are considered potential newly acquired spacers
since the two loci diverged from a common ancestor.

The rest of the alignment is used to calculate
'precent identity' [aln_len - (mismatches + gaps)) / alignment_len],
matches, mismatches, gaps, etc.

=head2 Trucations

Truncations due to incomplete assemblies could lead to false negatives.
Artificial truncations can be assessed (distance of array from end of scaffold)
and if the array is truncated, the alignment will only extend to the end
of the truncated array. 

=head2 Output

=head3 *_summary.txt

Summary of CRISPR array alignments.

=head3 *_scores.txt

This shows how each position in the alignment was scored 
(eg., match). Scoring scheme:

=over

=item * 'm' = match (same cluster_ID)

=item * 'x' = mismatch

=item * 'g' = gap

=back

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
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
use CLdb::query qw/table_exists
		   n_entries
		   join_query_opts
		   get_array_elem_pos
		   get_array_start_end
		   get_leader_info/;
use CLdb::genbank qw/get_loci_scaf/;
use CLdb::utilities qw/file_exists 
		       connect2db
		       get_seq_file_path/;
use CLdb::seq qw/revcomp/;
use CLdb::ArrayAlignWeighter;



#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $compare_all);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $cluster_cutoff = 1;
my $truncation = 500;
my $out_prefix = "psblNewSpacers";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "cutoff=f" => \$cluster_cutoff, 	
	   "query=s" => \$extra_query, 
	   "truncation=i" => \$truncation,
	   "all" => \$compare_all,
	   "out=s" => \$out_prefix,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "spacers");
table_exists($dbh, "spacer_clusters");


# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


# getting arrays of interest from database #
## locus_id : spacerID : (spacer_start, cluster_id, fasta_file, scaffold)
my %opts = (
        refine_sql => join("\n", $join_sql, $extra_query),
        cutoff => $cluster_cutoff,
        );
my $tbl_r = get_array_elem_pos($dbh, \%opts);


# truncation check:
my $trunc_loci_r; 
if ($truncation >= 0){
  # getting path to fasta files
  my $fastaDir = get_seq_file_path($database_file,
				   -seqType => 'fasta');

  # get scaffold lengths for each locus
  my $scafLens = get_loci_scaf($dbh, $fastaDir, [keys %$tbl_r], 
			       -len => 1, 
			       -seqType => 'fasta');

  # getting loci start-end for each locus
  my $arraySE = get_array_start_end($dbh, [keys %$tbl_r]);

  # mark each locus as truncated or not
  $trunc_loci_r = get_truncated($scafLens, $arraySE, $truncation);
}


# orient spacers by leader
## getting leader info for each locus_id
my $leader_info_r = get_leader_info($dbh, [keys %$tbl_r]);

## setting spacer position relative to leader, and thus order
my $spacer_order_r = get_spacer_order($tbl_r, $leader_info_r);

# comparing CRISPR arrays
my $array_cmp_r = compare_arrays($spacer_order_r, $tbl_r, 
				 -trunc => $trunc_loci_r,
				 -all => $compare_all);

# applying cutoffs?


# writing out tables
write_position_scores($array_cmp_r, $out_prefix);
write_summary($array_cmp_r, $out_prefix);



#--- disconnect ---#
$dbh->disconnect();
exit;



#--- subroutines ---#

=head2 filter_percentID

filtering out all CRISPR comparisons < percent ID

=cut

sub filter_percentID{
  my $array_cmp_r = shift or die "Provide array_cmp struct\n";

  foreach my $locus_id_i (keys %$array_cmp_r){
    foreach my $locus_id_j (keys %{$array_cmp_r->{$locus_id_i}}){
      my $percentID = exists $array_cmp_r->{$locus_id_i}{$locus_id_j}{percent_id}
	or die "KeyError: percent_id\n";
      
    }
  }
  
}

=head2 write_summary

Writing summary table of CRISPR array comparisons.

=cut

sub write_summary{
  my $array_cmp_r = shift or die "Provide array_cmp struct\n";
  my $out_prefix = shift or die "Provide output file prefix\n";

  # output file
  my $outFile = $out_prefix . "_summary.txt";
  open OUT, ">$outFile" or die $!;

  # header 
  my @fields = qw/aln_len_total first_match concord_matches matches mismatches 
		  gaps percent_id truncation/;
  print OUT join("\t", qw/locus_i locus_j/, @fields, 
		 qw/truncation_i truncation_j
		    psbl_new_spacers_i psbl_new_spacers_j/), "\n";

  # making 2d-array for sorting
  my @arr;
  foreach my $locus_id_i (keys %$array_cmp_r){
    foreach my $locus_id_j (keys %{$array_cmp_r->{$locus_id_i}}){
      exists $array_cmp_r->{$locus_id_i}{$locus_id_j}{percent_id}
	or die "KeyError: percent_id\n";
      push @arr, [$locus_id_i, $locus_id_j, 
		  $array_cmp_r->{$locus_id_i}{$locus_id_j}{percent_id}];     
    }
  }
  
  # writing after sorting by percent id
  foreach my $row (sort{$b->[2] <=> $a->[2]} @arr){    
    # check for keys
    map{ exists $array_cmp_r->{$row->[0]}{$row->[1]}{$_}
	   or die "KeyError: $_\n" } @fields;
    map{ exists $array_cmp_r->{$row->[0]}{$row->[1]}{$_}{truncation} 
	   or die "KeyError: $_ - truncation\n" } qw/i j/;

    # getting values
    my @vals = map{$array_cmp_r->{$row->[0]}{$row->[1]}{$_}} @fields;

    # write
    print OUT join("\t", $row->[0], $row->[1], 
		   @vals,
		   $array_cmp_r->{$row->[0]}{$row->[1]}{i}{truncation},
		   $array_cmp_r->{$row->[0]}{$row->[1]}{j}{truncation},
		   $array_cmp_r->{$row->[0]}{$row->[1]}{i}{psbl_new_spacers},
		   $array_cmp_r->{$row->[0]}{$row->[1]}{j}{psbl_new_spacers}
		  ), "\n";
    
  }
  close OUT or die $!;

  # status
  print STDERR "File written: $outFile\n";
}


=head2 write_position_scores

Writing a table of position scores.
Sorting by percent identity scores.

=cut

sub write_position_scores{
  my $array_cmp_r = shift or die "Provide array_cmp struct\n";
  my $out_prefix = shift or die "Provide output file prefix\n";
 
  # output files
  my $outFileW = $out_prefix . "_scoresW.txt";
  open OUTW, ">$outFileW" or die $!;
  my $outFileL = $out_prefix . "_scoresL.txt";
  open OUTL, ">$outFileL" or die $!;
  

  # headers
  print OUTW join("\t", qw/locus_i locus_j aln_percent_ID scores/), "\n";
  print OUTL join("\t", qw/locus_i locus_j aln_percent_ID abs_aln_pos
			   rel_aln_pos score_char score_val/), "\n";

  # making 2d-array for sorting
  my @arr;
  foreach my $locus_id_i (keys %$array_cmp_r){
    foreach my $locus_id_j (keys %{$array_cmp_r->{$locus_id_i}}){
      exists $array_cmp_r->{$locus_id_i}{$locus_id_j}{percent_id}
	or die "KeyError: percent_id\n";
      push @arr, [$locus_id_i, $locus_id_j, 
		  $array_cmp_r->{$locus_id_i}{$locus_id_j}{percent_id}];     
    }
  }
  
  # writing after sorting by percent id
  foreach my $row (sort{$b->[2] <=> $a->[2]} @arr){
    # position scores 
    my $pos_scores = exists $array_cmp_r->{$row->[0]}{$row->[1]}{position_score} ?
      $array_cmp_r->{$row->[0]}{$row->[1]}{position_score} :
	die "KeyError: position_score\n";
    # percent ID
    my $percentID = $array_cmp_r->{$row->[0]}{$row->[1]}{percent_id};


    # 'wide' table
    print OUTW join("\t", 
		    $row->[0], 
		    $row->[1], 
		    $percentID,
		    join("", @$pos_scores)
		   ), "\n";

    # 'long' table
    my $pos_scores_num = _scores2numeric($pos_scores);
    my $score_len = @$pos_scores;
    
    for my $i (1..$score_len){
      print OUTL join("\t", 
		      $row->[0], 
		      $row->[1], 
		      $percentID,
		      $i,
		      $i / $score_len,
		      $pos_scores->[$i-1],
		      $pos_scores_num->[$i-1]
		     ), "\n";      
    }
    
  }

  close OUTW or die $!;
  close OUTL or die $!;

  # status
  map{ print STDERR "File written: $_\n" } ($outFileW, $outFileL);
}


sub _scores2numeric{
  # scoring from character to numeric
  my $scores_r = shift or die $!;;

  # index: character => numeric
  my %index = ( 'm' => 1,
		'x' => 0,
		'g' => 'NA' );

  # iterating over alignment scores
  # converting
  my @num_scores =  map{  
    exists $index{$_} ? 
      $index{$_} :
	die "ERROR: illegal character in alignment scores -> '$_' \n";      
  } @$scores_r;

  return \@num_scores;
}



=head2 find_new_spacers

# comparisons of loci spacers
# using order hash: {locus_id : order : cluster_id}
## foreach locus_i
### foreach locus_j
#### alignment of CRISPR i,j arrays based on clusterIDs
#### determine end-gaps at start of alignment (psbl new spacers)
#### If either array a psbl truncation, trim end of alignment to shortest array
#### For rest of alignment: calc percent identity of alignment (aln_len - gaps / aln_len)
#### Save info in hash

=head3 IN


=head3 OUT

=cut

sub compare_arrays{
  my $spacer_order_r = shift or die "Provide spacer orderings\n";
  my $tbl_r = shift or die "Provide hashref output from get_array_elem_pos()\n";
  my %opts = @_;

  # unpack opts
  my $trunc_loci_r = $opts{-trunc};

  # just need lower triangle of pairwise
  my %array_cmp;
  my @locus_ids = keys %$spacer_order_r;

  # status
  my $nloci = scalar @locus_ids;
  printf STDERR "Number of loci to pairwise align: %i\n", $nloci;
  printf STDERR "Number of pairwise comparisons: %i\n", $nloci * ($nloci - 1) / 2;

  ## interating
  for my $i (0..$#locus_ids){
    # status
    printf STDERR "Processing array alignments to locus_id: %s\n", $locus_ids[$i];

    for my $j (0..$#locus_ids){
      next if $i >= $j;   # lower triangle

      # unpack
      my $locus_id_i = $locus_ids[$i];
      my $locus_id_j = $locus_ids[$j];
      my $locus_i = $spacer_order_r->{$locus_id_i}{cluster_ids};
      my $locus_j = $spacer_order_r->{$locus_id_j}{cluster_ids};

      # alignment of cluster_id arrays
      my $aln = ArrayAlignWeighter->new(left => $locus_i, right => $locus_j);

      # assessing alignment (eg., mismatches, matches, etc)
      my $aln_stats_r = _aln_stats($aln);
      
      # determining whether any psbl new spacers found
      ## must be mismatches in alignment to start
      next if $aln_stats_r->{first_match} <= 1 and not $opts{-all};

      # take into account truncations if exists
      ## if truncations will set aln_len_total to last alignment position of shortest truncated array
      ## matches, mismatches, and gaps changed to reflect 
      _apply_trunc($aln_stats_r, $locus_i, $locus_j, $trunc_loci_r)
	if defined $trunc_loci_r;

      # scoring and saving 
      ## psbl_new_spacers
      ## percent_id = (aln_len_total - (mismatches + gaps)) / aln_len_total
      _aln_score($aln_stats_r);

      # saving
      $array_cmp{$locus_id_i}{$locus_id_j} = $aln_stats_r;
    }
  }
  
  return \%array_cmp;
}


sub _aln_score{
  my $aln_stats_r = shift or die "Provide aln_stats\n";

  # sanity check
  map{ exists $aln_stats_r->{$_} or die "KeyError: $_\n" }
    qw/position_score first_match aln_len_total/;

  # unpack 
  my $pos_score_r = $aln_stats_r->{position_score};
  my $first_match = $aln_stats_r->{first_match};
  my $aln_len_total = $aln_stats_r->{aln_len_total};
  
  # scoreing
  ## initialize
  $aln_stats_r->{matches} = 0;
  $aln_stats_r->{mismatches} = 0;
  $aln_stats_r->{gaps} = 0;
  ## iterating over position score
  for my $i ($first_match..$aln_len_total){
    $i--;  # 0-index
    if($pos_score_r->[$i] eq 'm'){ $aln_stats_r->{matches}++; }
    elsif($pos_score_r->[$i] eq 'x'){ $aln_stats_r->{mismatches}++; }
    elsif($pos_score_r->[$i] eq 'g'){ $aln_stats_r->{gaps}++; }
    else{ die "Logic error\n"; }   
  }

  # calculating percent identity
  my $aln_len = $aln_len_total - $first_match + 1;
  $aln_stats_r->{percent_id} = (( $aln_len - 
				($aln_stats_r->{mismatches} + $aln_stats_r->{gaps})
			       ) / $aln_len) * 100;

    
  # getting number of concordinant matches beyond first (+ first)
  $aln_stats_r->{concord_matches} = 0;
  for my $i ($first_match..$aln_len_total){
    $i--;
    last if $pos_score_r->[$i] =~ /^[xg]$/;

    if ($pos_score_r->[$i] eq 'm'){ $aln_stats_r->{concord_matches}++; }
    else{ die "Logic error\n"; }
  }
}


sub _apply_trunc{
  my $aln_stats_r = shift or die "Provide aln_stats\n";
  my $locus_i = shift or die "Provide locus i id\n";
  my $locus_j = shift or die "Provide locus j id\n";
  my $trunc_loci_r = shift or die "Provide trunc_loci_r\n";

  # sanity check
  map{ exists $trunc_loci_r->{$_} or return 0 }
    ($locus_i, $locus_j);   # no truncation info for loci, then return 0

  # is either a truncation?
  if( $trunc_loci_r->{$locus_i} or $trunc_loci_r->{$locus_j} ){
    $aln_stats_r->{truncation} = 1;  # marked as truncation comparison

    # getting last spacer position of both
    my $last_sp_i = $aln_stats_r->{i}{last_spacer};
    my $last_sp_j = $aln_stats_r->{j}{last_spacer};
    
    # setting new last (truncated) position of alignment
    my $new_last = $aln_stats_r->{aln_len_total};
    if( $trunc_loci_r->{$locus_i} and $trunc_loci_r->{$locus_j} ){
      $aln_stats_r->{i}{truncation} = 1;
      $aln_stats_r->{j}{truncation} = 1;

      if( $last_sp_i < $last_sp_j){ $new_last = $last_sp_i; }
      else{ $new_last = $last_sp_j; }
    }
    elsif( $trunc_loci_r->{$locus_i} and not $trunc_loci_r->{$locus_j} ){
      $aln_stats_r->{i}{truncation} = 1;

      $new_last = $last_sp_i;
    }
    elsif( not $trunc_loci_r->{$locus_i} and $trunc_loci_r->{$locus_j} ){
      $aln_stats_r->{j}{truncation} = 1;

      $new_last = $last_sp_j;
    }
    else{ die "Logic error\n"; }
  
    # change aln_len_total to truncate just to 
    $aln_stats_r->{aln_len_total} = $new_last;

  }
  # no truncation
  else{
    return 1;
  }
}

sub _aln_stats{
  # processing alignment
  #### determine end-gaps at start of alignment (psbl new spacers)
  #### If either array a psbl truncation, trim end of alignment to shortest array
  #### For rest of alignment: calc percent identity of alignment 
  # positioning: 1-index
  my $aln = shift or die "Provide Array::Align object\n";
  my %opts = @_;
  
  my %ret = (
	     aln_len_total => 0,
	     first_match => 0,
	     last_match => 0,
	     position_score => [],
	     truncation => 0,
	     i => { psbl_new_spacers => 0,
		    last_spacer => 0,
		    len => 0,
		    truncation => 0},
	     j => { psbl_new_spacers => 0, 
		    last_spacer => 0,
		    len => 0,
		    truncation => 0}
	     );


  # interate over alignment
  for my $pair ($aln->pairwise){  # each column of the alignment
    # alignment running count in position
    $ret{aln_len_total}++;

    # match
    if (defined $pair->[0] and defined $pair->[1]
	and $pair->[0] == $pair->[1]){
      push @{$ret{position_score}}, 'm';  # m = match
      
      # first match?
      $ret{first_match} = $ret{aln_len_total}
	if not $ret{first_match};
		  
      # last match found
      $ret{last_match} = $ret{aln_len_total};      
    }
    # mismatch
    elsif(defined $pair->[0] and defined $pair->[1]
	  and $pair->[0] != $pair->[1]) {
      push @{$ret{position_score}}, 'x';  # x = mismatch
      
      # number of gaps for array at start of alignment
      $ret{i}{ngap_start}++ if not defined $pair->[0]
	and not $ret{first_match};	
      $ret{j}{ngap_start}++ if not defined $pair->[1]
	and not $ret{first_match};      
    }
    # gap
    elsif(not defined $pair->[0] or not defined $pair->[1]){
      push @{$ret{position_score}}, 'g';  # g = gap
    }
    else{
      die "Logic error\n";
    }	

    # stats for each array character (if present)
    if (defined $pair->[0]){  # charcter, no gap
      $ret{i}{len}++ if defined $pair->[0];
      $ret{i}{last_spacer} = $ret{aln_len_total};
      $ret{i}{psbl_new_spacers}++ unless $ret{first_match};
    }
    if (defined $pair->[1]){  # character, no gap
      $ret{j}{len}++ if defined $pair->[1];
      $ret{j}{last_spacer} = $ret{aln_len_total};
      $ret{j}{psbl_new_spacers}++ unless $ret{first_match};
    }

  }
  
  # sanity checks
  die "ERROR: position_score len != aln_len_total\n"
    unless scalar @{$ret{position_score}} == $ret{aln_len_total};

  #print Dumper %ret; exit;
  return \%ret;
}


sub _print_alignment{
  # simple printing of array alignment
  my $aln = shift or die "Provide Array::Align object\n";

  for my $pair ($aln->pairwise) {
    map{ $pair->[$_] = '' unless defined $pair->[$_]} 0..1;
    printf "%s - %s\n", @$pair;
  }
}



=head2 get_spacer_order

Setting spacer order in array relative to leader

=head3 IN

$tbl_r -- hashref --table returned by get_array_elem_pos()
$leader_info_r -- hashref -- leader table with locus_id as key

=head4 OUT

#Edited $tbl_r. "Order" category added.
#{locus_id : order : cluster_id}
{locus_id : [order-1, order-2, ..., order-n]

=cut

sub get_spacer_order{
  my $tbl_r = shift or die "Provide tbl produced by get_array_elem_pos()\n";
  my $leader_info_r = shift or die "Provide leader info as hashref\n";

  my %spacer_order;
  foreach my $locus_id (keys %$tbl_r){
    # check for leader info
    next unless exists $leader_info_r->{$locus_id};
    # unpack
    map{exists $leader_info_r->{$locus_id}{$_} or die "KeyError: $_\n"} qw/Leader_Start Leader_End/;
    my $leader_start = $leader_info_r->{$locus_id}{Leader_Start};
    my $leader_end = $leader_info_r->{$locus_id}{Leader_End};

    # setting leader start-end to + strand
    ($leader_start, $leader_end) = ($leader_end, $leader_start)
      if $leader_start > $leader_end;

    
    # setting order foreach spacer
    ## getting dist to leader
    my %leaderDist;
    foreach my $spacer_id (keys %{$tbl_r->{$locus_id}}){

      # unpack
      map{ exists $tbl_r->{$locus_id}{$spacer_id}{$_} or die "KeyError: $_\n" } 
	qw/Spacer_Start Spacer_End Cluster_ID/;
      my $spacer_start = $tbl_r->{$locus_id}{$spacer_id}{Spacer_Start};
      my $spacer_end = $tbl_r->{$locus_id}{$spacer_id}{Spacer_End};

      # setting spacer start-end to + strand
      ($spacer_start, $spacer_end) = ($spacer_end, $spacer_start)
	if $spacer_start > $spacer_end;

      # setting distance 
      ## leader_pos < spacer_pos
      if($leader_end <= $spacer_start){
	# dist from leader_end to spacer start
	$leaderDist{$spacer_id} = $spacer_start - $leader_end;
      }
      else{
	# dist from spacer_end to leader_start
	$leaderDist{$spacer_id} = $leader_start - $spacer_end;
      }
    }
    
    ## getting order: persisting as array
    foreach my $spacer_id (sort{$leaderDist{$a} <=> $leaderDist{$b}}
			   keys %leaderDist){
      my $cluster_id = exists $tbl_r->{$locus_id}{$spacer_id}{Cluster_ID} ?
	$tbl_r->{$locus_id}{$spacer_id}{Cluster_ID} : 
	  die "Cannot find 'Cluster_ID'\n";
	  
      push @{$spacer_order{$locus_id}{spacer_ids}}, $spacer_id;
      push @{$spacer_order{$locus_id}{cluster_ids}}, $cluster_id;
      #push @{$spacer_order{$locus_id}{spacer_start}}, $tbl_r->{$locus_id}{$spacer_id}{Spacer_Start};      
    }
  }

  #print Dumper %spacer_order; exit;
  return \%spacer_order;
}


=head2 get_truncated

Marking potentially truncated loci.

=head3 IN

$tbl_r -- hashref --table returned by get_array_elem_pos()
$scafLens -- hashref -- locus_id : '-len' : scaffold_len
$truncation -- int -- number bp that locus can be within to end of scaffold ot be considered truncated

=head4 OUT

locus_id : truncation

* truncation -- '1' = truncated

=cut

sub get_truncated{
  my $scafLens = shift or die "Provide hashref of scaffold lengths\n";
  my $arraySE = shift or die "Provide hash of array start-ends\n";
  my $truncation = shift or die "Provide truncation distinct (int)\n";

  my %trunc;
  foreach my $locus_id (keys %$scafLens){
    # sanity check
    exists $scafLens->{$locus_id}{-len} or die $!;
    exists $arraySE->{$locus_id} or die $!;

    # unpack
    my $scafLen = $scafLens->{$locus_id}{-len};
    my $array_start = $arraySE->{$locus_id}{array_start};
    my $array_end = $arraySE->{$locus_id}{array_end};
    
    # flip start-end to + strand if needed
    ($array_start, $array_end) = ($array_end, $array_start)
      if $array_start > $array_end;

    # determing whether truncated
    if($array_start < $truncation){
      $trunc{$locus_id} = 1;  
    }
    elsif($scafLen - $array_end < $truncation){
      $trunc{$locus_id} = 1;
    }
    else{
      $trunc{$locus_id} = 0;
    }
  }

  return \%trunc;
}
