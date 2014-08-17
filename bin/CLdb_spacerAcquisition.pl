#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_spacerAquisition.pl -- Identify instances of recent spacer acquisition

=head1 SYNOPSIS

CLdb_spacerAquisition.pl [flags] 

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

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -truncation  <int>

All loci ending within  '-truncation' bp away from scaffold end 
will be considered truncated, and only compare up to end of truncation. 
Negative values to turn option off. [500]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_spacerAquisition.pl

=head1 DESCRIPTION

Identify instances of recent spacer acquisition.

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

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



#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $cluster_cutoff = 1;
my $truncation = 500;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "cutoff=f" => \$cluster_cutoff, 		# clustering cutoff [1.00]
	   "query=s" => \$extra_query, 
	   "truncation=i" => \$truncation,
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


# finding newly acquired spacers
my $newAcqu_r = find_new_spacers($spacer_order_r, $trunc_loci_r, $tbl_r);




# pairwise comparison of arrays (parallel?)
## comparisons are directed: need lower & upper triangles
#my @locus_ids = keys %$tbl_r;
#foreach my $i (0..$#locus_ids){
#  foreach my $j (0..$#locus_ids){
# foreach my $locus_idi (keys %$tbl_r){
#   foreach my $locus_idj (keys %$tbl_r){
#     next if $locus_idi eq $locus_idj;

#     # unpack
#     #my $locus_i = $tbl_r->{$locus_ids[$i]};
#     #my $locus_j = $tbl_r->{$locus_ids[$j]};

#     # iterating over locus i 
#     my $first_match = 0;    # found first match yet?
#     foreach my $spacer_id (sort {$tbl_r->{$locus_idi}{$a}{Order} <=>
# 				   $tbl_r->{$locus_idi}{$b}{Order}}
# 			   keys %{$tbl_r->{$locus_idi}}){
#       # unpack
#       my $cluster_id = exists $tbl_r->{$locus_idi}{$spacer_id}{Cluster_ID} ?
# 	$tbl_r->{$locus_idi}{$spacer_id}{Cluster_ID} or
# 	  die "ERROR: cannot find 'Cluster_ID'\n";
      
#       # if clusterID in locus_j, where in order?
      


#     }

#   }
#}

# foreach array (array1):
#  foreach array (array2):
#   iterate over array1 positions:
#    match in array2? (1st match)
#    if match: is position for either > 0 [0-indexing]?
#    if no: next
#    if yes:
#     if no truncation check:
#      continue from 1st match to end of array
#     else:
#      check for end of array at end of scaffold (if yes, stop match check at end of array)
#     either:
#      how many matches for same-rank spacers?
#       number of matches must be >= to cutoff
#      how many consecutive matches from 1st match?
#       number of consecutive matches must be >= to cutoff
# output:
#  



#--- disconnect ---#
$dbh->disconnect();
exit;



#--- subroutines ---#

=head2 find_new_spacers

# comparisons of loci spacers
# using order hash: {locus_id : order : cluster_id}
## foreach locus_i
### foreach locus_j
#### scan locus_j starting at beginning of order locus_i
#### find spacer 1st match: flag
#### from 1st match: check 1-1 order matching; mark how many spacers prior to 1st match for locus_i
##### if no more spacers, see if truncation
      # if yes: stop (no more mismatches)
      # if no: mark all remaining spacer queries as mismatches
##### if no match, check next over with same query; mark as mismatch

=cut

sub align_arrays{
# aligning arrays based on ordering relative to spacer

}

sub find_new_spacers{
# comparing array alignments
}


sub find_new_spacers_OLD{
  my $spacer_order_r = shift or die "Provide spacer order\n";
  my $trunc_loci_r = shift or die "Provide trunc loci\n";
  my $tbl_r = shift or die "Provide output from get_array_elem_pos()\n";


  my %new_spacers;
  foreach my $locus_i (keys %$spacer_order_r){
    foreach my $locus_j (keys %$spacer_order_r){
      # no self-hits
      next if $locus_i eq $locus_j;

      # initialize
      $new_spacers{$locus_i}{$locus_j}{matches} = 0;
      $new_spacers{$locus_i}{$locus_j}{mismatches} = 0;
      $new_spacers{$locus_i}{$locus_j}{psbl_new_spacers} = 0;
      
      #-- first match --#
      ## iterating over locus_i spacers by order
      foreach my $order_i (sort{$a <=> $b} keys %{$spacer_order_r->{$locus_i}}){
	my $cluster_id_i = $spacer_order_r->{$locus_i}{$order_i};

	# looking for match in spacers
	foreach my $order_j (sort{$a <=> $b} keys %{$spacer_order_r->{$locus_j}}){
	  if ($spacer_order_r->{$locus_j}{$order_j} eq $cluster_id_i){
	    $new_spacers{$locus_i}{$locus_j}{first_match} = $order_i;
	    $new_spacers{$locus_i}{$locus_j}{first_match_order_i} = $order_i;
	    $new_spacers{$locus_i}{$locus_j}{first_match_order_j} = $order_j;
	    $new_spacers{$lcous_i}{$locus_j}{psbl_new_spacers} = $order_i - 1;
	  }
	}
	# end if found 1st match
	last if exists $new_spacers{$locus_i}{$locus_j}{first_match};
      }
 	
      #-- scoring from first match --#
      next unless eixsts $new_spacers{$locus_i}{$locus_j}{first_match};
      # unpack
      #my $first_match_order_i = $new_spacers{$locus_i}{$locus_j}{first_match_order_i};
      #my $first_match_order_j = $new_spacers{$locus_i}{$locus_j}{first_match_order_j};

      # iterating from first match
      #foreach my $i ($first_match_order_i..(


      ## iterating over locus_i spacers by order
      my $array_len = scalar keys %{$spacer_order_r->{$locus_i}};   # length of array_i
      foreach my $order_i (sort{$a <=> $b} keys %{$spacer_order_r->{$locus_i}}){

	# skipping up to first match
	my $cluster_id_i = $spacer_order_r->{$locus_i}{$order_i};
	next if $cluster_id_i <= $new_spacers{$locus_i}{$locus_j}{first_match};
	
	# looking for match in spacers beyond first match
	while(
	my $order_shift = 0;   # checking next in order with same one if necessary
	for my $i ($order_i..$array_len){
	  my $j = $i - $order_shift
	  # setting cluster id for current place in array_i
	  my $cluster_id_i = $spacer_order_r->{$locus_i}{$i};

	  # no more spacers in array_j
	  if (not exists $spacer_order_r->{$locus_j}{$i}){

	  }
	  # match in array_j
	  elsif ($spacer_order_r->{$locus_j}{$i} == $cluster_id_i){
	    
	  }
	  # mismatch in array_j
	  elsif ($spacer_order_r->{$locus_j}{$i} != $cluster_id_i){
	    
	  }

	}

	#foreach my $order_j (sort{$a <=> $b} keys %{$spacer_order_r->{$locus_j}}){
	#  # match
	#  if ($spacer_order_r->{$locus_j}{$order_j} eq $cluster_id_i){
	#  }
	#}	
      }
    }
  }
}


=head2 get_spacer_order

Setting spacer order in array relative to leader

=head3 IN

$tbl_r -- hashref --table returned by get_array_elem_pos()
$leader_info_r -- hashref -- leader table with locus_id as key

=head4 OUT

#Edited $tbl_r. "Order" category added.
{locus_id : order : cluster_id}

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
    
    ## getting order (1-indexing)
    my $order_cnt = 0;
    foreach my $spacer_id (sort{$leaderDist{$a} <=> $leaderDist{$b}}
			   keys %leaderDist){
      $order_cnt++;
      my $cluster_id = exists $tbl_r->{$locus_id}{$spacer_id}{Cluster_ID} ?
	$tbl_r->{$locus_id}{$spacer_id}{Cluster_ID} : 
	  die "Cannot find 'Cluster_ID'\n";
	  
      #$tbl_r->{$locus_id}{$spacer_id}{Order} = $order_cnt;
      $spacer_order{$locus_id}{$order_cnt} = $cluster_id;
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
