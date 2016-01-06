#!/usr/bin/env perl

=pod

=head1 NAME

make_dna_segs -- making dna_segs table for plotting

=head1 SYNOPSIS

make_dna_segs [flags] > dna_segs.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -ITEP  <char>

Get gene cluster info from ITEP. 2 arguments required: (ITEP_sqlite_file, cluster_runID).

=item -cutoff  <float>

Spacer clustering cutoff for spacer coloring (0.8 - 1). [1]

=item -locus_id  <char>

Refine query to specific a locus_id(s) (>1 argument allowed).

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -features  <char>

>=1 plotting feature to include. [spacer DR gene leader]

=item -verbose  <bool> 

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb -- plotLoci --perldoc -- make_dna_segs

=head1 DESCRIPTION

Make a basic dna_segs object needed for plotting.

=head2 Coloring genes by gene cluster (via ITEP)

Provide and ITEP sqlite file name and cluster run ID
to add cluster info to the table (used for coloring)

=head1 EXAMPLES

=head2 Plotting all loci classified as subtype 'I-A'

CLdb -- plotLoci -- make_dna_segs -d CLdb.sqlite -sub I-A 

=head2 Gene cluster info from ITEP

CLdb -- plotLoci -- make_dna_segs -d CLdb.sqlite -sub I-A  -I DATABASE.sqlite all_I_2.0_c_0.4_m_maxbit

=head2 No broken loci

CLdb --sql -q "AND loci.operon_status != 'broken'" -- plotLoci -- make_dna_segs -da CLdb.sqlite -sub I-A 

=head2 Just gene features for I-A subtype loci

CLdb -- plotLoci -- make_dna_segs -da CLdb.sqlite -sub I-A -f 'gene'

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../../lib";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $sqlite, @ITEP_sqlite, @feats);
my (@locus_id, @subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $spacer_cutoff = 1;
GetOptions(
	   "database=s" => \$sqlite,
	   "ITEP=s{,}" => \@ITEP_sqlite,
	   "locus_id=s{,}" => \@locus_id,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "cutoff=f" =>  \$spacer_cutoff,
	   "features=s{,}" => \@feats,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# checking CLdb #
file_exists($sqlite, "database");

# features
@feats = qw/spacer DR gene leader/
  unless @feats;

# checking ITEP input #
if(@ITEP_sqlite){
  die " ERROR: '-ITEP' must be 2 arguments (database_file, runID)!\n"
    unless scalar @ITEP_sqlite >= 2;
  die " ERROR: cannot find ITEP database file!\n"
    unless -e $ITEP_sqlite[0];
}
else{
  print STDERR "WARNING: no ITEP db provided!" .
    " All genes will have col value of '1'!\n";
}


#--- MAIN ---#
# connect 2 CLdb #
my $dbh = connect2db($sqlite);

# if cluster cutoff < 1; check for spacer_hclust entries #
die "ERROR: no entries in spacer_cluster table!",
  " Cannnot use a spacer clustering cutoff < 1!\n"
  unless n_entries($dbh, "spacer_clusters");

# connect to ITEP #
my $dbh_ITEP;
if(@ITEP_sqlite){
  my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
  $dbh_ITEP = DBI->connect("dbi:SQLite:dbname=$ITEP_sqlite[0]", '','', \%attr) 
    or die " Can't connect to $ITEP_sqlite[0]!\n";	
}

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@locus_id, "locus_id");
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting subtype #
my $subtypes_r = get_subtypes($dbh, $join_sql, $extra_query);

# getting spacer, DR, leader, gene from CLdb #
my %dna_segs; 
if( grep /^spacers*$/i, @feats ){
  my $spacer_clusters_r = get_spacer_info($dbh, 
					  \%dna_segs,
					  $join_sql, 
					  $extra_query, 
					  $spacer_cutoff);
}
if( grep /^(DRs*|direct.+repeats*)/i, @feats ){
  get_DR_info($dbh, 
	      \%dna_segs, 
	      $join_sql, 
	      $extra_query);
}
if( grep /^genes*$/, @feats ){
  get_gene_info($dbh, 
		\%dna_segs, 
		$join_sql, 
		$extra_query);
}
if( grep /^leaders*$/i, @feats){
  get_leader_info($dbh, 
		  \%dna_segs, 
		  $join_sql, 
		  $extra_query);
}

# getting gene cluster info if ITEP provided #
if (@ITEP_sqlite){
  my $gene_cluster_r = get_gene_cluster_info($dbh_ITEP, 
					     $ITEP_sqlite[0], 
					     $ITEP_sqlite[1], 
					     \%dna_segs);
}
else{	# ITEP note provided: gene cluster = 1 for all
  gene_cluster_same(\%dna_segs);
}

# mutli-loci / multi-subtype problem
## determining if multiple loci/subtypes per taxon_name ##
my ($multi_loci, $multi_subtype)  = check_multi(\%dna_segs, $subtypes_r);

## adding loci & subtypes to DNA_segs names ##
my $dna_segs_r = edit_dna_segs_taxon_name(\%dna_segs, 
					  $subtypes_r, 
					  $multi_loci, 
					  $multi_subtype);

# writing table #
write_dna_segs(\%dna_segs, $subtypes_r);

# disconnect #
$dbh->disconnect();
$dbh_ITEP->disconnect() if $dbh_ITEP;
exit;


### Subroutines
sub flip_genes_on_neg_strand{
# flipping start-end of genes on neg strand for plotting reasons #
# all gene features flipped for genes on - strand
  my ($dna_segs_r, $strand_r) = @_;
  
  foreach my $taxon (keys %$dna_segs_r){
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon}}){
      # only genes
      next unless exists $dna_segs_r->{$taxon}{$locus_id}{'gene'};	 
      next unless exists $strand_r->{$locus_id} && $strand_r->{$locus_id} == -1;
      
      # flipping s-e for all genes on neg CAS strand #
      foreach my $peg (keys %{$dna_segs_r->{$taxon}{$locus_id}{'gene'}}){
	# flipping all gene features
	my $x_r = $dna_segs_r->{$taxon}{$locus_id}{'gene'}{$peg};
	($$x_r[0], $$x_r[1]) = ($$x_r[1], $$x_r[0]); 		
	#print Dumper "flipped: $taxon, $locus_id, $peg";
      }	
    }
  }
  #exit;
}

sub get_cas_strand{
# determing strand of cas genes based on CAS_start & CAS_end
  my ($dbh, $dna_segs_r) = @_;
  
  # getting all loci #
  my @loci;
  map{ push @loci, keys %{$dna_segs_r->{$_}} } keys %$dna_segs_r;
	
  # querying db for CAS start-end #
  my $cmd = "SELECT CAS_start, CAS_end FROM loci WHERE locus_id = ?";
  my $sth = $dbh->prepare($cmd);
  
  my %strand;
  foreach my $locus_id (@loci){
    $sth->bind_param(1, $locus_id);
    $sth->execute() or confess $dbh->err;
    my $ret =$sth->fetchrow_arrayref();	
    die " ERROR: no CAS_start, CAS_end for locus: '$locus_id'!\n"
      unless $$ret[0];
    
    if($$ret[0] <= $$ret[1]){	# start < end; + strand
      $strand{$locus_id} = 1;
    }
    else{	                # end > start; - strand
      $strand{$locus_id} = -1;
    }
  }
  
  #print Dumper %strand; exit;
  return \%strand;
}


sub write_dna_segs{
# writing dna_segs precursor table #
## required columns ##
# name
# start
# end
# strand
# col
# lty
# lwd
# pch
# cex
# gene_type

## extra columns: ##
### taxon_name
### locus_ID
### feat_ID 

  my ($dna_segs, $subtypes_r) = @_;
  
  # header #
  print join("\t", qw/name start end strand col lty 
		      lwd pch cex gene_type taxon_name 
		      locus_id feat feat_id subtype dna_segs_id/), "\n";
  
  # body #
  foreach my $taxon_name (keys %$dna_segs){
    foreach my $locus_id (keys %{$dna_segs->{$taxon_name}}){
      # spacers #
      foreach my $id (keys %{$dna_segs->{$taxon_name}{$locus_id}{"spacer"}}){
	my $row = $dna_segs->{$taxon_name}{$locus_id}{"spacer"}{$id};
	
	# start-end, color #
	my $start = $$row[0];
	my $end = $$row[1];
	my $col = $$row[2];
	my $dna_segs_id = $$row[$#$row];
	
	print join("\t",  
		   $id,
		   $start,
		   $end,
		   1,			# strand
		   $col,
		   1, 0.5, 8, 1,	# plot formatting
		   "blocks", 		# end of required columns
		   $taxon_name,
		   $locus_id,
		   "spacer", 		# feature
		   $id,
		   $subtypes_r->{$taxon_name}{$locus_id},
		   $dna_segs_id
		  ), "\n";				
      }
      # DRs #
      foreach my $id (keys %{$dna_segs->{$taxon_name}{$locus_id}{"DR"}}){
	my $row = $dna_segs->{$taxon_name}{$locus_id}{"DR"}{$id};
	
	# start-end, color #
	my $start = $$row[0];
	my $end = $$row[1];
	my $dna_segs_id = $$row[$#$row];				
	
	print join("\t",  
		   $id,
		   $start,
		   $end,
		   1, 		        # strand 
		   1,			# col
		   1, 0.1, 8, 1,	# plot formatting
		   "blocks", 		# end of required columns
		   $taxon_name,
		   $locus_id,
		   "directrepeat", 	# feature
		   $id,
		   $subtypes_r->{$taxon_name}{$locus_id},
		   $dna_segs_id
		  ), "\n";				
      }
      # leaders #
      foreach my $id (keys %{$dna_segs->{$taxon_name}{$locus_id}{"leader"}}){
	my $row = $dna_segs->{$taxon_name}{$locus_id}{"leader"}{$id};
	
	# start-end, color #
	my $start = $$row[0];
	my $end = $$row[1];
	my $dna_segs_id = $$row[$#$row];				
	
	print join("\t",  
		   $id,			
		   $start,
		   $end,
		   1, 		        # strand 
		   1,			# col
		   1, 0.1, 8, 1,	# plot formatting
		   "blocks", 		# end of required columns
		   $taxon_name,
		   $locus_id,
		   "leader", 		# feature
		   $id,
		   $subtypes_r->{$taxon_name}{$locus_id},
		   $dna_segs_id
		  ), "\n";				
      }			
      # genes #
      foreach my $id (keys %{$dna_segs->{$taxon_name}{$locus_id}{"gene"}}){
	my $row = $dna_segs->{$taxon_name}{$locus_id}{"gene"}{$id};
	
	# start-end, color #			
	my $start = $$row[0];
	my $end = $$row[1];
	# getting strand before possible flip of start-end
	# needed for plotting genes in correct orientation
	my $strand = get_strand($start, $end); 	
	($start, $end) = flip_se($start, $end);	# all start-end to + strand
	
	my $alias = $$row[2];
	my $dna_segs_id = $$row[$#$row];
	
	
	# color #
	my $col = 1;
	$col = $$row[3]
	  if $$row[3];
	
	print join("\t",  
		   $alias,
		   $start,
		   $end,
		   $strand,				# strand 
		   $col,
		   1, 1, 8, 1,				# plot formatting
		   "arrows", 				# end of required columns
		   $taxon_name,
		   $locus_id,
		   "gene", 				# feature
		   $id,
		   $subtypes_r->{$taxon_name}{$locus_id},
		   $dna_segs_id			
		  ), "\n";				
      }				
    }
  }
}

sub flip_se{
  my ($start, $end) = @_;
  if($start <= $end){ return $start, $end; }
  else{ return $end, $start; }
}

sub edit_dna_segs_taxon_name{
  my ($dna_segs_r, $header_r, $multi_loci, $multi_subtype) = @_;
  
  foreach my $taxon_name (keys %$dna_segs_r){
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}}){
      # editing taxon_name in row #
      my $dna_segs_id = $taxon_name;
      $dna_segs_id = join("__", $taxon_name, 
			  $locus_id)
	if $multi_loci;
      $dna_segs_id = join("__", $dna_segs_id, 
			  $dna_segs_r->{$taxon_name}{$locus_id}{"subtype"})
	if $multi_subtype;
      
      foreach my $cat (keys %{$dna_segs_r->{$taxon_name}{$locus_id}}){  
	next if $cat eq 'subtype';
	foreach my $id (keys %{$dna_segs_r->{$taxon_name}{$locus_id}{$cat}}){
	  push(@{$dna_segs_r->{$taxon_name}{$locus_id}{$cat}{$id}}, $dna_segs_id);
	}                               
      }
    }
  }
  
  #print Dumper $dna_segs_r; exit;
}

sub check_multi{
  # checking for multiple entries per taxon #
  my ($dna_segs_r, $subtypes_r) = @_;  

  my $multi_loci = 0;				# mutliple loci per taxon_name
  my $multi_subtype = 0;			# multiple subtypes total 
  my %subtype_sum; 
  foreach my $taxon_name (keys %$dna_segs_r){
    $multi_loci = 1 if scalar keys %{$dna_segs_r->{$taxon_name}} > 1;
    
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}} ){
      # sanity check #
      die " ERROR: cannot find subtype for $taxon_name -> $locus_id!\n"
	unless exists $subtypes_r->{$taxon_name}{$locus_id};
      
      $subtype_sum{$subtypes_r->{$taxon_name}{$locus_id}}++;
      
      $dna_segs_r->{$taxon_name}{$locus_id}{'subtype'} = 
	$subtypes_r->{$taxon_name}{$locus_id};
    }
  }
  $multi_subtype = 1 if scalar keys %subtype_sum > 1;
  
  # status #
  print STDERR "...Found multiple loci for 1 or more taxa.",
    " Adding loci_ids to names in dna_segs table!\n"
    if $multi_loci;
  print STDERR "...Found multiple subtypes.",  
    " Adding subtype to names in dna_segs table!\n"
      if $multi_subtype;
  
  return $multi_loci, $multi_subtype;
}

sub get_strand{
  # start > end? #
  my ($start, $end) = @_;
  if($start <= $end){ return 1; }		# + strand
  else{ return -1; }					# - strand
}

sub gene_cluster_same{
  # gene cluster = 1 for all if ITEP not provided #
  my ($dna_segs_r) = @_;
  
  foreach my $taxon_id (keys %$dna_segs_r){
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon_id}}){
      foreach my $gene_id (keys %{$dna_segs_r->{$taxon_id}{$locus_id}{"gene"}}){
	push @{$dna_segs_r->{$taxon_id}{$locus_id}{"gene"}{$gene_id}}, 1;
      }
    }
  }
}

sub get_gene_cluster_info{
  # getting gene cluster info from ITEP #
  my ($dbh_ITEP, $ITEP_file, $runID, $dna_segs_r) = @_;
  
  my $query = "
SELECT clusterid
FROM clusters
WHERE runid = ?
AND geneid = ?
";
  $query =~ s/\n|\r/ /g;
  
  my $sth = $dbh_ITEP->prepare($query);
  
  my %gene_clusters;
  foreach my $taxon_id (keys %$dna_segs_r){
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon_id}}){
      foreach my $gene_id (keys %{$dna_segs_r->{$taxon_id}{$locus_id}{"gene"}}){
	$sth->bind_param(1, $runID);
	$sth->bind_param(2, $gene_id);
	$sth->execute() or confess $dbh->err;
	my $ret =$sth->fetchrow_arrayref();	
	die " ERROR: no matching entries for ITEP query!\n",
	  "  Query: '$query'\n"
	    unless $$ret[0];
	
	push @{$dna_segs_r->{$taxon_id}{$locus_id}{"gene"}{$gene_id}}, $$ret[0];
	
	$gene_clusters{$$ret[0]} = 1;
      }
    }
  }
  
  #print Dumper %gene_clusters; exit;
  return \%gene_clusters;
}

sub get_gene_info{
  # getting direct repeat info from CLdb #
  my ($dbh, $dna_segs_r, $join_sql, $extra_query) = @_;
  
  my $query = "
SELECT
loci.taxon_name,
loci.locus_id, 
genes.gene_id, 
genes.gene_start, 
genes.gene_end,
genes.gene_alias
FROM Loci, Genes
WHERE Loci.locus_id = Genes.locus_id
AND Genes.In_CAS = 'yes'
$join_sql
$extra_query
";
  $query =~ s/\n|\r/ /g;
  
  # status #
  print STDERR "$query\n" if $verbose;
  
  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  die " ERROR: no matching entries for CAS gene query!\n",
    "Query: '$query'" unless $$ret[0];
  
  # loading hash #
  foreach my $row (@$ret){
    $dna_segs_r->{$$row[0]}{$$row[1]}{"gene"}{$$row[2]} = [@$row[3..$#$row]];
  }
  
  #print Dumper %$dna_segs_r; exit;
  #return \%gene_ids;
}	

sub get_leader_info{
# getting leader info from CLdb (optional) #
## just using locus-ID for leader_ID
  my ($dbh, $dna_segs_r, $join_sql, $extra_query) = @_;
  
  my $query = "
SELECT 
loci.taxon_name,
loci.locus_id, 
leaders.locus_id,
leaders.leader_start, 
leaders.leader_end
FROM Loci, leaders
WHERE Loci.locus_id = leaders.locus_id
$join_sql
$extra_query
";
  $query =~ s/\n|\r/ /g;
  
  # status #
  print STDERR "$query\n" if $verbose;
  
  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  
  if($$ret[0]){
    foreach my $row (@$ret){
      $dna_segs_r->{$$row[0]}{$$row[1]}{"leader"}{$$row[2]} = [@$row[3..$#$row]];
    }	
  }
  else{
    print STDERR "WARNING: no matching entries in leaders table!",
      " Leaders not added to dna_segs table!\n";
  }  
  #print Dumper %$dna_segs_r; exit;
}

sub get_DR_info{
  # getting direct repeat info from CLdb #
  my ($dbh, $dna_segs_r, $join_sql, $extra_query) = @_;
  
  my $query = "
SELECT 
loci.taxon_name,
loci.locus_id, 
DRs.DR_id, 
DRs.DR_start, 
DRs.DR_end
FROM Loci, DRs
WHERE Loci.locus_id = DRs.locus_id
$join_sql
$extra_query
";
  $query =~ s/\n|\r/ /g;
  
  # status #
  print STDERR "$query\n" if $verbose;
  
  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  die " ERROR: no matching entries for DR query!\n"
    unless $$ret[0];
  
  foreach my $row (@$ret){
    $dna_segs_r->{$$row[0]}{$$row[1]}{"DR"}{$$row[2]} = [@$row[3..$#$row]];
  }
  
  #print Dumper %$dna_segs_r; exit;
}	

sub get_spacer_info{
# getting spacer info from CLdb #
  my ($dbh, $dna_segs_r, $join_sql, $extra_query, $spacer_cutoff) = @_;
  
  # checking for spacer_clusters #
  my $q = "SELECT count(*) FROM spacer_clusters";
  my $chk = $dbh->selectall_arrayref($q);
  die " ERROR! no entries in spacer_cluster table!  Run clusterArrayElements.pl prior to this script!\n"
    unless @$chk;
  
  # query of spacers #
  ## strand-agnostic clusters ##
  my $query = "
SELECT 
loci.taxon_name,
loci.locus_id, 
spacers.spacer_id, 
spacers.spacer_start, 
spacers.spacer_end, 
spacer_clusters.cluster_id
FROM Loci, Spacers, Spacer_clusters
WHERE Loci.locus_id = Spacers.locus_id
AND Spacers.locus_id = Spacer_clusters.locus_id
AND Spacers.spacer_id = Spacer_clusters.spacer_id
AND Spacer_clusters.cutoff = $spacer_cutoff
$join_sql
$extra_query
";

  $query =~ s/\n|\r/ /g;
  
  # status #
  warn "$query\n" if $verbose;
  
  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  die " ERROR: no matching entries for spacer query!\n"
    unless $$ret[0];
  
  my %spacer_clusters; 
  foreach my $row (@$ret){
    # sanity check #
    ## should only have 1 taxon_name->locus_id->spacer_id ##
    die " ERROR: multiple entries for taxon_name->$$row[0], locus_id->$$row[1], spacer_id->$$row[2]!\n"
      if exists $dna_segs_r->{$$row[0]}{$$row[1]}{"spacer"}{$$row[2]};
    
    # converting clusterID to just cluster number #
    $$row[5] =~ s/.+_//;
    
    # loading dna_segs_r #
    $dna_segs_r->{$$row[0]}{$$row[1]}{"spacer"}{$$row[2]} = [@$row[3..$#$row]];
    $spacer_clusters{$$row[5]}++;			# unique clusters 
  }
  
  #print Dumper %$dna_segs_r; exit;
  #print Dumper \%spacer_clusters; exit;
  #exit;
  return \%spacer_clusters;
}

sub get_subtypes{
  my ($dbh, $join_sql, $extra_query ) = @_;
  
  my $query = "
SELECT 
loci.taxon_name,
loci.locus_id,
loci.subtype
FROM Loci loci, Loci b
WHERE loci.locus_id = b.locus_id
$join_sql
$extra_query
GROUP BY loci.locus_id
";
  $query =~ s/\n|\r/ /g;
  
  
  # status #
  print STDERR "$query\n" if $verbose;
  
  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  die " ERROR: no matching entries for subtype query!\n",
    "  QUERY: '$query'\n"
      unless $$ret[0];
  
  my %subtypes; 
  # taxon_name => locus_id => subtype
  foreach my $row (@$ret){
    $subtypes{$$row[0]}{$$row[1]} = $$row[2]; 		
  }
  
  #print Dumper %subtypes; exit;
  return \%subtypes;
}




