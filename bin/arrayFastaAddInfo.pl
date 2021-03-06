#!/usr/bin/env perl

=pod

=head1 NAME

arrayFastaAddInfo.pl -- add CLdb info to fastas of array spacers or direct repeats

=head1 SYNOPSIS

arrayFastaAddInfo.pl [flags] < array.fasta > array_e.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <bool>

Add subtype? [FALSE]

=item -taxon_id  <bool>

Add taxon_id? [FALSE]

=item -taxon_name  <bool>

Add taxon_name? [FALSE]

=item -position  <bool>

Add start-stop position? [FALSE]

=item -order  <bool>

Add spacer-leader order? [FALSE]
# not implemented yet! #

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb --perldoc -- arrayFastaAddInfo.pl

=head1 DESCRIPTION

Add CLdb info to fastas of array elements (spacers/DRs).
Also, spacer/DR groups will be replicated (all spacers/DRs 
in the group will be written). 

=head1 EXAMPLES

=head2 Add element subtype and position 

arrayFastaAddInfo.pl -d CLdb.sqlite -sub -pos < elements.fna > elements_info.fna

=head2 Replicate spacer/DR groups (all elements in group)

arrayFastaAddInfo.pl -d CLdb.sqlite < grouped.fna > all.fna

=head2 Replicate spacer/DR groups & add element subtype 

arrayFastaAddInfo.pl -d CLdb.sqlite -sub < grouped.fna > all.fna

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
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
use CLdb::seq qw/
	read_fasta/;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my ($subtype_b, $taxon_id_b, $taxon_name_b, $pos_b, $order_b);
GetOptions(
	   "database=s" => \$database_file,
	   "subtype" => \$subtype_b,
	   "taxon_id" => \$taxon_id_b,
	   "taxon_name" => \$taxon_name_b,
	   "position" => \$pos_b,
	   "order" => \$order_b,
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

# loading fasta #
my $fasta_r = read_fasta();

# querying CLdb #
## adding to query ##
my $query = add_query_opts($subtype_b, $taxon_id_b, $taxon_name_b);

## querying ##
my $fasta2_r = get_info($dbh, $fasta_r, $query);

# writing edited fasta #
write_array_seq($fasta2_r);


#--- disconnect ---#
$dbh->disconnect();
exit;


### Subroutines
sub write_array_seq{
# writing arrays as fasta
  my ($fasta2_r) = @_;
  
  foreach my $seq (sort keys %$fasta2_r){
    foreach my $element (@{$fasta2_r->{$seq}{"data"}}){
      my $name = join("|", @$element);
      print join("\n", ">$name", $fasta2_r->{$seq}{"seq"}), "\n";
    }
  }
}

sub get_info{
  # getting necessary info from CLdb #
  ## query is sequence-dependent ##
  my ($dbh, $fasta_r, $query) = @_;
  
  my %res;
  foreach my $seq (keys %$fasta_r){
    my @seq = split /\|/, $seq;
    die "ERROR: '$seq' does not have 4 components separated by '|'\n"
      unless scalar @seq == 4;
    die "ERROR: '$seq[1]' must be 'spacer', 'DR' or 'leader'\n"
      unless $seq[1] =~ /^(spacer|DR|leader)$/i;
    
    my $cmd;
    if($seq[3] eq 'NA'){ # no spacer|DR grouping
      $cmd = query_noGroup(\@seq, $query,);
    }
    else{  # grouping
      $cmd = query_byGroup(\@seq, $query);
    }
    
    $cmd =~ s/[\n\t]+/ /g;
    $cmd =~ s/ +/ /g;

    # status
    print STDERR "$cmd\n" if $verbose;

    # query db
    my $ret = $dbh->selectall_arrayref($cmd);
    die "ERROR: no matching entries!\n"
      unless $$ret[0];
    
    # loading hash #
    $res{$seq}{"data"} = $ret;
    $res{$seq}{"seq"} = $fasta_r->{$seq};   
  }
  
  #print Dumper %res; exit;
  return \%res;
}

sub query_noGroup{
  my ($seq_r, $query) = @_;
  
  # table #  
  my $tbl_oi;
  if($seq_r->[1] =~ /spacer/i){
    $tbl_oi = "Spacers";
  }
  elsif($seq_r->[1] =~ /DR/i){
    $tbl_oi = "DRs";
  }
  elsif($seq_r->[1] =~ /leader/i){
    $tbl_oi = "Leaders";
  }
  else{ die "LOGIC ERROR: $!\n"; }

  my $prefix = $$seq_r[1];			# "spacer|DR|leader"
  my $cmd = "SELECT 
$tbl_oi.Locus_ID,
'$prefix'";
  $cmd .= $tbl_oi =~ /Spacers|DRs/ ? ",$tbl_oi.$prefix\_ID,'NA'" :
    ",'NA','NA'";
  
  # adding other query options #
  $cmd .= ", $query" if $query;
  if($pos_b){
    $cmd .= ", loci.Scaffold";
    $cmd .= ", $tbl_oi.$prefix\_start";
    $cmd .= ", $tbl_oi.$prefix\_end";
  }	
  
  # where statement #
  $cmd .= "
FROM $tbl_oi, loci
WHERE loci.locus_id = $tbl_oi.locus_id
AND loci.locus_id = $$seq_r[0]";

  
  # selecting by element id
  $cmd .= " AND $tbl_oi.$prefix\_ID = $$seq_r[2]" 
    if $tbl_oi =~ /Spacers|DRs/;
     
  #print Dumper $cmd; exit;
  return $cmd;
}

sub query_byGroup{
  my ($seq_r, $query) = @_;
  
  # table #
  my $tbl_oi = "Spacers";
  $tbl_oi = "DRs" if $$seq_r[1] =~ /DR/i;
  
  my $prefix = $$seq_r[1];			# "spacer|DR"
  my $cmd = "SELECT 
$tbl_oi.Locus_ID,
'$prefix',
$tbl_oi.$prefix\_ID,
$prefix\_clusters.cluster_id";
  
  # adding other query options #
  $cmd .= ", $query" if $query;
  if($pos_b){
    $cmd .= ", loci.Scaffold";
    $cmd .= ", $tbl_oi.$prefix\_start";
    $cmd .= ", $tbl_oi.$prefix\_end";
  }
  
  # where statement #
  $cmd .= "
FROM $tbl_oi, loci, $prefix\_clusters
WHERE loci.locus_id = $tbl_oi.locus_id 
AND $prefix\_clusters.locus_ID = $tbl_oi.locus_ID
AND $prefix\_clusters.$prefix\_ID = $tbl_oi.$prefix\_ID
AND $prefix\_clusters.cluster_id = '$$seq_r[3]'";
  
  return $cmd;
}

sub add_query_opts{
  my ($subtype_b, $taxon_id_b, $taxon_name_b) = @_;
  
  my @query;
  push @query, "loci.subtype" if $subtype_b;
  push @query, "loci.taxon_id" if $taxon_id_b;
  push @query,"loci.taxon_name" if $taxon_name_b;
  
  return join(",", @query);
}








