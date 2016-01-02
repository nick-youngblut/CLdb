#!/usr/bin/env perl

=pod

=head1 NAME

getGeneInfo.pl -- get gene & loci (& ITEP cluster) info for >=1 gene (fig.peg format)

=head1 SYNOPSIS

getGeneInfo.pl [flags] < tbl_withGeneIds.txt > tbl_of_gene_info.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -column  <int>

Column containing gene_IDs (fig.peg format).

=item -ITEP  <char>

ITEP database (if you want clusterIDs).

=item -runID  <char>

ITEP cluster runID (for selecting cluster IDs).

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb_perldoc getGeneInfo.pl

=head1 DESCRIPTION

Get info on gene IDs (ITEP FIG.PEG format).
Gene info is joined with loci info on the gene,
and the results are written to STDOUT.

ITEP cluster run IDs will be in the last column
if -ITEP and -runid flags are provided.

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
use List::Util qw/max/;

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
	revcomp/;



#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $ITEP_db, $extra_query, $runid);
my $column = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "ITEP=s" => \$ITEP_db,
	   "runid=s" => \$runid,
	   "column=i" => \$column,
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
$column--;
file_exists($database_file, "database");
file_exists($ITEP_db, "ITEP_database") if defined $ITEP_db;
die "ERROR: provide an ITEP cluster runID (-r flag) along with the ITEP database file name!\n"
  if defined $ITEP_db && ! defined $runid;

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);
my $ITEP_dbh = connect2db($ITEP_db) if defined $ITEP_db;

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "genes");

# loading gene IDs
my $genes_r = load_gene_ids($column);

# querying CLdb
geneInfo($genes_r, $dbh);

# querying ITEP
ITEP_geneInfo($genes_r, $ITEP_dbh, $runid) if defined $ITEP_db;

# writing out table
write_table($genes_r);

#--- disconnect ---#
$dbh->disconnect();
$ITEP_dbh->disconnect() if defined $ITEP_dbh;
exit;


#--- Subroutines ---#
sub write_table{
  my ($genes_r) = @_;

  foreach my $id (keys %$genes_r){
    #map{ $_ = '' unless defined $_ } @{$genes_r->{$id}{CLdb}};
    #map{ $_ = '' unless defined $_ } @{$genes_r->{$id}{ITEP}};
    print join("\t", $id, 
	       @{$genes_r->{$id}{CLdb}},
	       @{$genes_r->{$id}{ITEP}}), "\n";
  }
}

sub ITEP_geneInfo{
# getting clusterID from gene id
  my ($genes_r, $ITEP_dbh, $runid) = @_;

  # query
  my $q = "
SELECT clusterid
FROM clusters
WHERE geneid = ?
AND runid='$runid'
";

  $q =~ s/\n/ /g;
  my $sth = $ITEP_dbh->prepare($q);

  foreach my $id (keys %$genes_r){
    $sth->bind_param(1, $id);
    $sth->execute();
    while(my $r = $sth->fetchrow_arrayref()){
      die "ERROR: >1 entry for gene \"$id\"\n"
	if exists $genes_r->{$id}{ITEP};
      my @r = @$r;
      map{ $_ = 'NULL' unless $_ } @r;
      $genes_r->{$id}{ITEP} = \@r;
    }
    $genes_r->{$id}{ITEP} = ['NULL'] 
      unless exists $genes_r->{$id}{ITEP};
  }

  #print Dumper %$genes_r; exit;
}

sub geneInfo{
# querying CLdb for gene info; genes-loci table join
  my ($genes_r, $dbh) = @_;

  # query
  my $q = "
SELECT g.*, l.*
FROM genes g, loci l
WHERE g.locus_id=l.locus_id
AND g.gene_id = ?
";
  $q =~ s/\n/ /g;

  my $sth = $dbh->prepare($q);

  my @Nrow = (0);
  foreach my $i (keys %$genes_r){
    $sth->bind_param(1, $i);
    $sth->execute();
    while(my $r = $sth->fetchrow_arrayref()){
      die "ERROR: >1 entry for gene \"$i\"\n"
	if exists $genes_r->{$i}{CLdb};
      my @r = @$r;
      map{ $_ = 'NULL' unless $_ } @r;
      $genes_r->{$i}{CLdb} = \@r;
      push @Nrow, scalar @r;
    }
  }

  # adding '' for entries with nothing
  my $max_row = max @Nrow;
  foreach my $i (keys %$genes_r){
    my $nrow = 0;
    $nrow = scalar @{$genes_r->{$i}{CLdb}} 
      if exists $genes_r->{$i}{CLdb};
    die "LOGIC ERROR $!\n" unless $nrow <= $max_row;

    for(my $ii=($nrow+1);$ii<=$max_row;$ii++){
      push @{$genes_r->{$i}{CLdb}}, 'NULL';
    }
  }

  #print Dumper $genes_r; exit;
}

sub load_gene_ids{
  my ($column) = @_;

  my %genes;
  while(<>){
    chomp;
    next if /^\s*$/;
    my @l = split /\t/;
    die "ERROR: line $. does not contain column ", $column+1, "\n"
      unless defined $l[$column];

    $genes{$l[$column]}{tmp} = 1; # unique geneIDs
  }

  #print Dumper keys %genes; exit;
  return \%genes;
}



