#!/usr/bin/env perl

use strict;
use warnings;
use Data::Dumper;
use DBI;

my $usage = <<HERE;
Usage:
\t$0 CLdbFile spacerBlast taxon_name [query_regex]

Description:
Getting CLdb spacer clusters for taxon
and then pulling those (self-hits)
out of a spacer blast table file.

WARNING:
Assuming spacer blast query (qseqid) at 
least starts as "NA|spacer|NA|cluster_ID"

HERE

# I/O check
die $usage if scalar @ARGV != 3;


# querying CLdb
## connect 2 db
my $db_file = $ARGV[0];
my %attr = (RaiseError => 0, PrintError=>1, AutoCommit=>0 );
my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr)
  or die " Can't connect to $db_file!\n";


## query
my $sql = <<HERE;
SELECT sc.cluster_ID
FROM loci l, spacer_clusters sc
WHERE l.locus_id = sc.locus_id
AND l.taxon_name == "$ARGV[2]"
AND sc.cutoff == 1
HERE

my $ret = $dbh->selectall_arrayref($sql);
die "ERROR: no entries found for query: '$sql'\n"
  unless scalar @$ret > 0;
printf STDERR "Number of spacer clusters found: %i\n", scalar @$ret;

$dbh->disconnect();

# convert to dict for lookup
my %ret;
map{ $ret{$_->[0]} = 1 } @$ret;
$ret = ();


### getting blast hits from blast file
open IN, $ARGV[1] or die $!;
my $found_cnt = 0;
while(<IN>){
  chomp;
  next if /^\s*$/;
  my @l = split /\t/;

  my @query = split /\|/, $l[0];
  die "ERROR: line $. does not have a recognized clusterID in qseqid\n"
    unless defined $query[3];

  # checking for self hit
  ## sseqid
  next unless $l[1] =~ /^$ARGV[2]/;

  ## qseqid
  #print $query[3], "\n";
  if (exists $ret{$query[3]}){
    print $_,"\n";
    $found_cnt++;
  }
}
close IN;


# status
printf STDERR "Number of blast hits written: %i\n", $found_cnt;
