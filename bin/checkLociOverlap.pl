#!/usr/bin/env perl

=pod

=head1 NAME

checkLociOverlap.pl -- checking for loci position overlap

=head1 SYNOPSIS

checkLociOverlap.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb --perldoc -- checkLociOverlap.pl

=head1 DESCRIPTION

Check for overlap of loci based on locus_start-locus_end
position info in the CLdb loci table. Overlap 
may be caused by erroneous entries in the database.

Overlapping entries will be written to STDOUT.

=head1 EXAMPLES

=head2 Usage:

CLdb -- checkLociOverlap.pl -d CLdb.sqlite 

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
use Set::IntervalTree;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
	

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	  );


#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);
my $column_list_r = list_columns("loci", 1);


# getting input loci table #
my ($loci_r) = get_loci_table($dbh, $column_list_r);

# making interval trees #
make_interval_trees($loci_r, $column_list_r);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub make_interval_trees{
# making interval trees and checking for overlap at the same time #
# only comparing arrays in the same scaffold in the name taxon #
  my ($loci_r, $column_list_r) = @_;
  
  # status #
  print STDERR "...Checking for loci overlap on the same scaffold\n";
  
  # checking for overlap #
  my $overlap_cnt = 0;
  foreach my $taxon_name (keys %$loci_r){
    foreach my $scaf (keys %{$loci_r->{$taxon_name}}){
      my $tree = Set::IntervalTree->new;
      
      foreach my $locus_id (keys %{$loci_r->{$taxon_name}{$scaf}} ){		
	# getting loci start-end info #
	my $locus_start = ${$loci_r->{$taxon_name}{$scaf}{$locus_id}}[$column_list_r->{"locus_start"}];
	my $locus_end = ${$loci_r->{$taxon_name}{$scaf}{$locus_id}}[$column_list_r->{"locus_end"}];
	
	# checking for overlap w/ already added start-end #
	my $results = $tree->fetch($locus_start, $locus_end);
	if(@$results){
	  $overlap_cnt++;
	  foreach (@$results){
	    print join("\t", "Overlapping loci",
		       @{$loci_r->{$taxon_name}{$scaf}{$_}}), "\n";
	  }
	  print join("\t", "Overlapping loci",
		     @{$loci_r->{$taxon_name}{$scaf}{$locus_id}}), "\n";				
	}
	
	# adding to interval tree #
	$tree->insert($locus_id, $locus_start, $locus_end);
      }
    }
  }
  
  print STDERR "...No overlapping loci found!\n" unless $overlap_cnt;
}

sub get_loci_table{
  # getting whole loci table #
  my ($dbh, $column_list_r) = @_;
  
  # getting loci table #
  my $cmd = join(" ", "SELECT", 
		 join(",", sort{$column_list_r->{$a}<=> $column_list_r->{$b}} keys %$column_list_r), 
		 "from loci");
  my $loci_r = $dbh->selectall_arrayref($cmd);
  
  die " ERROR: no entries in loci table!\n" unless @$loci_r;
  
  # hash: taxon_name=>locus_id=>scaffold=>entry # 
  my %loci;
  foreach my $row (@$loci_r){
    # uninitialized to blank #
    map{$_ = "" unless $_ }  @$row;
    
    # loading hash #
    my $taxon_name = $$row[$column_list_r->{"taxon_name"}];
    my $locus_id = $$row[$column_list_r->{"locus_id"}];
    my $scaffold = $$row[$column_list_r->{"scaffold"}];
    $loci{$taxon_name}{$scaffold}{$locus_id} = $row;
  }
  #print Dumper %loci; exit;
  return \%loci;
}

sub check_for_loci_table{
  my ($table_list_r) = @_;
  die " ERROR: loci table not found in database!\n"
    unless grep(/^loci$/i, @$table_list_r);
}

sub list_tables{
  my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
  return [keys %$all];
}

sub list_columns{
  my ($column_list, $silent_ret) = @_;
  my $all = $dbh->selectall_arrayref("pragma table_info($column_list)");
  #print Dumper $$all[0]; exit;
  my %tmp;
  for my $i (0..$#$all){ 
    $$all[$i][1] =~ tr/A-Z/a-z/;		# lower case for matching
    $tmp{$$all[$i][1]} = $i; 
  }
  if($silent_ret){ return \%tmp; }
  else{  print "Columns:\n", join(",\n", keys %tmp), "\n\n";  exit; }
}


