#!/usr/bin/env perl

=pod

=head1 NAME

writeTable.pl -- write out >=1 table from a CLdb

=head1 SYNOPSIS

writeTable.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -tables <char>

>=1 table to write out. 

=item -list  <int> 

=over

=item '1' = list table names and number of entries per table

=item '2' = list table name and fields per table

=back

=item -prefix  <char>

Output file prefix. [""]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc writeTable.pl

=head1 DESCRIPTION

List or write out one or more tables 
in CLdb. Output tables are tab-delimited.
Table headers are included. A warning is
given if a table has 0 entries.

=head1 EXAMPLES

=head2 Write out all tables

writeTable.pl -d CLdb.sqlite 

=head2 Write out just 'Loci' table

writeTable.pl -d CLdb.sqlite -t loci

=head2 Write out 'Spacers' and 'DRs' table

writeTable.pl -d CLdb.sqlite -t spacers drs

=head2 List table names and number of entries

writeTable.pl -d CLdb.sqlite -l

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

my ($verbose, $database_file, $listOpt, $prefix);
my (@tables);
GetOptions(
	   "database=s" => \$database_file,
	   "table=s{,}" => \@tables,
	   "prefix=s" => \$prefix,
	   "list=i" => \$listOpt,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "-list must be '1' or '2'\n" if defined $listOpt
  and $listOpt != 1 and $listOpt != 2;

file_exists($database_file, "database");


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);


# getting all tables if list of tables not provided
if(@tables){
  # checking for user-defined tables of interest #
  map{ table_exists($dbh, $_) or die "ERROR: $_ does not exist!\n" } @tables;
}
else{
  # using all tables
  @tables = listTables($dbh) unless @tables;
}


# listing tables & number_of_entries or tables & fields #
## tables
if (defined $listOpt and $listOpt == 1){
  # getting number of entries per table
  list_table_entries($dbh, \@tables);

  $dbh->disconnect();
  exit;
}
## fields
if (defined $listOpt and $listOpt == 2){
  # getting number of entries per table
  list_table_fields($dbh, \@tables);

  $dbh->disconnect();
  exit;
}


	
# writing tables #
write_tables($dbh, 
	     \@tables, 
	     $prefix);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines

sub write_tables{
# writing out each table to a file
  my $dbh = shift or die "ERROR: provide dbh\n";
  my $tables_r = shift or die "ERROR: provide arrayref of tables\n";
  my $prefix = shift; # or die "ERROR: provide output file prefix\n";
  $prefix = "" unless defined $prefix;

  
  foreach my $table (sort @$tables_r){
    # skipping if not entries for table
    unless( n_entries($dbh, $table) > 0){    
      print STDERR " WARNING: no entries in $table table! Skipping!\n";
      next;
    }
    
    # sql fed to database #
    my $outfile = "$table.txt";
    $outfile = join("_", $prefix, "$table.txt") if $prefix;
    open OUT, ">$outfile" or die $!;
    
   
    # sql 
    my $sql = "SELECT * from $table";
    my $sth = $dbh->prepare($sql) or die $dbh->err;

    # header
    print OUT join("\t", @{$sth->{NAME}}), "\n";
    
    # body
    $sth->execute();
    my $ret = $sth->fetchall_arrayref();
    foreach my $row (@$ret){
      map{ $_ = 'NULL' unless defined $_ } @$row;
      print OUT join("\t", @$row), "\n";
    }

    # finish up
    close OUT or die $!;
    print STDERR "...$table table written: $outfile\n";
  }
}


sub list_table_fields{
  # listing fields in certain tables
  my $dbh = shift or die "Provide dbh\n";
  my $tables_r = shift or die "Provide arrayref of table names\n";

  # header
  print join("\t", qw/Table Field/), "\n";

  # body
  foreach my $table (sort @$tables_r){
    my $sql = "SELECT * from $table limit 1";

    my $sth = $dbh->prepare($sql) or die $dbh->err;
    my $names_r = $sth->{NAME};
    
    foreach my $name (@$names_r){
      print join("\t", $table, $name), "\n";
    }
  }
}



sub list_table_entries{
# listing the number entries per table
  my $dbh = shift or die "Provide dbh\n";
  my $tables_r = shift or die "Provide arrayref of table names\n";

  # header
  print join("\t", qw/Table Number_entries/), "\n";

  # body
  foreach my $table (sort @$tables_r){
    my $cmd = "SELECT count(*) from $table";
    my $ret = $dbh->selectrow_arrayref($cmd) or die $dbh->err;

    $ret->[0] = 0 if not defined $ret->[0];
    print join("\t", $table, $ret->[0]), "\n";
  }
}


sub listTables{
# returns array of table names
  my $dbh = shift or die "Provide dbh\n";

  my %tbls;
  foreach my $tbl (grep /"main"\./, $dbh->tables()){
    $tbl =~ s/("main"\."|")//g;
    $tbls{$tbl} = 1;
  }

  return keys %tbls;
}




