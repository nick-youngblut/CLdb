#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_writeTable.pl -- write out one or more tables in CLdb

=head1 SYNOPSIS

CLdb_writeTable.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -table  <char>

>=1 table to write out. 

=item -list  <bool>

Just list table names and number of entries in each table. [FALSE]

=item -prefix  <char>

Output file prefix. [""]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_writeTable.pl

=head1 DESCRIPTION

List or write out one or more tables 
in CLdb. Output tables are tab-delimited.
Table headers are included. A warning is
given if a table has 0 entries.

=head1 EXAMPLES

=head2 Write out all tables

CLdb_writeTable.pl -d CLdb.sqlite 

=head2 Write out just 'Loci' table

CLdb_writeTable.pl -d CLdb.sqlite -t loci

=head2 Write out 'Spacers' and 'DRs' table

CLdb_writeTable.pl -d CLdb.sqlite -t spacers drs

=head2 List table names and number of entries

CLdb_writeTable.pl -d CLdb.sqlite -l

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

my ($verbose, $database_file, $list_tables, $prefix);
my (@tables);
GetOptions(
	   "database=s" => \$database_file,
	   "table=s{,}" => \@tables,
	   "prefix=s" => \$prefix,
	   "list" => \$list_tables,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking for tables of interest #
my $table_list_r = list_tables($dbh);
$table_list_r = entry_count($dbh, $table_list_r);

# listing tables #
if ($list_tables){
	list_table_entries($table_list_r);
	$dbh->disconnect();
	exit;
	}
	
# writing tables #
write_tables($database_file, \@tables, $table_list_r, $prefix);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub list_table_entries{
	my ($table_list_r) = @_;
	
	print join("\t", qw/Table Number_entries/), "\n";
	foreach my $table (sort keys %$table_list_r){
		print join("\t", $table, $table_list_r->{$table}), "\n";
		}
	}

sub write_tables{
	my ($database_file, $tables_r, $table_list_r, $prefix) = @_;
	
	foreach my $table (keys %$table_list_r){
		next if @$tables_r && ! grep /^$table$/i, @$tables_r;	# just tables selected (if -t provided)

		# no entry warning #
		unless( $table_list_r->{$table} > 0){
			print STDERR " WARNING: no entries in $table table! Skipping!\n";
			next;
			}
			
		# sql fed to database #
		my $outfile = "$table.txt";
		$outfile = join("_", $prefix, "$table.txt") if $prefix;
		
		my $sql = "
.header on
.sep \"\\t\"
.output $outfile
SELECT * from $table;";
		
		open PIPE, "| sqlite3 $database_file" or die $!;
		print PIPE $sql;
		close PIPE;

		print STDERR "...$table table written: $table.txt\n";
		}

	}

sub entry_count{
	my ($dbh, $table_list_r) = @_;
	my %table;
	foreach my $table (@$table_list_r){
		my $q = "SELECT count(*) FROM $table";
		my $res = $dbh->selectall_arrayref($q);
		$table =~ tr/A-Z/a-z/;
		$table{$table} = $$res[0][0];
		}
		#print Dumper %table; exit;
	return \%table;
	}

sub list_tables{
	my $dbh = shift;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}




