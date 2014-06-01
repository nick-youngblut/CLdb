#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_loadArrays.pl -- loading direct repeats & spacers into CLdb

=head1 SYNOPSIS

CLdb_loadArrays.pl [flags] 

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

This help message.

=back

=head2 For more information:

perldoc CLdb_loadArrays.pl

=head1 DESCRIPTION

Parse array files (CRISPRFinder format) and load
the spacers and direct repeats into the CLdb database.

Array file names are obtained from the loci table in the
CRISPR database.

Array files must be in the array directory in $CLdb_HOME

Start-stop is based on CRISPRFinder orientation, which is
always + strand.

=head1 EXAMPLES

=head2 Usage:

CLdb_loadArrays.pl -d CLdb.sqlite 

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
use CLdb::utilities qw/
	file_exists 
	connect2db
	lineBreaks2unix
	get_file_path/;
use CLdb::query qw/
	table_exists/;
use CLdb::load qw/
	load_db_table/;
use CLdb::load::loadArray qw/
			      make_headers
			      parse_array_file
			    /;


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
my $db_path = get_file_path($database_file);
my $array_path = File::Spec->catdir($db_path, "array");
$array_path = File::Spec->rel2abs($array_path);
die "ERROR: cannot find '$array_path'\n" unless -d $array_path;

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# getting array file names from database #
my $array_r = get_array_file_names($dbh);

# loading array tables #
my %arrays;
foreach my $array_file (@$array_r){
  parse_array_file($array_path, $$array_file[1], $$array_file[0],\%arrays);
}
	
my ($dr_header_r, $sp_header_r) = make_headers();

# adding spacers/DR to database #
## DR entries ##
load_db_table($dbh, 'DRs', $dr_header_r, $arrays{'DR'});
## spacer entries ##
load_db_table($dbh, 'spacers', $sp_header_r, $arrays{'spacer'});

# disconnect #
$dbh->disconnect();
exit;


### Subroutines 
sub get_array_file_names{
# querying genbank names from sqlite loci table #
  my ($dbh) = @_;
  
  my $cmd = "SELECT locus_id, array_file from loci where array_file is not null";
  my $names_r = $dbh->selectall_arrayref($cmd);
  
  die " ERROR: no array file names found!\n" unless $names_r;
  #print Dumper $names_r; exit;
  return $names_r;
}



