#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $spacer_bool);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_bool,
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# getting arrays of interest from database #
my $arrays_r = get_arrays($dbh, $spacer_bool, $extra_query);

# writing fasta #
write_arrays_fasta($arrays_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_arrays_fasta{
# writing arrays as fasta
	my ($arrays_r) = @_;
	
	foreach my $locus_id (keys %$arrays_r){
		foreach my $x_id (keys %{$arrays_r->{$locus_id}}){
			print join("\n", ">lci.$locus_id\__$x_id", $arrays_r->{$locus_id}{$x_id}), "\n";
			}
		}
	}

sub get_arrays{
	my ($dbh, $spacer_bool, $extra_query) = @_;
	
	# make query #
	my $query;
	if($spacer_bool){		# direct repeat
		$query = "SELECT Locus_ID, Repeat_ID, Repeat_sequence FROM directrepeats";
		}
	else{					# spacer
		$query = "SELECT Locus_ID, Spacer_ID, Spacer_sequence FROM spacers";
		}
	$query = join(" ", $query, $extra_query);
	
	#print Dumper $query; exit;
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	
	my %arrays;
	foreach my $row (@$ret){
		$arrays{$$row[0]}{$$row[1]} = $$row[2];
		}
	
	#	print Dumper %arrays; exit;
	return \%arrays;
	}


#my $query = "SELECT FROM locus_id";

__END__

=pod

=head1 NAME

Cdb_array2fasta.pl -- write CRISPR array spacers or direct repeats to fasta

=head1 SYNOPSIS

Cdb_array2fasta.pl [options] > array.fasta

=head2 options

=over

=item -d 	sqlite3 database (required).

=item -r 	Get direct repeats instead of spacers.

=item -q 	Extra sql to refine which sequences are returned.

=item -h	This help message

=back

=head2 For more information:

perldoc Cdb_array2fasta.pl

=head1 DESCRIPTION

Get spacer or direct repeat sequences from the CRISPR database
and write them to a fasta.

By default, all spacers or direct repeats (if '-r') will be written.
The '-q' flag can be used to refine the query to certain sequences (see examples).

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

Cdb_array2fasta.pl -data CRISPR.sqlite 

=head2 Write all direct repeats to a fasta:

Cdb_array2fasta.pl -data CRISPR.sqlite -r

=head2 Refine spacer sequence query:

Cdb_array2fasta.pl -data CRISPR.sqlite -q "where LOCUS_ID=1" 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CRISPR_db/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

