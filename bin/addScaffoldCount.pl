#!/usr/bin/env perl

=pod

=head1 NAME

addScaffoldCount.pl -- adding number of scaffolds in unmerged genbanks to loci table

=head1 SYNOPSIS

addScaffoldCount.pl [flags] 

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

CLdb_perldoc addScaffoldCount.pl

=head1 DESCRIPTION

Add the number of scaffolds in the orginal, unmerged genbank file 
for each organism. 

Scaffolds are counted from 'source' in the genbank files provided
in the loci table.

Genbank files must be in $HOME/genbank/

=head1 EXAMPLES

=head2 Usage:

addScaffoldCount.pl -data CRISPR.sqlite

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
	n_entries
	join_query_opts/;
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
my $genbank_path = path_by_database($database_file);
$genbank_path = File::Spec->rel2abs($genbank_path);

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);


# getting genbank file names from database #
my $genbank_r = get_genbank_names($dbh);

# getting number of scaffold in each genbank #
my %scaf_cnt;
foreach my $genbank (@$genbank_r){
	count_scaffolds($genbank_path, $$genbank[0], \%scaf_cnt);
	}

# updating db #
update_db($dbh, \%scaf_cnt);

# disconnect #
$dbh->disconnect();
exit;

### Subroutines
sub update_db{
# updating scaffold info #
	my ($dbh, $scaf_cnt_r) = @_;
	
	my $cmd = "Update Loci SET scaffold_count = ? where genbank_file = ?";
	my $sql = $dbh->prepare($cmd);
	
	foreach my $genbank (keys %$scaf_cnt_r){
		
		$sql->execute( ($scaf_cnt_r->{$genbank}, $genbank) );
					
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", $genbank, $scaf_cnt_r->{$genbank}), "\n";
			}
		}
	$dbh->commit;	
	
	print STDERR "...Scaffold information added to loci table\n" unless $verbose;
	}

sub count_scaffolds{
# counting scaffolds of genbank #
	my ($genbank_path, $genbank, $scaf_cnt_r) = @_;
	
	die " ERROR: $genbank_path/$genbank not found!\n"
		unless -e "$genbank_path/$genbank";
		
	open IN, "$genbank_path/$genbank" or die $!;
	while(<IN>){
		$scaf_cnt_r->{$genbank}++ if /^ *source +/;
		}
	close IN;
	}

sub get_genbank_names{
# querying genbank names from sqlite loci table #
	my ($dbh) = @_;
	
	my $cmd = "SELECT distinct genbank_file from loci";
	my $names_r = $dbh->selectall_arrayref($cmd);
	
	die " ERROR: no genbank names found!\n" unless $names_r;
	
		#print Dumper $names_r; exit;
	return $names_r;
	}

sub path_by_database{
	my ($database_file) = @_;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "genbank");
	}


