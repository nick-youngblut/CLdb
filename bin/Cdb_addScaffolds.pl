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

my ($verbose, $database_file);
my $genbank_path = "./genbank/";
GetOptions(
	   "database=s" => \$database_file,
	   "path=s" => \$genbank_path, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$genbank_path = File::Spec->rel2abs($genbank_path);

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

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
	
	my $cmd = "Update Loci SET scaffolds = ? where genbank = ?";
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
		$scaf_cnt_r->{$genbank}++ if /^\s*source\s+/;
		}
	close IN;
	}

sub get_genbank_names{
# querying genbank names from sqlite loci table #
	my ($dbh) = @_;
	
	my $cmd = "SELECT genbank from loci";
	my $names_r = $dbh->selectall_arrayref($cmd);
	
	die " ERROR: no genbank names found!\n" unless $names_r;
		#print Dumper $names_r; exit;
	return $names_r;
	}


__END__

=pod

=head1 NAME

Cdb_addScaffolds.pl -- adding number of scaffolds in unmerged genbanks to loci table

=head1 SYNOPSIS

Cdb_addScaffolds.pl [options] 

=head2 options

=over

=item -d 	sqlite3 database (required).

=item -p 	Path to the genbank files. [./genbank]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc Cdb_addScaffolds.pl

=head1 DESCRIPTION

Add the number of scaffolds in the orginal, unmerged genbank file 
for each organism. 

Scaffolds are counted from 'source' in the genbank files provided
in the loci table.

=head1 EXAMPLES

=head2 Usage:

Cdb_addScaffolds.pl -data CRISPR.sqlite

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CRISPR_db/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

