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
GetOptions(
	   "database=s" => \$database_file,
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

# database metadata #
my $table_list_r = list_tables($dbh);
check_for_genes_table($table_list_r);
#get_genes_table_info($dbh);

# loading genes table #
my ($genes_r, $header_r) = get_gene_table();

# updating / loading_db #
load_new_entries($dbh, $genes_r, $header_r);

# disconnect to db #
$dbh->disconnect();
exit;


### Subroutines
sub load_new_entries{
# adding new entries to db #
## conflicts should be replaced (ie. updated ##
	my ($dbh, $genes_r, $header_r) = @_;

	# making locus_id = NULL #
	my @keys = keys %$header_r;
	@keys = grep(!/^sequence$/i, @keys);		# removed from table; not needed
	my @values = @$header_r{@keys};

	# loading entries #
	my $cmd = "INSERT INTO genes(" . join(",", @keys) . ") VALUES(?".",?"x$#keys . ")";
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $locus_id (keys %$genes_r){
		foreach my $gene_id (keys %{$genes_r->{$locus_id}}){			
			#	print Dumper @{$genes_r->{$locus_id}{$gene_id}}[@values]; exit;
			$sql->execute( @{$genes_r->{$locus_id}{$gene_id}}[@values] );	
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @{$genes_r->{$locus_id}{$gene_id}}[@values]), "\n";
				}
			else{ $cnt++; }
			}
		}
	$dbh->commit;

	print STDERR "...Number of entries added/updated in gene table: $cnt\n"
		unless $verbose;
	}


sub get_gene_table{
	my %genes;
	my %header;
	while(<>){
		chomp;
		next if /^\s*$/;

		if($. == 1){ # loading header
			tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
			my @line = split /\t/;
			for my $i (0..$#line){
				$header{$line[$i]} = $i;		# column_name => index
				}
			}
		else{
			my @line = split /\t/;
			$line[$header{"locus_id"}] =~ s/^cli\.//;
			$genes{$line[$header{"locus_id"}]}{$line[$header{"gene_id"}]} = \@line;
			}
		}
		#print Dumper %header; exit;
		#print Dumper %genes; exit; 
	return (\%genes, \%header);;
	}


sub check_for_genes_table{
	my ($table_list_r) = @_;
	die " ERROR: Genes table not found in database!\n"
		unless grep(/^Genes$/i, @$table_list_r);
	}

sub list_tables{
	my ($dbh) = @_;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}

__END__

=pod

=head1 NAME

Cldb_loadGenes.pl -- adding/updating gene table in CRISPR database

=head1 SYNOPSIS

Cldb_loadGenes.pl [flags] < gene_table.txt

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

perldoc Cldb_loadGenes.pl

=head1 DESCRIPTION

Add/update genes entries in Genes table within the
specified CRISPR database.

A gene table made by CLdb_getGenesInLoci.pl can be
piped directly into CLdb_loadGenes.pl, or 
the aliases (or other values) can be currated (manually)
first. 

=head1 EXAMPLES

=head2 Basic usage:

Cldb_loadGenes.pl -d CRISPR.sqlite < genes_table.txt

=head2 Piping from CLdb_getGenesInLoci.pl 

CLdb_getGenesInLoci.pl -d CRISPR.sqlite | Cldb_loadGenes.pl -d CRISPR.sqlite

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

