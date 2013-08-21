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


my ($verbose, $CLdb_sqlite, @ITEP_sqlite);
my $spacer_cutoff = 1;
GetOptions(
	   "database=s" => \$CLdb_sqlite,
	   "ITEP=s{,}" => \@ITEP_sqlite,
	   "cutoff=f" =>  \$spacer_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# checking CLdb #
die " ERROR: provide a database file name!\n"
	unless $CLdb_sqlite;
die " ERROR: cannot find CLdb database file!\n"
	unless -e $CLdb_sqlite;
	
# checking ITEP input #
if(@ITEP_sqlite){
	die " ERROR: '-ITEP' must be 2 arguments (database_file, runID)!\n"
		unless scalar @ITEP_sqlite >= 2;
	die " ERROR: cannot find ITEP database file!\n"
		unless -e $ITEP_sqlite[0];
	}

### MAIN
# connect 2 CLdb #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$CLdb_sqlite", '','', \%attr) 
	or die " Can't connect to $CLdb_sqlite!\n";

# connect to ITEP #
my $dbh_ITEP;
if(@ITEP_sqlite){
	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	$dbh_ITEP = DBI->connect("dbi:SQLite:dbname=$ITEP_sqlite[0]", '','', \%attr) 
		or die " Can't connect to $ITEP_sqlite[0]!\n";	
	}


# loading dna_segs table #
my $dna_segs_r = load_dna_segs();

# querying ITEP #
## for querying ITEP, need:
	# runID
	# geneIDs (by dna_segs_id)

# querying CLdb #
## for querying spacer_hclust, need: ##
	# spacerIDs (by dna_segs_id)
## for querying spacer_pairwise_blast, need: ##
	# spacerIDs (by dna_segs_id)
	
		

## joining query options (for table join) ##
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


# disconnect #
$dbh->disconnect();
$dbh_ITEP->disconnect() if $dbh_ITEP;
exit;


### Subroutines
sub load_dna_segs{
# loading dna_segs from stdin #

	my @dna_segs_order;
	my %dna_segs;
	my %header;
	my %dna_segs_ids;
	my %name_index;
	while(<>){
		chomp;
		next if /^\s*$/;

		# header #
		if($.==1){
			# indexing header #
			my @tmp = split /\t/;
			for my $i (0..$#tmp){
				$header{$tmp[$i]} = $i;
				}
			
			# checking for existence of columns #
			my @check = qw/taxon_name locus_id subtype dna_segs_id/;
			map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
				unless exists $header{$_} } @check;
			}
		else{
			my @line = split /\t/;
			my $taxon_name = $line[$header{"taxon_name"}];
			my $locus_id = $line[$header{"locus_id"}];
			my $subtype = $line[$header{"subtype"}];
			my $dna_seg_id = $line[$header{"dna_segs_id"}];
			
			# order of loci#
			push( @dna_segs_order, $taxon_name)
				unless exists $dna_segs{$taxon_name};
			
			# dna_segs  #
			push( @{$dna_segs{$taxon_name}{$locus_id}{"entries"}}, \@line );
			$dna_segs{$taxon_name}{$locus_id}{"subtype"} = $subtype;
			
			# dna_segs_id index #
			$dna_segs_ids{$locus_id} = $dna_seg_id;
			
			# name index #
			$name_index{$dna_seg_id} = $taxon_name;
			}
		}
		
		#print Dumper @dna_segs_order;
		#print Dumper %dna_segs; exit;
		#print Dumper %dna_segs_ids; exit;
		#print Dumper %name_index; exit;
	return \%dna_segs, \@dna_segs_order, \%header, \%dna_segs_ids ,\%name_index;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}



__END__

=pod

=head1 NAME

CLdb_dna_segs_make.pl -- making dna_segs table for plotting

=head1 SYNOPSIS

CLdb_dna_segs_make.pl [flags] > dna_segs.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -ITEP

Get gene cluster info from ITEP. 2 arguments required: (ITEP_sqlite_file, cluster_runID).

=item -cutoff

Spacer clustering cutoff for spacer coloring (0.8 - 1). [1]

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine which sequences are returned.

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_dna_segs_make.pl

=head1 DESCRIPTION

Make a basic dna_segs object needed for plotting.

=head2 Coloring genes by gene cluster

Provide and ITEP sqlite file name and cluster run ID
to add cluster info to the table (used for coloring)

=head1 EXAMPLES

=head2 Plotting all loci classified as subtype 'I-A'

CLdb_dna_segs_make.pl -d CLdb.sqlite -sub I-A 

=head2 Gene cluster info from ITEP

CLdb_dna_segs_make.pl -d CLdb.sqlite -sub I-A  -I DATABASE.sqlite all_I_2.0_c_0.4_m_maxbit

=head2 No broken loci

CLdb_dna_segs_make.pl -da CLdb.sqlite -sub I-A -q "AND loci.operon_status != 'broken'"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

