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


my ($verbose, $database_file, $iupac_consensus);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "IUPAC" => \$iupac_consensus, 		# IUPAC instead of threshold 
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

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting arrays of interest from database #
my $arrays_r = get_consensus($dbh, $iupac_consensus, $extra_query, $join_sql);


# writing fasta #
write_DR_fasta($arrays_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_DR_fasta{
# writing arrays as fasta
	my ($arrays_r) = @_;
	
	foreach my $locus_id (keys %$arrays_r){
		print join("\n", ">cli.$locus_id", $arrays_r->{$locus_id}), "\n";
		}
	}
	
sub get_consensus{
	my ($dbh, $iupac_consensus, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query;
	if($iupac_consensus){
		$query = "SELECT DR_consensus.Locus_ID, DR_consensus.Consensus_Sequence_IUPAC
FROM DR_consensus, Loci WHERE DR_consensus.locus_id = Loci.locus_id $join_sql";
		}
	else{			# consensus threshold
		$query = "SELECT DR_consensus.Locus_ID, DR_consensus.Consensus_Sequence_Threshold
FROM DR_consensus, Loci WHERE DR_consensus.locus_id = Loci.locus_id $join_sql";
		}

	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %DRs;
	foreach my $row (@$ret){
		$DRs{$$row[0]} = $$row[1];
		}
	
		#print Dumper %DRs; exit;
	return \%DRs;
	}
	
sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND Loci.$cat IN (", join(", ", @$vals_r), ")");
	}



__END__

=pod

=head1 NAME

CLdb_DRconsensus2fasta.pl -- write direct repeat consensus sequences to fasta

=head1 SYNOPSIS

CLdb_DRconsensus2fasta.pl [flags] > DR_consensus.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -IUPAC

Get consensus sequences with IUPAC nucleotide ambiguity codes
instead of consensus by threshold? [FALSE]

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

perldoc CLdb_DRconsensus2fasta.pl

=head1 DESCRIPTION

Get direct repeat (DR) consenesus sequences. Output options:

=head3 Threshold consensus (default)

The consensus residue must have a frequency of >= than '-cutoff',
otherwise a 'N' will be used.

=head3 IUPAC consensus

Base ambiguity letters are used.


DR consensus must be calculated prior!

=head1 EXAMPLES

=head2 Write all DR consensus sequences:

CLdb_DRconsensus2fasta.pl -d CLdb.sqlite 

=head2 Write all DR consensus sequences (IUPAC):

CLdb_DRconsensus2fasta.pl -d CLdb.sqlite -I

=head2 Refine query to a specific subtype & 2 taxon_id's

CLdb_DRconsensus2fasta.pl -da CLdb.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

