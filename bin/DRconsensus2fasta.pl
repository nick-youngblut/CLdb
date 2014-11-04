#!/usr/bin/env perl

=pod

=head1 NAME

DRconsensus2fasta.pl -- write direct repeat consensus sequences to fasta

=head1 SYNOPSIS

DRconsensus2fasta.pl [flags] > DR_consensus.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -IUPAC  <bool>

Get consensus sequences with IUPAC nucleotide ambiguity codes
instead of consensus by threshold? [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -verbose  <bool> 

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb_perldoc DRconsensus2fasta.pl

=head1 DESCRIPTION

Get direct repeat (DR) consenesus sequences. Output options:

=head3 Threshold consensus (default)

The consensus residue must have a frequency of >= than '-cutoff',
otherwise a 'N' will be used.

=head3 IUPAC consensus

Base ambiguity letters are used.

Output sequence names: ">locus_ID|subtype"

DR consensus must be calculated prior!

=head1 EXAMPLES

=head2 Write all DR consensus sequences:

DRconsensus2fasta.pl -d CLdb.sqlite 

=head2 Write all DR consensus sequences (IUPAC):

DRconsensus2fasta.pl -d CLdb.sqlite -I

=head2 Refine query to a specific subtype & 2 taxon_id's

DRconsensus2fasta.pl -da CLdb.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

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
	join_query_opts
	get_array_seq/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;


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

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

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
		print join("\n", 
					join('|', ">$locus_id", $arrays_r->{$locus_id}{'subtype'}),
					$arrays_r->{$locus_id}{'seq'}), "\n";
		}
	}
	
sub get_consensus{
	my ($dbh, $iupac_consensus, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query;
	if($iupac_consensus){
		$query = "SELECT DR_consensus.Locus_ID, DR_consensus.Consensus_Sequence_IUPAC,
				loci.subtype
			FROM DR_consensus, Loci WHERE DR_consensus.locus_id = Loci.locus_id $join_sql";
		}
	else{			# consensus threshold
		$query = "SELECT DR_consensus.Locus_ID, DR_consensus.Consensus_Sequence_Threshold,
				loci.subtype
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
		$DRs{$$row[0]}{'seq'} = $$row[1];
		$DRs{$$row[0]}{'subtype'} = $$row[2];		
		}
	
		#print Dumper %DRs; exit;
	return \%DRs;
	}
	
sub join_query_opts_OLD{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND Loci.$cat IN (", join(", ", @$vals_r), ")");
	}


