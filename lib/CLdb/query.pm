package CLdb::query;

# module use #
## core ##
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use DBI;

## CLdb ##
use CLdb::seq qw/
	revcomp/;


# export #
use base 'Exporter';
our @EXPORT_OK = qw/
table_exists
n_entries
join_query_opts
get_array_seq
get_leader_pos
get_arrays_seq_byLeader
list_columns
/;

	
=head1 NAME

CLdb::query

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for querying CLdb

=head1 EXPORT_OK

table_exists
n_entries
join_query_opts
get_array_seq

=cut

sub list_columns{
#-- description --#
# listing all columns in specified table (table must exist) 
# lists column names
#-- input --#
# $dbh = DBI object
# $tbl = sql table of interest
# $silent_ret = no verbose & exit; return column names
	my ($dbh, $tbl, $silent_ret) = @_;
	my $all = $dbh->selectall_arrayref("pragma table_info($tbl)");

	my %tmp;
	foreach (@$all){ 
		$$_[1] =~ tr/A-Z/a-z/;		# lower case for matching
		$tmp{$$_[1]} = 1; 
		}
		
	if(defined $silent_ret){ return \%tmp; }
	else{  print "Columns:\n", join(",\n", keys %tmp), "\n\n";  exit; }
	}

sub get_arrays_seq_byLeader{
# getting array sequences and orienting by leaders #
	#my ($dbh, $spacer_bool, $extra_query, $join_sql) = @_;
	my ($dbh, $opts_r) = @_;
	
	# checking for opts #
	map{ confess "ERROR: cannot find option: '$_'" 
		unless exists $opts_r->{$_} } qw/spacer_DR_b extra_query join_sql/;

	# getting table info #
	my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
	($tbl_oi, $tbl_prefix) = ("DRs","DR") if defined $opts_r->{"spacer_DR_b"};

	
	my $query = "SELECT
$tbl_oi.Locus_ID, 
'$tbl_prefix', 
$tbl_oi.$tbl_prefix\_ID,
$tbl_oi.$tbl_prefix\_sequence,
loci.array_end,
leaders.leader_end
FROM loci, leaders, $tbl_oi
WHERE loci.locus_id = $tbl_oi.locus_id
AND loci.locus_id = leaders.locus_id
AND leaders.locus_id IS NOT NULL
";	

	$query .= "GROUP BY $tbl_oi.$tbl_prefix\_group" if $opts_r->{"by_group"};
	$query =~ s/\n/ /g;
	
	#print Dumper $query; exit;
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	confess " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	# making fasta #
	my @arrays;
	foreach my $row (@$ret){
		
		# loading fasta #
		confess " ERROR: not spacer/repeat group found!\n\tWas CLdb_groupArrayElements.pl run?\n"
			unless $$row[1]; 
		
		# rev-comp if leader end > array end #
		if ($$row[6] > $$row[5]){
			$$row[4] = revcomp($$row[4]);
			}
	
		@$row = @$row[0..4];

		# 'NA' if needed #
		if(defined $$row[3] && defined $opts_r->{"by_group"}){		# group found 
			$$row[0] = 'NA';
			$$row[2] = 'NA';
			}
		else{
			$$row[3] = 'NA' if ! defined $$row[3];
			confess "ERROR: no group ID found for an entry. Cannot group sequences!\n"
				if defined $opts_r->{"by_group"};
			}
		}
		#print Dumper @$ret; exit;
	return $ret;
	}

sub get_array_seq{
# getting spacer or DR sequence from either table #
#-- input --#
# $dbh = DBI database object
# $opts_r = hash of options 
#-- options --#
# spacer_DR_b = spacer or DR [spacer]
# extra_query = extra sql
# by_cluster	= spacer/DR group (undef|1)? [undef]
# 	cutoff = sequenceID cutoff [1]
# 	strand_spec = strand specific (0|1)? [0]
# join_sql = "AND" statements 
	
	my ($dbh, $opts_r) = @_;
	
	# checking for opts #
	map{ confess "ERROR: cannot find option: '$_'" 
		unless exists $opts_r->{$_} } qw/spacer_DR_b extra_query join_sql/;
	
	# getting table info #
	my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
	($tbl_oi, $tbl_prefix) = ("DRs","DR") if $opts_r->{"spacer_DR_b"};		# DR instead of spacer

# columns:
#	locus_ID
#	spacer|DR
#	spacer/DR_ID
#	start
#	end
# 	clusterID
# 	sequence

my $query;
if(defined $opts_r->{"by_cluster"}){		# by group
	# by group default options #
	$opts_r->{"cutoff"} = 1 unless exists $opts_r->{'cutoff'}; 	
	$opts_r->{"strand_spec"} = 1 unless exists $opts_r->{'strand_spec'}; 
	
	$query = "SELECT 
'NA',
'$tbl_prefix', 
'NA',
'NA',
'NA',
$tbl_prefix\_clusters.cluster_ID,
$tbl_oi.$tbl_prefix\_sequence
FROM $tbl_oi, $tbl_prefix\_clusters, loci 
WHERE loci.locus_id = $tbl_oi.locus_id 
AND $tbl_prefix\_clusters.$tbl_prefix\_ID = $tbl_oi.$tbl_prefix\_ID
AND $tbl_prefix\_clusters.cutoff = $opts_r->{'cutoff'}
AND $tbl_prefix\_clusters.strand_spec = $opts_r->{'cutoff'}
GROUP BY $tbl_prefix\_clusters.cluster_ID
ORDER BY $tbl_prefix\_clusters.cluster_ID
$opts_r->{'join_sql'}";
	}
else{		# not by clusters
	$query = "SELECT
$tbl_oi.Locus_ID, 
'$tbl_prefix', 
$tbl_oi.$tbl_prefix\_ID,
$tbl_oi.$tbl_prefix\_start,
$tbl_oi.$tbl_prefix\_end,
'NA',
$tbl_oi.$tbl_prefix\_sequence
FROM $tbl_oi, loci 
WHERE loci.locus_id = $tbl_oi.locus_id $opts_r->{'join_sql'}";
	}
	

	$query =~ s/[\n\t]+/ /g;
	$query = join(" ", $query, $opts_r->{"extra_query"});
	
		#print Dumper $query;
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	confess "ERROR: no matching entries!\n"
		unless $$ret[0];

		#print Dumper $ret; exit;

	# 'NA' for ID & locus if cluster defined #
	#foreach my $row (@$ret){
		# @$row = locus,spacer/DR,ID,groupID,seq
	#	if(defined $$row[3] && defined $opts_r->{"by_cluster"}){		# cluster found 
	#		$$row[0] = 'NA';
	#		$$row[2] = 'NA';
	#		}
	#	else{
	#		$$row[3] = 'NA' if ! defined $$row[3];
	#		confess "ERROR: no group ID found for an entry. Cannot group sequences!\n"
	#			if defined $opts_r->{"by_cluster"};
	#		}
	#	}

		#print Dumper @$ret; exit;
	return $ret;
	}

sub join_query_opts{
# joining query options for with 'AND' 
	my ($vals_r, $cat) = @_;
	
	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}

sub table_exists {
# checking for the existence of a table #
	my ($dbh, $table) = @_;
	confess "ERROR: Provide a DBI database object!\n" if ! defined $dbh;
	confess "ERROR: Provide a CLdb table name!\n" if ! defined $table;
	
	my $res = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", "tbl_name"); 
	
	confess "ERROR: '$table' table not found in CLdb!\n" 
		unless grep(/^$table$/i, keys %$res);
	}

sub n_entries {
# getting number of entries in a table #
	my ($dbh, $table) = @_;
	confess "ERROR: Provide a DBI database object!\n" if ! defined $dbh;
	confess "ERROR: Provide a CLdb table name!\n" if ! defined $table;
	
	my $q = "SELECT count(*) FROM $table";
	my $res = $dbh->selectall_arrayref($q);

	return $$res[0][0];
	}

sub get_leader_pos{
# getting spacer or DR sequence from either table #
#-- input --#
# $dbh = DBI database object
# $opts_r = hash of options 
#-- options --#
# extra_query = extra sql
# join_sql = "AND" statements 
	
	my ($dbh, $opts_r) = @_;
	
	my $query = "SELECT leaders.locus_ID, leaders.leader_start, leaders.leader_end
		FROM leaders, loci WHERE loci.locus_id = leaders.locus_id $opts_r->{'join_sql'}";
	
	$query =~ s/[\n\t]+/ /g;
	$query = join(" ", $query, $opts_r->{"extra_query"});
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	confess " ERROR: no matching entries!\n"
		unless $$ret[0];

	# making hash of sequences #
	my %leaders;
	foreach my $row (@$ret){
		$leaders{$$row[0]}{"start"} = $$row[1];
		$leaders{$$row[0]}{"end"} = $$row[2];
		}

		#print Dumper %leaders; exit;
	return \%leaders;	
	}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::query


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=CRISPR_db>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/CRISPR_db>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/CRISPR_db>

=item * Search CPAN

L<http://search.cpan.org/dist/CRISPR_db/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2013 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.

=cut

1; 
