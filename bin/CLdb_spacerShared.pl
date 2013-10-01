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


my ($verbose, $database_file, $spacer_bool, $by_group, $leader_bool);
my ($subtype, $taxon_id, $taxon_name, $locus_id);
my $extra_query = "";
my $cutoff = 1;				# cd-hit cutoff of 1
GetOptions(
	   "database=s" => \$database_file,
	   "subtype" => \$subtype,
	   "id" => \$taxon_id,
	   "name" => \$taxon_name,
	   "locus" => \$locus_id,				# by locus_id
	   "cutoff=f" => \$cutoff,
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

# listing existing tables #
my $table_list_r = list_tables($dbh);
$table_list_r = entry_count($dbh, $table_list_r);

# making group by sql
my $group_by_r = make_group_by_sql($subtype, $taxon_id, $taxon_name, $locus_id);


# getting arrays of interest from database #
my ($arrays_r, $groups_r);
if( exists $table_list_r->{"spacer_hclust"} && $table_list_r->{"spacer_hclust"} > 0){		
	($arrays_r, $groups_r) = get_arrays_join_hclust($dbh, $extra_query, $group_by_r);
	}
else{
	($arrays_r, $groups_r) = get_arrays_join($dbh, $extra_query, $group_by_r);
	}

# write shared matrix (cluster ~ group)#
write_shared_matrix($arrays_r, $groups_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_shared_matrix{
# writing matrix of shared spacers (based on spacer groups/clusters) #
	my ($arrays_r, $groups_r) = @_;		# cluster_ID => grouping_ID => count
	
	my @groups = sort{$a cmp $b} keys %$groups_r;
	
	# header #
	print join("\t", "Spacer_cluster", @groups), "\n"; 
	
	# body #
	foreach my $clust (sort{$a<=>$b} keys %$arrays_r){
		my @row = ("$clust");
		foreach my $group ( @groups ){
			if(exists $arrays_r->{$clust}{$group}){
				push @row, $arrays_r->{$clust}{$group};
				}
			else{
				push @row, 0;				# count of zero for cluster => group
				}
			}
		print join("\t", @row), "\n";		# writing line of matrix
		}
	}

sub get_arrays_join_hclust{
	my ($dbh, $extra_query, $group_by_r) = @_;
	
	my $sel = join(",", qw/c.cluster_ID count(c.cluster_ID)/, @$group_by_r);
	my $group_by = join(",", @$group_by_r, "c.Cluster_ID");
	
	# make query #
	my $q = "SELECT 
$sel
FROM loci a, spacers b, spacer_hclust c
WHERE a.locus_id = b.locus_id
AND b.locus_id = c.locus_id
AND b.spacer_id = c.spacer_ID
AND c.cutoff = $cutoff
$extra_query
GROUP BY $group_by";
	$q =~ s/\n/ /g;
	
	#print Dumper $q; exit;
	
	# status #
	print STDERR "$q\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($q);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %arrays;
	my %groups;
	foreach my $row (@$ret){
		die " ERROR: not matching entries!\n"
			unless $$row[0]; 
		my $id;
		if(! @$group_by_r){			# just total for each cluster
			$id = "Total";
			}
		elsif(! $$row[2]){			# group_by ID
			$id = "NULL";
			}
		else{
			$id = join("__", @$row[2..$#$row]);			
			}
		$arrays{$$row[0]}{$id} = $$row[1];
		$groups{$id} = 1;
		}
	
		#print Dumper %arrays; exit;
	return \%arrays, \%groups;
	}

sub get_arrays_join{
	my ($dbh, $extra_query, $group_by_r) = @_;
	
	my $sel = join(",", qw/b.spacer_group count(b.spacer_group)/, @$group_by_r);
	my $group_by = join(",", @$group_by_r, "b.spacer_group");
	
	# make query #
	my $q = "SELECT 
$sel
FROM loci a, spacers b
WHERE a.locus_id = b.locus_id
$extra_query
GROUP BY $group_by";
	$q =~ s/\n/ /g;
	
	#print Dumper $q; exit;
	
	# status #
	print STDERR "$q\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($q);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %arrays;
	my %groups;
	foreach my $row (@$ret){
		die " ERROR: not matching entries!\n"
			unless $$row[0]; 
		my $id;
		if(! @$group_by_r){			# just total for each cluster
			$id = "Total";
			}
		elsif(! $$row[2]){			# group_by ID
			$id = "NULL";
			}
		else{
			$id = join("__", @$row[2..$#$row]);			
			}
		$arrays{$$row[0]}{$id} = $$row[1];
		$groups{$id} = 1;
		}
	
		#print Dumper %arrays; exit;
	return \%arrays, \%groups;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}

sub make_group_by_sql{
	my ($subtype, $taxon_id, $taxon_name, $locus_id) = @_;
	my @group_by;
	push @group_by, "a.subtype" if $subtype;
	push @group_by, "a.taxon_id" if $taxon_id;
	push @group_by,"a.taxon_name" if $taxon_name;
	push @group_by, "a.locus_id" if $locus_id;
	return \@group_by;
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


__END__

=pod

=head1 NAME

CLdb_spacerShared.pl -- write a matrix of spacers shared among taxa, subtypes, and/or loci

=head1 SYNOPSIS

CLdb_spacerShared.pl [flags] > shared.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <bool>

Group count data by subtype? [FALSE]

=item -id  <bool> 

Group count data by taxon_id? [FALSE]

=item -name  <bool>

Group count data by taxon_name? [FALSE]

=item -locus  <bool>

Group count data by locus_id? [FALSE] 

=item -cutoff  <float>

Which Spacer/DR clustering cutoffs to summarize (>= 1 argument)? [1]

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerShared.pl

=head1 DESCRIPTION

How similar are CRISPRs in spacer content?

Get a table showing the number of a particular
spacer (by default, exact same sequence, change with
'-cutoff') among subtypes, taxa, or both.

By default, total counts of each spacer cluster are written
(no parsing by group).

WARNING: '-cutoff' can only be used if CLdb_hclusterArrays.pl
has been run first.

=head1 EXAMPLES

=head2 Write table of spacers shared among subtypes

CLdb_spacerShared.pl -d CLdb.sqlite -s > subtype_shared.txt

=head2 Write table of spacers shared among taxa (by taxon_name)

CLdb_spacerShared.pl -d CLdb.sqlite -n > taxa_shared.txt

=head2 Write table of spacers shared among subtypes/taxa

CLdb_spacerShared.pl -d CLdb.sqlite -s -n > subtype_taxa_shared.txt

=head2 Write table of spacers shared among loci (names include taxon_name & subtype)

CLdb_spacerShared.pl -d CLdb.sqlite -l -s -n > loci_subtype-taxa_shared.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

