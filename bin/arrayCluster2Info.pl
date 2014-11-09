#!/usr/bin/env perl

=pod

=head1 NAME

arrayCluster2Info.pl -- get Info on spacer or repeat groups

=head1 SYNOPSIS

arrayCluster2Info.pl [flags] < groups.txt > array_info.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -repeat  <bool>

Get direct repeats instead of spacers. [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb_perldoc arrayCluster2Info.pl

=head1 DESCRIPTION

Get info for spacer or direct-repeat clusters.

For instance, getting info on
all individual spacers in each spacer cluster
that was BLASTed against >=1 genome. 

The input is any file with spacer/DR cluster IDs (e.g. '1_1.00_72').

A table is returned with the following columns:

=over 

=item * Spacer/Repeat_cluster

=item * Spacer/Repeat_ID

=item * Locus_ID

=item * Taxon_name

=item * Taxon_id

=item * Subtype 

=back 

=head1 EXAMPLES

=head2 Get all info for Spacer_group 1

echo 1 | arrayCluster2Info.pl -d CLdb.sqlite < groups.txt > array_info.txt

=head2 Get all info for a list of repeat groups

arrayCluster2Info.pl -d CLdb.sqlite -r < repeat_groups.txt > array_info.txt

=head2 Get just subtype 'I-B' info for a list of spacer groups > array_info.txt

arrayCluster2Info.pl -d CLdb.sqlite -sub I-B < spacer_groups.txt > array_info.txt

=head2 Get info for table of BLAST hits where spacer groups were the query

arrayCluster2Info.pl -d CLdb.sqlite < spacer_groups_blastn.txt > array_info.txt

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


my ($verbose, $database_file, $spacer_bool);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $group_column = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_bool,			 # spacers or repeats? [spacers]
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "column=i" => \$group_column, 		# column containing groupID
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# loading groups #
my $groups_r = load_groups($group_column);

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting arrays of interest from database #
my $arrays_r = get_arrays_join($dbh, $groups_r, $spacer_bool, $extra_query, $join_sql);


# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub load_groups{
# loading groups from STDIN #
	my ($group_column) = @_;
	
	my @groups;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		if( /[01]_[0-9.]+_[0-9]+/ ){		# group ID
			s/.+([01]_[0-9.]+_[0-9]+).+/$1/;
			push @groups, $_;
			}
		
		#my @line = split /\t/;
		#$line[$group_column - 1] =~ s/^group\.//i;
		#die " ERROR: group column naming not reconcnized!\n"
		#	unless $line[$group_column - 1] =~ /^\d+$/;
		#push @groups, $line[$group_column - 1];
		}
	die " ERROR: provide at least 1 group!\n"
		unless @groups;
	
		#print Dumper @groups; exit;
	return \@groups;
	}

sub get_arrays_join{
	my ($dbh, $groups_r, $spacer_bool, $extra_query, $join_sql) = @_;
	
	# make query #
	my $cat;
	if($spacer_bool){	$cat = 'DR'; }		# DR
	else{ $cat = 'spacer'; }
	my $cats = $cat . 's';
	
	my $query = "SELECT b.cluster_ID, b.$cat\_ID, a.Locus_ID, a.Taxon_name, a.Taxon_id, a.Subtype 
		FROM loci a, $cat\_clusters b
		WHERE a.locus_id = b.locus_id
		AND b.cluster_ID = ?
		$join_sql";
	$query = join(" ", $query, $extra_query);
	$query =~ s/\s+/ /g;

	print STDERR $query, "\n" if $verbose;

	# preparing query #
	my $sql = $dbh->prepare($query);
	
	# querying #
	foreach my $group (@$groups_r){
		$sql->execute($group);
		while(my @row = $sql->fetchrow_array()){
			map{$_ = "" unless $_} @row;
			print join("\t", @row), "\n";
			}
		}

	}
	
sub join_query_opts_OLD{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND a.$cat IN (", join(", ", @$vals_r), ")");
	}



