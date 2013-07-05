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
		
		my @line = split /\t/;
		$line[$group_column - 1] =~ s/^group//i;
		die " ERROR: group column naming not reconcnized!\n"
			unless $line[$group_column - 1] =~ /^\d+$/;
		push @groups, $line[$group_column - 1];
		}
	die " ERROR: provide at least 1 group!\n"
		unless @groups;
	
		#print Dumper @groups; exit;
	return \@groups;
	}

	
sub get_arrays_join{
	my ($dbh, $groups_r, $spacer_bool, $extra_query, $join_sql) = @_;
	
	# make query #
	my $group_str = join(",", @$groups_r);
	my $query;
	if($spacer_bool){		# direct repeat
		$query = "SELECT b.repeat_group, b.repeat_ID, a.Locus_ID, a.Taxon_name, a.Taxon_id, a.Subtype FROM loci a, directrepeats b WHERE a.locus_id = b.locus_id $join_sql AND Repeat_group = ?";
		}
	else{					# spacer
		$query = "SELECT b.spacer_group, b.spacer_ID, a.Locus_ID, a.Taxon_name, a.Taxon_id, a.Subtype FROM loci a, spacers b WHERE a.locus_id = b.locus_id $join_sql AND Spacer_group = ?";
		}
	$query = join(" ", $query, $extra_query);

	print STDERR $query, "\n" if $verbose;

	# preparing query #
	my $sql = $dbh->prepare($query);
	
	# querying #
	foreach my $group (@$groups_r){
		$sql->execute($group);
		while(my @row = $sql->fetchrow_array()){
			$row[0] =~ s/^/Group/;
			$row[2] =~ s/^/cli./;
			map{$_ = "" unless $_} @row;
			print join("\t", @row), "\n";
			}
		}

	}
	
sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND a.$cat IN (", join(", ", @$vals_r), ")");
	}



__END__

=pod

=head1 NAME

CLdb_arrayGroup2Info.pl -- get Info on spacer or repeat groups

=head1 SYNOPSIS

CLdb_arrayGroup2Info.pl [flags] < groups.txt > array.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -repeat

Get direct repeats instead of spacers.

=item -column

The group column in the input table (index by 1). [1]

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

perldoc CLdb_arrayGroup2Info.pl

=head1 DESCRIPTION

Get info for spacer or direct-repeat groups.
For instance, getting all spacers for spacer groups
that blasted to a genome. 

The input is a list of group IDs (integers or 'GroupINTEGER').

A table is returned with the following columns:

=over 

=item * Spacer/Repeat_group

=item * Spacer/Repeat_ID

=item * Locus_ID

=item * Taxon_name

=item * Taxon_id

=item * Subtype 

=back 

=head1 EXAMPLES

=head2 Get all info for Spacer_group 1

echo 1 | CLdb_arrayGroup2Info.pl -d CRISPR.sqlite 

=head2 Get all info for a list of repeat groups

CLdb_arrayGroup2Info.pl -d CRISPR.sqlite -r < repeat_groups.txt

=head2 Get just subtype 'I-B' info for a list of spacer groups

CLdb_arrayGroup2Info.pl -d CRISPR.sqlite -sub I-B < spacer_groups.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

