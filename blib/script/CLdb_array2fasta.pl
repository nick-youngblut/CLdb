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


my ($verbose, $database_file, $spacer_bool, $by_group);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_bool,			 # spacers or repeats?
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "group" => \$by_group,
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

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting arrays of interest from database #
#if($join_sql){		# table-join query
my $arrays_r = get_arrays_join($dbh, $spacer_bool, $extra_query, $join_sql);
#	}
#else{				# simple query of all spacers/DR
#	$arrays_r = get_arrays($dbh, $spacer_bool, $extra_query);
#	}

# writing fasta #
write_arrays_fasta($arrays_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_arrays_fasta{
# writing arrays as fasta
	my ($arrays_r) = @_;
	
	foreach my $locus_id (keys %$arrays_r){
		foreach my $x_id (keys %{$arrays_r->{$locus_id}}){
			if($by_group){
				print join("\n", ">Group$x_id", $arrays_r->{$locus_id}{$x_id}), "\n";
				}
			else{
				print join("\n", ">cli.$locus_id\__$x_id", $arrays_r->{$locus_id}{$x_id}), "\n";
				}
			}
		}
	}

sub get_arrays{
	my ($dbh, $spacer_bool, $extra_query) = @_;
	
	# make query #
	my $query;
	if($spacer_bool){		# direct repeat
		if($by_group){
			$query = "SELECT Locus_ID, repeat_group, Repeat_sequence FROM directrepeats GROUP BY repeat_group";
			}
		else{
			$query = "SELECT Locus_ID, Repeat_ID, Repeat_sequence FROM directrepeats";
			}
		}
	else{					# spacer
		if($by_group){
			$query .= "SELECT Locus_ID, spacer_group, Spacer_sequence FROM spacers GROUP BY spacer_group";		
			}
		else{
			$query = "SELECT Locus_ID, Spacer_ID, Spacer_sequence FROM spacers";
			}
		}
	$query = join(" ", $query, $extra_query);
	print STDERR $query, "\n" if $verbose;
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	
	my %arrays;
	foreach my $row (@$ret){
		die " ERROR: not spacer/repeat group found!\n\tWas CLdb_groupArrayElements.pl run?\n"
			unless $$row[1]; 
		$arrays{$$row[0]}{$$row[1]} = $$row[2];
		}
	
	#	print Dumper %arrays; exit;
	return \%arrays;
	}
	
sub get_arrays_join{
	my ($dbh, $spacer_bool, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query;
	if($spacer_bool){		# direct repeat
		if($by_group){
			$query = "SELECT directrepeats.Locus_ID, directrepeats.Repeat_group, directrepeats.Repeat_sequence FROM directrepeats, loci WHERE loci.locus_id = directrepeats.locus_id $join_sql GROUP BY directrepeats.repeat_group";
			}
		else{
			$query = "SELECT directrepeats.Locus_ID, directrepeats.Repeat_ID, directrepeats.Repeat_sequence FROM directrepeats, loci WHERE directrepeats.locus_id = loci.locus_id $join_sql";
			}
		}
	else{					# spacer
		if($by_group){
			$query = "SELECT spacers.Locus_ID, spacers.Spacer_group, spacers.spacer_sequence FROM spacers, loci WHERE spacers.locus_id = loci.locus_id $join_sql GROUP BY spacers.spacer_group";
			}
		else{
			$query = "SELECT spacers.Locus_ID, spacers.spacer_ID, spacers.spacer_sequence FROM spacers, loci WHERE spacers.locus_id = loci.locus_id $join_sql";
			}
		}
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %arrays;
	foreach my $row (@$ret){
		$arrays{$$row[0]}{$$row[1]} = $$row[2];
		}
	
	#	print Dumper %arrays; exit;
	return \%arrays;
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

CLdb_array2fasta.pl -- write CRISPR array spacers or direct repeats to fasta

=head1 SYNOPSIS

CLdb_array2fasta.pl [flags] > array.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -repeat

Get direct repeats instead of spacers.

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -group

Get array elements de-replicated by group (ie. all unique sequences).

=item -query

Extra sql to refine which sequences are returned.

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_array2fasta.pl

=head1 DESCRIPTION

Get spacer or direct repeat sequences from the CRISPR database
and write them to a fasta.

By default, all spacers or direct repeats (if '-r') will be written.
The '-q' flag can be used to refine the query to certain sequences (see examples).

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite 

=head2 Write all direct repeats to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite -r

=head2 Write all unique spacers

CLdb_array2fasta.pl -d CLdb.sqlite -g

=head2 Refine spacer sequence query:

CLdb_array2fasta.pl -d CLdb.sqlite -q "AND loci.Locus_ID=1" 

=head2 Refine spacer query to a specific subtype & 2 taxon_id's

CLdb_array2fasta.pl -d CLdb.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

