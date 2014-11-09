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


my ($verbose, $database_file, $gap_rm);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "query=s" => \$extra_query,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name, 
	   "gaps" => \$gap_rm,
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


# getting leaderss of interest from database #
my $leaders_r;
if($join_sql){
	$leaders_r = get_leaders_join($dbh, $extra_query, $join_sql);
	}
else{
	$leaders_r = get_leaders($dbh, $extra_query);
	}

# writing fasta #
write_leaders_fasta($leaders_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_leaders_fasta{
# writing arrays as fasta
	my ($leaders_r) = @_;
	
	foreach my $locus_id (keys %$leaders_r){
		print join("\n", ">cli.$locus_id", $leaders_r->{$locus_id}), "\n";
		}
	}

sub get_leaders{
	my ($dbh, $extra_query) = @_;
	
	# make query #
	my $query = "SELECT Locus_ID, Leader_sequence FROM Leaders";
	$query = join(" ", $query, $extra_query);
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %leaders;
	foreach my $row (@$ret){
		$$row[1] =~ s/-//g if $gap_rm; 		# removing gaps
		$leaders{$$row[0]}= $$row[1];
		}
	
	#	print Dumper %leaders; exit;
	return \%leaders;
	}
	
sub get_leaders_join{
	my ($dbh, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query = "SELECT a.Locus_ID, a.Leader_sequence FROM Leaders a, Loci b WHERE a.locus_id = b.locus_id $join_sql";
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %leaders;
	foreach my $row (@$ret){
		$$row[1] =~ s/-//g if $gap_rm; 		# removing gaps
		$leaders{$$row[0]} = $$row[1];
		}
	
	#	print Dumper %leaders; exit;
	return \%leaders;
	}
	
sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND b.$cat IN (", join(", ", @$vals_r), ")");
	}



#my $query = "SELECT FROM locus_id";

__END__

=pod

=head1 NAME

CLdb_leader2fasta.pl -- write CRISPR leader sequences in fasta format

=head1 SYNOPSIS

CLdb_leader2fasta.pl [flags] > leaders.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -gap  <bool>

Remove gaps from the sequences? [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_leader2fasta.pl

=head1 DESCRIPTION

Get leader sequences from the CRISPR database
and write them in fasta format.

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

CLdb_leader2fasta.pl -data CRISPR.sqlite 

=head2 Write all direct repeats to a fasta:

CLdb_leader2fasta.pl -data CRISPR.sqlite -r

=head2 Refine spacer sequence query:

CLdb_leader2fasta.pl -data CRISPR.sqlite -q "where LOCUS_ID=1" 

=head2 Refine spacer query to a specific subtype & 2 taxon_id's

CLdb_leader2fasta.pl -da CRISPR.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

