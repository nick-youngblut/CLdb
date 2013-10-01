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
my (@subtype, @taxon_id, @taxon_name);
my $kill;
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "kill=s" => \$kill,
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

# joining query 
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# checking for mass deletion #
check_deletion($join_sql, $extra_query) unless $kill && $kill =~ /^y/i;

# database metadata #
my $loci_r = get_loci_entries($dbh, $join_sql, $extra_query);

# writing out which loci to be deleted #
write_selected_loci($loci_r) unless $verbose;

# deleting all associated entries #
## Spacers *
## DirectRepeats *
## DirectRepeatConsensus *
## LeaderSeqs *
## Draft *
## Genes *
## directrepeat_hclust
## spacer_hclust
## blast_hits
## blast_subject
## spacer_pairwise_blast
	# '*' = simple join on locus_id

# listing tables to delete #
my $tables_r = list_tables($dbh);

# need to join; get unique IDs and delete on unique IDs #
## hclust tables ##
delete_hclust($dbh, $loci_r, "spacers", "spacer_id", "spacer_hclust") 
	if exists $tables_r->{"spacers_hclust"};
delete_hclust($dbh, $loci_r, "DRs", "DR_id", "DR_hclust")
	if exists $tables_r->{"DR_hclust"};

## blast tables ##
delete_spacer_pairwise_blast($dbh, $loci_r) if exists $tables_r->{"spacer_pairwise_blast"};


# deleting tables without needing to join; locus_id present #
simple_delete($dbh, $loci_r, "spacers") if exists $tables_r->{"spacers"};
simple_delete($dbh, $loci_r, "drs") if exists $tables_r->{"drs"};
simple_delete($dbh, $loci_r, "dr_consensus") if exists $tables_r->{"dr_consensus"};
simple_delete($dbh, $loci_r, "leaders") if exists $tables_r->{"leaders"};
#simple_delete($dbh, $loci_r, "draft") if exists $tables_r->{"draft"};
simple_delete($dbh, $loci_r, "genes") if exists $tables_r->{"genes"};
simple_delete($dbh, $loci_r, "loci") if exists $tables_r->{"loci"};
simple_delete($dbh, $loci_r, "spacer_hclust") if exists $tables_r->{"spacer_hclust"};


$dbh->commit();
exit;

### Subroutines
sub list_tables{
	my $dbh = shift;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	map{delete $all->{$_}; tr/A-Z/a-z/; $all->{$_} = 1} keys %$all;
		#print Dumper %$all; exit;
	return $all;
	}

sub delete_spacer_pairwise_blast{
	my ($dbh, $loci_r) = @_;
	
	# deleting blast hits by locus_id#
	my $cmd = "DELETE FROM spacer_pairwise_blast where query_locus_ID = ? OR subject_locus_id = ?";
	my $sql = $dbh->prepare($cmd);

	my $delete_cnt = 0;
	foreach my $locus (@$loci_r){
		$sql->bind_param(1,$$locus[0]);
		$sql->bind_param(2,$$locus[0]);		
		$sql->execute();
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for locus_id: $$locus[0]\n",
			}
		else{ $delete_cnt++; }
		}
	
	# status #
	if($delete_cnt){
		print STDERR "...Entries deleted from: spacer_pairwise_blast\n"
			unless $verbose;
		}
	else{
		print STDERR "...No entries deleted from: spacer_pairwise_blast!\n"
			unless $verbose;
		}	
	}

sub delete_hclust{
# deleting spacer hclust entries #
	my ($dbh, $loci_r, $table, $id, $del_table) = @_;
	
	# getting spacer/repeat_IDs #
	my $cmd = "SELECT $table.$id FROM $table, loci
WHERE $table.locus_id = loci.locus_id
AND loci.locus_id = ?";
	$cmd =~ s/[\n\r]/ /g;
	
	my $sql = $dbh->prepare($cmd);
	
	my @ret;
	foreach my $row (@$loci_r){
		$sql->bind_param(1, $$row[0]);
		$sql->execute();
		push(@ret, $sql->fetchall_arrayref());
		}
	
	# deleting spacer_IDs from spacer_hclust #
	$cmd = "DELETE FROM $del_table where $id = ?";
	$sql = $dbh->prepare($cmd);

	my $delete_cnt = 0;
	foreach my $locus (@ret){
		foreach my $row (@$locus){
			$sql->bind_param(1,$$row[0]);
			$sql->execute();

			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr for $id: $$row[0]\n",
				}
			else{ $delete_cnt++; }
			}
		}
	
	# status #
	if($delete_cnt){
		print STDERR "...Entries deleted from: $del_table\n"
			unless $verbose;
		}
	else{
		print STDERR "...No entries deleted from: $del_table!\n"
			unless $verbose;
		}
	}

sub simple_delete{
# deleting based on locus_id
	my ($dbh, $loci_r, $table) = @_;
	
	my $cmd = "DELETE FROM $table where locus_id = ?";
	my $sql = $dbh->prepare($cmd);
	
	my $delete_cnt = 0;
	foreach my $row (@$loci_r){
		$sql->bind_param(1, $$row[0]);
		$sql->execute();
	
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for locus_ID: $$row[0]\n",
			}
		else{ $delete_cnt++; }
		}

	# status #
	if($delete_cnt){
		print STDERR "...Entries deleted from: $table\n"
			unless $verbose;
		}
	else{
		print STDERR "...No entries deleted from: $table!\n"
			unless $verbose;
		}
	}

sub write_selected_loci{
# writing selected loci to delete #
	my ($loci_r) = @_;
	
	print STDERR "### Loci to be deleted ###\n";
	print STDERR join("\t", qw/Locus_ID Taxon_name Taxon_ID/), "\n";
	map{print STDERR join("\t", join(".", "cli", $$_[0]),  @$_[1..$#$_]), "\n"} @$loci_r;
	print STDERR "\n";
	}

sub get_loci_entries{
# getting loci entries that match  query #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	my $cmd = "SELECT locus_id, taxon_name, taxon_id from loci where locus_id=locus_id $join_sql $extra_query";
	my $ret = $dbh->selectall_arrayref($cmd);

		#print Dumper $ret; exit;
	return $ret;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND $cat IN (", join(", ", @$vals_r), ")");
	}

sub check_deletion{
# deleting all of a taxon, subtype, or some other match #
	my ($join_sql, $extra_query) = @_;
	
	unless($join_sql || $extra_query){
		die " ERROR: no specified query provided! All entries will be deleted! Use '-k Y' if you are sure!\n";
		}
	}



__END__

=pod

=head1 NAME

CLdb_delete.pl -- deleting loci-associated entries throughout database

=head1 SYNOPSIS

CLdb_delete.pl [options]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine query.

=item -kill  <bool>

Will override the script's attempts to prevent the deletion of all loci. [FALSE]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_delete.pl

=head1 DESCRIPTION

Select loci that should be deleted.
All entries in all tables that
are associated with the selected loci.

=head1 EXAMPLES

=head2 Delete all loci in a genome 

CLdb_delete.pl -d CLdb.sqlite -taxon_id 6666666.403

=head2 Delete locus_ids 1 through 4 

CLdb_delete.pl -d CLdb.sqlite -q "AND locus_id IN (1,2,3,4)

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

