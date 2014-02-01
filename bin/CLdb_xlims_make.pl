#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_xlims_make.pl -- making xlims table for plotting

=head1 SYNOPSIS

CLdb_xlims_make.pl [flags] > xlims.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -locus_id  <char>

Refine query to specific a locus_id(s) (>1 argument allowed).

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

perldoc CLdb_xlims_make.pl

=head1 DESCRIPTION

Make a basic xlims table needed for plotting.

=head1 EXAMPLES

=head2 Plotting all loci classified as subtype 'I-A'

CLdb_xlims_make.pl -d CLdb.sqlite -sub I-A 

=head2 No broken loci

CLdb_xlims_make.pl -da CLdb.sqlite -sub I-A -q "AND loci.operon_status != 'broken'"

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
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $CLdb_sqlite, @ITEP_sqlite);
my (@locus_id, @subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $spacer_cutoff = 1;
my $xlim_out = "xlims.txt";
GetOptions(
	   "database=s" => \$CLdb_sqlite,
	   "ITEP=s{,}" => \@ITEP_sqlite,
	   "locus_id=s{,}" => \@locus_id,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "cutoff=f" =>  \$spacer_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# checking CLdb #
file_exists($CLdb_sqlite, "database");

#--- MAIN ---#
# connect 2 CLdb #
my $dbh = connect2db($CLdb_sqlite);

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@locus_id, "locus_id");
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting loci start-end #
my $xlims_r = load_xlims($dbh, $join_sql, $extra_query);

# mutli-loci / multi-subtype problem
## determining if multiple loci/subtypes per taxon_name ##
my ($multi_loci, $multi_subtype)  = check_multi($xlims_r);

## adding loci & subtypes to DNA_segs names ##
edit_xlims_taxon_name($xlims_r, $multi_loci, $multi_subtype);

# writing table #
write_xlims($xlims_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_xlims{
	my ($xlims_r) = @_;
	
	# header #
	print join("\t", qw/start end taxon_name locus_id subtype dna_segs_id/), "\n";
	
	# body #
	foreach my $taxon_name (keys %$xlims_r){
		foreach my $locus_id (keys %{$xlims_r->{$taxon_name}}){
			print join("\t", @{$xlims_r->{$taxon_name}{$locus_id}{"entry"}}), "\n";
			}
		}
	}

sub edit_xlims_taxon_name{
	my ($xlims_r, $multi_loci, $multi_subtype) = @_;
	
	foreach my $taxon_name (keys %$xlims_r){
		foreach my $locus_id (keys %{$xlims_r->{$taxon_name}}){
			# editing taxon_name in row #
			my $new_taxon_name = $taxon_name;
			$new_taxon_name = join("__", $taxon_name, $locus_id)
					if $multi_loci;
			$new_taxon_name = join("__", $new_taxon_name, $xlims_r->{$taxon_name}{$locus_id}{"subtype"})
					if $multi_subtype;
			
			push(@{$xlims_r->{$taxon_name}{$locus_id}{"entry"}}, $new_taxon_name);
			}
		}
		#print Dumper %$xlims_r; exit;
	}

sub check_multi{
# checking for multiple entries per taxon #
	my ($xlims_r) = @_;
	
	my $multi_loci = 0;				# mutliple loci per taxon_name
	my $multi_subtype = 0;			# multiple subtypes total 
	my %subtype_sum; 
	foreach my $taxon_name (keys %$xlims_r){
		$multi_loci = 1 if scalar keys %{$xlims_r->{$taxon_name}} > 1;
			
		foreach my $locus_id (keys %{$xlims_r->{$taxon_name}} ){
			# sanity check #
			die " ERROR: cannot find subtype for $taxon_name -> $locus_id!\n"
				unless exists  $xlims_r->{$taxon_name}{$locus_id}{"subtype"};
			
			$subtype_sum{ $xlims_r->{$taxon_name}{$locus_id}{"subtype"} }++;
			}
		}
	$multi_subtype = 1 if scalar keys %subtype_sum > 1;

	# status #
	print STDERR "...Found multiple loci for 1 or more taxa. Adding leaves to the tree! Adding loci_ids to leaves & xlims table!\n"
		if $multi_loci;
	print STDERR "...Found multiple subtypes. Adding subtype to names in tree & xlims table!\n"
		if $multi_subtype;
		
	return $multi_loci, $multi_subtype;
	}

sub load_xlims{
	my ($dbh, $join_sql, $extra_query) = @_;
	
# same table join #
	my $query = "
SELECT 
loci.locus_start,
loci.locus_end,
loci.taxon_name,
loci.locus_id,
loci.subtype
FROM Loci Loci, Loci b
WHERE Loci.locus_id = b.locus_id
$join_sql
$extra_query
GROUP BY loci.locus_id
";
	$query =~ s/\n|\r/ /g;
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	# loading hash #
	my %xlims;
	foreach my $row (@$ret){
		$xlims{$$row[2]}{$$row[3]}{"entry"} = $row;
		$xlims{$$row[2]}{$$row[3]}{"subtype"} = $$row[4];
		}
	
		#print Dumper %xlims; exit;
	return \%xlims;
	}

sub join_query_opts_OLD{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}


