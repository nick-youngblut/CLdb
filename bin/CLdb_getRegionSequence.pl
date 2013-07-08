#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Bio::SeqIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $region_oi = "locus";
GetOptions(
	   "database=s" => \$database_file,
	   "region=s" => \$region_oi,				# region of interest
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
my $genbank_path = path_by_database($database_file);
$genbank_path = File::Spec->rel2abs($genbank_path);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# joining query options (for table join) #
my $join_sql = "";
$join_sql = join_query_opts(\@subtype, "subtype", $join_sql);
$join_sql = join_query_opts(\@taxon_id, "taxon_id", $join_sql);
$join_sql = join_query_opts(\@taxon_name, "taxon_name", $join_sql);

# region_oi #
my $region_oi_r = region_oi_check($region_oi);

# getting regions of interest from database #
my $regions_r = get_regions($dbh, $region_oi_r, $extra_query, $join_sql);

# parsing genbanks to get sequences #
parse_genbanks($regions_r, $genbank_path);

# writing fasta #
#write_arrays_fasta($arrays_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub parse_genbanks{
# parsing region sequences out of genbanks #
	my ($regions_r, $genbank_path) = @_;

	# output of genbank regions #
	my $out = Bio::SeqIO->new(-format => "fasta", -fh => \*STDOUT);

	# parsing each genbank #
	foreach my $genbank (keys %$regions_r){
		my $infile = "$genbank_path/$genbank";
		die " ERROR: $infile not found!\n" unless -e $infile;
		
		my $seqio = Bio::SeqIO->new(-format => "genbank", -file =>$infile)->next_seq();

		foreach my $locus (keys %{$regions_r->{$genbank}}){
			my $r_start = ${$regions_r->{$genbank}{$locus}}[0];
			my $r_end = ${$regions_r->{$genbank}{$locus}}[1];
			($r_start, $r_end) = flip_se($r_start, $r_end);
			
			my $region_seqio = $seqio->trunc($r_start, $r_end);	

			$region_seqio->display_id("cli.$locus");
			$region_seqio->desc(" $r_start\__$r_end");			
			
			$out->write_seq($region_seqio);
			}
		}
	}
	
sub flip_se{
	my ($start, $end) = @_;
	if($start <= $end){ return $start, $end; }
	else{ return $end, $start; }
	}

sub region_oi_check{
	my ($region_oi) = @_;
	if($region_oi =~ /locus/i){
		return [qw/locus_start locus_end/];
		}
	elsif($region_oi =~ /operon/i){
		return [qw/operon_start operon_end/];
		}
	elsif($region_oi =~ /crispr_array/i){
		return [qw/CRISPR_array_start CRISPR_array_end/];
		}
	else{ die " ERROR: region of interest not found in table!\n"; }
	}

sub write_arrays_fasta{
# writing arrays as fasta
	my ($arrays_r) = @_;
	
	foreach my $locus_id (keys %$arrays_r){
		foreach my $x_id (keys %{$arrays_r->{$locus_id}}){
			print join("\n", ">cli.$locus_id\__$x_id", $arrays_r->{$locus_id}{$x_id}), "\n";
			}
		}
	}

sub get_regions{
	my ($dbh, $region_oi_r, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query = join(" ", "SELECT locus_id, genbank, ", 
					join(",", @$region_oi_r),
					"FROM loci",
					$join_sql);
	$query = join(" ", $query, $extra_query);
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	
	my %regions;
	foreach my $row (@$ret){
		$regions{$$row[1]}{$$row[0]} = [@$row[2..$#$row]];		# genbank=>locus
		}
	
		#print Dumper %regions; exit;
	return \%regions;
	}
	
sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat, $join_sql) = @_;

	return $join_sql unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;

	if($join_sql){
		return join("", $join_sql, " AND $cat IN (", join(", ", @$vals_r), ")");
		}
	else{
		return join("", $join_sql, " WHERE $cat IN (", join(", ", @$vals_r), ")");
		}
	}

sub path_by_database{
	my ($database_file) = @_;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "genbank");
	}


__END__

=pod

=head1 NAME

CLdb_getRegionSequence.pl -- get nucleotide sequence of CRISPR loci regions

=head1 SYNOPSIS

CLdb_getRegionSequence.pl [flags] > regions.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -region

The CRISPR loci region of interest ('locus', 'operon', or 'CRISPR_array'). [locus] 

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

perldoc CLdb_getRegionSequence.pl

=head1 DESCRIPTION

Get nucleotide sequences for CRISPR locus regions of
interest. Output in fasta format. Regions of interest
can be the whole locus, the operon, or the CRISPR array.

Genbank files must be in $CLdb_HOME/genbank/

=head1 EXAMPLES

=head2 Write all locus regions:

CLdb_getRegionSequence.pl -d CRISPR.sqlite 

=head2 Write all CRISPR array regions:

CLdb_getRegionSequence.pl -data CRISPR.sqlite -r crispr_array

=head2 Refine region sequence query:

CLdb_getRegionSequence.pl -data CRISPR.sqlite -q "where LOCUS_ID=1" 

=head2 Refine region query to a specific subtype & 2 taxon_id's

CLdb_getRegionSequence.pl -da CRISPR.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

