#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use List::Util qw/max/;
use Bio::AlignIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $path);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
my $consensus_cutoff = 51;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "path=s" => \$path,
	   "cutoff=i" => \$consensus_cutoff,  
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# making an consensus temp directory #
my $dir = make_consensus_dir($database_file, $path); 


# getting loci in DR table #
my $DR_loci_r = get_DR_loci($dbh);

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting arrays of interest from database #
my %arrays;
if($join_sql){		# table-join query
	get_arrays_join(\%arrays, $dbh, $DR_loci_r, $extra_query, $join_sql);
	}
else{				# simple query of all spacers/DR
	get_arrays(\%arrays, $dbh, $DR_loci_r, $extra_query);
	}


# getting consensus of sequences #
my %DR_con;
foreach my $locus (keys %arrays){
	my $outfile = write_DR_fasta($arrays{$locus}, $dir, $locus);
	call_mafft($outfile, \%DR_con, $locus, $consensus_cutoff);
	}

# loading consensus sequences to db #
add_entries($dbh, \%DR_con);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub add_entries{
# adding entries to DR consensus table #
	my ($dbh, $DR_con_r) = @_;
	
	my $cmd = "INSERT INTO DR_consensus(Locus_ID, Consensus_Sequence_IUPAC, Consensus_Sequence_Threshold) values (?,?,?)";
	
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $locus_id (keys %$DR_con_r){
		$sql->execute( $locus_id, $DR_con_r->{$locus_id}{"IUPAC"}, $DR_con_r->{$locus_id}{"Threshold"} );
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for cli.$locus_id\n";
			}
		else{ $cnt++; }
		}
	$dbh->commit;

	print STDERR "...$cnt entries added/updated in database\n" unless $verbose;
	}

sub call_mafft{
	my ($outfile, $DR_con_r, $locus, $consensus_cutoff) = @_;
	
	my $alnIO = Bio::AlignIO->new(
				-file => "mafft --quiet $outfile |",
				-format => "fasta" );
	
	while(my $aln = $alnIO->next_aln){
		$DR_con_r->{$locus}{"IUPAC"} = $aln->consensus_iupac();
		$DR_con_r->{$locus}{"Threshold"} = $aln->consensus_string($consensus_cutoff);		
		$DR_con_r->{$locus}{"Threshold"} =~ s/\?/N/g;			# ambiguous letter = 'N'
		last;	
		}

		#print Dumper %$DR_con_r; exit;
	}


sub write_DR_fasta{
# writing arrays as fasta
	my ($DR_r, $dir, $locus) = @_;

	my $outfile = "$dir/cli$locus\_DR.fna";
	open OUT, ">$outfile" or die $!;
	foreach my $row (@$DR_r){
		print OUT join("\n", ">$$row[0]", $$row[1]), "\n";
		}
	close OUT;
	
	return $outfile;
	}

sub get_DR_loci{
# getting DR loci from directrepeats table #
	my ($dbh) = @_;

	my $query = "SELECT distinct(locus_id) FROM DRs";	
	my $ret = $dbh->selectall_arrayref($query);

	die " ERROR: no direct repeat entries found in database (DRs table)!\n"
		unless @$ret;
	
	return $ret;
	}

sub get_arrays{
	my ($arrays_r, $dbh, $DR_loci_r, $extra_query) = @_;
	
	# make query #
	my $query = "SELECT DR_ID, DR_sequence FROM DRs WHERE Locus_id = ?";
	$query = join(" ", $query, $extra_query);
	
	my $sql = $dbh->prepare($query);

	foreach my $locus (@$DR_loci_r){
		$locus = $$locus[0];
		$sql->execute($locus);
		my $ret = $sql->fetchall_arrayref();
		
		die " ERROR: no direct repeats entries for $locus!\n"
			unless @$ret;
		
		$arrays_r->{$locus} = $ret;	
		}
		
		#print Dumper %$arrays_r; exit;
	}
	
sub get_arrays_join{
	my ($arrays_r, $dbh, $DR_loci_r, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query = "SELECT a.DR_ID, a.DR_sequence FROM DRs a left outer join loci b ON a.locus_id = b.locus_id AND b.locus_id = ? $join_sql";
	$query = join(" ", $query, $extra_query);
	
	my $sql = $dbh->prepare($query);

	# query db #
	foreach my $locus (@$DR_loci_r){
		$locus = $$locus[0];
		$sql->execute($locus);
		my $ret = $sql->fetchall_arrayref();
		
		die " ERROR: no direct repeats entries for $locus!\n"
			unless @$ret;
		
		$arrays_r->{$locus} = $ret;	
		}
		
		#print Dumper %$arrays_r; exit;
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND b.$cat IN (", join(", ", @$vals_r), ")");
	}

sub make_consensus_dir{
# making a directory for grouping files	#
	my ($database_file, $path) = @_;
	my $dir;
	if($path){
		$dir = $path . "consensus/";
		}
	else{
		my @parts = File::Spec->splitpath($database_file);
		$dir = $parts[1] . "consensus/";
		}
	
	mkdir $dir unless -d $dir;
	
	return $dir;
	}


__END__

=pod

=head1 NAME

CLdb_loadDRConsensus.pl -- determining DR consensus sequences & adding to database

=head1 SYNOPSIS

CLdb_loadDRConsensus.pl [flags]

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

sql to refine query (see EXAMPLES).

=item -cutoff  <int>

Cutoff (threshold) using a residue in the consensus (%). [51]

=item -path  <char>

Directory where mafft alignments are performed. [$CLdb_HOME/consensus/]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_loadDRConsensus.pl

=head1 DESCRIPTION

Calculate the consensus sequence (IUPAC and
threshold) for all direct repeats 
in a CRISPR array.

For each array, DR sequences are aligned and
then the consensus sequence is calculated with bioperl.

=head3 IUPAC consensus

Base ambiguity letters are used.

=head3 Threshold consensus

The consensus residue must have a frequency of >= than '-cutoff',
otherwise a 'N' will be used.

=head1 EXAMPLES

=head2 Calculate & upload DB consensus sequences for all CRISPR arrays

CLdb_loadDRConsensus.pl -d CRISPR.sqlite 

=head2 Calculate & upload DB consensus sequences for just locus_id '1'

CLdb_loadDRConsensus.pl -d CRISPR.sqlite -q "where LOCUS_ID=1" 

=head2 Calculate & upload DB consensus sequences for specific subtype & 2 taxon_id's

CLdb_loadDRConsensus.pl -d CRISPR.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489
 
=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

