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
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "path=s" => \$path, 	   
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
	call_mafft($outfile, \%DR_con, $locus);
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
	
	my $cmd = "INSERT INTO DirectRepeatConsensus(Locus_ID, Repeat_Consensus_Sequence) values (?,?)";
	
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $locus_id (keys %$DR_con_r){
		$sql->execute( $locus_id, $DR_con_r->{$locus_id} );
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", $locus_id, $DR_con_r->{$locus_id} ), "\n";
			}
		else{ $cnt++; }
		}
	$dbh->commit;

	print STDERR "...$cnt entries added/updated in database\n" unless $verbose;
	}

sub get_consensus{
	my ($DR_con_r, $locus, $DR_fasta_r) = @_;

	# counting number of each nucleotide at each position #
	my %by_pos;	
	foreach my $DR (keys %$DR_fasta_r){
		$DR_fasta_r->{$DR} =~ tr/A-Z/a-z/;			# lower case
		my @line = split //, $DR_fasta_r->{$DR};
		for my $i (0..$#line){
			$by_pos{$i}{$line[$i]}++;				# position, nucleotide
			}
		}
		
	# determining majority at each position #
	my $con = "";
	foreach my $pos (sort {$a<=>$b} keys %by_pos){
		my $max_val = max (values %{$by_pos{$pos}});

		#print Dumper $max_val; exit;
		}
	
		
	#print Dumper %by_pos; exit;
	}

sub call_mafft{
	my ($outfile, $DR_con_r, $locus) = @_;
	
	my $alnIO = Bio::AlignIO->new(
				-file => "mafft --quiet $outfile |",
				-format => "fasta" );
	
	while(my $aln = $alnIO->next_aln){
		$DR_con_r->{$locus} = $aln->consensus_iupac();
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

	my $query = "SELECT distinct(locus_id) FROM DirectRepeats";	
	my $ret = $dbh->selectall_arrayref($query);

	die " ERROR: no direct repeat entries found in database!\n"
		unless @$ret;
	
	return $ret;
	}

sub get_arrays{
	my ($arrays_r, $dbh, $DR_loci_r, $extra_query) = @_;
	
	# make query #
	my $query = "SELECT Repeat_ID, Repeat_sequence FROM directrepeats WHERE Locus_id = ?";
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
	my $query = "SELECT a.Repeat_ID, a.Repeat_sequence FROM directrepeats a left outer join loci b ON a.locus_id = b.locus_id AND b.locus_id = ? $join_sql";
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

CLdb_loadDRConsensuss.pl -- adding/updating DR consensus sequences

=head1 SYNOPSIS

CLdb_loadDRConsensuss.pl [flags] leader.fasta leader_aligned.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -t 	Number of bp (& gaps) to trim of end of alignment. [0]

=item -g 	Leave gaps in leader sequence entered into the DB? 

=item -v	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_loadDRConsensuss.pl

=head1 DESCRIPTION

After determining the actual leader region by aligning
approximate leader regions written by CLdb_getLeaderRegions.pl,
automatically trim and load leader sequences & start-end
information into the CRISPR database.

Trimming will be done from the leader region end most
distant from the CRISPR array. Reversing & 
reverse-complementation of sequences during the alignment
will be accounted for (that's why both the aligned and
'raw' sequenced must be provided).

Not all sequences in the aligned fasta need to be in
the 'raw' leader fasta (eg. if both ends of the array
were written because the side containing the leader
region could not be determined from direct-repeat
degeneracy.


=head1 EXAMPLES

=head2 Triming off the 50bp of unconserved alignment 

CLdb_loadDRConsensuss.pl -d ../CRISPR.sqlite test_leader_Ib.fna test_leader_Ib_aln.fna -t 50

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

