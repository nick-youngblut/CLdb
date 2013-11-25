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
my $range = 150;
my $gap_cutoff = 50;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "range=i" => \$range,
	   "gap=i" => \$gap_cutoff,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;

### MAIN
# path #
my $db_path = get_database_path($database_file);
my $fasta_dir = make_fasta_dir($db_path);

# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting array start-end, scaffold, genbank, & fasta info #
my $array_r = get_array_info($dbh, $join_sql);

# determining whether array is at the end of a contig (end of scaffold or next to gap) #
## getting scaffold start-end & gap position info ##
my $contig_pos_r = get_contig_pos($array_r, $db_path, $fasta_dir);

## compare array start-end to contig_pos info ##
compare_start_end($array_r, $contig_pos_r, $range);


# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub compare_start_end{
# finding broken arrays and writing #
# broken/intact 
# scaffold_end/contig_end 
# bp from end 
	my ($array_r, $contig_pos_r, $range) = @_;

	# header #
	print join("\t", qw/status break_ID bp_from_break locus_id taxon_name taxon_id subtype locus_start locus_end array_start array_end genbank_file fasta_file scaffold/), "\n";

	# body #	
	foreach my $entry (@$array_r){
		
		# if no array start-end; NA #
		unless($$entry[6] && $$entry[7]){
			$$entry[6] = "" unless $$entry[6];
			$$entry[7] = "" unless $$entry[7];
			print join("\t", "missing", "NA", "NA", @$entry), "\n";
			next;
			}

		# setting to + strand #
		my ($array_start, $array_end) = flip_se(@$entry[6..7]);
		
		# at beginning of contig #
		if($array_start < $range ){		
			print join("\t", "broken", "scaffold_end", $array_start, @$entry), "\n";
			}
		# no fasta found; skipping #
		elsif(! exists $contig_pos_r->{$$entry[9]}){
			print STDERR " WARNING: no fasta for taxon_name->$$entry[1] taxon_id->$$entry[2]. Cannot determine whether array is at the end of a contig!\n";
			print join("\t", "NA", "NA", "NA", @$entry), "\n";			
			}
		# no scaffold! shouldn't happen #
		elsif(! exists $contig_pos_r->{$$entry[9]}{$$entry[10]} ){
			die " ERROR: $$entry[10] not found in $$entry[9]! Is the scaffold name in the loci table wrong?\n";
			}
		# if scaffold exists; check scaffold end & gap start-end #
		else{		
			# checking scaffold start-end #
			my $scaf_len = $contig_pos_r->{$$entry[9]}{$$entry[10]}{"length"};
			if($scaf_len - $array_end < $range){
				print join("\t", "broken", "scaffold_end", $scaf_len - $array_end, @$entry), "\n";
				next;
				}
			
			# checking contig start-end #
			foreach my $gap (@{$contig_pos_r->{$$entry[9]}{$$entry[10]}{"gap"}}){		# gap start-end
				if($$gap[0] - $array_end > 0 && $$gap[0] - $array_end < $range){
					print join("\t", "broken", "contig_end", $$gap[0] - $array_end, @$entry), "\n";	
					next;
					}
				elsif($array_start - $$gap[1] > 0 && $array_start - $$gap[1] < $range){
					print join("\t", "broken", "contig_end", $array_start - $$gap[1], @$entry), "\n";	
					next;
					}
				}
			
			# if passed other checks; not broken #
			print join("\t", "intact", "NA", "NA", @$entry), "\n";	
			}
		}
	}

sub get_contig_pos{
	my ($array_r, $db_path, $fasta_dir) = @_;
	
	my %contig_pos;
	my %loci2update;
	foreach my $entry (@$array_r){
			#print Dumper $entry; exit;

		# skipping genomes already loaded #
		next if $$entry[9] && exists $contig_pos{$$entry[9]};
		
		# getting fasta if needed #
		if(! $$entry[9] || ! -e "$db_path/fasta/$$entry[9]"){		# if no fasta file
			print STDERR " WARNING: no fasta for taxon_name->$$entry[1] taxon_id->$$entry[2]! Trying to extract sequence from genbank...\n";
			unless($$entry[8]){
				print STDERR " WARNING: no 'genbank_file' value for taxon_name->$$entry[1] taxon_id->$$entry[2]! Skipping!\n";
				next;
				}
			$$entry[9] = genbank2fasta_extract($$entry[8], $fasta_dir);
			$loci2update{$$entry[8]} = $$entry[9];			# genbank => fasta
			}

		# loading fasta #
		print STDERR "...loading $db_path/fasta/$$entry[9]\n" unless $verbose;
		open IN, "$db_path/fasta/$$entry[9]" or die $!;
		my $tmp = "";
		while(<IN>){
			# getting scaffold length #
			chomp;
			next if /^\s*$/;
			if(/^>/){
				($tmp = $_) =~ s/^>//;
				}
			else{
				$contig_pos{$$entry[9]}{$tmp}{"length"} += length $_;
             	while ( $_ =~ /(N+|-+)/gi){             
                    push @{$contig_pos{$$entry[9]}{$tmp}{"gap"}}, 
                    	[ $contig_pos{$$entry[9]}{$tmp}{"length"} + $-[0], 
                    	  $contig_pos{$$entry[9]}{$tmp}{"length"} + $+[0] ]
                    	unless $+[0] - $-[0] < $gap_cutoff;
                    }
				}
		
			}
		close IN;
		
			#print Dumper %contig_pos; exit;
		}
	
	# updating loci table #
	update_loci($dbh, \%loci2update) if %loci2update;
	
		#print Dumper %contig_pos; exit;
	return \%contig_pos;
	}

sub update_loci{
# updating loci w/ fasta file info for newly written fasta #
	my ($dbh, $loci_r) = @_;
		
	# status#
	print STDERR "...updating Loci table with newly made fasta files\n";
	
	# sql #
	my $q = "UPDATE loci SET fasta_file=? WHERE genbank_file=?";
	my $sth = $dbh->prepare($q);
	
	foreach my $genbank (keys %$loci_r){

		$sth->bind_param(1, $loci_r->{$genbank} );				# fasta_file (just file name)
		$sth->bind_param(2, $genbank );				# array_file
		$sth->execute( );
		
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", $genbank, $loci_r->{$genbank}), "\n";
			}
		}
	$dbh->commit();
	
	print STDERR "...updates committed!\n";
	}

sub genbank2fasta_extract{
	my ($genbank_file, $fasta_dir) = @_;
	# sanity check #
	die " ERROR: cannot find '$db_path/genbank/$genbank_file!'\n"
		unless -e "$db_path/genbank/$genbank_file";
		
	my $seqio = Bio::SeqIO->new(-file => "$db_path/genbank/$genbank_file", -format => "genbank");
	
	# output #
	my @parts = File::Spec->splitpath($genbank_file);
	$parts[2] =~ s/\.[^.]+$|$/.fasta/;
	my $fasta_out = "$fasta_dir/$parts[2]";
	open OUT, ">$fasta_out" or die $!;
	
	# writing fasta #
	my $seq_cnt = 0;
	while(my $seqo = $seqio->next_seq){
		$seq_cnt++;
		
		# seqID #
		my $scafID = $seqo->display_id;
		print OUT ">$scafID\n";
			for my $feato (grep { $_->primary_tag eq 'source' } $seqo->get_SeqFeatures){
			print OUT $feato->seq->seq, "\n";
			}
		}
	close OUT;
	
	# if genome seq found and fasta written, return fasta #
	if($seq_cnt == 0){
		print STDERR " WARNING: no genome sequnece found in Genbank file: $genbank_file!\nSkipping BLAST!\n";
		unlink $fasta_out;
		return 0;
		}
	else{ 
		print STDERR "...fasta file extracted from $genbank_file written: $fasta_out\n";
		return $parts[2]; 			# just fasta name
		}
	}

sub check_loci_range{
# checking for arrays that have loci range beyond set range #
	my ($array_r, $range) = @_;
	
	my @new_array;
	foreach my $entry (@$array_r){
		# setting to + strand #
		my ($locus_start, $locus_end) = flip_se(@$entry[4..5]);
		my ($array_start, $array_end) = flip_se(@$entry[6..7]);
		
		if($array_start - $locus_start < $range &&		# if locus ends at end of array, keep
			$locus_end - $array_end < $range){
			push @new_array, $entry;
			}
		}
	return \@new_array;
	}
	
sub flip_se{
	my ($start, $end) = @_;
	if($start <= $end){ return $start, $end; }
	else{ return $end, $start; }
	}

sub get_array_info{
# getting info on arrays of interest #
## arrays must have a start-end ##
	my ($dbh, $join_sql) = @_;
	
	my $q = "
SELECT 
locus_id,
taxon_name,
taxon_id,
subtype,
locus_start,
locus_end,
array_start,
array_end,
genbank_file,
fasta_file,
scaffold
from loci
where locus_id IS NOT NULL
$join_sql
";
	$q =~ s/\n/ /g;

	my $ret = $dbh->selectall_arrayref($q);
	
	die " ERROR: no matching entries!\n"
		unless @$ret;
		
	foreach (@$ret){
		map{$_ = "" unless $_} @$_[1..3,10];
		}
		
		#print Dumper @$ret; exit;
	return $ret;
	}
	
sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND $cat IN (", join(", ", @$vals_r), ")");
	}

sub get_database_path{
	my $database_file = shift;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return $parts[1];
	}

sub make_fasta_dir{
	my $db_path = shift;
	
	my $dir = join("/", $db_path, "fasta");
	mkdir $dir unless -d $dir;

	return $dir;
	}

sub check_for_loci_table{
	my ($table_list_r) = @_;
	die " ERROR: loci table not found in database!\n"
		unless grep(/^loci$/i, @$table_list_r);
	}

sub list_tables{
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}


__END__

=pod

=head1 NAME

CLdb_findBrokenArrays.pl -- arrays broken by assembly errors

=head1 SYNOPSIS

CLdb_findBrokenArrays.pl [flags] > summary.txt

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

Extra sql to refine which sequences are returned.

=item -range  <int>

How close a contig end must be to a CRISPR array end to call it broken. [150]

=item -gap  <int>

The size cutoff for using a gap. [50]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_findBrokenArrays.pl

=head1 DESCRIPTION

Determine which CRISPR loci appear to
be broken due to breaks in the genome assembly.

The output is tab-delimited table with some
breakage summary columns appended to the beginning
of some loci table columns.

=head2 summary columns:

=over 

=item 'status'

The array is either intact, broken (end of contig), or missing

=item break_ID

Is the array next to the end of the entire scaffold or next to a gap inside the scaffold?

=item bp_from_break

Distance from the scaffold end or gap

=back

=head1 EXAMPLES

=head2 Find broken arrays among all loci

CLdb_findBrokenArrays.pl -d CLdb.sqlite 

=head2 Find broken arrays in Subtype I-A

CLdb_findBrokenArrays.pl -d CLdb.sqlite -sub I-A

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

