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

my ($verbose, $database_file, $genbank_path, $quiet);
my $query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "genbank=s" => \$genbank_path,
	   "sql=s" => \$query,
	   "verbose" => \$verbose,				# TRUE
	   "quiet" => \$quiet, 					# turn off warnings
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
my $db_path = get_database_path($database_file);
$genbank_path = "$db_path/genbank" unless $genbank_path;
$genbank_path = File::Spec->rel2abs($genbank_path);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# getting loci start - end from loci table (& genbank file names)
my $loci_se_r = get_loci_start_end($dbh, $query);

# getting CDS info for CDS features in loci regions #
my $loci_tbl_r = call_genbank_get_region($loci_se_r, $genbank_path);

# determine if CDS are in operons #
check_in_operon($loci_se_r, $loci_tbl_r);

# write locus-gene table #
write_loci_tbl($loci_tbl_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_loci_tbl{
# writing out loci table for manual editting of aliases #
	my ($loci_tbl_r) = @_;

	print join("\t", qw/Locus_ID Gene_Id Gene_start Gene_end Gene_length__AA In_Operon Gene_Alias Sequence/), "\n";
	foreach my $loci (keys %$loci_tbl_r){
		foreach my $feature (sort{$a<=>$b} keys %{$loci_tbl_r->{$loci}}){

			unless (${$loci_tbl_r->{$loci}{$feature}}[2] =~ /fig\|.+peg/){		# db_xref must be fig|.+peg
				print STDERR " WARNING: Locus$loci -> $feature does not have a FIG-PEG ID in a db_xref tag!\n"
					unless $quiet;
				next;
				}
			print join("\t", "lci.$loci", 
				${$loci_tbl_r->{$loci}{$feature}}[2], 			# Gene_id (db_xref)
				@{$loci_tbl_r->{$loci}{$feature}}[(0..1)], 		# start, end 
				length ${$loci_tbl_r->{$loci}{$feature}}[4], 	# length (AA)
				${$loci_tbl_r->{$loci}{$feature}}[5], 			# In Operon	
				${$loci_tbl_r->{$loci}{$feature}}[3],
				${$loci_tbl_r->{$loci}{$feature}}[4]
				), "\n";
			}
		}
	}

sub check_in_operon{
# checing whether CDS features are in the designated operon locations #
	my ($loci_se_r, $loci_tbl_r) = @_;
	
	foreach my $locus (keys %$loci_tbl_r){
		foreach my $feature (keys %{$loci_tbl_r->{$locus}}){
			die " LOGIC ERROR: $!\n" unless 
				exists $loci_se_r->{$locus};
			
			# defining start - end #
			## f = feature; o = operion; c = crispr array ##
			my $f_start = ${$loci_tbl_r->{$locus}{$feature}}[0];
			my $f_end = ${$loci_tbl_r->{$locus}{$feature}}[1];
			
			my $o_start = ${$loci_se_r->{$locus}}[3];
			my $o_end = ${$loci_se_r->{$locus}}[4];			
			
			my $c_start = ${$loci_se_r->{$locus}}[5];
			my $c_end = ${$loci_se_r->{$locus}}[6];
			
	#			print Dumper $f_start, $f_end, $o_start, $o_end, $c_start, $c_end; 

			# making all on the + strand #
			($f_start, $f_end) = set_to_pos_strand($f_start, $f_end);
			($o_start, $o_end) = set_to_pos_strand($o_start, $o_end);
			($c_start, $c_end) = set_to_pos_strand($c_start, $c_end);			
			
			# determining location #
			## gene must fall in operon, without falling into crispr array ##
			if($f_start >= $o_start && $f_end <= $o_end){	# gene in operon
				if( ($f_start < $c_start && $f_end < $c_end) ||
					($f_start > $c_start && $f_end > $c_end) ){		# not in crispr array
					push(@{$loci_tbl_r->{$locus}{$feature}}, 1);
					}
				else{ push(@{$loci_tbl_r->{$locus}{$feature}}, 1); }
				}
			else{ push(@{$loci_tbl_r->{$locus}{$feature}}, 1); }
			}
		}
		#print Dumper %$loci_tbl_r; exit;
	}

sub set_to_pos_strand{
# setting all start-end so start is <= end #
	my ($start, $end) = @_;
	if ($start > $end){ return $end, $start; }
	else{ return $start, $end; }
	}

sub get_database_path{
	my $database_file = shift;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return $parts[1];
	}

sub call_genbank_get_region{
# calling genbank_get_region to get CDS info from genbank files #
	my ($loci_se_r, $genbank_path) = @_;
	
	my %loci_tbl;
	foreach my $locus (keys %$loci_se_r){
		my $genbank_file = join("/", $genbank_path, ${$loci_se_r->{$locus}}[2]);
		die " ERROR: $genbank_file not found!\n"
			unless -e $genbank_file;
		
		my $cmd; 
		my $start = ${$loci_se_r->{$locus}}[0];
		my $end = ${$loci_se_r->{$locus}}[1];
		if($start <= $end){
			$cmd = "genbank_get_region.pl -r $start $end < $genbank_file |";
			}
		else{
			$cmd = "genbank_get_region.pl -r $start $end < $genbank_file |";			
			}
			
		print STDERR "$cmd\n" unless $verbose;
		open PIPE, $cmd or die $!;
		while(<PIPE>){
			chomp;
			next if $. == 1;
			
			my @line = split /\t/;
			$loci_tbl{$locus}{$line[0]} = [@line[1..$#line]]; 	# locusID=>feature_num = feature
			}
		close PIPE;
		}
		#print Dumper %loci_tbl; exit;
	return \%loci_tbl;
	}

sub get_loci_start_end{
# getting locus start and locus end from loci table in db #
	my ($dbh, $query) = @_;
	
	my %loci_se;
	my $cmd = "SELECT Locus_id, Locus_start, Locus_end, Genbank, Operon_start, Operon_end, CRISPR_Array_Start, CRISPR_Array_End from loci $query";
	my $loci_se_r = $dbh->selectall_arrayref($cmd);	
	die " ERROR: no locus start-end values found!\n" unless $loci_se_r;	

	foreach my $row (@$loci_se_r){
		$loci_se{$$row[0]} = [@$row[1..$#$row]];
		}
		#print Dumper %loci_se; exit;
	return \%loci_se;
	}



__END__

=pod

=head1 NAME

CLdb_getGenesInLoci.pl -- getting all genes in CRISPR Loci; output is a table of gene information

=head1 SYNOPSIS

CLdb_getGenesInLoci.pl [flags] 

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -s 	sql to refine loci table query (see EXAMPLES). 

=item -q	Turn off all warnings.

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_getGenesInLoci.pl

=head1 DESCRIPTION

Get all CDS features from >= genbanks that fall
into CRISPR loci. 

The output can be piped directly into
CLdb_loadGenes.pl or the aliases (or other values)
can be edited first.

Determining whether genes fall into the designated
operon for the locus (and not in the CRISPR array)
is done by determining if genes fall into the operon 
range while falling outside of the CRISPR array range.

=head2 WARNING:

CDS features in genbanks must have db_xref tags with 
the fig|peg ID (example: "fig|6666666.40253.peg.2362");
otherwise, the CDS will not be written to the output table.

=head2 Requires:

genbank_get_region.pl

=head1 EXAMPLES

=head2 Basic usage:

CLdb_getGenesInLoci.pl -d CRISPR.sqlite > gene_info.txt

=head2 Direct loading to CRISPR database

CLdb_getGenesInLoci.pl -d CRISPR.sqlite | CLdb_loadGenes.pl -d CRISPR.sqlite

=head2 Refining query to just 1 subtype

CLdb_getGenesInLoci.pl -d CRISPR.sqlite -s "WHERE subtype='I-B'"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

