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

my ($verbose, $database_file, $leader_path);
my $extra_query = "";
my $genbank_path = "";
my $length = 400;
GetOptions(
	   "database=s" => \$database_file,
	   "query=s" => \$extra_query, 
	   "path=s" => \$leader_path,
	   "genbank=s" => \$genbank_path,
	   "length=i" => \$length,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);
$genbank_path = path_by_database($database_file) unless $genbank_path;


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# making leader directory #
$leader_path = make_leader_dir($leader_path, $database_file);

# getting array start-end from loci table #
my $array_se_r = get_array_se($dbh, $extra_query);

# determinig leader end of array based on DR degeneracy #
## pulling out direct repeat sequences for each loci ##
my $leader_loc_r = get_DR_seq($dbh, $array_se_r);

## getting sequences from genbanks ##
get_leader_seq($array_se_r, $leader_loc_r, $genbank_path, $length);

# writing fasta #
#write_arrays_fasta($arrays_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub path_by_database{
	my ($database_file) = @_;
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "genbank");
	}

sub get_leader_seq{
# getting leader seqs from genbanks #
	my ($array_se_r, $leader_loc_r, $genbank_path, $length) = @_;

	# getting all leader regions for each genbank #
	my %genbank_loci_index;
	foreach my $locus (keys %$array_se_r){
		push(@{$genbank_loci_index{ ${$array_se_r->{$locus}}[3] }}, $locus);
		}
	
	# opening genbanks & getting sequences #
	foreach my $genbank_file (keys %genbank_loci_index){
		my $seqio = Bio::SeqIO->new(-format => "genbank", -file => "$genbank_path/$genbank_file")->next_seq();
		foreach my $locus (@{$genbank_loci_index{ $genbank_file }} ){
			
			my $strand = shift @{$leader_loc_r->{$locus}};
			foreach my $loc (@{$leader_loc_r->{$locus}}){
#				print Dumper $loc; exit;
				
				# determining leader region start-end
				my $leader_start;
				my $leader_end;
			
				if( $loc eq "start" ){
					$leader_start = ${$array_se_r->{$locus}}[1] - $length;
					$leader_end = ${$array_se_r->{$locus}}[1]; 
					}
				}
			#elsif( $leader_loc eq "end" ){
			#	
			#	}
			#else{
				
			#	}

								#$leader_seqio = $seqio->trunc($region_start, $region_end);
			#my $leader_seqio;
			#if(${$array_se_r->{$locus}}[4] eq "+"){			# positive strand					
			#	$leader_seqio = $seqio->trunc($leader_start, $leader_end);	
			#	}
			}
		}
	
	}

sub get_DR_seq{
	my ($dbh, $array_se_r) = @_;
	
	my $cmd = "SELECT repeat_id, repeat_start, repeat_end, repeat_group from DirectRepeats where Locus_ID = ?";
	my $sql = $dbh->prepare($cmd);
	
	my %leader_loc;
	foreach my $locus (keys %$array_se_r){
		$sql->execute($locus);
		my $ret = $sql->fetchall_arrayref();
		$leader_loc{$locus} = determine_leader($ret, $array_se_r->{$locus});
		}
	
		#print Dumper %leader_loc; exit;
	return \%leader_loc;
	}

sub determine_leader{
# determing leader region #
	my ($ret, $array_r) = @_;
		#print Dumper $array_r; exit;	
	
	# getting halfway point in array #
	my $array_half = ($$array_r[2] +  $$array_r[1]) / 2;
		#print Dumper $array_half; exit;
	
	# counting groups for each half of the array #
	my %group_cnt;
	foreach my $DR (@$ret){		# each direct repeat
		# $DR = (repeat, start, end, group_ID)
		if($$DR[1] < $array_half && $$DR[2] < $array_half){
			$group_cnt{1}{$$DR[3]}++;
			}
		elsif($$DR[1] > $array_half && $$DR[2] > $array_half){
			$group_cnt{2}{$$DR[3]}++;
			}
		}
		
	my $first_cnt = scalar keys %{$group_cnt{1}};
	my $second_cnt = scalar keys %{$group_cnt{2}};

	# strand #
	my $strand;
	if($$array_r[1] <= $$array_r[2]){ $strand = "+"; }			# + strand
	else{ $strand = "-"; }
		
	# start or end? #		
	if($first_cnt > $second_cnt){
		return [$strand, "start"];		# start & stop must be inverted if negative strand
		}
	elsif($first_cnt < $second_cnt){
		return [$strand, "end"];		
		}
	else{
		return [$strand, "start", "end"];
		}
	}

sub get_array_se{
# getting the array start-end from loci table #
	my ($dbh, $extra_query) = @_;
	
	my $cmd = "SELECT Locus_ID, crispr_array_start, crispr_array_end, genbank FROM loci $extra_query";
	my $ret = $dbh->selectall_arrayref($cmd);

	my %array_se;
	foreach my $row (@$ret){
		$array_se{$$row[0]} = $row;
		}
		#print Dumper %array_se; exit;
	return \%array_se; 
	}

sub make_leader_dir{
# making a directory for leader region output #
	my ($leader_path, $database_file) = @_;
	
	unless ($leader_path){
		my @parts = File::Spec->splitpath($database_file);
		$leader_path = join("/", $parts[1], "leader");
		}
	mkdir $leader_path unless -d $leader_path;
	return $leader_path;
	}


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

