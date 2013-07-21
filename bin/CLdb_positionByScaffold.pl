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
use Set::IntervalTree;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
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

# loading genbanks #
my $gen_list_r = get_genbank_list();
my $gen_io_r = load_genbank_io($gen_list_r);

# making merged => unmerged location index #
foreach my $unmerged (keys %$gen_list_r){
	# getting scaffold lengths & sequence #
	my $unmerged_r = parse_unmerged_genbank($gen_list_r, $gen_io_r, $unmerged);
	my $merged_r = parse_merged_genbank($gen_list_r, $gen_io_r, $gen_list_r->{$unmerged});

	# finding merged-unmerged contig 1-to-1 assocations in contig length #
	my ($scaf_index_r, $needs_blasting_r) = find_same_length($gen_list_r, $unmerged_r, $merged_r);
	
	# determine if sequences are the same (either identical or rev-comp)
	
	# making an interval tree relating absolute location to scaffold-based location #
	my $itree = make_loc_itree($scaf_index_r, $merged_r, $unmerged_r);
	
	# getting values from database #
	my $loci_r = get_loci_oi($dbh, $gen_list_r->{$unmerged});
	
	}

# disconnect #
$dbh->disconnect();
exit;



### Subroutines
sub get_loci_oi{
# getting loci with genbank of interest #
	my ($dbh, $merged_genbank) = @_;

	#my $cmd = "SELECT * from loci where 
	}

sub make_loc_itree{
# making a location itree #
## input: genome-wide loc; output: unmerged scaf loc ##
	my ($scaf_index_r, $merged_r, $unmerged_r) = @_;
	
	my $itree = Set::IntervalTree->new();

	foreach my $merged_contig (keys %$scaf_index_r){
		foreach my $unmerged_contig (@{$scaf_index_r->{$merged_contig}}){
			$itree->insert($unmerged_r->{$unmerged_contig},
						$merged_r->{$merged_contig}{"start"},
						$merged_r->{$merged_contig}{"end"});
			}
		}

	#print_itree($itree, 1, 10000);
	#my $res_r = $itree->fetch(200, 300);
	#print Dumper $res_r; 
	
	return $itree;
	}


sub print_itree{
	my ($itree, $start, $end) = @_;
	
	for my $i ($start..$end){
		my $res = $itree->fetch($i, $i);
		
		print join("\t", "pos:$i", join("::", @$res)), "\n";
		}
	}

sub find_same_length{
# foreach merged-unmerged, find scaffolds that match in length #
## finding 1-to-1 associations in length ##
	my ($gen_list_r, $unmerged_r, $merged_r) = @_;
	
	my %scaf_index;
	foreach my $mcontig (keys %$merged_r){				# merged => unmerged
		foreach my $ucontig (keys %$unmerged_r){
			push( @{$scaf_index{$mcontig}}, $ucontig) if
				$merged_r->{$mcontig}{"length"} == $unmerged_r->{$ucontig}{"length"};
			}
		unless (exists $scaf_index{$mcontig} ){
			die " ERROR: Contig$mcontig did not have a match in the merged genbank!\n";
			}
		}		

	# making list of taxa that need blasting #
	my %needs_blasting;
	foreach my $contig (keys %scaf_index){
		if(scalar @{$scaf_index{$contig}} > 1){
			push (@{$needs_blasting{"merged"}}, $contig);
			foreach ( @{$scaf_index{$contig}} ){
				push (@{$needs_blasting{"unmerged"}}, $_);
				}
			delete $scaf_index{$contig};
			}
		elsif( scalar @{$scaf_index{$contig}} < 1){
			die " LOGIC ERROR: $!\n";
			}
		}
	
		#print Dumper %needs_blasting;
		#print Dumper %scaf_index; exit;	
	return \%scaf_index, \%needs_blasting;
	}

sub parse_unmerged_genbank{
# pulling out start-stop & sequence for all 'source' features #
	my ($gen_list_r, $gen_io_r, $unmerged) = @_;
	
	my %unmerged;			
	print STDERR "...loading: $unmerged\n" unless $verbose;
	while(my $seqo = $gen_io_r->{$unmerged}->next_seq){
		my @sources = grep{ $_->primary_tag eq "source" } $seqo->get_SeqFeatures;
		foreach my $source (@sources){
			my $contig_id = $source->seq->display_id; 
			$unmerged{$contig_id}{"start"} = $source->location->start;
			$unmerged{$contig_id}{"end"} = $source->location->end;
			$unmerged{$contig_id}{"length"} = abs($unmerged{$contig_id}{"end"} -
															$unmerged{$contig_id}{"start"});
			$unmerged{$contig_id}{"seq"} = $source->seq->seq;
			$unmerged{$contig_id}{"id"} = $contig_id;
			}
		}

	
		#print Dumper %unmerged; exit;
	return \%unmerged;
	}

sub parse_merged_genbank{
# pulling out start-stop & sequence for all 'source' features #
## start-end is absolute positioning ##
	my ($gen_list_r, $gen_io_r, $merged) = @_;
	
	my %merged;
	my $source_cnt = 0;
	print STDERR "...loading: $merged\n" unless $verbose;
	while(my $seqo = $gen_io_r->{$merged}->next_seq){
		my @sources = grep{ $_->primary_tag eq "source" } $seqo->get_SeqFeatures;
		foreach my $source (@sources){
			$merged{$source_cnt}{"start"} = $source->location->start;
			$merged{$source_cnt}{"end"} = $source->location->end;
			$merged{$source_cnt}{"length"} = abs($merged{$source_cnt}{"end"} -
														$merged{$source_cnt}{"start"});
			$merged{$source_cnt}{"seq"} = $source->seq->seq;
			$source_cnt++;
			}		
		}
		#print Dumper %merged; exit;
	return \%merged;
	}

sub load_genbank_io{
# parsing the merged genbank into scaffold lengths & sequence #
	my ($gen_list_r) = @_;
	
	my %gen_io;
	foreach my $file (%$gen_list_r){
		$gen_io{$file} = Bio::SeqIO->new(-format => "genbank", -file => $file);
		}
	
	return \%gen_io;
	}

sub get_genbank_list{
# getting list of genbanks #
## 2 column: unmerged, merged ##
	
	my %list;
	while(<>){
		chomp;
		next if /^\s*$/;
		
		my @tmp = split /\t/;
		die " ERROR: the genbank list must have 2 columns!\n"
			unless scalar @tmp == 2;
		
		$list{$tmp[0]} = $tmp[1];
		}
	
		#print Dumper %list; exit;
	return \%list;
	}


__END__

=pod

=head1 NAME

CLdb_makeDB.pl -- Initial DB construction

=head1 SYNOPSIS

CLdb_makeDB.pl [options] [DATABASE_name]

=head2 options

=over

=item -r 	Replace existing database.

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_makeDB.pl

=head1 DESCRIPTION

Make all of the CRISPR_db tables.

=head1 EXAMPLES

=head2 Naming database 'CRISPR_db1'

CLdb_makeDB.pl CRISPR_db1

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

