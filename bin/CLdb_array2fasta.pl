#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_array2fasta.pl -- write CRISPR array spacers or direct repeats to fasta

=head1 SYNOPSIS

CLdb_array2fasta.pl [flags] > array.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -repeat  <bool>

Get direct repeats instead of spacers. [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -group  <bool>

Get array elements de-replicated by group (ie. all unique sequences). [FALSE]

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_array2fasta.pl

=head1 DESCRIPTION

Get spacer or direct repeat sequences from the CRISPR database
and write them to a fasta.

By default, all spacers or direct repeats (if '-r') will be written.
The query can be refined with many of the flags. 

=head2 Output

If not using grouping: ">locus_ID|spacer/DR|elementID|groupID"

=head2 WARNING

Using '-leader' will only write out sequences from arrays with
identified leaders!

=head1 EXAMPLES

=head2 Write all spacers to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite 

=head2 Write all direct repeats to a fasta:

CLdb_array2fasta.pl -d CLdb.sqlite -r

=head2 Write all unique spacers

CLdb_array2fasta.pl -d CLdb.sqlite -g

=head2 Refine spacer sequence query:

CLdb_array2fasta.pl -d CLdb.sqlite -q "AND loci.Locus_ID=1" 

=head2 Refine spacer query to a specific subtype & 2 taxon_id's

CLdb_array2fasta.pl -d CLdb.sqlite -sub I-B -taxon_id 6666666.4038 6666666.40489

=head2 Orienting sequence by leaders

CLdb_array2fasta.pl -d CLdb.sqlite -l 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
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
	join_query_opts
	get_array_seq
	get_arrays_seq_byLeader/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
use CLdb::seq qw/
	revcomp/;



#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my ($spacer_DR_b, $by_group, $DR_x, $strand_b);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_DR_b,			 # spacers or repeats? [spacers]
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "group" => \$by_group,
	   "x=f" => \$DR_x,						# add DR extension to spacers?
	   "strand" => \$strand_b, 				# orient by loci table array_start, array_end strand? [TRUE]
	   "query=s" => \$extra_query, 
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "spacers") unless $spacer_DR_b;
table_exists($dbh, "DRs") unless $spacer_DR_b;


# determining table of interest
my $tbl_oi = "spacers";
$tbl_oi = "DRs" if $spacer_DR_b;
table_exists($dbh, $tbl_oi); 
die "ERROR: no entries in $tbl_oi table!\n" 
	unless n_entries($dbh, $tbl_oi);


# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");


# getting arrays of interest from database #
## defining options for query subroutines ##
my %opts = (
	extra_query => $extra_query,
	join_sql => $join_sql,
	by_group => $by_group,
	spacer_DR_b => $spacer_DR_b
	);

## querying CLdb ##
my $arrays_r = get_array_seq($dbh,\%opts);

## adding DR extension to spacers if needed ##
#add_DR_x($dbh, $arrays_r, $DR_x) if ! $spacer_DR_b && defined $DR_x;

## orienting by loci table strand ##
my $spacer_strand_r = get_spacer_strand($dbh, $arrays_r)
	if ! $spacer_DR_b & ! $strand_b;

# writing fasta #
write_array_seq($arrays_r, \%opts);


#--- disconnect ---#
$dbh->disconnect();
exit;


### Subroutines
sub add_DR_x{
	my ($dbh, $arrays_r, $DR_x) = @_;
	
	# queries #
	my $cmd1 = "SELECT DR_sequence from DRs where locus_ID = ? 
					AND DR_end <= ?
					AND DR_end >= ? - 1";
	my $sth1 = $dbh->prepare($cmd1);

	my $cmd2 = "SELECT DR_sequence from DRs where locus_ID = ? 
					AND DR_start <= ? + 1
					AND DR_start >= ?";
	my $sth2 = $dbh->prepare($cmd2);

	foreach my $entry (@$arrays_r){
			#print Dumper $entry; exit;
		$sth1->bind_param(1, $$entry[0]);
		$sth1->bind_param(2, $$entry[5]);
		$sth1->bind_param(3, $$entry[5]);
		$sth1->execute();
		my $ret1 = $sth1->fetchrow_arrayref();
		die "ERROR: no DR 3' matches for locus_id->$$entry[0]!\n"
			unless defined $ret1;
		my $DR_seq_5p = $$ret1[0];
		
		$sth2->bind_param(1, $$entry[0]);
		$sth2->bind_param(2, $$entry[6]);
		$sth2->bind_param(3, $$entry[6]);		
		$sth2->execute();
		my $ret2 = $sth2->fetchrow_arrayref();
		die "ERROR: no DR 5' matches for locus_id->$$entry[0]!\n"
			unless defined $ret2;
		my $DR_seq_3p = $$ret2[0];

		# lower case (for blastn masking) #
		map{ tr/A-Z/a-z/ } ($DR_seq_5p, $DR_seq_3p);
		
		# trimming DR #
		if($DR_x <= 1){						# assumed to take a fraction
			my $len = length $DR_seq_5p;
			$DR_seq_5p = substr($DR_seq_5p, 
						$len - int($len * $DR_x),			
						$len); 				# 3' end 
			$DR_seq_3p = substr($DR_seq_3p, 
						0,			
						int($len * $DR_x)); 				# 3' end 
			}
		
		# appending sequences to array 
		$$entry[4] = join("", $DR_seq_5p, $$entry[4], $DR_seq_3p);
		}
	
		#print Dumper @$arrays_r; exit;
	}

sub get_spacer_strand{
	my ($dbh, $arrays_r) = @_;
	
	# getting loci of array elements #	
	my %loci;
	map {$loci{$$_[0]} = 1} @$arrays_r;
	
	# getting strand of loci #
	my $cmd = "SELECT array_start, array_end from loci where locus_ID = ?";
	my $sth = $dbh->prepare($cmd);
	
	# querying & getting strand #
	foreach my $locus_id (keys %loci){
		$sth->bind_param(1, $locus_id);
		$sth->execute();
		my $ret = $sth->fetchrow_arrayref();

		die "ERROR: not matches for locus_id->$locus_id!\n"
			unless defined $ret;
		
		# flipping orientation if needed #
		if($$ret[0] > $$ret[1]){
			$loci{$locus_id} = '-';
			print STDERR "...loci '$locus_id' array is designated - strand. Reversing the orientation of all spacers in the array\n";
			}
		else{ $loci{$locus_id} = '+'; }
		}
	
	return \%loci;
	}

sub write_array_seq{
# writing arrays as fasta
	my ($arrays_r, $opts_r) = @_;
	
	foreach (@$arrays_r){
		print join("\n", 
			join("|", ">$$_[0]", @$_[1..3]), 
			$$_[4]), "\n";
		}
	}






