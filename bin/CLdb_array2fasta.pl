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

=item -leader  <bool>

Orient array sequences to leader (leader on 5')? [FALSE]

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

If not using grouping: ">cli.[locus_ID]__[start]-[end]__[array_order]"

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

=head2 Ordering by leaders 

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


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $spacer_DR_b, $by_group, $leader_b);
my (@subtype, @taxon_id, @taxon_name);
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "repeat" => \$spacer_DR_b,			 # spacers or repeats?
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "group" => \$by_group,
	   "leader" => \$leader_b,				# orient by leader?
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
my $arrays_r;
if($leader_b){ $arrays_r = get_arrays_seq_byLeader($dbh,\%opts); }
else{ $arrays_r = get_array_seq($dbh,\%opts); }

# writing fasta #
write_array_seq($arrays_r, \%opts);


#--- disconnect ---#
$dbh->disconnect();
exit;


### Subroutines
sub write_array_seq{
# writing arrays as fasta
	my ($arrays_r, $opts_r) = @_;
	
	foreach my $seq_id (sort keys %$arrays_r){
		# grouping #
		if(defined $opts_r->{"by_group"}){
			foreach my $start ( keys %{$arrays_r->{$seq_id}} ){	
				if(defined $opts_r->{"spacer_DR_b"}){
					print join("\n", ">DR-Group.$seq_id", $arrays_r->{$seq_id}{$start}{"seq"}), "\n";
					}
				else{
					print join("\n", ">Spacer-Group.$seq_id", $arrays_r->{$seq_id}{$start}{"seq"}), "\n";
					}
				}
			}
		# no grouping #
		else{
			my $order = 0;
			foreach my $start (sort{$a<=>$b} keys %{$arrays_r->{$seq_id}}){	
				$order++;
				print join("\n", 
						join("__", ">cli.$seq_id", 
							join("-", $start, $arrays_r->{$seq_id}{$start}{"stop"}),
							$order),
						$arrays_r->{$seq_id}{$start}{"seq"}), "\n";
				}
			}
		}
	}






