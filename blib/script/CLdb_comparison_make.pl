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


my ($verbose, $CLdb_sqlite, @ITEP_sqlite);
my $spacer_cutoff = 1;
my $spacer_opt = "cluster";
GetOptions(
	   "database=s" => \$CLdb_sqlite,
	   "ITEP=s{,}" => \@ITEP_sqlite,
	   "cutoff=f" =>  \$spacer_cutoff,
	   "spacer=s" => \$spacer_opt, 		# 'blast' or 'cluster'
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# checking CLdb #
die " ERROR: provide a database file name!\n"
	unless $CLdb_sqlite;
die " ERROR: cannot find CLdb database file!\n"
	unless -e $CLdb_sqlite;
	
# checking ITEP input #
if(@ITEP_sqlite){
	die " ERROR: '-ITEP' must be 2 arguments (database_file, runID)!\n"
		unless scalar @ITEP_sqlite >= 2;
	die " ERROR: cannot find ITEP database file!\n"
		unless -e $ITEP_sqlite[0];
	}

### MAIN
# connect 2 CLdb #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$CLdb_sqlite", '','', \%attr) 
	or die " Can't connect to $CLdb_sqlite!\n";

# connect to ITEP #
my $dbh_ITEP;
if(@ITEP_sqlite){
	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	$dbh_ITEP = DBI->connect("dbi:SQLite:dbname=$ITEP_sqlite[0]", '','', \%attr) 
		or die " Can't connect to $ITEP_sqlite[0]!\n";	
	}


# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r) = load_dna_segs();
my %compare;

# querying CLdb #
## for querying spacer_hclust, need: ##
if( exists $dna_segs_r->{"spacer"} ){
	if($spacer_opt =~ /cluster/i){
		get_spacer_hclust($dbh, $dna_segs_r, $dna_segs_order_r, \%compare, $spacer_cutoff);
		}
	elsif($spacer_opt =~ /blast/i){
	
		}
	}
	
## for querying spacer_pairwise_blast, need: ##
	# spacerIDs (by dna_segs_id)


# querying ITEP #
## for querying ITEP, need:
	# runID
	# geneIDs (by dna_segs_id)


# writing out compare table #
write_compare(\%compare, $dna_segs_order_r);

# disconnect #
$dbh->disconnect();
$dbh_ITEP->disconnect() if $dbh_ITEP;
exit;


### Subroutines
sub write_compare{
# writing out compare table #
	my ($compare_r, $dna_segs_order_r) = @_;
	
	# header #
	print join("\t", qw/start1 end1 start2 end2 col dna_segs_id1 dna_segs_id2 feat feat_id/), "\n";
	
	# body #
	foreach my $dna_seg1 (@$dna_segs_order_r){
		foreach my $feat (keys %$compare_r){
			foreach my $dna_seg2 (keys %{$compare_r->{$feat}{$dna_seg1}}){
				foreach my $feat_id (keys %{$compare_r->{$feat}{$dna_seg1}{$dna_seg2}}){
					print join("\t", @{$compare_r->{$feat}{$dna_seg1}{$dna_seg2}{$feat_id}}, 
							$dna_seg1, $dna_seg2, $feat, $feat_id), "\n";
					}
				}
			}
		}
	}

sub get_spacer_hclust{
# querying CLdb for spacer_hclust info #
	my ($dbh, $dna_segs_r, $dna_segs_order_r, $compare_r, $spacer_cutoff) = @_;
	
	# preparing query #
	my $query = "
SELECT b.spacer_start, b.spacer_end 
FROM spacer_hclust a, spacers b 
WHERE a.locus_id=b.locus_id 
AND  a.spacer_id=? 
AND b.spacer_id=?
AND a.cutoff = ?
AND a.locus_id IN (?, ?);
";
	$query =~ s/\r|\n/ /g;

	my $sth = $dbh->prepare($query);
	
	# querying each spacer ID for matches against adjacent locus #
	for my $i (0..($#$dna_segs_order_r-1)){
		my $dna_seg_id1 = $$dna_segs_order_r[$i];
		my $dna_seg_id2 = $$dna_segs_order_r[$i+1];
		# getting loci to compare #
		my ($locus_id1, $locus_id2);
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}}){
			$locus_id1 = $locus_id;
			}
		foreach my $locus_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id2}}){
			$locus_id2 = $locus_id;
			}
		foreach my $feat_id (keys %{$dna_segs_r->{"spacer"}{$dna_seg_id1}{$locus_id1}}){
			#print Dumper $feat_id, $spacer_cutoff, $locus_id1, $locus_id2; exit;
			$sth->bind_param(1, $feat_id);
			$sth->bind_param(2, $feat_id);			
			$sth->bind_param(3, $spacer_cutoff);
			$sth->bind_param(4, $locus_id1);
			$sth->bind_param(5, $locus_id2);
			$sth->execute();
			my $res = $sth->fetchall_arrayref();
			
			# skipping any spacers not in both loci #
			next if @$res && scalar @$res != 2;
				#print Dumper $res if @$res && scalar @$res != 2;
			
			# loading hash #
			$compare_r->{"spacer"}{$dna_seg_id1}{$dna_seg_id2}{$feat_id} = [$$res[0][0], $$res[0][1], $$res[1][0], $$res[1][1], 100];
			}
		}
		#print Dumper %$compare_r; exit;
	}

sub load_dna_segs{
# loading dna_segs from stdin #

	my @dna_segs_order;
	my %dna_segs;
	my %header;
	my %name_index;
	while(<>){
		chomp;
		next if /^\s*$/;

		# header #
		if($.==1){
			# indexing header #
			my @tmp = split /\t/;
			for my $i (0..$#tmp){
				$header{$tmp[$i]} = $i;
				}
			
			# checking for existence of columns #
			my @check = qw/locus_id dna_segs_id feat feat_id/;
			map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
				unless exists $header{$_} } @check;
			}
		else{
			my @line = split /\t/;
			my $locus_id = $line[$header{"locus_id"}];
			my $dna_seg_id = $line[$header{"dna_segs_id"}];
			my $feat = $line[$header{"feat"}];
			my $feat_id = $line[$header{"feat_id"}];						
			
			# order of loci#
			push( @dna_segs_order, $dna_seg_id)
				unless exists $name_index{$dna_seg_id};
			$name_index{$dna_seg_id} = 1;
			
			# dna_segs #
			die " ERROR: $feat -> $dna_seg_id -> $locus_id -> $feat_id is not unique!\n"
				if exists $dna_segs{$feat}{$dna_seg_id}{$locus_id}{$feat_id};
			$dna_segs{$feat}{$dna_seg_id}{$locus_id}{$feat_id} = \@line;
			
			}
		}
		
		#print Dumper @dna_segs_order; exit;
		#print Dumper %dna_segs; exit;
	return \%dna_segs, \@dna_segs_order, \%header;
	}


__END__

=pod

=head1 NAME

CLdb_dna_segs_make.pl -- making dna_segs table for plotting

=head1 SYNOPSIS

CLdb_dna_segs_make.pl [flags] > dna_segs.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -ITEP

Get gene cluster info from ITEP. 2 arguments required: (ITEP_sqlite_file, cluster_runID).

=item -cutoff

Spacer clustering cutoff for spacer coloring (0.8 - 1). [1]

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

perldoc CLdb_dna_segs_make.pl

=head1 DESCRIPTION

Make a basic dna_segs object needed for plotting.

=head2 Coloring genes by gene cluster

Provide and ITEP sqlite file name and cluster run ID
to add cluster info to the table (used for coloring)

=head1 EXAMPLES

=head2 Plotting all loci classified as subtype 'I-A'

CLdb_dna_segs_make.pl -d CLdb.sqlite -sub I-A 

=head2 Gene cluster info from ITEP

CLdb_dna_segs_make.pl -d CLdb.sqlite -sub I-A  -I DATABASE.sqlite all_I_2.0_c_0.4_m_maxbit

=head2 No broken loci

CLdb_dna_segs_make.pl -da CLdb.sqlite -sub I-A -q "AND loci.operon_status != 'broken'"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

