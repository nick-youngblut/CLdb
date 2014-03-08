#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_clusterArrayElements.pl -- cluster spacers and DRs at various sequence identity cutoffs

=head1 SYNOPSIS

CLdb_clusterArrayElements.pl [flags] 

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -spacer  <bool>

Cluster spacers. [FALSE]

=item -repeat  <bool>

Cluster direct repeats. [FALSE]

=item -cluster  <float>

CD-HIT-EST cluster cutoff range. [0.80 1.00 0.01]

=item -path  <char>

Directory where intermediate files are written. [$CLdb_HOME/clustering/]

=item -threads <int>

Number of threads used by CD-HIT-EST. [1]

=item -verbose  <char>

Verbose output. [TRUE]

=item -help  <char>

This help message

=back

=head2 For more information:

perldoc CLdb_clusterArrayElements.pl

=head1 DESCRIPTION

Pseudo-hierarchical clustering of spacers and/or
direct repeats in the CRISPR database.
CD-HIT-EST is used to cluster the array elements
at different similarity cutoffs (default: 0.8 to 1 by 0.01 steps).

The spacer/DR cluster IDs are added to
the spacer_cluster & DR_cluster tables in CLdb.

Temporary spacer and DR fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/clustering/' by default.

Sequences must be the same length to be in the same cluster
(cd-hit-est -s 1). 

=head2 Array elements will be clustered in 2 ways:

=over 

=item Strand-by-leader: 

elements must be on the same strand (+/+ or -/-) with strand defined
by leader position. If no leader defined, either strand can be used.

=item Strand-agnostric: 

elements must be different strands (+/+ or -/- or +/- or -/+) and still
fall into the same cluster

=back

=head2 Requires:

cd-hit-est

=head1 EXAMPLES

=head2 Clustering spacers

CLdb_clusterArrayElements.pl -d CLdb.sqlite -s

=head2 Clustering spacers & DRs

CLdb_clusterArrayElements.pl -d CLdb.sqlite -s -r

=head2 Clustering spacers & DRs (only at 0.9 & 1 cutoffs)

CLdb_clusterArrayElements.pl -d CLdb.sqlite -s -r -c 0.9 1 0.1


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use List::Util qw/max/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
use CLdb::seq qw/
	revcomp/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $spacer_bool, $dr_bool, $path, @cluster_cut);
my $cluster = 1;
my $threads = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f{,}" => \@cluster_cut,
	   "path=s" => \$path,
	   "threads=i" => \$threads,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
print STDERR " WARNING: neither '-s' nor '-r' flags used. No clustering will be performed!\n"
	unless $spacer_bool || $dr_bool;

# path #
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;
my $curdir = File::Spec->rel2abs(File::Spec->curdir());

# cluster #
my $cluster_cut_r = check_cluster(\@cluster_cut);

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# making a clustering directory #
my $dir = make_cluster_dir($database_file, $path); 

# making fastas of array sequences #
## status ##
print STDERR "Getting array element sequences. Orienting sequences by leader or loci_start-end\n";
## spacers ##
my (%aliases, $spacer_fasta, $dr_fasta);
if($spacer_bool){
	my %opts = (
		extra_query => "",
		join_sql => "",
		spacer_DR_b => 0,
		leader => 1,
		 );
	my $arrays_r = get_array_seq_preCluster($dbh,\%opts); 

	chdir $dir or die $!;
	write_array_seq($arrays_r, \%aliases, 'spacers');
    chdir $curdir or die $!;
	$spacer_fasta = "spacers.fna";
	}
if($dr_bool){
	my %opts = (
		extra_query => "",
		join_sql => "",
		spacer_DR_b => 1,
		leader => 1
		);
	my $arrays_r = get_array_seq_preCluster($dbh,\%opts); 
	
	chdir $dir or die $!;
	write_array_seq($arrays_r, \%aliases, 'DRs');
	chdir $curdir or die $!;
	$dr_fasta = "DRs.fna";
	}

# clustering at each cutoff #
## strand-specific clustering (sequence orientation based on leader or loci orientation ##
chdir $dir or die $!;
my $spacer_cluster_r = cluster_cutoffs($cluster_cut_r, $spacer_fasta, \%aliases, 'spacers', $threads)
		if $spacer_bool;
my $DR_cluster_r = cluster_cutoffs($cluster_cut_r, $dr_fasta, \%aliases, 'DRs',  $threads)
		if $dr_bool;		
chdir $curdir or die $!;

## updating db ##
add_entry($dbh, $spacer_cluster_r, "spacer") if $spacer_bool;
add_entry($dbh, $DR_cluster_r, "DR") if $dr_bool;


# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub add_entry{
# add_entry to *_cluster #
	my ($dbh, $clust_r, $cat) = @_;

	# sql #
	my $cmd = "INSERT into $cat\_clusters(locus_id, $cat\_id, cutoff, cluster_id, Rep_sequence) values(?,?,?,?,?)";
	my $sql = $dbh->prepare($cmd);
	
	my $entry_cnt = 0;
	foreach my $cutoff (keys %$clust_r){
		foreach my $entry ( @{$clust_r->{$cutoff}}){
			
			$sql->execute($$entry[0], $$entry[2], $cutoff, $$entry[4], $$entry[5]);
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$entry ), "\n";
				}
			else{ $entry_cnt++; }
			}
		}
	$dbh->commit;

	print STDERR "...$entry_cnt $cat cluster entries added/updated\n" unless $verbose;
	}

sub cluster_cutoffs{
# clustering at each cutoff #
	my ($cluster_cut_r, $fasta, $aliases_r, $cat, $threads) = @_;
	
	# status #
	print STDERR "Clustering array sequences...\n";
	
	my %cluster;
	my %cluster_seq;
	for (my $i=$$cluster_cut_r[0]; 
			$i<=$$cluster_cut_r[1]; 		# rounding issues
			$i+=$$cluster_cut_r[2]){
		$i = sprintf('%.2f', $i);
		print STDERR "...Clustering $cat at cutoff: $i\n" unless $verbose;
	
		## call cd-hit-est ##
		my $cdhit_out = call_cdhit($fasta, $i, $threads);

		## parse cd-hit-est cluster representative sequence output ##
		parse_cdhit_seq(\%cluster_seq, $cdhit_out);

		## parse cd-hit-est cluster output ##
		$cluster{$i} = parse_cdhit_clstr($cdhit_out, $i, $aliases_r, $cat, \%cluster_seq);
		}
		
		#print Dumper %cluster; exit;
	return \%cluster;
	}

sub check_cluster{
# setting defaults for clustering if not provided #
	my $cluster_r = shift;
	
	$$cluster_r[0] = 0.8 unless defined $$cluster_r[0];
	$$cluster_r[1] = 1.0 unless defined $$cluster_r[1];
	$$cluster_r[2] = 0.01 unless defined $$cluster_r[2];
	
	return $cluster_r;
	}

sub parse_cdhit_clstr{
# parsing the cdhit results #
	my ($cdhit_out, $cluster, $aliases_r, $cat, $cluster_seq_r) = @_;

	open IN, "$cdhit_out.clstr" or die $!;
	
	my @clusters;			
	my $cluster_id;
	my %clust_chk;
	my %clust_cnt;
	while(<IN>){
		chomp;
		if(/^>/){
			($cluster_id = $_) =~ s/>.+ //g;
			next;
			}
			
		my @line = split /[\t ]+/;		# 1=number; 2=length; 3=name; 4=[at|*]; 5=%id
		
		# checking cluster ID #
		if($line[4]){
			
			$line[4] =~ s/[\+\-%\/]*//g;
			$line[4] = sprintf '%.3f', $line[4];
				
			die " ERROR: cluster '$cluster_id' has SeqID < cutoff of '$cluster'\n"
				if $line[4] / 100 < sprintf '%.3f', $cluster;
			}
			
		# getting ID (aliase -> id) #
		$line[2] =~ s/^>|\.{3}$//g;	
		die "ERROR: cannot find '$line[2]' in alias list!\n"
			unless exists $aliases_r->{$cat}{$line[2]};
		my @entry = @{$aliases_r->{$cat}{$line[2]}};
		
		# editing clusterID #
		my $cluster_id2 = $cluster_id;
		$cluster_id2 = join("_", $cluster, $cluster_id);

		# assigning cluster ID
		$entry[4] = $cluster_id2;
		
		# assigning rep sequence
		die "ERROR: cannot find rep-seq for clusterID: $cluster_id\n"
			unless exists $cluster_seq_r->{$cluster_id};
		$entry[5] = $cluster_seq_r->{$cluster_id};
		
		# loading & checking #
		push @clusters, \@entry;
		$clust_chk{$line[2]} = 1;
		
		$clust_cnt{$cluster_id} = 1;
		}
	close IN;

	# checking to see if all clusters are accounted for #
	my @not_found;
	foreach (keys %{$aliases_r->{$cat}}){
		push @not_found, $_ unless exists $clust_chk{$_};
		}
		
	if(@not_found){
		print STDERR "ERROR: some sequences do not have a cluster:\n";
		print STDERR join(",\n", @not_found), "\n";
		exit(1);
		}
	
	print STDERR "Number of clusters for $cat at cutoff '$cluster': ", scalar keys %clust_cnt, "\n";
	
		#print Dumper @clusters; 
	return \@clusters;
	}
 
sub parse_cdhit_seq{
# parsing the representative fasta sequence file from cdhit
	my ($clusters_seq_r, $cdhit_out) = @_;
	
	open IN, $cdhit_out or die $!;
	my $seq_cnt = -1;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		if(/^>/){
			$seq_cnt++;
			$clusters_seq_r->{$seq_cnt} = "";
			}
		else{
			$clusters_seq_r->{$seq_cnt} .= $_;
			}
		}
	close IN;
	
		#print Dumper %$clusters_seq_r; exit;
	}

sub call_cdhit{
# calling cd-hit-est on spacers/DRs #
	my ($fasta, $cluster, $threads) = @_;
	
	(my $cdhit_out = $fasta) =~ s/\.fna/.txt/;
	my $cmd = "cd-hit-est -i $fasta -o $cdhit_out -c $cluster -n 8 -s 1 -T $threads -r 0";
	`$cmd`;
	if ($? == -1) {
    	die "ERROR: failed to execute: $!\n";
		}
	elsif ($? & 127) {
    	die "ERROR: child died with signal %d, %s coredump\n",
        	($? & 127),  ($? & 128) ? 'with' : 'without';
		}
	else {
    	#die "ERROR: child exited with value %d\n", $? >> 8;
		}

	return $cdhit_out;
	}

sub make_cluster_dir{
# making a directory for grouping files	#
	my ($database_file, $path) = @_;
	my $dir;
	if($path){
		$dir = $path . "clustering/";
		}
	else{
		my @parts = File::Spec->splitpath($database_file);
		$dir = $parts[1] . "clustering/";
		}
	
	mkdir $dir unless -d $dir;
	
	return $dir;
	}

sub write_array_seq{
# writing arrays as fasta
	my ($arrays_r, $aliases_r, $cat) = @_;

	open OUT, ">$cat.fna" or die $!;

	my $alias_cnt = 0;
	foreach my $row (@$arrays_r){	
		$alias_cnt++;
			
		$aliases_r->{$cat}{$alias_cnt} = $row;
		print OUT join("\n", ">$alias_cnt", $$row[$#$row]), "\n";	# sequence & alias
		}
	close OUT;
	
		#print Dumper %$aliases_r; exit;
	}
	
sub get_array_seq_preCluster{
# getting array sequences and orienting by leader (if selected & if leader present) or by loci orientation
  my ($dbh, $opts_r) = @_;
  
  # checking for opts #
  map{ die "ERROR: cannot find option: '$_'" 
	 unless exists $opts_r->{$_} } qw/spacer_DR_b extra_query join_sql/;

  # getting table info (spacer|DR)#
  my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
  ($tbl_oi, $tbl_prefix) = ("DRs","DR") if defined $opts_r->{"spacer_DR_b"};

  my $query = "SELECT
b.Locus_ID, 
'$tbl_prefix', 
b.$tbl_prefix\_ID,
b.$tbl_prefix\_sequence,
a.leader_end,
a.array_start,
a.array_end,
a.locus_start,
a.locus_end
FROM 
	(SELECT * 
	FROM loci 
	LEFT OUTER JOIN leaders ON loci.locus_id=leaders.locus_id) a,
$tbl_oi b
WHERE a.locus_id = b.locus_id";

	$query =~ s/\n/ /g;
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	#print Dumper %$opts_r; exit;
	
	# making fasta #
	my @arrays;
	foreach my $row (@$ret){
		
		# if leader: revcomp if leader_end > max(array_end, array_start)
		if($opts_r->{leader} && defined $$row[4] && $$row[4] ne ""){
			$$row[3] = revcomp($$row[3]) if $$row[4] > max(@$row[5..6]);
			}
		else{ # if no leader: recomp if locus_start > locus_end
			$$row[3] = revcomp($$row[3]) if $$row[7] >= $$row[8];
			}
		
		@$row = @$row[0..3];
		}
		
		#print Dumper @$ret; exit;
	return $ret;
}

