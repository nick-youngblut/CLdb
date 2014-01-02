#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_hclusterArrayElements.pl -- group spacers and DRs by 100% sequence ID; add to database

=head1 SYNOPSIS

CLdb_hclusterArrayElements.pl [flags] 

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

CD-HIT-EST cluster cutoff range. [0.8 1 0.01]

=item -path  <char>

Directory where intermediate files are written. [$CLdb_HOME/grouping/]

=item -threads <int>

Number of threads used by CD-HIT-EST. [1]

=item -verbose  <char>

Verbose output. [TRUE]

=item -help  <char>

This help message

=back

=head2 For more information:

perldoc CLdb_hclusterArrayElements.pl

=head1 DESCRIPTION

Pseudo-hierarchical clustering of spacers and/or
direct repeats in the CRISPR database.
CD-HIT-EST is used to cluster the array elements
at different similarity cutoffs (default: 0.8 to 1 by 0.01 steps).

The spacer/DR cluster IDs are added to
the spacer_hclust & directrepeat_hclust tables in CLdb.

Temporary spacer and DR fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/grouping/' by default.

Sequences must be the same length to be in the same group
(cd-hit-est -s 1). 

Array elements only clustered if same strand (+/+ or -/-)!

=head2 Requires:

cd-hit-est, CLdb_array2fasta.pl

=head1 EXAMPLES

=head2 Clustering spacers

CLdb_hclusterArrayElements.pl -d CRISPR.sqlite -s

=head2 Clustering spacers & DRs

CLdb_hclusterArrayElements.pl -d CRISPR.sqlite -s -r

=head2 Clustering spacers & DRs (only at 0.9 & 1 cutoffs)

CLdb_hclusterArrayElements.pl -d CRISPR.sqlite -s -r -c 0.9 1 0.1


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

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	get_array_seq/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $spacer_bool, $dr_bool, $path);
my $cluster = 1;
my $threads = 1;
my @cluster = (0.8, 1.0, 0.01);
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f{,}" => \@cluster,
	   "path=s" => \$path,
	   "threads=i" => \$threads,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
print STDERR " WARNING: neither '-s' nor '-r' flags used. No grouping will be performed!\n"
	unless $spacer_bool || $dr_bool;

# path #
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;
my $curdir = File::Spec->rel2abs(File::Spec->curdir());

# cluster #
my $cluster_r = check_cluster(\@cluster);


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# making a grouping directory #
my $dir = make_group_dir($database_file, $path); 

# making fastas of array sequences #
## spacers ##
my (%aliases, $spacer_fasta, $dr_fasta);
if($spacer_bool){
	my %opts = (
		extra_query => "",
		join_sql => "",
		spacer_DR_b => 0,
		);
	my $arrays_r = get_array_seq($dbh,\%opts); 
	
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
		);
	my $arrays_r = get_array_seq($dbh,\%opts); 
	chdir $dir or die $!;
	write_array_seq($arrays_r, \%aliases, 'DRs');
	chdir $curdir or die $!;
	$dr_fasta = "DRs.fna";
	}
	
# 'hcluster' #
chdir $dir or die $!;
my %spacer_hclust;
my %DR_hclust;
for (my $i=$$cluster_r[0]; $i<=$$cluster_r[1] + $$cluster_r[2]; $i+=$$cluster_r[2]){
	print STDERR "...clustering at cutoff $i\n" unless $verbose;
	
	## call cd-hit-est ##
	my $spacer_cdhit_out = call_cdhit($spacer_fasta, $i, $threads) if $spacer_bool;
	my $dr_cdhit_out = call_cdhit($dr_fasta, $i, $threads) if $dr_bool;

	## parse cd-hit-est output ##
	$spacer_hclust{$i} = parse_cdhit($spacer_cdhit_out, $i, \%aliases, 'spacers') if $spacer_bool;
	$DR_hclust{$i} = parse_cdhit($dr_cdhit_out, $i, \%aliases, 'DRs') if $dr_bool;
	}
		
# updating db #
add_entry($dbh, \%spacer_hclust, "spacer") if $spacer_bool;
add_entry($dbh, \%DR_hclust, "DR") if $dr_bool;


# disconnect #
$dbh->disconnect();
exit;

### Subroutines
sub check_cluster{
	my $cluster_r = shift;
	
	$$cluster_r[1] = 1 unless $$cluster_r[1];
	$$cluster_r[2] = 0.1 unless $$cluster_r[2];
	
	return $cluster_r;
	}

sub add_entry{
# add_entry to *_hcluster #
	my ($dbh, $clust_r, $cat) = @_;
	
	my $cmd = "INSERT into $cat\_hclust(locus_id, $cat\_id, cutoff, cluster_id) values(?,?,?,?)";
	my $sql = $dbh->prepare($cmd);
	
	my $entry_cnt = 0;
	foreach my $cutoff (keys %$clust_r){
		foreach my $entry ( @{$clust_r->{$cutoff}}){
			$sql->execute($$entry[0], $$entry[2], $cutoff, $$entry[3] );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$entry ), "\n";
				}
			else{ $entry_cnt++; }
			}
		}
	$dbh->commit;

	print STDERR "...$entry_cnt $cat cluster entries added/updated\n" unless $verbose;
	}
 
sub parse_cdhit{
# parsing the cdhit results #
	my ($cdhit_out, $cluster, $aliases_r, $cat) = @_;

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
		
		${$aliases_r->{$cat}{$line[2]}}[3] = $cluster_id;
		push @clusters, $aliases_r->{$cat}{$line[2]};
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
	
	print STDERR "Number of clusters for $cat at cutoff '$cluster': \t", scalar keys %clust_cnt, "\n";
	
		#print Dumper @clusters; exit;
	return \@clusters;
	}

sub call_cdhit{
# calling cd-hit-est on spacers #
	my ($fasta, $cluster, $threads) = @_;
	
	(my $cdhit_out = $fasta) =~ s/\.fna/.txt/;
	my $cmd = "cd-hit-est -i $fasta -o $cdhit_out -c $cluster -n 8 -s 1 -T $threads -r 0";
	if($verbose){ system("$cmd"); }
	else{ `$cmd`; }

	return $cdhit_out;
	}

sub make_group_dir{
# making a directory for grouping files	#
	my ($database_file, $path) = @_;
	my $dir;
	if($path){
		$dir = $path . "grouping/";
		}
	else{
		my @parts = File::Spec->splitpath($database_file);
		$dir = $parts[1] . "grouping/";
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
		print OUT join("\n", ">$alias_cnt", $$row[4]), "\n";
		}
	close OUT;
	
		#print Dumper %$aliases_r; exit;
	}

