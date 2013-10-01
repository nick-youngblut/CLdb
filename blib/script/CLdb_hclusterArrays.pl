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


my ($verbose, $database_file, $spacer_bool, $dr_bool, $path);
my $cluster = 1;
my @cluster = (0.8, 1.0, 0.01);
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f{,}" => \@cluster,
	   "path=s" => \$path,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# required input #
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
print STDERR " WARNING: neither '-s' nor '-r' flags used. No grouping will be performed!\n"
	unless $spacer_bool || $dr_bool;

# path #
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;

# cluster #
my $cluster_r = check_cluster(\@cluster);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# making a grouping directory #
my $dir = make_group_dir($database_file, $path); 

# make fasta 
my $spacer_fasta = call_array2fasta($dir, $database_file, "") if $spacer_bool;
my $dr_fasta = call_array2fasta($dir, $database_file, "-r") if $dr_bool;

# 'hcluster' #
my %spacer_hclust;
my %DR_hclust;
for (my $i=$$cluster_r[0]; $i<=$$cluster_r[1] + $$cluster_r[2]; $i+=$$cluster_r[2]){
	print STDERR "...clustering at cutoff $i\n" unless $verbose;
	
	## call cd-hit-est ##
	my $spacer_cdhit_out = call_cdhit($spacer_fasta, $i) if $spacer_bool;
	my $dr_cdhit_out = call_cdhit($dr_fasta, $i) if $dr_bool;

	## parse cd-hit-est output ##
	$spacer_hclust{$i} = parse_cdhit($spacer_cdhit_out, $i) if $spacer_bool;
	$DR_hclust{$i} = parse_cdhit($dr_cdhit_out, $i) if $dr_bool;
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
	
	foreach my $cutoff (keys %$clust_r){
		foreach my $locus_id (keys %{$clust_r->{$cutoff}}){
			foreach my $x_id ( keys %{$clust_r->{$cutoff}{$locus_id}} ){	# spacer/DR_id => cluster_id		
				$sql->execute($locus_id, $x_id, $cutoff, $clust_r->{$cutoff}{$locus_id}{$x_id} );
				if($DBI::err){
					print STDERR "ERROR: $DBI::errstr in: ", join("\t", $x_id, $cutoff, $clust_r->{$cutoff}{$locus_id}{$x_id} ), "\n";
					}
				}
			}
		}
	$dbh->commit;

	print STDERR "...$cat clusters added/updated\n" unless $verbose;
	}
 
sub parse_cdhit{
# parsing the cdhit results #
	my ($cdhit_out, $cluster) = @_;

	open IN, "$cdhit_out.clstr" or die $!;
	
	my %clusters;			# file=>spacer/dr=>cluster_ID
	my $cluster_id;
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
			}
			
		# getting file an spacer number #
		$line[2] =~ s/^>|\.\.\.$//g;			# name 
		my @name_parts = split /__/, $line[2];
		$name_parts[0] =~ s/^cli\.//i;
		
		# length #
		$line[1] =~ s/nt//g;
	
		# loading hashes #
		$clusters{$name_parts[0]}{$name_parts[1]} = $cluster_id + 1;			# groups start at 1
		}
	close IN;

		#print Dumper %clusters; exit;
	return \%clusters;
	}

sub call_cdhit{
# calling cd-hit-est on spacers #
	my ($fasta, $cluster) = @_;
	
	(my $cdhit_out = $fasta) =~ s/\.fna/.txt/;
	my $cmd = "cd-hit-est -i $fasta -o $cdhit_out -c $cluster -n 8 -s 1";
	if($verbose){ system("$cmd"); }
	else{ `$cmd`; }

	return $cdhit_out;
	}

sub call_array2fasta{
# calling array2fasta to make fasta sequences of array elements #
	my ($dir, $database_file, $flag) = @_;
	
	my $outfile;
	if($flag){ $outfile = "$dir/direct_repeats.fna"; }
	else{ $outfile = "$dir/spacers.fna"; }
	`CLdb_array2fasta.pl -d $database_file $flag > $outfile`;
	
	return $outfile;
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


__END__

=pod

=head1 NAME

CLdb_hclustArrays.pl -- group spacers and DRs by 100% sequence ID; add to database

=head1 SYNOPSIS

CLdb_hclustArrays.pl [flags] 

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

=item -verbose  <char>

Verbose output. [TRUE]

=item -help  <char>

This help message

=back

=head2 For more information:

perldoc CLdb_hclustArrays.pl

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

=head2 Requires:

cd-hit-est, CLdb_array2fasta.pl

=head1 EXAMPLES

=head2 Clustering spacers

CLdb_hclustArrays.pl -d CRISPR.sqlite -s

=head2 Clustering spacers & DRs

CLdb_hclustArrays.pl -d CRISPR.sqlite -s -r

=head2 Clustering spacers & DRs (only at 0.9 & 1 cutoffs)

CLdb_hclustArrays.pl -d CRISPR.sqlite -s -r -c 0.9 1 0.1


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

