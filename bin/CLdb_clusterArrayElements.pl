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

=item -length  <float> 

Length difference cutoff used by CD-HIT-EST. [1]

=item -threads  <int>

Number of threads used by CD-HIT-EST. [1]

=item -clear  <bool>

Clear the cluster table(s) of all entries prior to clustering? [TRUE]

=item -verbose  <char>

Verbose output. [TRUE]

=item -help  <char>

This help message

=back

=head2 For more information:

perldoc CLdb_clusterArrayElements.pl

=head1 DESCRIPTION

Clustering of spacers and/or
direct repeats in the CRISPR database.
CD-HIT-EST is used to cluster the array elements
at different similarity cutoffs (default: 0.8 to 1 by 0.01 steps).

=head2 Method

=over

=item Spacers or DRs selected from CLdb & written to fasta

=item Fasta is clustered with cd-hit-est at various clustering cutoffs.
This is strand-specific. By default sequences must be the same length 
(edit with -length).

=item cd-hit-est output is parsed and the clusterIDs for each spacer|DR
at each clustering cutoff is entered into CLdb.

=back

Temporary spacer and DR fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/clustering/' by default.

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
		    n_entries
		    get_array_seq_preCluster
		  /;
use CLdb::utilities qw/
			file_exists 
			connect2db/;
use CLdb::seq qw/
		  revcomp
		  write_fasta
		/;
use CLdb::cluster qw/
		      cluster_cutoffs
		      insert_cluster
		      clear_cluster_table
		    /;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $path, @cluster_cut);
my ($spacer_bool, $dr_bool, $clear_entries);
my $cluster = 1;
my $threads = 1;
my $length = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f{,}" => \@cluster_cut,
	   "path=s" => \$path,
	   "clear" => \$clear_entries,  # deleting all entries currently in CLdb
	   "threads=i" => \$threads,
	   "length=f" => \$length,
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

# deleting any existing entries if requested
unless($clear_entries){
  clear_cluster_table($dbh, 'spacer_clusters') if $spacer_bool;
  clear_cluster_table($dbh, 'DR_clusters') if $dr_bool;  
}

# making a clustering directory #
my $dir = make_cluster_dir($database_file, $path); 

# making fastas of array sequences #
## status ##
print STDERR "Getting array element sequences. Orienting sequences by array_sense_strand\n";
## spacers ##
my ($spacer_fasta, $dr_fasta);
my (%aliases, $invariant); #debug
my %fasta_files;
if($spacer_bool){
  my %opts = (
	      refine_sql => "",
	      spacer_DR_b => 0,
	      verbose => $verbose
	     );
  my $fasta_r = get_array_seq_preCluster($dbh,
					  { refine_sql => "",
					    spacer_DR_b => 0,
					    verbose => $verbose
					  } ); 
  $fasta_files{spacer} = "spacers.fna";
  write_fasta(fasta => $fasta_r, dir => $dir, file => ">$dir/$fasta_files{spacer}");
}
if($dr_bool){  ## DRs 
  my $fasta_r = get_array_seq_preCluster($dbh, 
					  { refine_sql => "",
					    spacer_DR_b => 1,
					    verbose => $verbose
					    } ); 
  $fasta_files{DR} = "DRs.fna";
  write_fasta(fasta => $fasta_r, dir => $dir, file => ">$dir/$fasta_files{DR}");
}


# clustering at each cutoff #
## strand-specific clustering (sequence orientation based array_sense_strand) ##
chdir $dir or die $!;  # moving to 'clustering' directory for clustering step
my $spacer_cluster_r;
if($spacer_bool){
  print STDERR "\nClustering spacer sequences...\n";
  $spacer_cluster_r = cluster_cutoffs( cluster_cutoff => $cluster_cut_r, 
				       fasta_file => $fasta_files{spacer}, 
				       element => 'spacers', 
				       threads => $threads,
				       length => $length,
				       verbose => $verbose
				     );
}
my $DR_cluster_r;
if($dr_bool){
  print STDERR "\nClustering DR sequences...\n";
  $DR_cluster_r = cluster_cutoffs( cluster_cutoff => $cluster_cut_r, 
				   fasta_file => $fasta_files{DR},
				   element => 'DRs',  
 				   threads => $threads,
				   length => $length,
				   verbose => $verbose
				 );
}
chdir $curdir or die $!;

## updating db ##
print STDERR "\nInserting/updating entries in CLdb...\n";
insert_cluster(dbh => $dbh, 
	       cluster => $spacer_cluster_r, 
	       element => "spacer",
	       verbose => $verbose 
	      ) if $spacer_bool;
insert_cluster(dbh => $dbh, 
	       cluster => $DR_cluster_r, 
	       element =>"DR",
	       verbose => $verbose
	      ) if $dr_bool;


# disconnect #
$dbh->disconnect();
exit;


#--- Subroutines ---#
sub check_cluster{
# setting defaults for clustering if not provided #
  my $cluster_r = shift;
  
  $$cluster_r[0] = 0.8 unless defined $$cluster_r[0];
  $$cluster_r[1] = 1.0 unless defined $$cluster_r[1];
  $$cluster_r[2] = 0.01 unless defined $$cluster_r[2];
  
  return $cluster_r;
}


sub make_cluster_dir{
  # making a directory for grouping files
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

