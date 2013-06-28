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
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f" => \$cluster,
	   "path=s" => \$path,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;

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

# call cd-hit-est #
my $spacer_cdhit_out = call_cdhit($spacer_fasta, $cluster) if $spacer_bool;
my $dr_cdhit_out = call_cdhit($dr_fasta, $cluster) if $dr_bool;

# parse cd-hit-est output #
my $spacer_clust_r = parse_cdhit($spacer_cdhit_out, $cluster) if $spacer_bool;
my $dr_clust_r = parse_cdhit($dr_cdhit_out, $cluster) if $dr_bool;

# updating db #
update_db($dbh, $spacer_clust_r, "spacers") if $spacer_bool;
update_db($dbh, $dr_clust_r, "DR") if $dr_bool;

# disconnect #
$dbh->disconnect();
exit;

### Subroutines
sub update_db{
# update group column #
	my ($dbh, $clust_r, $cat) = @_;

	my $cmd;
	if($cat eq "spacers"){
		$cmd = "UPDATE spacers SET spacer_group = ? WHERE locus_id = ? and spacer_id = ?";
		}
	elsif($cat eq "DR"){
		$cmd = "UPDATE directrepeats SET repeat_group = ? WHERE locus_id = ? and repeat_id = ?";
		}
	else{ die " LOGIC ERROR: $!\n"; }
	
	my $sql = $dbh->prepare($cmd);
	
	foreach my $locus_id (keys %$clust_r){
		foreach my $x_id (keys %{$clust_r->{$locus_id}}){
			$sql->execute( ($clust_r->{$locus_id}{$x_id}, $locus_id, $x_id) );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", $clust_r->{$locus_id}{$x_id}, $locus_id, $x_id), "\n";
				}
			}
		}
	$dbh->commit;
	print STDERR "...$cat groups added\n" unless $verbose;
	}
 
sub parse_cdhit{
# parsing the cdhit results #
	my ($cdhit_out, $cluster) = @_;

	open IN, "$cdhit_out.clstr" or die $!;
	
	my %clusters;			# file=>spacer/dr=>cluster_ID
	my $cluster_id;
	while(<IN>){
			#print Dumper $_; next;
		chomp;
		if(/^>/){
			($cluster_id = $_) =~ s/>.+ //g;
			next;
			}
			
		my @line = split /[\t ]+/;		# 1=number; 2=length; 3=name; 4=[at|*]; 5=%id
		
		# checking cluster ID #
		if($line[4]){
			$line[4] =~ s/[\+\-%\/]*//g;
			die " ERROR: cluster '$cluster_id' has SeqID < cutoff of '$cluster'\n"
				if $line[4] / 100 < $cluster;
			}
			
		# getting file an spacer number #
		$line[2] =~ s/^>|\.\.\.$//g;			# name 
		my @name_parts = split /__/, $line[2];
		$name_parts[0] =~ s/^lci\.//i;
		
		# length #
		$line[1] =~ s/nt//g;
	
		# loading hashes #
		$clusters{$name_parts[0]}{$name_parts[1]} = $cluster_id;
		}
	close IN;

		#print Dumper %clusters; exit;
	return \%clusters;
	}

sub call_cdhit{
# calling cd-hit-est on spacers #
	my ($fasta, $cluster) = @_;
	
	(my $cdhit_out = $fasta) =~ s/\.fna/.txt/;
	my $cmd = "cd-hit-est -i $fasta -o $cdhit_out -c $cluster -n 8";
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

CLdb_groupArrayElements.pl -- group spacers and DRs by 100% sequence ID; add to database

=head1 SYNOPSIS

CLdb_groupArrayElements.pl [flags] 

=head2 Required flags

=over

=item -d 	CRISPR database.

=back

=head2 Optional flags

=over

=item -s 	Cluster spacers

=item -d 	Cluster direct repeats

=item -c 	CD-HIT-EST cluster cutoff. [1]

=item -p 	Directory where intermediate files are written. [./grouping/]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_groupArrayElements.pl

=head1 DESCRIPTION

Group the spacers and/or direct repeats in the CRISPR
database using CD-HIT-EST and add the group ID of
each spacer/DR to the CRISPR database.

Spacer and DR fasta files and CD-HIT-EST files
are written to './grouping/' by default.

=head2 Requires:

cd-hit-est, CLdb_array2fasta.pl

=head1 EXAMPLES

=head2 Grouping spacers

CLdb_groupArrayElements.pl -d CRISPR.sqlite -s

=head2 Grouping spacers & DRs

CLdb_groupArrayElements.pl -d CRISPR.sqlite -s -r

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CRISPR_db/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

