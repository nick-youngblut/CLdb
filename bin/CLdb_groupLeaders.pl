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
my $leader_fasta = call_leader2fasta($dir, $database_file);

# call cd-hit-est #
my $leader_cdhit_out = call_cdhit($leader_fasta, $cluster);

# parse cd-hit-est output #
my $leader_clust_r = parse_cdhit($leader_cdhit_out, $cluster);

# updating db #
update_db($dbh, $leader_clust_r);

# disconnect #
$dbh->disconnect();
exit;

### Subroutines
sub update_db{
# update group column #
	my ($dbh, $clust_r) = @_;

	my $cmd = "UPDATE Leaders SET leader_group = ? WHERE locus_id = ?";
	
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $locus_id (keys %$clust_r){
		(my $locus_tmp = $locus_id) =~ s/cli\.//;
		$sql->execute( $clust_r->{$locus_id}, $locus_tmp );
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", $clust_r->{$locus_id}, $locus_id), "\n";
			}
		else{ $cnt++; }
		}
	$dbh->commit;
	print STDERR "...$cnt groups added\n" unless $verbose;
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
		$clusters{$name_parts[0]} = $cluster_id;
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

sub call_leader2fasta{
# calling leader2fasta to make fasta sequences of array elements #
	my ($dir, $database_file) = @_;
	
	my $outfile = "$dir/leaders.fna"; 
	`CLdb_leader2fasta.pl -d $database_file -g > $outfile`;
	
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

CLdb_groupLeaders.pl -- group leader sequences at 100% sequence ID; add to database

=head1 SYNOPSIS

CLdb_groupLeaders.pl [flags] 

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -s 	Cluster spacers

=item -c 	CD-HIT-EST cluster cutoff. [1]

=item -p 	Directory where intermediate files are written. [$CLdb_HOME/grouping/]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_groupLeaders.pl

=head1 DESCRIPTION

Group the leader sequences in the CRISPR
database using CD-HIT-EST and add the group ID of
each leader sequence to the CRISPR database.

All leader sequences are grouped by default.

Leader sequence fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/grouping/' by default.

=head2 Requires:

cd-hit-est, CLdb_leader2fasta.pl

=head1 EXAMPLES

=head2 Grouping all leaders

CLdb_groupLeaders.pl -d CRISPR.sqlite -s

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

