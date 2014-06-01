#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_clusterLeaders.pl -- cluster leader sequences & adding to CLdb

=head1 SYNOPSIS

CLdb_clusterLeaders.pl [flags] 

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -cluster  <float>

CD-HIT-EST cluster cutoff. [1]

=item -path  <char>

Directory where intermediate files are written. [$CLdb_HOME/grouping/]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc CLdb_clusterLeaders.pl

=head1 DESCRIPTION

Cluster the leader sequences in the CRISPR
database using CD-HIT-EST and add the group ID of
each leader sequence to the CRISPR database.

All leader sequences are clustered by default.

Sequences must be the same length to be in the same cluster
(cd-hit-est -s 1).

Leader sequence fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/grouping/' by default.

=head2 Requires:

cd-hit-est, CLdb_leader2fasta.pl

=head1 EXAMPLES

=head2 Grouping all leaders

CLdb_clusterLeaders.pl -d CRISPR.sqlite -s

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
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

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

#--- I/O error & defaults ---#
file_exists($database_file, "database");
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# making a grouping directory #
my $dir = make_group_dir($database_file, $path); 

# make fasta 
my $leaders_r = get_leaders($dbh);
my ($name_index_r, $outfile) = write_leaders($leaders_r, $dir);

# call cd-hit-est #
my $leader_cdhit_out = call_cdhit($outfile, $cluster);

# parse cd-hit-est output #
my $leader_clust_r = parse_cdhit($leader_cdhit_out, $cluster, $name_index_r);

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
	print STDERR "...$cnt cluster IDs added\n" unless $verbose;
	}
 
sub parse_cdhit{
# parsing the cdhit results #
	my ($cdhit_out, $cluster, $name_index_r) = @_;

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
		## converting cluster_id (index) to locus_ID ##
		$line[2] =~ s/^>|\.\.\.$//g;			# name 
		
		die "ERROR: cluster_id: $cluster_id not found in name index!"
			unless exists $name_index_r->{$line[2]};
		$line[2] = $name_index_r->{$line[2]};
		
			
		# loading hashes #
		$clusters{$line[2]} = $cluster_id;
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

sub write_leaders{
# writing leaders w/ numeric index  (cdhit issue) #
	my ($leaders_r, $dir) = @_;
	
	my $outfile = "$dir/leaders.fna";
	open OUT, ">$outfile" or die $!;

	my %name_index;
	my $seq_cnt = 0;
	foreach my $locus (keys %$leaders_r){
		$seq_cnt++;
		print OUT join("\n", ">$seq_cnt", 
						$leaders_r->{$locus}{"Leader_sequence"}), "\n";
						
		$name_index{$seq_cnt} = $locus;
		}
	close OUT;
	
		#print Dumper %name_index; exit;
	return \%name_index, $outfile;
	}

sub get_leaders{
	my ($dbh, $extra_query) = @_;
	
	# make query #
	my $query = "SELECT Locus_ID, Leader_start, Leader_end, Leader_sequence FROM Leaders";
	$query .= " $extra_query" if defined $extra_query;
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	my %leaders;
	foreach my $row (@$ret){
		$$row[3] =~ s/-//g; 		# removing gaps
		$leaders{$$row[0]}{'Leader_start'} = $$row[1];
		$leaders{$$row[0]}{'Leader_end'} = $$row[2];
		$leaders{$$row[0]}{'Leader_sequence'} = $$row[3];
		}
	
		#print Dumper %leaders; exit;
	return \%leaders;
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

sub call_leader2fasta{
# calling leader2fasta to make fasta sequences of array elements (no gaps) #
	my ($dir, $database_file) = @_;
	
	my $outfile = "$dir/leaders.fna"; 
	`CLdb_leader2fasta.pl -d $database_file -g > $outfile`;
	
	return $outfile;
	}

