#!/usr/bin/env perl

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


my ($verbose, $database_file, $spacer_bool, $dr_bool, $path, $debug);
my $cluster = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "spacer" => \$spacer_bool,
	   "repeat" => \$dr_bool,
	   "cluster=f" => \$cluster,
	   "path=s" => \$path,
	   "verbose" => \$verbose,
	   "z" => \$debug,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
print STDERR " WARNING: neither '-s' nor '-r' flags used. No grouping will be performed!\n"
	unless $spacer_bool || $dr_bool;
$database_file = File::Spec->rel2abs($database_file);
$path = File::Spec->rel2abs($path) if $path;


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "spacers") if $spacer_bool;
table_exists($dbh, "DRs") if $dr_bool;

# making a grouping directory; getting current directory #
my $dir = make_group_dir($database_file, $path); 
my $curdir = File::Spec->rel2abs(File::Spec->curdir());

# making fastas of array sequences #
## spacers ##
my %aliases;
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
	}

# call cd-hit-est #
chdir $dir or die $!;
my $spacer_cdhit_out = call_cdhit("spacers.fna", $cluster) if $spacer_bool;
my $dr_cdhit_out = call_cdhit("DRs.fna", $cluster) if $dr_bool;

# parse cd-hit-est output #
my $spacer_clust_r = parse_cdhit($spacer_cdhit_out, $cluster, \%aliases, 'spacer') if $spacer_bool;
my $dr_clust_r = parse_cdhit($dr_cdhit_out, $cluster, \%aliases, 'DR') if $dr_bool;

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
		$cmd = "UPDATE DRs SET DR_group = ? WHERE locus_id = ? and DR_id = ?";
		}
	else{ die " LOGIC ERROR: $!\n"; }
	
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $row (@$clust_r){
		# @$row = locus,spacer/DR,ID,groupID,seq
			#print Dumper @$row; exit;
		map{ undef $_ if $_ eq 'NA' } @$row;
		$sql->execute( $$row[3], $$row[0], $$row[2] );
		
		if($DBI::err){
			my $entry = join("|", @$row);
			print STDERR "ERROR: $DBI::errstr for $entry\n";
			}
		else{ $cnt++; }
		}
		
	$dbh->commit;
	print STDERR "...$cat group info added to $cnt entries\n" unless $verbose;
	}
 
sub parse_cdhit{
# parsing the cdhit results #
	my ($cdhit_out, $cluster, $aliases_r, $cat) = @_;

	open IN, "$cdhit_out.clstr" or die $!;
	
	my @clusters;			
	my $cluster_id;
	my %clust_chk;
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
			die " ERROR: cluster '$cluster_id' has SeqID < cutoff of '$cluster'\n"
				if $line[4] / 100 < $cluster;
			}
			
		# getting ID (aliase -> id) #
		$line[2] =~ s/^>|\.{3}$//g;	
		die "ERROR: cannot find $line[2] in alias list!\n"
			unless exists $aliases_r->{$line[2]};
		
		${$aliases_r->{$line[2]}}[3] = $cluster_id;
		push @clusters, $aliases_r->{$line[2]};
		$clust_chk{$line[2]} = 1;
		}
	close IN;

	# checking to see if all clusters are accounted for #
	my @not_found;
	foreach (keys %$aliases_r){
		push @not_found, $_ unless exists $clust_chk{$_};
		}
	if(@not_found){
		print STDERR "ERROR: some sequences do not have a cluster:\n";
		print STDERR join(",\n", @not_found), "\n";
		exit(1);
		}
		
		#exit;
	return \@clusters;
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
	if($debug){ `perl $FindBin::RealBin/CLdb_array2fasta.pl -d $database_file $flag > $outfile`; }
	else{ `CLdb_array2fasta.pl -d $database_file $flag > $outfile`; }
	
	
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
	
sub write_array_seq{
# writing arrays as fasta
	my ($arrays_r, $aliases_r, $cat) = @_;

	open OUT, ">$cat.fna" or die $!;

	my $alias_cnt = 0;
	foreach my $row (@$arrays_r){	
		$alias_cnt++;
			
		$aliases_r->{$alias_cnt} = $row;
		print OUT join("\n", ">$alias_cnt", $$row[$#$row]), "\n";
		}
	close OUT;
	
		#print Dumper %$aliases_r; exit;
	}
	


__END__

=pod

=head1 NAME

CLdb_groupArrayElements.pl -- group spacers and DRs by 100% sequence ID; add to database

=head1 SYNOPSIS

CLdb_groupArrayElements.pl [flags] 

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

CD-HIT-EST cluster cutoff. [1]

=item -path  <char>

Directory where intermediate files are written. [$CLdb_HOME/grouping/]

=item -verbose  <bool>	

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_groupArrayElements.pl

=head1 DESCRIPTION

Group the spacers and/or direct repeats in the CRISPR
database using CD-HIT-EST and add the group ID of
each spacer/DR to the CRISPR database. 

Spacer and DR fasta files and CD-HIT-EST files
are written to '$CLdb_HOME/grouping/' by default.

Sequences must be the same length to be in the same group
(cd-hit-est -s 1).

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

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

