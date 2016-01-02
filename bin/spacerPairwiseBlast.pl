#!/usr/bin/env perl

=pod

=head1 NAME

spacerPairwiseBlast.pl -- pairwise BLASTn of spacers and loading into CLdb

=head1 SYNOPSIS

spacerPairwiseBlast.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -blast  <float>

BLASTn parameters (besides required flags). [-evalue 0.00001]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb --perldoc -- spacerPairwiseBlast.pl

=head1 DESCRIPTION

Perform pairwise BLASTn-short of all (or some) spacers
in CLdb. Only the top blast hit is kept ('top' defined by
percent-identity, then length). 

Blast hits are loaded into the 'spacer_pairwise_blast' table in CLdb. 
Entries overlapping (same query & subject locus_ID and spacer_ID)
will be overwritten.

Temporary blast files will be written in $HOME/spacer_blast/

=head1 EXAMPLES

=head2 Pairwise spacer BLASTn-short of all spacers

CLdb -- spacerPairwiseBlast -d CLdb.sqlite

=head2 Pairwise spacer BLASTn-short for 1 subtype

CLdb -- spacerPairwiseBlast -da CLdb.sqlite -subtype I-B

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	get_file_path/;
use CLdb::seq qw/
	revcomp/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $spacer_bool, $by_group);
my (@subtype, @taxon_id, @taxon_name);
my @subject_in;
my $debug;
my $blast_params = "-evalue 0.00001";
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "blast=s" => \$blast_params,
	   "verbose" => \$verbose,
	   "x" => \$debug, 
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

my $db_path = get_file_path($database_file);


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# making blast directory #
my $blast_dir = make_blast_dir($db_path);

# getting spacers for blasting #
my $spacer_fasta = call_array2fasta($database_file, \@subtype, \@taxon_id, \@taxon_name, $extra_query, "spacer", $blast_dir);

# making blastdb #
make_blastdb($spacer_fasta);

# pairwise blast #
my $blast_res_r = pairwise_spacer_blast($spacer_fasta, $blast_params, $blast_dir);

# loading DB #
add_entry($dbh, $blast_res_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub add_entry{
# adding gentry for blast subject #
	my ($dbh, $blast_res_r) = @_;
	
	# cmd #
	my @cols = qw/Query_locus_ID Query_spacer_ID Subject_locus_ID Subject_spacer_ID pident length mismatch gapopen qstart qend sstart send evalue bitscore/;
	my $cmd = join(" ", "INSERT INTO spacer_pairwise_blast (", 
				join(",", @cols), ") values (",
				join(",", ("?") x scalar @cols), ")");
	
	print STDERR $cmd, "\n" if $verbose;
	
	my $sql = $dbh->prepare($cmd);
	
	# loading each entry #
	my $entry_cnt = 0;
	foreach my $query (keys %$blast_res_r){
		foreach my $subject (keys %{$blast_res_r->{$query}}){
			my @query = split /__/, $query;
			my @subject = split /__/, $subject;
			push(my @total, @query, @subject, @{$blast_res_r->{$query}{$subject}});
			die " ERROR: not enough values in entry!\n"
				unless scalar @cols == scalar @total;
			$sql->execute( @total );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", $query, $subject), "\n";
				}
			else{ $entry_cnt++; }
			}
		}
	
	$dbh->commit;
	
	print STDERR "...Number of entries added/updated to spacer_pairwise_blast: $entry_cnt\n";
	}

sub pairwise_spacer_blast{
	my ($spacer_fasta, $blast_params, $blast_dir) = @_;

	my $cmd = "blastn -task 'blastn-short' -outfmt 6 -query $spacer_fasta -db $spacer_fasta $blast_params";
	print STDERR "$cmd\n" if $verbose;
	open PIPE, "$cmd |" or die $!;

	my %blast_res;
	while(<PIPE>){
		chomp;
		my @line = split /\t/;
		
		# taking best hit (percentID, then length) #
		if(exists $blast_res{$line[0]}{$line[1]}){
			$blast_res{$line[0]}{$line[1]} = [@line[2..$#line]]
				if $line[2] >= ${$blast_res{$line[0]}{$line[1]}}[0] &&		# percent_ID
					$line[3] >= ${$blast_res{$line[0]}{$line[1]}}[1];		# length
			}
		else{
			$blast_res{$line[0]}{$line[1]} = [@line[2..$#line]];
			}
		}
	close PIPE;
		#print Dumper %blast_res; exit;
	return \%blast_res;
	}

sub make_blastdb{
# making blastdb for spacers #
	my ($spacer_fasta) = @_;
	
	my $cmd = "makeblastdb -dbtype nucl -in $spacer_fasta";
	`$cmd`;
	}

sub call_array2fasta{
# calling array2fasta to get spacers to blast #
# using grouped spacers #
	my ($database_file, $subtype, $taxon_id, $taxon_name, $extra_query, $element, $blast_dir) = @_;
	
	my $out = "$blast_dir/spacers.fna";
	
	# adding flags #
	my $cmd = "array2fasta.pl -database $database_file";	
	$cmd = append2cmd($cmd, "-subtype", $subtype) if @$subtype;
	$cmd = append2cmd($cmd, "-taxon_id", $taxon_id) if @$taxon_id;
	$cmd = append2cmd($cmd, "-taxon_name", $taxon_name) if @$taxon_name;
	$cmd = join(" ", $cmd, "-extra_query", "\"$extra_query\"") if $extra_query;
	$cmd = join(" ", $cmd, "-r") if $element eq "repeat";
	$cmd = join(" ", $cmd,  ">$out");
		
		#print Dumper $cmd; exit;
	`$cmd`;
	
	return $out;
	}

sub append2cmd{
	my ($cmd, $flag, $params) = @_;
	map{ $_=~ s/(.+)/"$1"/ } @$params;
	$cmd = join(" ", $cmd, "$flag", @$params);
	return $cmd;
	}

sub make_blast_dir{
	my $db_path = shift;
	
	my $dir = join("/", $db_path, "spacer_blast");
	mkdir $dir unless -d $dir;

	return $dir;
	}

sub get_database_path{
	my $database_file = shift;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return $parts[1];
	}


