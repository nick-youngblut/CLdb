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


my ($verbose, $database_file, $spacer_bool, $by_group);
my (@subtype, @taxon_id, @taxon_name);
my @subject_in;
my $blast_params = "-evalue 0.00001";
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "blast=s" => \$blast_params,
	   "subject=s{,}" => \@subject_in,		# blast subject
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: provide >= 1 blast subject fasta file\n"
	unless @subject_in;
	
my $db_path = get_database_path($database_file);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# making blast directory #
my $blast_dir = make_blast_dir($db_path);

# getting array commands #
my $spacer_fasta = call_CLdb_array2fasta($database_file, \@subtype, \@taxon_id, \@taxon_name, $extra_query, "spacer", $blast_dir);
my $DR_fasta = call_CLdb_array2fasta($database_file, \@subtype, \@taxon_id, \@taxon_name, $extra_query, "repeat", $blast_dir);

# blasting #
foreach my $subject (@subject_in){
	my $blast_db = make_blast_db($subject, $blast_dir);
	spacer_blast($spacer_fasta, $blast_params, $blast_db, $blast_dir);
	}


### Subroutines
sub spacer_blast{
	my ($spacer_fasta, $blast_params, $blast_db, $blast_dir) = @_;
	
	my $cmd = "blastn -task 'blastn-short' -outfmt 6 -db $blast_db -query $spacer_fasta > $blast_dir/spacer_blast.txt $blast_params";
		#print Dumper $cmd; exit;
	system("$cmd");
	}

sub make_blast_db{
# making blastdb for subject #
	my ($subject, $blast_dir) = @_;
	
	my @parts = File::Spec->splitpath(File::Spec->rel2abs($subject));
	(my $out = $parts[2]) =~ s/\.[^.]+$|$/_blast.db/;
	$out = join("/", $blast_dir, $out);
	
	my $cmd = "makeblastdb -dbtype nucl -in $subject -out $out";
	`$cmd`;
	
	return $out;
	}

sub call_CLdb_array2fasta{
# calling CLdb_array2fasta to get spacers to blast #
# using grouped spacers #
	my ($database_file, $subtype, $taxon_id, $taxon_name, $extra_query, $element, $blast_dir) = @_;
	
	my $out = "$blast_dir/spacers.fna";
	
	# adding flags #
	my $cmd = "CLdb_array2fasta.pl -database $database_file -g";	
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


__END__

=pod

=head1 NAME

template.pl -- script template

=head1 SYNOPSIS

template.pl [options] < input > output

=head2 options

=over

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc template.pl

=head1 DESCRIPTION

The flow of execution is roughly:
   1) Step 1
   2) Step 2
   3) Step 3

=head1 EXAMPLES

=head2 Usage method 1

template.pl <read1.fastq> <read2.fastq> <output directory or basename>

=head2 Usage method 2

template.pl <library file> <output directory or basename>

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

