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
my $range = 30;		# spacer-DR blast hit overlap (bp)
my $Ncpu = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "blast=s" => \$blast_params,
	   "subject=s{,}" => \@subject_in,		# blast subject
	   "range=i" => \$range,
	   "cpu=i" => \$Ncpu,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: provide a blast_subject file or 1 subject (& taxon_id) (& taxon_name)\n"
	unless @subject_in;
	
my $db_path = get_database_path($database_file);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# getting all taxon_name/taxon_id from loci table #
my $taxon_name_id_r = get_taxon_name_id($dbh);

# loading subject #
my $sub_in_r = load_subject(\@subject_in);

# making blast directory #
my $blast_dir = make_blast_dir($db_path);

# getting array commands #
my $spacer_fasta = call_CLdb_array2fasta($database_file, \@subtype, \@taxon_id, \@taxon_name, $extra_query, "spacer", $blast_dir);
my $DR_fasta = call_CLdb_array2fasta($database_file, \@subtype, \@taxon_id, \@taxon_name, $extra_query, "repeat", $blast_dir);

# blasting #
foreach my $subject (@$sub_in_r){
	print STDERR "### BLAST subject: \"$$subject[0]\" ###\n";
	
	# checking for existing taxon_name / taxon_in in CLdb #
	#next if check_taxon_name_id($taxon_name_id_r, \@subject_in);
	
	# change line breaks #
	conv_line_breaks($$subject[0]);
	
	# make blast db #
	my $blast_db = make_blast_db($subject, $blast_dir);
	
	# blasting #
	my $spacer_blast_out = spacer_blast($spacer_fasta, $blast_params, $blast_db, $blast_dir);
	my $DR_blast_out = DR_blast($DR_fasta, $blast_params, $blast_db, $blast_dir);

	# filtering blast #
	my $filt_blast_out = call_CLdb_spacerBlastDRFilter($blast_dir, $spacer_blast_out, $DR_blast_out, $range);
	
	# loading blast into db #
	call_CLdb_loadBlastHits($filt_blast_out, $database_file, $subject);
	}

### Subroutines
sub check_taxon_name_id{
	my ($taxon_name_id_r, $subject) = @_;
	
	unless ($$subject[1] || $$subject[2]){
		print STDERR " WARNING: no Taxon_ID or Taxon_Name provided for $$subject[0]! Skipping!\n";
		next;
		}
		
	unless (exists $taxon_name_id_r->{"taxon_id"}{$$subject[1]} ||
			exists $taxon_name_id_r->{"taxon_name"}{$$subject[2]} ){
		if($$subject[2]){		# taxon_name
			print STDERR " WARNING: '$$subject[2]' does not exist in CLdb Loci table! Skipping! Use '-v' to get a list of taxon_name values in the Loci table\n";	
			print join("\n", "### taxon_name values in loci table ###", keys %{$taxon_name_id_r->{"taxon_name"}}), "\n" if $verbose;
			return 1;
			}
		elsif($$subject[1]){	# taxon_id
			print STDERR " WARNING: '$$subject[1]' does not exist in CLdb Loci table! Skipping! Use '-v' to get a list of taxon_name values in the Loci table\n";
			print join("\n", "### taxon_id values in loci table ###", keys %{$taxon_name_id_r->{"taxon_id"}}), "\n" if $verbose;
			return 1;
			}
		else{ die " LOGIC ERROR: $!\n"; }
		}
	return 0;
	}

sub get_taxon_name_id{
	my ($dbh) = @_;
	my $q = "SELECT taxon_name, taxon_id from loci";
	my $res = $dbh->selectall_arrayref($q);
	
	my %taxon_name_id;
	foreach my $row (@$res){
		$taxon_name_id{"taxon_name"}{$$row[0]} = 1 if $$row[0];
		$taxon_name_id{"taxon_id"}{$$row[1]} = 1 if $$row[1];
		}
	
		#print Dumper %taxon_name_id; exit;
	return \%taxon_name_id;
	}

sub conv_line_breaks{
	my ($subject_fasta) = @_;
	die " ERROR: $subject_fasta not found!\n" unless -e $subject_fasta;
	
	# windows to unix #
	my $cmd = "perl -pi -e 's/\\r\$//' $subject_fasta";
	`$cmd`;

	# mac to unix #
	$cmd = "perl -pi -e 's/\\r//' $subject_fasta";
	`$cmd`;	
	}

sub call_CLdb_loadBlastHits{
	my ($filt_blast_out, $database_file, $subject) = @_;

	my $cmd = "perl ~/perl/projects/CLdb/bin/CLdb_loadBlastHits.pl -database $database_file -subject $$subject[0] < $filt_blast_out";
	#my $cmd = "CLdb_loadBlastHits.pl -database $database_file -subject $$subject[0] < $filt_blast_out";
	$cmd = join(" ", $cmd, "-taxon_id", $$subject[1]) if $$subject[1];
	$cmd = join(" ", $cmd, "-taxon_name", $$subject[2]) if $$subject[2];
	system($cmd);
	}

sub call_CLdb_spacerBlastDRFilter{
# calling CLdb_spacerBlastDRFilter for adding spacer blasts to database #
	my ($blast_dir, $spacer_blast_out, $DR_blast_out, $range) = @_;
	
	(my $out = $spacer_blast_out) =~ s/\.[^.]+$|$/_filt.txt/;
	my $cmd = "CLdb_spacerBlastDRFilter.pl -a -r $range $spacer_blast_out $DR_blast_out > $out";
		#print Dumper $cmd; 
	`$cmd`;
	
	return $out;
	}

sub DR_blast{
	my ($DR_fasta, $blast_params, $blast_db, $blast_dir) = @_;
	
	my $out = "$blast_dir/DR_blast.txt";
		#my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen";
	my $cmd = "blastn -task 'blastn-short' -outfmt 6 -db $blast_db -query $DR_fasta -num_threads $Ncpu > $out $blast_params";
		#print Dumper $cmd; exit;
	`$cmd`;
	
	return $out;
	}

sub spacer_blast{
	my ($spacer_fasta, $blast_params, $blast_db, $blast_dir) = @_;

	my $out = "$blast_dir/spacer_blast.txt";
		#my $outfmt = "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen";
	my $cmd = "blastn -task 'blastn-short' -outfmt 6 -db $blast_db -query $spacer_fasta -num_threads $Ncpu > $out $blast_params";
		#print Dumper $cmd; exit;
	`$cmd`;
	
	return $out;
	}

sub make_blast_db{
# making blastdb for subject #
	my ($subject, $blast_dir) = @_;
	
	die " ERROR: $$subject[0] not found!\n" unless -e $$subject[0];
	
	my @parts = File::Spec->splitpath(File::Spec->rel2abs($$subject[0]));
	(my $out = $parts[2]) =~ s/\.[^.]+$|$/_blast.db/;
	$out = join("/", $blast_dir, $out);
	
	my $cmd = "makeblastdb -dbtype nucl -in $$subject[0] -out $out";
	`$cmd`;
	
	return $out;
	}

sub load_subject{
# loading subject, either a file or 1-3 arguments #
# return @@: subject_fasta_file, taxon_id, taxon_name #
	my ($subject_in_r) = @_;
	
	unless($$subject_in_r[1]){
		my @tmp = $$subject_in_r[0];
		$subject_in_r = \@tmp;
		}
	
	if(scalar @$subject_in_r == 1 && -e $$subject_in_r[0]){		# possibly file
		open IN, $$subject_in_r[0] or die $!;
		my @subjects;
		
		while(<IN>){
			chomp;
			if($.==1 && /^>/){		# fasta; no subject file 
				return [$subject_in_r];
				}
			
			next if /^\s*$/;
			my @line = split /\t/;
			push @subjects, \@line;
			}
		return \@subjects;
		}
	else{
		return [$subject_in_r];
		}

	}

sub call_CLdb_array2fasta{
# calling CLdb_array2fasta to get spacers to blast #
# using grouped spacers #
	my ($database_file, $subtype, $taxon_id, $taxon_name, $extra_query, $element, $blast_dir) = @_;
	
	my $out = "$blast_dir/$element.fna";
	
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

CLdb_spacerBlast.pl -- wrapper for spacer blasting

=head1 SYNOPSIS

CLdb_spacerBlast.pl [flags]

=head2 Required flags

=over

=item -database

CLdb database.

=item -subject

Either subject file or 1-3 arguments (see DESCRIPTION).

=back

=head2 Optional flags

=over

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine which sequences are returned.

=item -blast

BLASTn parameters (besides required flags). [-evalue 0.00001]

=item -range

Range allowable between spacer & DR blast hit (bp). [30]

=item -cpu

Use for '-num_threads' parameter in BLASTn. [1]

=item -v	Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlast.pl

=head1 DESCRIPTION

A wrapper around other CLdb scripts for blasting all or
a subset of spacer groups.

=head2 The script's procedure is:

=over

=item * Select spacer & direct repeat (DR) groups (can refine to particular taxa & subtypes)

=item * BLASTn-short of spacer & DR groups against provide subjects (e.g. genomes)

=item * Determining which spacer blast hits are hitting CRISPR arrays (if DRs hit adjacent)

=item * Adding the blast hits and subjects (sequences with a hit) to CLdb

=back

=head2 '-subject' flag

Provide either a subject fasta file and a Taxon_ID and/or a Taxon_Name,
or provide a tab-delimited file (3 columns) with subject fasta files and
a Taxon_IDs and/or Taxon_Names (columns: fasta_file, Taxon_id, Taxon_Name).

Example1 "-subject ecoli.fna 666666.452 escherichia_coli"

Example2 "-subject ecoli.fna '' escherichia_coli"

The Taxon_IDs and Taxon_names can be used, for example, to see if CRISPR
spacers are hitting other places in the same genome (and not in other CRISPR arrays).

=head1 EXAMPLES

=head2 Spacer blast of subtype I-B I-C spacers against Ecoli

CLdb_spacerBlast.pl -da CLdb.sqlite -subtype I-B I-C -subject ecoli.fna 666666.452 Escherichia_coli

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

