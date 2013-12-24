#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use File::Path;
use File::Copy;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $quiet);
GetOptions(
	   "database=s" => \$database_file,
	   "quiet" => \$quiet,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;

exit;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);
my $column_list_r = list_columns("loci", 1);

#print Dumper $column_list_r; exit;

# getting input loci table #
my ($loci_r, $header_r) = load_loci_table($column_list_r);

# checking for required headers #
check_headers($header_r);

# adding scaffold values if not found #
add_scaffold_ID($loci_r, $header_r);

# copying array files & genbanks if not in ./genbank & ./array #
exit;
make_genbank_array_dirs($loci_r, $header_r);

# striping off paths from file columns #
	#remove_paths($loci_r, $header_r);

# updating / loading_db #
load_new_entries($dbh, $loci_r->{"new_entry"}, $header_r);
update_db($dbh, $loci_r, $header_r);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub make_genbank_array_dirs{
	my ($loci_r, $header_r) = @_;
	
	# current directory #
	my $dir = File::Spec->rel2abs(File::Spec->curdir());

	print Dumper $loci_r; exit;

	#print STDERR "...Removing any paths in file names\n" unless $verbose;
	
	my @cp_warning;
	foreach my $entry_type (keys %$loci_r){
		foreach my $row (@{$loci_r->{$entry_type}}){
			# genbank file #
			if( $$row[$header_r->{"genbank_file"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"genbank_file"}] );

				# making genbank dir; copying file #
				if(File::Spec->rel2abs($parts[1]) ne "$dir/genbank"){
					mkdir "$dir/genbank" unless -d "$dir/genbank";
					unless(-e "$dir/genbank/$parts[2]"){
						die " ERROR: cannot find ", $$row[$header_r->{"genbank_file"}], "\n"
							unless -e $$row[$header_r->{"genbank_file"}];
							
						copy($$row[$header_r->{"genbank_file"}], "$dir/genbank/$parts[2]") or die $!;
						print STDERR "...Copied ", $$row[$header_r->{"genbank_file"}], 
							" to $dir/genbank/$parts[2]\n" unless $quiet;
						}
					}
				# stripping path from genbank value #
				$$row[$header_r->{"genbank_file"}] = $parts[2];
				# sanity check #
				my $file_chk = join("/", $dir, "genbank", $$row[$header_r->{"genbank_file"}]);
				die " ERROR: cannot find ", $file_chk, "\n"
					unless -e $file_chk;
				} 
			# array file #
			if( exists $header_r->{"array_file"} && $$row[$header_r->{"array_file"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"array_file"}] );

				# making genbank dir; copying file #
				if(File::Spec->rel2abs($parts[1]) ne "$dir/array"){
					mkdir "$dir/array" unless -d "$dir/array";
					unless(-e "$dir/array/$parts[2]"){
						die " ERROR: cannot find ", $$row[$header_r->{"array_file"}], "\n"
							unless -e $$row[$header_r->{"array_file"}];
							
						copy($$row[$header_r->{"array_file"}], "$dir/array/$parts[2]") or die $!;
						print STDERR "...Copied ", $$row[$header_r->{"array_file"}], 
							" to $dir/array/$parts[2]\n" unless $quiet;
						}
					}
				# stripping path from array value #
				$$row[$header_r->{"array_file"}] = $parts[2];
				# sanity check #
				my $file_chk = join("/", $dir, "array", $$row[$header_r->{"array_file"}]);
				die " ERROR: cannot find ", $file_chk, "\n"
					unless -e $file_chk;
				}
			# fasta file #
			if( exists $header_r->{"fasta_file"} && $$row[$header_r->{"fasta_file"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"fasta_file"}] );

				# making genbank dir; copying file #
				if(File::Spec->rel2abs($parts[1]) ne "$dir/fasta"){
					mkdir "$dir/fasta" unless -d "$dir/fasta";
					unless(-e "$dir/fasta/$parts[2]"){
						die " ERROR: cannot find ", $$row[$header_r->{"fasta_file"}], "\n"
							unless -e $$row[$header_r->{"fasta_file"}];
							
						copy($$row[$header_r->{"fasta_file"}], "$dir/fasta/$parts[2]") or die $!;
						print STDERR "...Copied ", $$row[$header_r->{"fasta_file"}], 
							" to $dir/fasta/$parts[2]\n" unless $quiet;
						}
					}
				# stripping path from array value #
				$$row[$header_r->{"fasta_file"}] = $parts[2];
				# sanity check #
				my $file_chk = join("/", $dir, "fasta", $$row[$header_r->{"fasta_file"}]);
				die " ERROR: cannot find ", $file_chk, "\n"
					unless -e $file_chk;
				}				
			}
		}
	}

sub remove_paths{
	my ($loci_r, $header_r) = @_;
	
	print STDERR "...Removing any paths in file names\n" unless $verbose;
	
	foreach my $entry_type (keys %$loci_r){
		foreach my $row (@{$loci_r->{$entry_type}}){
			if( $$row[$header_r->{"genbank"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"genbank"}] );
				$$row[$header_r->{"genbank"}] = $parts[2];
				} 
			if( $$row[$header_r->{"array_file"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"array_file"}] );
				$$row[$header_r->{"array_file"}] = $parts[2];			
				}
			if( $$row[$header_r->{"fasta_file"}] ){
				my @parts = File::Spec->splitpath( $$row[$header_r->{"fasta_file"}] );
				$$row[$header_r->{"fasta_file"}] = $parts[2];			
				}				
			}
		}
	}

sub update_db{
# updating any entries with CLI identifiers #
	my ($dbh, $loci_r, $header_r) = @_;
	
	# entries ordered #
	my @keys = keys %$header_r;
	@keys = grep(!/^locus_id$/i, @keys);
	my @values = @$header_r{@keys};

	# setting up update for all columns #
	my @set;
	foreach my $key (@keys){
		(my $tmp = $key) =~ s/$/ = ?/;
		push @set, $tmp;
		}

	# preparing sql #	
	my $cmd = join(" ", "UPDATE Loci SET", join(",", @set), "where locus_id = ?");
	my $sql = $dbh->prepare($cmd);

	# updating #
	foreach my $entry (keys %$loci_r){
		next if $entry eq "new_entry";
		my $row = $loci_r->{$entry};
		
		$sql->execute( (@$row[@values], $$row[$header_r->{"locus_id"}]) );
		
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$row[@values]), "\n";
			}
		}
	
	$dbh->commit;	
	
	print STDERR "...Number of entries updated in loci table (locus_ID provided in loci table): ", (scalar keys %$loci_r) -1, "\n"
		unless $verbose;
	}

sub load_new_entries{
# adding new entries to db #
	my ($dbh, $loci_new_r, $header_r) = @_;

	# making locus_id = NULL #
	my @keys = keys %$header_r;
	@keys = grep(!/^locus_id$/i, @keys);
	my @values = @$header_r{@keys};

	# loading entries #
	my $cmd = "INSERT INTO loci(" . join(",", @keys) . ") VALUES(?".",?"x$#keys . ")";
	my $sql = $dbh->prepare($cmd);
	foreach my $row (@$loci_new_r){
		# "" values changed to undef = null in sqlite #
		my @null = @$row;
		map{ undef $_ if $_ eq ""} @null;
		
		# loading #
		$sql->execute( @null[@values] );	
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: '", join("\t", @$row[@values]), "'\n";
			}
		}
	$dbh->commit;

	print STDERR "...Number of entries added in loci table (no locus_ID provided in loci table): ", scalar @$loci_new_r, "\n"
		unless $verbose;
	}

sub add_scaffold_ID{
# adding scaffold ID if not found; needed to make sure entries are unique #
	my ($loci_r, $header_r) = @_;
	
	$header_r->{"scaffold"} = scalar keys %$header_r unless exists $header_r->{"scaffold"};
	
	my $scaf_add_bool = 0;
	foreach my $entry_type (keys %$loci_r){
		foreach my $row (@{$loci_r->{$entry_type}}){
			unless( $$row[$header_r->{"scaffold"}] ){
				$$row[$header_r->{"scaffold"}] = "CLDB__ONE_CHROMOSOME";
				$scaf_add_bool = 1;
				}
			}
		}

	# status #
	print STDERR "...Did not find 'scaffold' values for some entries. Adding scaffold names as: 'CLDB__ONE_CHROMOSOME'\n"
		if $scaf_add_bool and ! $verbose;
	}

sub check_headers{
	my ($header_r) = @_;
	my @req = qw/taxon_id taxon_name locus_start locus_end operon_status array_status genbank_file array_file author/;
	
	# checking for required headers not found #
	my (@not_found, @found);
	foreach (@req){
		if (exists $header_r->{$_}){
			push @found, $_;
			}
		else{
			push @not_found, $_;
			}
		}
	
	if(@not_found){
		print STDERR "ERROR: Required columns not found in loci table!\n\n";
		print STDERR "### Required columns not found (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort @not_found), "\n";
		print STDERR "\n### Correct headers found in the loci table (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort @found), "\n";
		exit;
		}
	
	}

sub load_loci_table{
# loading tab-delimited table from STDIN #
	my $column_list_r = shift;
		
	# checking line breaks #
	my @table = <>;
	my @tmp;
	map{ s/\r$//; s/\r/\n/g; s/\n+/\n/g; push @tmp, split /\n/;  } @table;
	@table = @tmp;

	# loading into a hash #
	my %loci;
	my %header;
	my $line_cnt = 0;
	foreach (@table){
		chomp;
		$line_cnt++;
		next if /^\s*$/;

		if($line_cnt == 1){ 					# loading header
			tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
			tr/ /_/;
			my @line = split /\t/;
			for my $i (0..$#line){
				next unless exists $column_list_r->{$line[$i]};		# only including columns in loci DB table
				next if exists $header{$line[$i]}; 					# only using 1st column found with a particular name
				$header{$line[$i]} = $i;		# column_name => index
				}
			}
		else{
			my @line = split /\t/;
			
			# locus ID #
			if( exists $header{"locus_id"} && $line[$header{"locus_id"}] ){		# if updating lci
				$loci{$line[$header{"locus_id"}]} = \@line;
				}
			else{
				push (@{$loci{"new_entry"}}, \@line);
				}
			}
		}
		#print Dumper %loci; exit; 
		#print Dumper %header; exit;
	# sanity check #
	die " ERROR: entries found in loci table!\n" unless %loci;	
	
	return (\%loci, \%header);
	}

sub check_for_loci_table{
	my ($table_list_r) = @_;
	die " ERROR: loci table not found in database!\n"
		unless grep(/^loci$/i, @$table_list_r);
	}

sub list_tables{
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}

sub list_columns{
	my ($column_list, $silent_ret) = @_;
	my $all = $dbh->selectall_arrayref("pragma table_info($column_list)");
		#print Dumper $$all[0]; exit;
	my %tmp;
	foreach (@$all){ 
		$$_[1] =~ tr/A-Z/a-z/;		# lower case for matching
		$tmp{$$_[1]} = 1; 
		}
	if($silent_ret){ return \%tmp; }
	else{  print "Columns:\n", join(",\n", keys %tmp), "\n\n";  exit; }
	}


__END__

=pod

=head1 NAME

CLdb_loadLoci.pl -- adding/updating loci entries in to CRISPR_db

=head1 SYNOPSIS

CLdb_loadLoci.pl [flags] < loci_table.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_loadLoci.pl

=head1 DESCRIPTION

Load and/or update CRISPR loci table entries.

Loci entries will be added if no Locus_id is provided;
otherwise, the entrie will be updated.

=head2 WARNING!

The loci table must be in tab-delimited format.

Loci entries must be unique by 'taxon_name', 'taxon_id', 'scaffold', 'locus_start', 'locus_end'.
Only the 1st of duplicated entries will be kept.

A scaffold name will be added to any entries lacking one. 

File names in the loci table should match files in the 
$CLdb_HOME/genbank/ & $CLdb_HOME/array/ directories.

Extra columns (not in CLdb Loci table) can exist in the 
input loci table. They just won't be added to CLdb.

=head1 EXAMPLES

=head2 Usage:

CLdb_loadLoci.pl -d CLdb.sqlite < loci.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

