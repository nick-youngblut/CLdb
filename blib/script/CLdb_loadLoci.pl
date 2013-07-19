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

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);

# getting input loci table #
my ($loci_r, $header_r) = load_loci_table();

# checking for required headers #
check_headers($header_r);

# striping off paths from file columns #
remove_paths($loci_r, $header_r);

# updating / loading_db #
load_new_entries($dbh, $loci_r->{"new_entry"}, $header_r);
update_db($dbh, $loci_r, $header_r);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
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
	
	print STDERR "...Number of entries updated in loci table: ", (scalar keys %$loci_r) -1, "\n"
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

	print STDERR "...Number of entries added to loci table: ", scalar @$loci_new_r, "\n"
		unless $verbose;
	}

sub check_headers{
	my ($header_r) = @_;
	my @req = qw/taxon_id taxon_name locus_start locus_end operon_status crispr_array_status genbank author/;
	
	my @not_found;
	foreach (@req){
		unless (exists $header_r->{$_}){
			push @not_found, $_;
			}
		}
	
	if(@not_found){
		print STDERR "ERROR: Required columns not found in loci table!\n\n";
		print STDERR "### Required columns not found (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort @not_found), "\n";
		print STDERR "\n### Headers found in the loci table (capitalization does not matter) ###\n";
		print STDERR join(",\n", sort keys %$header_r), "\n";
		exit;
		}
	
	}

sub load_loci_table{
	# checking line breaks #
	my @table = <>;
	map{ s/\r$//; s/\r/\n/g; push @table, split /\n/;  } @table;
	shift @table;

	# loading into a hash #
	my %loci;
	my %header;
	my $line_cnt = 0;
	foreach (@table){
		chomp;
		$line_cnt++;
		next if /^\s*$/;

		if($line_cnt == 1){ # loading header
			tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
			my @line = split /\t/;
			for my $i (0..$#line){
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
	return (\%loci, \%header);;
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


__END__

=pod

=head1 NAME

CLdb_loadLoci.pl -- adding/updating loci entries in to CRISPR_db

=head1 SYNOPSIS

CLdb_loadLoci.pl [flags] < loci_table.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_loadLoci.pl

=head1 DESCRIPTION

Load and/or update CRISPR loci table entries.

Loci entries will be added if no Locus_id is provided;
otherwise, the entrie will be updated.

=head2 WARNING!

The loci table must be in tab-delimited format.

File names in the loci table should match files in the 
$CLdb_HOME/genbank/ & $CLdb_HOME/array/ directories.

=head1 EXAMPLES

=head2 Usage:

CLdb_loadLoci.pl -d CRISPR.sqlite < loci.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

