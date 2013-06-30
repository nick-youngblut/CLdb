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
my $trim_length = 0;
GetOptions(
	   "database=s" => \$database_file,
	   "trim=i" => \$trim_length,
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

# load alignment #
my $fasta_r = load_fasta();

# getting leader info from fasta names #
my $leader_info_r = parse_names($fasta_r); 


# determining side trimmed #
my $side_r = which_side($leader_info_r);


# updating / loading_db #
#load_new_entries($dbh, $loci_r->{"new_entry"}, $header_r);
#update_db($dbh, $loci_r, $header_r);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub which_side{
# which side is farthest from array? need to account for rev-comp during alignment #
# mafft provides rev-comp name 
	my ($leader_info_r) = @_;

	foreach my $leader (@$leader_info_r){
		
		}
	}

sub parse_names{
# parsing leader names to get leader start-stop info #
	my ($fasta_r) = @_;

	my %leader_info;
	foreach my $name (keys %$fasta_r){
		my @parts = split /  |__/, $name;
		die " ERROR: leader sequence name not formatted correctly!\n"
			unless scalar @parts == 4;
		#$parts[0] =~ s/>cli\.//;
		$leader_info{$name} = \@parts;
		}
		print Dumper %leader_info; exit;
	return \%leader_info;
	}

sub load_fasta{
# loading fasta #
	my (%fasta, $tmpkey);
	while(<>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
		#print Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta


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
		$sql->execute( @$row[@values] );	
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$row[@values]), "\n";
			}
		}
	$dbh->commit;

	print STDERR "...Number of entries added to loci table: ", scalar @$loci_new_r, "\n"
		unless $verbose;
	}

sub load_loci_table{
	my %loci;
	my %header;
	while(<>){
		chomp;
		next if /^\s*$/;

		if($. == 1){ # loading header
			tr/A-Z/a-z/; 						# all to lower case (caps don't matter)
			my @line = split /\t/;
			for my $i (0..$#line){
				$header{$line[$i]} = $i;		# column_name => index
				}
			}
		else{
			my @line = split /\t/;
			if( exists $header{"locus_id"} && $line[$header{"locus_id"}] ){		# if updating lci
				$loci{$line[$header{"locus_id"}]} = \@line;
				}
			else{
				push (@{$loci{"new_entry"}}, \@line);
				}
			}
		}
		#print Dumper %loci; exit; 
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

CLdb_loadLeader.pl -- adding/updating leader entries in to CRISPR_db

=head1 SYNOPSIS

CLdb_loadLeader.pl [options] < loci_table.txt

=head2 options

=over

=item -d 	CLdb database (required).

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_loadLeader.pl

=head1 DESCRIPTION

Load and/or update CRISPR loci table entries.

Loci entries will be added if no Locus_id is provided;
otherwise, the entrie will be updated.

=head1 EXAMPLES

=head2 Usage:

CLdb_loadLeader.pl -d CRISPR.sqlite < loci.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

