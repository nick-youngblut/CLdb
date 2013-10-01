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
my $array_path = path_by_database($database_file);
$array_path = File::Spec->rel2abs($array_path);


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# getting array file names from database #
my $array_r = get_array_file_names($dbh);

# loading array tables #
my %arrays;
foreach my $array_file (@$array_r){
	load_array_table($array_path, $$array_file[1], $$array_file[0],\%arrays);
	}

# adding spacers/DR to database #
add_entries($dbh, \%arrays, "spacer");
add_entries($dbh, \%arrays, "DR");

# disconnect #
$dbh->disconnect();
exit;


### Subroutines 
sub add_entries{
# adding entries to spacer & dr tables #
	my ($dbh, $arrays_r, $cat) = @_;
	
	
	my $cmd;
	if($cat eq "spacer"){
		$cmd = "INSERT INTO Spacers(Locus_ID, Spacer_ID, Spacer_start, Spacer_end, Spacer_sequence) values (?,?,?,?,?)";
		}
	elsif($cat eq "DR"){
		$cmd = "INSERT INTO DRs(Locus_ID, DR_ID, DR_start, DR_end, DR_sequence) values (?,?,?,?,?)";
		}
	else{ die " LOGIC ERROR: $!\n"; }
	
	my $sql = $dbh->prepare($cmd);
	
	foreach my $locus_id (keys %$arrays_r){
		foreach my $x_id (keys %{$arrays_r->{$locus_id}{$cat}}){
			$sql->execute( ($locus_id, $x_id, @{$arrays_r->{$locus_id}{$cat}{$x_id}}) );
			if($DBI::err){
				print STDERR "ERROR: $DBI::errstr in: ", join("\t", $locus_id, $x_id, @{$arrays_r->{$locus_id}{$cat}{$x_id}}), "\n";
				}
			}
		}
	$dbh->commit;

	print STDERR "...$cat entries added to database\n" unless $verbose;
	}

sub load_array_table{
# loading array table #
	my ($array_path, $array_file, $locus_id, $arrays_r) = @_;
	
	die " ERROR: '$array_path/$array_file' not found!\n"
		unless -e "$array_path/$array_file";
	open IN, "$array_path/$array_file" or die $!;
	
	my $seq;
	my $repeat_cnt = 0;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		s/^\s+//;
		my @line = split /\t/;
		map{$_ =~ s/\s+//g} @line;
		
		# sanity check #
		die " ERROR: table not formatted correctly (>=4 columns required)\n"
			if scalar(@line) < 4;
		die " ERROR: $array_file not formatted correctely!\n"
			unless $line[0] =~ /^\d+$/;
		die " ERROR: $array_file not formatted correctely!\n"
			unless $line[1] =~ /^[A-Z-]+$/;
		
		# repeat count #
		$repeat_cnt++;
		
		# getting positions #
		my $pos_r = determine_positions(\@line, 1, $seq);
	
		# loading hash #
		$arrays_r->{$locus_id}{"DR"}{"dri.$repeat_cnt"} = [ $$pos_r[0], $$pos_r[1], $line[1] ];		# repeat=>ID=>start=>end => cluster (always1)
		
		# loading spacer #
		if($$pos_r[2] && $$pos_r[3]){
			$arrays_r->{$locus_id}{"spacer"}{"sid.$repeat_cnt"} = [ $$pos_r[2], $$pos_r[3], $line[2] ];
			}
	
		# catting array sequence #
		$seq .= join("", @line[1..2]);				# total array sequence
		$seq =~ s/ +//g;		
		}
		
	close IN;

		#print Dumper %array; exit;	
	}

sub determine_positions{
	my ($line_r, $location, $seq) = @_;
	
	my @pos;
	my $seq_length = length($seq);
	$seq_length = 0 if ! $seq_length;
	
	$$line_r[0] = $seq_length + 1 unless $location; 				# start = sequence length
		
	$pos[0] = $$line_r[0] if $$line_r[1];													# repeat start
	#$pos[1] = $$line_r[0] + length($$line_r[1]) - 1 if $$line_r[1];						# repeat end
	$pos[1] = $$line_r[0] + length($$line_r[1]) if $$line_r[1];								# repeat end
	$pos[2] = $$line_r[0] + length($$line_r[1]) if $$line_r[2]; 							# spacer start
	#$pos[3] = $$line_r[0] + length($$line_r[1]) + length($$line_r[2]) -2 if $$line_r[2];	# spacer end
	$pos[3] = $$line_r[0] + length($$line_r[1]) + length($$line_r[2]) if $$line_r[2];		# spacer end
	
		die "HERE" unless $$line_r[0];
	
	return \@pos;
	}

sub get_array_file_names{
# querying genbank names from sqlite loci table #
	my ($dbh) = @_;
	
	my $cmd = "SELECT locus_id, array_file from loci where array_file is not null";
	my $names_r = $dbh->selectall_arrayref($cmd);
	
	die " ERROR: no array file names found!\n" unless $names_r;
		#print Dumper $names_r; exit;
	return $names_r;
	}

sub path_by_database{
	my ($database_file) = @_;
	$database_file = File::Spec->rel2abs($database_file);
	my @parts = File::Spec->splitpath($database_file);
	return join("/", $parts[1], "array");
	}

__END__

=pod

=head1 NAME

CLdb_loadArrays.pl -- loading direct repeats & spacers into CRISPR database

=head1 SYNOPSIS

CLdb_loadArrays.pl [flags] 

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

This help message.

=back

=head2 For more information:

perldoc CLdb_loadArrays.pl

=head1 DESCRIPTION

Parse array files (CRISPRFinder format) and load
the spacers and direct repeats into the CLdb database.

Array file names are obtained from the loci table in the
CRISPR database.

Array files must be in $CLdb_HOME/array/

=head1 EXAMPLES

=head2 Usage:

CLdb_loadArrays.pl -d CLdb.sqlite 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

