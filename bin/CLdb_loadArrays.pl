#!/usr/bin/env perl

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
use CLdb::utilities qw/
	file_exists 
	connect2db
	lineBreaks2unix
	get_file_path/;
use CLdb::query qw/
	table_exists/;
use CLdb::load qw/
	load_db_table/;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
my $db_path = get_file_path($database_file);
my $array_path = "$db_path/array/";
$array_path = File::Spec->rel2abs($array_path);
die "ERROR: cannot find '$array_path'\n" unless -d $array_path;

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# getting array file names from database #
my $array_r = get_array_file_names($dbh);

# loading array tables #
my %arrays;
foreach my $array_file (@$array_r){
	get_array_table($array_path, $$array_file[1], $$array_file[0],\%arrays);
	}
	
my ($dr_header_r, $sp_header_r) = make_headers();

# adding spacers/DR to database #
## DR entries ##
load_db_table($dbh, 'DRs', $dr_header_r, $arrays{'DR'});
## spacer entries ##
load_db_table($dbh, 'spacers', $sp_header_r, $arrays{'spacer'});

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

sub make_headers{
# headers for CLdb entry loading #
	my %dr_header = (
		locus_id => 1,
		dr_id => 2,
		dr_start => 3,
		dr_end => 4,
		dr_sequence => 5
		);
	my %sp_header = (
		locus_id => 1,
		spacer_id => 2,
		spacer_start => 3,
		spacer_end => 4,
		spacer_sequence => 5
		);

		#print Dumper %dr_header; exit;
	return \%dr_header, \%sp_header;
	}

sub get_array_table{
# getting array table info #
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
		my $uID = join("_", $locus_id, "dri.$repeat_cnt");
		$arrays_r->{"DR"}{$uID}{"locus_id"} = $locus_id;
		$arrays_r->{"DR"}{$uID}{"dr_id"} = "dr.$repeat_cnt";
		$arrays_r->{"DR"}{$uID}{"dr_start"} = $$pos_r[0];
		$arrays_r->{"DR"}{$uID}{"dr_end"} = $$pos_r[1];
		$arrays_r->{"DR"}{$uID}{"dr_sequence"} = $line[1];
		
		# loading spacer #
		if($$pos_r[2] && $$pos_r[3]){
			my $uID = join("_", $locus_id, "$repeat_cnt");
			$arrays_r->{"spacer"}{$uID}{"locus_id"} = $locus_id;
			$arrays_r->{"spacer"}{$uID}{"spacer_id"} = "sp.$repeat_cnt";
			$arrays_r->{"spacer"}{$uID}{"spacer_start"} = $$pos_r[2];
			$arrays_r->{"spacer"}{$uID}{"spacer_end"} = $$pos_r[3];
			$arrays_r->{"spacer"}{$uID}{"spacer_sequence"} = $line[2];			
			}
	
		# catting array sequence #
		$seq .= join("", @line[1..2]);				# total array sequence
		$seq =~ s/ +//g;		
		}
		
	close IN;

		#print Dumper %$arrays_r; exit;	
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



