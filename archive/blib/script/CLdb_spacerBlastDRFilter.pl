#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Set::IntervalTree;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $array_bool, $database_file);
my $range = 30;
my $DR_cnt = 1;
my $full_length = 0.7;
GetOptions(
	   "database=s" => \$database_file,
	   "range=i" => \$range,			# range beyond hit to consider overlap
	   "full=f" => \$full_length, 		# DRs hit length must be x fraction of total query length
	   "DR=i" => \$DR_cnt,				# number of sides that a DR must hit (adjacent to spacer)
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

my $itrees_r = make_DR_itrees($dbh, $range);
filter_spacer_blast($dbh, $itrees_r);

### Subroutines
sub filter_spacer_blast{
	my ($dbh, $itrees_r) = @_;
	
	# selecting spacers #
	# selecting DRs from DB #
	my $q = "SELECT 
S_taxon_id, 
S_taxon_name, 
S_accession, 
sseqid,
sstart, 
send, 
qstart, 
qend, 
qlen, 
array_hit, 
blast_id
FROM blast_hits
WHERE spacer_DR='Spacer'";
	$q =~ s/\n/ /g;
	
	my $res = $dbh->selectall_arrayref($q);
	die " ERROR: no spacer blast hits found! Nothing to filter!\n" unless @$res;
	
	# updating sql #
	my $u = "UPDATE blast_hits
SET array_hit=?
WHERE blast_id=?";
	$u =~ s/\n/ /g;
	my $sth = $dbh->prepare($u);
	
	my %summary = ("proto" => 0, "array" => 0, "update" => 0);
	foreach my $line (@$res){
	
		map{$_ = "" unless $_} @$line[0..2];
		my $genome_id = join("_", @$line[0..2]);
		
		# check for existence of genome/scaffold #
		my @LR = (0,0); 				# DR hit to left & right of spacer?
		if( exists $itrees_r->{$genome_id}{$$line[3]} ){

			# querying itree #
			my $res;
			my ($qstart, $qend);
			if( $$line[4] <= $$line[5] ){ 
				$qstart = $$line[4];
				$qend = $$line[5];
				}
			else{
				$qstart = $$line[5];
				$qend = $$line[4];
				}
			$res = $itrees_r->{$genome_id}{$$line[3]}->fetch($qstart, $qend);
		
			# updating 'array_hit' value # spacer hit adjacent to DR hits? #

			foreach my $line (@$res){
					#print Dumper "$$line[0], $qstart, $$line[1], $qend";
				if($$line[0] <= $qstart){ $LR[0]++; }
				elsif( $$line[1] >= $qend){ $LR[1]++; }
				}
			}
			
		if($LR[0] > 0 && $LR[1] > 0){		# updating array_hit to 1
			$sth->bind_param(1, "yes");
			$summary{"array"}++;
			}
		elsif($DR_cnt == 1 && ($LR[0] > 0 || $LR[1] > 0) ){
			$sth->bind_param(1, "yes");
			$summary{"array"}++;		
			}
		else{								# updating array_hit to 0
			$sth->bind_param(1, "no");
			$summary{"proto"}++;			# protospacer hit
			}
		$sth->bind_param(2, $$line[10]);
		
		$sth->execute();
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join(",", @$line), "\n";
			}
		else{ $summary{"update"}++; }
		}
	$dbh->commit();
		
	# status #
	print STDERR "...Number of spacer blast hits identified as hitting an array: ", $summary{"array"}, "\n";
	print STDERR "...Number of spacer blast hits identified as hitting a protospacer: ", $summary{"proto"}, "\n";
	print STDERR "...Number of spacer blast hits updated in CLdb blast_hits table: ", $summary{"update"}, "\n";		
	}


sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}

sub make_DR_itrees{
	my ($dbh, $range) = @_;
	
	# selecting DRs from DB #
	my $q = "SELECT 
S_taxon_id, 
S_taxon_name, 
S_accession, 
sseqid, 
sstart, 
send,
qstart, 
qend, 
qlen,
Group_ID
FROM blast_hits
WHERE spacer_DR='DR'";
	$q =~ s/\n/ /g;

	my $res = $dbh->selectall_arrayref($q);
	die " ERROR: no DR blast hits found! Cannot filter spacer blast hits!\n" unless @$res;	
	
	my %itrees;
	my %summary = ("short" => 0, "added" => 0);
	foreach my $line (@$res){
	
		# checking for full length hits #
		if(abs($$line[7] - $$line[6] + 1) / $$line[8] < $full_length){		# if short
			$summary{"short"}++;			
			next;
			}

		# genome_id: initializing itree #
		map{$_ = "" unless $_} @$line[0..2];
		my $genome_id = join("_", @$line[0..2]);
		unless( exists $itrees{$genome_id}{$$line[3]} ){  
			$itrees{$genome_id}{$$line[3]} = Set::IntervalTree->new();
			}
		
		# loading itree #
		if($$line[4] <= $$line[5]){		# start <= end
			my $xstart =  $$line[4] - $range -1;
			my $xend =  $$line[5] + $range + 1;
			$itrees{$genome_id}{$$line[3]}->insert([$xstart,$xend], $xstart, $xend);
			}
		else{
			my $xstart =  $$line[5] - $range -1;
			my $xend =  $$line[4] + $range + 1;
			$itrees{$genome_id}{$$line[3]}->insert([$xstart,$xend], $xstart, $xend);
			}
		$summary{"added"}++;
		}

	# status #
	print STDERR "...Number of DRs hits selected: ", scalar @$res, "\n";
	print STDERR "...Number of DRs hits < length cutoff: ", $summary{"short"}, "\n";
	print STDERR "...Number of DRs hits used for spacer filtering: ", $summary{"added"}, "\n";

		#print Dumper %itrees; exit;
	return \%itrees;
	}


__END__

=pod

=head1 NAME

CLdb_spacerBlastDRFilter.pl -- which spacer blast hits hit CRISPR arrays or protospacers?

=head1 SYNOPSIS

CLdb_spacerBlastDRFilter.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -range  <int>

Range allowable between spacer & DR blast hit (bp). [30]

=item -DR  <int>

Number of adjacent DR hits to a spacer to ID the spacer hit as an array hit. [2]

=item -full_length  <float>

Length cutoff for DR hit (DR_hit / DR_seq_length), (>=). [0.7]

=item -verbose  <bool>

Verbose output.

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlastDRFilter.pl

=head1 DESCRIPTION

Filter out spacer blast hits to CRISPR arrays by
removing all spacer blast hits that hit adjacent
to direct repeat blast hits.

All DR hits (or just full length hits) will be
used for the filtering process.

The 'array_hit' column in the CLdb blast_hits
table will be updated with either 'yes' or 'no'
if the spacer hit is determined to hit an
array or protospacer, respectively.

=head1 EXAMPLES

=head2 Basic Usage:

CLdb_spacerBlastDRFilter.pl -d CLdb.sqlite

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

