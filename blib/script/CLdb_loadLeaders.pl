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

my ($verbose, $database_file, $rm_gap);
my $trim_length = 0;
GetOptions(
	   "database=s" => \$database_file,
	   "trim=i" => \$trim_length,
	   "gap" => \$rm_gap,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: provide the original leader fasta file!\n"
	unless $ARGV[0];
die " ERROR: provide the aligned leader fasta file!\n"
	unless $ARGV[1];
map{die " ERROR: $_ not found!\n" unless -e $_} @ARGV[0..1];


### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);

# load alignment #
my $fasta_raw_r = load_fasta($ARGV[0]);
my $fasta_aln_r = load_fasta($ARGV[1]);

# determining orientation of aligned sequence #
my $ori_r = get_orientation($fasta_raw_r, $fasta_aln_r);

# getting array start-end #
##my $array_se_r = get_array_se($dbh, $ori_r);

# determining trim start-end #
trim_se($fasta_aln_r, $ori_r, $trim_length, $rm_gap);

# updating / loading_db #
load_leader($dbh, $fasta_aln_r);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub load_leader{
# loading leader info into db #
	my ($dbh, $fasta_aln_r) = @_;
	
	my $cmd = "INSERT INTO LeaderSeqs(Locus_ID, Leader_start, Leader_end, Leader_sequence) values (?,?,?,?)";
	my $sql = $dbh->prepare($cmd);
	
	my $cnt = 0;
	foreach my $locus (keys %$fasta_aln_r){
		(my $tmp_locus = $locus) =~ s/cli\.//;
		$sql->execute($tmp_locus, int ${$fasta_aln_r->{$locus}{"meta"}}[2], 
							int ${$fasta_aln_r->{$locus}{"meta"}}[3], 
							$fasta_aln_r->{$locus}{"seq"});
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for $locus\n";
			}
		else{ $cnt++; }
		}
	$dbh->commit;

	print STDERR "...$cnt entries added/updated in database\n";
	}

sub trim_se{
# modifying raw fasta to relect alignment #
## trimming sequence and modifying start-end locations accordingly ##
## trimming from leader end most distant from array ##
	my ($fasta_aln_r, $ori_r, $trim_length, $rm_gap) = @_;
	
		#print Dumper $fasta_aln_r; exit;
	
	foreach my $locus (keys %$fasta_aln_r){
		die " LOGIC ERROR: $locus not found in ori hash!\n"
			unless exists $ori_r->{$locus};
		
		# getting seq info #
		my $seq_len = length $fasta_aln_r->{$locus}{"seq"};
		
		# trimming sequences #
		if($ori_r->{$locus} eq "norm"){
			if(${$fasta_aln_r->{$locus}{"meta"}}[1] eq "start"){		# trimming from 5' end (+ strand)
				# trim seq #
				$fasta_aln_r->{$locus}{"seq"} = substr( $fasta_aln_r->{$locus}{"seq"}, $trim_length, $seq_len - $trim_length);
				# changing start position (moving up toward 3') #
				my $ngap = count_gaps($fasta_aln_r->{$locus}{"seq"});
				${$fasta_aln_r->{$locus}{"meta"}}[2] += ($trim_length + $ngap);
				}
			elsif( ${$fasta_aln_r->{$locus}{"meta"}}[1] eq "end" ){		# trimming from 3' end (+ strand)
				$fasta_aln_r->{$locus}{"seq"} = substr($fasta_aln_r->{$locus}{"seq"}, 0, $seq_len - $trim_length);
				# changing end position (moving back toward 5') #
				my $ngap = count_gaps($fasta_aln_r->{$locus}{"seq"});
				${$fasta_aln_r->{$locus}{"meta"}}[3] -= ($trim_length + $ngap);				
				}
			else{ die $!; }
			}
		elsif($ori_r->{$locus} eq "rev" || $ori_r->{$locus} eq "rev-comp"){
			if(${$fasta_aln_r->{$locus}{"meta"}}[1] eq "start"){		# trimming from 3' end (+ strand)
				$fasta_aln_r->{$locus}{"seq"} = substr($fasta_aln_r->{$locus}{"seq"}, 0, $seq_len - $trim_length);						
				# changing end position (moving back toward 5') #
				my $ngap = count_gaps($fasta_aln_r->{$locus}{"seq"});
				${$fasta_aln_r->{$locus}{"meta"}}[3] -= ($trim_length + $ngap);			
				}
			elsif( ${$fasta_aln_r->{$locus}{"meta"}}[1] eq "end" ){		# trimming from 5' end (+ strand)
				$fasta_aln_r->{$locus}{"seq"} = substr($fasta_aln_r->{$locus}{"seq"}, $trim_length, $seq_len - $trim_length);
				# changing start position (moving up toward 3') #
				my $ngap = count_gaps($fasta_aln_r->{$locus}{"seq"});
				${$fasta_aln_r->{$locus}{"meta"}}[2] += ($trim_length + $ngap);
				}		
			else{ die $!; }
			}
		else{ die $!; }
		
		# removing gaps if specified #
		$fasta_aln_r->{$locus}{"seq"} =~ s/-//g if $rm_gap;
		}
		
		#print Dumper %$fasta_aln_r; exit;
	}

sub count_gaps{
	my $seq = shift;
	my $ngap = 0;
	$ngap++ while $seq =~ /-/g;
	return $ngap;
	}

sub get_array_se{
# getting the array start-end from loci table #
	my ($dbh, $ori_r) = @_;
	
	my $cmd = "SELECT crispr_array_start, crispr_array_end FROM loci where locus_id = ?";
	my $sql = $dbh->prepare($cmd);

	my %array_se;
	foreach my $locus (keys %$ori_r){
		(my $tmp = $locus) =~ s/cli\.//;
		$sql->execute($tmp);
		my $ret = $sql->fetchall_arrayref();
		die " ERROR: no array start-end info in Loci table for $locus!\n"
			unless $$ret[0];
		
		$array_se{$locus} = $ret;
		}

		#print Dumper %array_se; exit;
	return \%array_se; 
	}
	
sub get_orientation{
# which side is farthest from array? need to account for rev-comp during alignment #
# mafft provides rev-comp name 
	my ($fasta_raw_r, $fasta_aln_r) = @_;

	my %ori;
	foreach my $aln (keys %$fasta_aln_r){
		die " ERROR: \"$aln\" not found in raw leader sequence file!\n"
			unless exists $fasta_raw_r->{$aln}; 
		# making permutations of aligned seq #
		my $perms_r = make_perms($fasta_aln_r->{$aln}{"seq"});
		
		if($fasta_raw_r->{$aln}{"seq"} =~ /$$perms_r[0]/i){		# same orientation
			print STDERR "$aln: normal orientation\n" if $verbose;
			$ori{$aln} = "norm";
			}
		elsif($fasta_raw_r->{$aln}{"seq"} =~ /$$perms_r[1]/i){	# reversed in alignment
			print STDERR "$aln: reverse orientation\n" if $verbose;
			$ori{$aln} = "rev";
			}
		elsif($fasta_raw_r->{$aln}{"seq"} =~ /$$perms_r[2]/i){	# rev-comp in alignment
			print STDERR "$aln: rev-comp orientation\n" if $verbose;		
			$ori{$aln} = "rev-comp";
			}
		else{ die " ERROR: could not match aligned sequence of $aln with raw sequence!\n"; }
		}
	return \%ori;
	}
	
sub make_perms{
# removing gaps; rev; rev-comp of sequence #
	my ($seq) = @_;
	(my $no_gap = $seq) =~ s/-//g;
	my $rev = reverse $no_gap;
	(my $rev_comp = $rev) =~ tr/ACGTN?acgtn/TGCAN?tgcan/;
	
	return [$no_gap, $rev, $rev_comp];
	}

sub load_fasta{
# loading fasta #
	my $fasta_in = shift;
	open IN, $fasta_in or die $!;	
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		
 		if(/>.+/){
 			my @line = split /  |__/;
 			$line[0] =~ s/.+(cli\.\d+).*/$1/;
 			die " ERROR: seq name: \"$_\" does not have identifiablle cli ID!\n"
 				unless $line[0] =~ /^cli\.\d+$/;
 			$fasta{$line[0]}{"meta"} = \@line;
 			
 			$tmpkey = $line[0];	# changing key
 			}
 		else{$fasta{$tmpkey}{"seq"} .= $_; }
		}
	close IN;
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

CLdb_loadLeaders.pl -- adding/updating leader entries in to CRISPR_db

=head1 SYNOPSIS

CLdb_loadLeaders.pl [flags] leader.fasta leader_aligned.fasta

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -t 	Number of bp (& gaps) to trim of end of alignment. [0]

=item -g 	Leave gaps in leader sequence entered into the DB? 

=item -v	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_loadLeaders.pl

=head1 DESCRIPTION

After determining the actual leader region by aligning
approximate leader regions written by CLdb_getLeaderRegions.pl,
automatically trim and load leader sequences & start-end
information into the CRISPR database.

Trimming will be done from the leader region end most
distant from the CRISPR array. Reversing & 
reverse-complementation of sequences during the alignment
will be accounted for (that's why both the aligned and
'raw' sequenced must be provided).

Not all sequences in the aligned fasta need to be in
the 'raw' leader fasta (eg. if both ends of the array
were written because the side containing the leader
region could not be determined from direct-repeat
degeneracy.


=head1 EXAMPLES

=head2 Triming off the 50bp of unconserved alignment 

CLdb_loadLeaders.pl -d ../CRISPR.sqlite test_leader_Ib.fna test_leader_Ib_aln.fna -t 50

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

