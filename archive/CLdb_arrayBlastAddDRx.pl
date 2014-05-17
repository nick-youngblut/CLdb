#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddDRx.pl -- adding direct repeat extension to spacer queries in blast table

=head1 SYNOPSIS

CLdb_arrayBlastAddDRx.pl [flags] < blast_results.txt > blast_results_DRx.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -x 	length of direct repeat extension to spacer queries (<= 1 is fraction of DR length). [0.5]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddDRx.pl

=head1 DESCRIPTION

Use to get the full crDNA sequence.
Three new fields will be added: 'query_seq_full_DRx', 'DRx_5p', & 'DRx_3p'

=head2 WARNING

All spacer groups used for querying must be replicated
to the individual spacers (use CLdb_arrayBlastAddInfo.pl),
which is needed to get the proper DR sequences adjacent
to each spacer. 

ALSO, the full length query sequence ('query_seq_full' field)
must be added to the blast table prior to
running this script (use CLdb_BlastAddFullQuery.pl)! 

=head1 EXAMPLES

=head2 Basic usage:

CLdb_arrayBlastAddDRx.pl -d CLdb.sqlite < spacer_blast.txt > spacer_blast_DRx.txt

=head2 DRx of 10 nucleotides

CLdb_arrayBlastAddDRx.pl -d CLdb.sqlite -x 10 < spacer_blast.txt > spacer_blast_DRx.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
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
use CLdb::query qw/
	table_exists
	n_entries/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
use CLdb::seq qw/
	revcomp/;
use CLdb::blast qw/
	read_blast_file
	write_blast_file/;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my $DR_x = 0.5;
GetOptions(
	   "database=s" => \$database_file,
	   "x=f" => \$DR_x,						# add DR extension to spacers?
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# checking table existence #
table_exists($dbh, "loci"); 
table_exists($dbh, "spacers");
table_exists($dbh, "DRs");

# loading file #
my ($lines_r) = read_blast_file();

# querying CLdb to get DR extensions #
my $strand_r = get_locus_strand($dbh, $lines_r);
add_DR_x($dbh, $lines_r, $strand_r, $DR_x);

# writing edited fasta #
write_blast_file($lines_r);


#--- disconnect ---#
$dbh->disconnect();
exit;


### Subroutines
sub write_blast_file_OLD2{
	my ($lines_r) = @_;
	
	foreach my $query ( sort keys %$lines_r ){
		foreach my $db ( sort keys %{$lines_r->{$query}} ){
			foreach my $blast ( keys %{$lines_r->{$query}{$db}} ){
				print $blast, "\n";
				print $query, "\n";
				print $db, "\n";
				print $lines_r->{$query}{$db}{$blast}{'fields'}, "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'fields'};
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
				# printing all hits #
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'hits'}}), "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'hits'};
				}
			}
		}
	}

sub write_blast_file_OLD{
	my ($index_r, $lines_r) = @_;
	
	foreach my $query (sort keys %$lines_r){
		(my $q = $query) =~ s/^# Query:\s([^|]+\|(spacer|DR)\|[^|]+\|[^|]+).*/$1/;
		foreach my $db (sort keys %{$lines_r->{$query}}){
			foreach my $blast (keys %{$lines_r->{$query}{$db}}){
				if(exists $index_r->{$q}){
					# duplicating blast hits #
					foreach my $new_query (@{$index_r->{$q}}){
						print $blast, "\n";
						print "# Query: $new_query\n";						
						print $db, "\n";
						print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
						foreach my $l (@{$lines_r->{$query}{$db}{$blast}{'hits'}}){
							my @l = split /\t/, $l;
							print join("\t", $new_query, @l[1..$#l]), "\n";
							}
						}				
					}
				else{
					print $blast, "\n";
					print $query, "\n";
					print $db, "\n";
					print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
					print join("\n", @{$lines_r->{$query}{$db}{$blast}{'hits'}}), "\n";
					}
				}
			}
		}
	
	
	}

sub add_DR_x{
	my ($dbh, $lines_r, $strand_r, $DR_x) = @_;
	
	# query for spacer start-end #
	my $cmd_spacer = "SELECT spacer_start, spacer_end FROM spacers WHERE locus_id = ? AND spacer_id = ?";
	my $sth_spacer = $dbh->prepare($cmd_spacer);
	
	# queries to get DR sequence for each spacer #
	## allowing for spacer-DR position overlap or adjacency ##
	my $cmd1 = "SELECT DR_sequence from DRs where locus_ID = ? 
					AND DR_end <= ?
					AND DR_end >= ? - 1";
	my $sth1 = $dbh->prepare($cmd1);

	my $cmd2 = "SELECT DR_sequence from DRs where locus_ID = ? 
					AND DR_start <= ? + 1
					AND DR_start >= ?";
	my $sth2 = $dbh->prepare($cmd2);


	foreach my $query (sort keys %$lines_r){
		
		# query #
		(my $name = $query) =~ s/# Query: //i;
		my @n = split /\|/, $name;			# locus_id = $l[0]
		die "ERROR: '$name' is not a spacer hit (should have 'spacer' in 2nd field)\n"
			unless $n[1] =~ /^spacer$/i;

		# getting spacer start-end (+ strand in spacers table) #
		$sth_spacer->bind_param(1, $n[0]);		# locus_id
		$sth_spacer->bind_param(2, $n[2]);		# spacer_id

		$sth_spacer->execute();
		my $ret = $sth_spacer->fetchrow_arrayref();
		die "ERROR: no spacer matches for locus_id->$n[0], spacer_id->$n[2]!\n"
			unless defined $ret;

		# getting 5' DR extension (+ strand orientation) #
		$sth1->bind_param(1, $n[0]);
		$sth1->bind_param(2, $$ret[0]);
		$sth1->bind_param(3, $$ret[0]);
		$sth1->execute();
		my $ret1 = $sth1->fetchrow_arrayref();
		die "ERROR: no DR 5' matches for locus_id->$n[0], spacer_start->$$ret[0]!\n"
			unless defined $ret1;
		my $DR_seq_5p = $$ret1[0];

		# getting 3' DR extension (+ strand orientation) #
		$sth2->bind_param(1, $n[0]);
		$sth2->bind_param(2, $$ret[1]);
		$sth2->bind_param(3, $$ret[1]);
		$sth2->execute();
		my $ret2 = $sth2->fetchrow_arrayref();
		die "ERROR: no DR 3' matches for locus_id->$n[0], spacer_start->$$ret[1]!\n"
			unless defined $ret2;
		my $DR_seq_3p = $$ret2[0];

		# lower case #
		map{ tr/A-Z/a-z/ } ($DR_seq_5p, $DR_seq_3p);

		# trimming DR #
		$DR_seq_5p = trim_DR($DR_seq_5p, $DR_x, '5p');
		$DR_seq_3p = trim_DR($DR_seq_3p, $DR_x, '3p');
		


		# inverting if spacer is - strand #
		## not needed! array file will begin inverted relative to others if on other strand in genome #
		
		#die "ERROR: cannot find strand for $n[0]\n"
		#	unless exists $strand_r->{$n[0]};
		#if ($strand_r->{$n[0]} eq '-'){
		#	$DR_seq_5p = revcomp $DR_seq_5p;
		#	$DR_seq_3p = revcomp $DR_seq_3p;
		#	($DR_seq_5p, $DR_seq_3p) = ($DR_seq_3p, $DR_seq_5p);
		#	}
			

		# adding to 'query_seq_full' #
		foreach my $db (sort keys %{$lines_r->{$query}}){
			foreach my $blast (keys %{$lines_r->{$query}{$db}}){
				# checking for 'query_seq_full' #
				unless(exists $lines_r->{$query}{$db}{$blast}{'fields_sep'}{'query_seq_full'}){
					warn "WARNING: 'query_seq_full' not found for query->$query, db->$db, blast->$blast. Skipping!\n";
					next;
					}
				my $query_seq_full_i = 	$lines_r->{$query}{$db}{$blast}{'fields_sep'}{'query_seq_full'};
					
				# adding 'query_seq_full_DRx' to fields #
				$lines_r->{$query}{$db}{$blast}{'fields_sep'}{'query_seq_full_DRx'} = 
					scalar keys %{$lines_r->{$query}{$db}{$blast}{'fields_sep'}};
				$lines_r->{$query}{$db}{$blast}{'fields'} .= ', query_seq_full_DRx';
			
					
				# adding 'DRx_5p' to fields #
				$lines_r->{$query}{$db}{$blast}{'fields_sep'}{'DRx_5p'} = 
					scalar keys %{$lines_r->{$query}{$db}{$blast}{'fields_sep'}};
				$lines_r->{$query}{$db}{$blast}{'fields'} .= ', DRx_5p';

				# adding 'DRx_3p' to fields #
				$lines_r->{$query}{$db}{$blast}{'fields_sep'}{'DRx_3p'} = 
					scalar keys %{$lines_r->{$query}{$db}{$blast}{'fields_sep'}};
				$lines_r->{$query}{$db}{$blast}{'fields'} .= ', DRx_3p';

				# adding values to each hit #
				foreach my $hit ( @{$lines_r->{$query}{$db}{$blast}{'hits'}} ){
					my @h = split /\t/, $hit;
					push @h, join("", $DR_seq_5p, $h[$query_seq_full_i], $DR_seq_3p),
								$DR_seq_5p, $DR_seq_3p;
					$hit = join("\t", @h);
					}
				}
			}
		}

		#print Dumper $lines_r; exit;
	}
	
sub trim_DR{
# trimming DR to just amount of extension needed #
## can trim by fraction or N-nucleotides ##
	my ($DR_seq, $DR_x, $cat) = @_;
	
	my $len = length($DR_seq);
	my $ext_len;
	if($DR_x <= 1){						# assumed to take a fraction
		$ext_len = int($len * $DR_x);		# fraction of total length
		}
	else{
		$ext_len = int($DR_x);				# N-nucleotides 
		}
	
	if($cat eq '5p'){
		return substr($DR_seq, $len - $ext_len, $len); 				# 3' end 
		}
	elsif($cat eq '3p'){
		return substr($DR_seq, 0, $ext_len); 				# 3' end 
		}
	else{
		die "LOGIC ERROR: $!\n";
		}
	}
	
sub get_locus_strand{
# getting spacer strand according to how array_start array_end in loci table is coded #
	my ($dbh, $lines_r) = @_;
	
	my $cmd = "SELECT array_start, array_end from loci where locus_id=?";
	my $sth = $dbh->prepare($cmd);
	
	my %locus_strand;
	foreach my $query ( keys %$lines_r ){
		# checking for formating of query #
		(my $name = $query) =~ s/# Query: //i;
		die "ERROR: '$name' not formatted as 'locus_id|spacer/DR|spacerID|spacer_group'\n"
			unless $query =~ /^[^|]+\|/;
		my @l = split /\|/, $name;
		
		# querying CLdb #
		$sth->bind_param(1, $l[0]);
		$sth->execute();
		my $ret = $sth->fetchrow_arrayref();
		
		if(defined $ret){
			# determining strand #
			if($$ret[0] <= $$ret[1]){
				$locus_strand{$l[0]} = '+';
				}
			else{
				$locus_strand{$l[0]} = '-';
				}
			}
		else{
			die "ERROR: not matching entries for locusID: '$l[0]'\n";
			}
		}		
		
		#print Dumper %locus_strand; exit;
	return \%locus_strand;
	}

sub read_blast_file_OLD{
# reading in each blast entry & extracting names and line numbers #
	my %lines;
	my $blast;
	my $query;
	my $db;
	while(<>){
		chomp;
		
		if(/^# BLAST/i){
			$blast = $_;
			}
		elsif(/^# Query/i){
			$query = $_;
			}
		elsif(/^# Database/i){
			$db = $_;
			}	
		elsif(/^# /){
			push @{$lines{$query}{$db}{$blast}{'comments'}}, $_;
			}
		else{
			push @{$lines{$query}{$db}{$blast}{'hits'}}, $_;
			}
			
		}		
		#print Dumper %lines; exit;	
	return \%lines;
	}








