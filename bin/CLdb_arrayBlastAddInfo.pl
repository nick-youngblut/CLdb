#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddInfo.pl -- add CLdb info to spacer or DR IDs ('locusID|spacer/DR|spacer/DR_ID|groupID') in Blast file

=head1 SYNOPSIS

CLdb_arrayBlastAddInfo.pl [flags] < blast_results.txt > blast_results_info.txt

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <bool>

Add subtype? [FALSE]

=item -taxon_id  <bool>

Add taxon_id? [FALSE]

=item -taxon_name  <bool>

Add taxon_name? [FALSE]

=item -position  <bool>

Add start-stop position? [FALSE]

=item -order  <bool>

Add spacer-leader order? [FALSE]
# not implemented yet! #

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddInfo.pl

=head1 DESCRIPTION

Add CLdb info to spacer or DR IDs (locusID|spacer/DR|spacer/DR_ID|groupID).
Also, all blast hit entries for a spacer/DR group will be
duplicated for each individual spacer/DR in the spacer/DR group (replication).

=head1 EXAMPLES

=head2 Add spacer subtype and position (BLAST hits will be replicated for each individual element in the group)

CLdb_arrayBlastAddInfo.pl -d CLdb.sqlite -sub -pos < spacer-group_blast.txt > spacers_ALL_blast.txt

=head2 Add DR subtype and position (BLAST hits will be replicated for each individual element in the group)

CLdb_arrayBlastAddInfo.pl -d CLdb.sqlite -sub -pos < DR-group_blast.txt > DR_ALL_blast.txt

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
	read_fasta/;
use CLdb::blast qw/
	read_blast_file/;


#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my ($subtype_b, $taxon_id_b, $taxon_name_b, $pos_b, $order_b);
GetOptions(
	   "database=s" => \$database_file,
	   "subtype" => \$subtype_b,
	   "taxon_id" => \$taxon_id_b,
	   "taxon_name" => \$taxon_name_b,
	   "position" => \$pos_b,
	   "order" => \$order_b,
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
#table_exists($dbh, "leaders") if $leader_b;

# loading file #
my ($lines_r) = read_blast_file();

# querying CLdb #
## adding to query ##
my $query = add_query_opts($subtype_b, $taxon_id_b, $taxon_name_b);

## querying ##
my $index_r = get_CLdb_info($dbh, $lines_r, $query);

# writing edited fasta #
write_blast_file($index_r, $lines_r);


#--- disconnect ---#
$dbh->disconnect();
exit;


### Subroutines
sub write_blast_file{
	my ($index_r, $lines_r) = @_;
	
	foreach my $query (sort keys %$lines_r){
		(my $q = $query) =~ s/^# Query:\s([^|]+\|(spacer|DR)\|[^|]+\|[^|]+).*/$1/;
		foreach my $db (sort keys %{$lines_r->{$query}}){
			foreach my $blast (keys %{$lines_r->{$query}{$db}}){
				if(exists $index_r->{$q}){		# spacer group
					# duplicating blast hits #
					foreach my $new_query (@{$index_r->{$q}}){
						print $blast, "\n";
						print "# Query: $new_query\n";						
						print $db, "\n";
						print $lines_r->{$query}{$db}{$blast}{'fields'}, "\n"
							if exists $lines_r->{$query}{$db}{$blast}{'fields'};
						print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
						# printing each hit #
						foreach my $l (@{$lines_r->{$query}{$db}{$blast}{'hits'}}){
							my @l = split /\t/, $l;
							
							# changing query ID column #
							die "ERROR: cannot find 'query id' for '$l'\n"
								unless exists $lines_r->{$query}{$db}{$blast}{'fields_sep'}{'query id'};
							my $query_id_i = $lines_r->{$query}{$db}{$blast}{'fields_sep'}{'query id'};
							$l[$query_id_i] = $new_query;	
							print join("\t", @l), "\n";
							}
						}				
					}
				else{		# individual spacer
					print $blast, "\n";
					print $query, "\n";
					print $db, "\n";
					print $lines_r->{$query}{$db}{$blast}{'fields'}, "\n"
						if exists $lines_r->{$query}{$db}{$blast}{'fields'};
					print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
					# printing all hits #
					print join("\n", @{$lines_r->{$query}{$db}{$blast}{'hits'}}), "\n";
					}
				}
			}
		}
	
	
	}

sub get_CLdb_info{
# getting necessary info from CLdb #
## query is sequence-dependent ##
	my ($dbh, $lines_r, $query) = @_;
	
	my %index;
	foreach my $l (keys %$lines_r){

		next unless ($l =~ /^# Query:\s[^|]+\|(spacer|DR)\|[^|]+\|[^|]+.*/i);

		# processing line w/ name #
		(my $name = $l) =~ s/^# Query:\s([^|]+\|(spacer|DR)\|[^|]+\|[^|]+).*/$1/;
		my @name = split /\|/, $name;
		
		# sanity check #
		die "ERROR: '$name' does not have 4 components separated by '|'\n"
			unless scalar @name == 4;
		
			#print Dumper @name;
		
		# table #			
		my $tbl_oi = "Spacers";
		$tbl_oi = "DRs" if $name[1] =~ /DR/i;
	
		my $prefix = $name[1];			# "spacer|DR"
		my $cmd = "SELECT 
			$tbl_oi.Locus_ID,
			'$prefix',
			$tbl_oi.$prefix\_ID,
			$tbl_oi.$prefix\_group";
		$cmd =~ s/\t+/ /g;

		# adding other query options #
		$cmd .= ", $query" if $query;
		if($pos_b){
			$cmd .= ", loci.Scaffold";
			$cmd .= ", $tbl_oi.$prefix\_start";
			$cmd .= ", $tbl_oi.$prefix\_end";
			}	
			
		# WHERE statement #
		$cmd .= "
			FROM $tbl_oi, loci 
			WHERE loci.locus_id = $tbl_oi.locus_id";
		$cmd =~ s/\t+//g;
	
		# group or element ID #
		if( $name[2] eq "NA" ){
			$cmd .= " AND $tbl_oi.$prefix\_group == $name[3]";
			}
		else{
			$cmd .= " AND $tbl_oi.locus_ID == '$name[0]'";				# locus id
			$cmd .= " AND $tbl_oi.$prefix\_ID == $name[2]";		# elementID
			}
	
		$cmd =~ s/[\n\t]+/ /g;
		$cmd =~ s/ +/ /g;
	
		my $ret = $dbh->selectall_arrayref($cmd);
		die "ERROR: no matching entries!\n"
			unless $$ret[0];
	
		# loading hash #
		foreach my $x (@$ret){
			my $new_name = join("|", @$x);
			push @{$index{$name}}, $new_name;
			}
		}
		#print Dumper %index; exit;
	return \%index;
	}
	
sub add_query_opts{
	my ($subtype_b, $taxon_id_b, $taxon_name_b) = @_;
	
	my @query;
	push @query, "loci.subtype" if $subtype_b;
	push @query, "loci.taxon_id" if $taxon_id_b;
	push @query,"loci.taxon_name" if $taxon_name_b;
	
	return join(",", @query);
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








