#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlast2Align.pl -- align query (crDNA) vs protospacer (full length)

=head1 SYNOPSIS

CLdb_arrayBlast2Align.pl [flags] < spacer_blast_proto-DRx.txt

=head2 Optional flags

=over

=item -directory  <char>

Output directory. [./arrayBlastAlign/]

=item -query  <bool>

Parse output by blast query? [FALSE]

=item -database  <bool>

Parse output by blast database? [FALSE]

=item -v  <bool>

Verbose output. [TRUE]

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlast2Align.pl

=head1 DESCRIPTION

Make an alignment of spacers (crDNA) and protospacers (full length).

=head3 CANNOT run this script until:

The full length spacer with direct-repeat extensions must be added to
the blast table via CLdb_arrayBlastAddFullQuery.pl & CLdb_arrayBlastAddDRx.pl

The full length protospacer & extensions must be added to the blast table
via CLdb_arrayBlastAddProto.pl 

=head2 Output

If no -directory provided, output written to ./arrayBlastAlign/

The output can be parsed into separate alignment fastas using 
-query and -database.

=head1 EXAMPLES

=head2 Basic usage

$ CLdb_arrayBlast2Align.pl < spacer_blast_proto-DRx.txt

=head2 Output alignment file for each spacer query

$ CLdb_arrayBlast2Align.pl -query < spacer_blast_proto-DRx.txt

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
use File::Path qw/rmtree/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db
	get_file_path/;
use CLdb::blast qw/
	read_blast_file/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $out_dir, $by_query, $by_db);
GetOptions(
	   "directory=s" => \$out_dir,
	   "query" => \$by_query, 		# output by query
	   "database" => \$by_db, 		# output by db
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
unless( defined $out_dir ){
	my $curdir = File::Spec->rel2abs(File::Spec->curdir());
	$out_dir = join("/", $curdir, "arrayBlastAlign");
	}
mkdir $out_dir unless -d $out_dir;			# making output directory 

#--- MAIN ---#
my $blast_r = read_blast_file();

blast_2_align($blast_r, $out_dir);

#--- Subroutines ---#
sub blast_2_align{
# writing a fasta for each query #
	my ($blast_r, $out_dir) = @_;
	
		#print Dumper $blast_r; exit;
	
	open OUT, ">$out_dir/all.fna" or die $!
		if ! defined $by_query && ! defined $by_db;		# output all together	
	
	foreach my $query ( keys %$blast_r ){
		
		# checking for formating of query #
		(my $name = $query) =~ s/# Query: //i;
		die "ERROR: '$name' not formatted as 'locus_id|spacer/DR|spacerID|spacer_group'\n"
			unless $query =~ /^[^|]+\|/;
		$name =~ s/\|/-/g;
		
		# output file #
		open OUT, ">$out_dir/$name" or die $!
			if $by_query && ! defined $by_db;
		
		foreach my $db ( keys %{$blast_r->{$query}} ){
			
			my @parts = File::Spec->splitpath($db);
			
			# output file #
			open OUT, ">$out_dir/$name\__$parts[2]" or die $!
				if $by_query && $by_db;
			
			foreach my $blast ( keys %{$blast_r->{$query}{$db}} ){
				
				# checking for 'proto_seq_fullx' & 'query_seq_full_DRx' #
				die "ERROR: cannot find 'proto_seq_fullx' & 'query_seq_full_DRx' for $query -> $db -> $blast\n" 
					unless exists $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'proto_seq_fullx'}
						&& exists $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'query_seq_full_DRx'};
				
				foreach my $hit ( @{$blast_r->{$query}{$db}{$blast}{'hits'}} ){
					
					# indices of for output #
					my $query_id_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'query id'};
					my $query_seq_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'query_seq_full_DRx'};
					
					my $proto_seq_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'proto_seq_fullx'};
					my $proto_start_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'s_start_full'};
					my $proto_end_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'s_end_full'};
					my $proto_id_index = $blast_r->{$query}{$db}{$blast}{'fields_sep'}{'subject id'};
					
					my @row = split /\t/, $hit;
				
					# trimming query & proto to same length #
					my ($query_seq, $proto_seq) = trim_seqs($row[$query_seq_index], $row[$proto_seq_index]);
				
					# proto ID #
					(my $db_e = $db) =~ s/.+\///;
					my $proto_name = join("|", $db_e, $row[$proto_id_index], 
											$row[$proto_start_index],
											$row[$proto_end_index]);
					
					# writing alignment #
					print OUT join("\n", ">$row[$query_id_index]", $query_seq), "\n";
					print OUT join("\n", ">$proto_name", $proto_seq), "\n";
					}
				}
			close OUT if $by_query && $by_db;
			}
		close OUT if $by_query && ! defined $by_db;
		}
	close OUT if ! defined $by_query && ! defined $by_db;		# output all together	
	}

sub trim_seqs{
# trimming sequences so that query & proto are same length #
	my ($query_seq, $proto_seq) = @_;
	
	# determine extension lengths #
	sub ext_len{
		my $seq = shift;
		my @seq = split /[A-Z]+/, $seq;				# just lower-case extensions remaining
		map{ $seq[$_] = "" unless defined $seq[$_] } 0..1;
		my @x_len = (length $seq[0], length $seq[1] );
		return @x_len;
		}
	my @q_x_lens = ext_len($query_seq);
	my @p_x_lens = ext_len($proto_seq);
	
	# trimming extensions if needed #
	## left side ##
	if($q_x_lens[0] > $p_x_lens[0]){
		my $extra_len = $q_x_lens[0] - $p_x_lens[0];
		$query_seq = substr($query_seq, $extra_len, length($query_seq) - $extra_len);
		}
	elsif($q_x_lens[0] < $p_x_lens[0]){
		my $extra_len = $p_x_lens[0] - $q_x_lens[0];
		$proto_seq = substr($proto_seq, $extra_len, length($proto_seq) - $extra_len);
		}
	
	## right side ##
	if($q_x_lens[1] > $p_x_lens[1]){
		my $extra_len = $q_x_lens[1] - $p_x_lens[1];
		$query_seq = substr($query_seq, 0, length($query_seq) - $extra_len );
		}
	elsif($q_x_lens[1] < $p_x_lens[1]){
		my $extra_len = $p_x_lens[1] - $q_x_lens[1];
		$proto_seq = substr($proto_seq, 0, length($proto_seq) - $extra_len );
		}	
	
	return $query_seq, $proto_seq;
	}

