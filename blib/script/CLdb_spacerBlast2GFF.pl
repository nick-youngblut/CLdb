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

my ($verbose, $database_file);
my (@staxon_id, @staxon_name, @sacc); 				# blast subject query refinement
my $extra_query = "";
my $len_cutoff = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "sacc=s{,}" => \@sacc,
	   "length=f" => \$len_cutoff,			# blast hit must be full length of query
	   "query=s" => \$extra_query,
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

# joining query options (for table join) #
my $join_sqls = join_query_opts_or(\@staxon_id, \@staxon_name, \@sacc);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sqls, $extra_query);

# filter spacer hits to just full length hits if needed #
#$blast_hits_r = apply_hit_len_cutoff($blast_hits_r, $spacer_len_r, $hit_len_cutoff) 
#	unless $hit_len_cutoff == 1;

# convert blast data to GFF3 format #
blast2gff3($blast_hits_r);


### Subroutines
sub blast2gff3{
# blast results to gff3 # 

# blast table 2 GFF (index) #
# seqid = spacer_group_ID
# source = 'CLdb_spacer_blast'
# feature = 'region'
# start = 'hit start'
# end = 'hit end'
# score = e-value
# strand = strand
# phase = 0
# attributes
	# ID = taxon_ID
	# name = taxon_name
	# alias = locus_id
	# note = N-mismatches

	my ($blast_hits_r) = @_;
	
	#print Dumper $blast_hits_r; exit;
	foreach my $hit (@$blast_hits_r){
		# applying length cutoff #
		next if abs($$hit[15] - $$hit[14] + 1) / $$hit[16] < $len_cutoff;
		
		# getting strand #
		my $strand;
		
		if($$hit[9] <= $$hit[10]){ $strand = "+"; }
		else{ $strand = "-"; }
		
		# query_id #
		map{$_ = "" unless $_} @$hit[1..4];
		(my $locus_sid = join("__", $$hit[1], $$hit[4])) =~ s/^/cli./;
		
		# writing table #
		print join("\t",
			$$hit[8],					# subject scaffold
			"CLdb_spacer_blast",		# source
			"region",
			$$hit[9],			# hit start (sstart)
			$$hit[10],			# hit end	(send)
			$$hit[11],			# evalue
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$hit[3]\"",			# query taxon_name
				"Name=\"$$hit[2]\"",		# query taxon_id
				"Alias=\"$locus_sid\"",			# locus_id
				"Note='mismatch:$$hit[12]  pident:$$hit[13]  spacer_group:$$hit[0]'"
				)
			), "\n";
		}
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sqls, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
blast_hits.group_id, 
loci.locus_id,
loci.taxon_name,
loci.taxon_id,
spacers.spacer_id,
blast_hits.S_Taxon_ID, 
blast_hits.S_Taxon_name,
blast_hits.S_accession,
blast_hits.sseqid,
blast_hits.sstart, 
blast_hits.send, 
blast_hits.evalue,
blast_hits.mismatch,
blast_hits.pident,
blast_hits.qstart,
blast_hits.qend,
blast_hits.qlen
FROM Loci, Spacers, Blast_hits
WHERE Loci.locus_id == Spacers.locus_id
AND blast_hits.group_id == spacers.spacer_group
AND blast_hits.array_hit == 'no'
AND blast_hits.spacer_DR == 'Spacer'
$join_sqls";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
		#print Dumper @$ret; exit;
	return $ret;
	}

sub join_query_opts_or{
	my ($staxon_id_r, $staxon_name_r, $sacc_r) = @_;
	
	return "" unless @$staxon_id_r || @$staxon_name_r || @$sacc_r;
	
	# adding quotes #
	map{ s/"*(.+)"*/'$1'/ } @$staxon_id_r;
	map{ s/"*(.+)"*/'$1'/ } @$staxon_name_r;	
	map{ s/"*(.+)"*/'$1'/ } @$sacc_r;	
	
	if(@$staxon_id_r && @$staxon_name_r){
		return join("", " AND (blast_hits.s_taxon_id IN (", join(", ", @$staxon_id_r),
						") OR blast_hits.s_taxon_name IN (", join(", ", @$staxon_name_r),
						"))");
		}
	elsif(@$staxon_id_r){
		return join("", " AND blast_hits.s_taxon_id IN (", join(", ", @$staxon_id_r), ")");
		}
	elsif(@$staxon_name_r){
		return join("", " AND blast_hits.s_taxon_name IN (", join(", ", @$staxon_name_r), ")");	
		}
	elsif(@$sacc_r){
		return join("", " AND blast_hits.s_accession IN (", join(", ", @$sacc_r), ")");		
		}
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}




__END__

=pod

=head1 NAME

CLdb_spacerBlast2GFF.pl -- make GFF3 file from spacer blasts

=head1 SYNOPSIS

CLdb_spacerBlast2GFF.pl [flags] > spacer_hits.gff

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 options

=over

=item -staxon_id  <char>

Refine query to specific a subject taxon_id(s) (>1 argument allowed).

=item -staxon_name  <char>

Refine query to specific a subject taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine the query.

=item -length  <float>

Length cutoff for blast hit (>=; fraction of spacer length). [1]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_spacerBlast2GFF.pl

=head1 DESCRIPTION

Convert spacer blasts in CLdb (see CLdb_spacerBlast.pl)
to GFF3 file format for viewing spacer blast hits
on a subject genome. 

=head2 Output columns

=over

=item seqid = scaffold (subject_genome)

=item source = 'CLdb_spacer_blast'

=item feature = 'region'

=item start = 'hit start'

=item end = 'hit end'

=item score = e-value

=item strand = strand

=item phase = 0

=back 

=head3 'Attributes' column

=over

=item ID = taxon_ID (spacer)

=item name = taxon_name (spacer)

=item alias = locus_id (spacer)

=item note = mismatches; percentID; spacer_group

=back

=head2 WARNING!

Spacer blasting must be done prior!

'-staxon_name' or '-staxon_id' should be used
to limit hits to just 1 subject genome.

=head1 EXAMPLES

=head2 GFF3 of all spacers that hit E.coli

CLdb_spacerBlast2GFF.pl -d CLdb.sqlite -staxon_name "E.coli" > ecoli_hits.gff

=head2 GFF3 of all spacers that hit FIG 2209.27

CLdb_spacerBlast2GFF.pl -d CLdb.sqlite -staxon_id 2209.27 > 2209.27_hits.gff

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

