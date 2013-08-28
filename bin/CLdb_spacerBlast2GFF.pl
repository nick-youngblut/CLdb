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
my (@subtype, @taxon_id, @taxon_name);
my (@staxon_id, @staxon_name); 				# blast subject query refinement
my $extra_query = "";
my $hit_len_cutoff = 0.9;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "length=f" => \$hit_len_cutoff,			# blast hit must be full length of query
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
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");
my $join_sqls .= join_query_opts_or(\@staxon_id, \@staxon_name);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sqls, $extra_query);

# getting spacer group length #
my $spacer_len_r = get_spacer_group_length($dbh, $join_sql, $extra_query);

# filter spacer hits to just full length hits if needed #
$blast_hits_r = apply_hit_len_cutoff($blast_hits_r, $spacer_len_r, $hit_len_cutoff) 
	unless $hit_len_cutoff == 1;

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
		# getting strand #
		my $strand;
		if($$hit[4] <= $$hit[5]){ $strand = "+"; }
		else{ $strand = "-"; }
		
		# writing table #
		print join("\t",
			$$hit[10],					# subject scaffold
			"CLdb_spacer_blast",		# source
			"region",
			$$hit[4],			# hit start
			$$hit[5],			# hit end
			$$hit[6],			# evalue
			$strand,
			".",				# phase
			join(";", 
				"ID=\"$$hit[1]\"",		# query taxon_name
				"Name=\"$$hit[2]\"",	# query taxon_id
				"Alias=$$hit[9]",		# locus_id
				"Note='mismatch:$$hit[7]  pident:$$hit[8]  spacer_group:$$hit[0]'"
				)
			), "\n";
		}
	}

sub apply_hit_len_cutoff{
# hit must be so much of total query length #
	my ($blast_hits_r, $spacer_len_r, $hit_len_cutoff) = @_;

	my @filt_hits;
	foreach my $hit (sort @$blast_hits_r){	
		# start - end #
		my $hit_start = $$hit[4];
		my $hit_end = $$hit[5];
		
		# length cutoff #
		die " ERROR: $$hit[0] not found in spacer length index!\n Did you mean to use '-staxon' instead of '-taxon'?"
			unless exists $spacer_len_r->{$$hit[0]};
		push (@filt_hits, $hit) if abs($hit_end - $hit_start + 1)  /  $spacer_len_r->{$$hit[0]}
			>= $hit_len_cutoff;
			
		}
	
		#print Dumper @filt_hits; exit;
	return \@filt_hits;
	}

sub get_spacer_group_length{
# getting spacer group length for removing spacers that do not meet length cutoff #
	my ($dbh, $join_sql, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT spacers.spacer_group, spacers.spacer_start, spacers.spacer_end
from Loci, Spacers
WHERE Loci.locus_id == Spacers.locus_id
$join_sql";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	
	# making hash (spacer_group => spacer_length) #
	my %spacer_length;
	foreach my $row (@$ret){
		$spacer_length{$$row[0]} = abs($$row[2] - $$row[1]);
		}
	
	return \%spacer_length;
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sqls, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
blast_hits.Spacer_group, 
loci.Taxon_ID, 
loci.Taxon_name, 
blast_hits.Subject, 
blast_hits.sstart, 
blast_hits.send, 
blast_hits.evalue,
blast_hits.mismatch,
blast_hits.pident,
loci.locus_id,
blast_hits.subject
FROM Loci, Spacers, Blast_hits
WHERE Loci.locus_id == Spacers.locus_id
AND blast_hits.spacer_group == spacers.spacer_group
AND blast_hits.CRISPR_array == 'no'
$join_sqls";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	return $ret;
	}

sub join_query_opts_or{
	my ($staxon_id_r, $staxon_name_r) = @_;
	
	return "" unless @$staxon_id_r || @$staxon_name_r;
	
	# adding quotes #
	map{ s/"*(.+)"*/'$1'/ } @$staxon_id_r;
	map{ s/"*(.+)"*/'$1'/ } @$staxon_name_r;	
	
	if(@$staxon_id_r && @$staxon_name_r){
		return join("", " AND (blast_hits.taxon_id IN (", join(", ", @$staxon_id_r),
						") OR blast_hits.taxon_name IN (", join(", ", @$staxon_name_r),
						"))");
		}
	elsif(@$staxon_id_r){
		return join("", " AND blast_hits.taxon_id IN (", join(", ", @$staxon_id_r), ")");
		}
	elsif(@$staxon_name_r){
		return join("", " AND blast_hits.taxon_name IN (", join(", ", @$staxon_name_r), ")");	
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

=item -d 	CLdb database.

=back

=head2 options

=over

=item -subtype

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query

Extra sql to refine the query.

=item -staxon_id

BLAST subject taxon_id (>1 argument allowed)

=item -staxon_name

BLAST subject taxon_name (>1 argument allowed)

=item -length

Length cutoff for blast hit (>=; fraction of spacer length). [0.9]

=item -v	Verbose output. [FALSE]

=item -h	This help message

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

