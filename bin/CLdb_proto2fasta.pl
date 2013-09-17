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
my (@subtype, @taxon_id, @taxon_name);		# query refinement
my (@staxon_id, @staxon_name, @sacc); 		# blast subject query refinement
my $extra_query = "";
my $len_cutoff = 1;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query, 
	   "staxon_id=s{,}" => \@staxon_id,
	   "staxon_name=s{,}" => \@staxon_name,
	   "saccession=s{,}" => \@sacc,
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
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");
$join_sql .= join_query_opts_or(\@staxon_id, \@staxon_name, \@sacc);

# getting blast hits #
my $blast_hits_r = get_blast_hits($dbh, $join_sql, $extra_query);

# binning mismatches #
write_fasta($blast_hits_r);


### Subroutines
sub write_fasta{
	my ($blast_hits_r) = @_;
		
	# body #
	foreach my $entry (sort{$a->[0]<=>$b->[0]} @$blast_hits_r){
		# applying len cutoff #
		next if abs($$entry[15] - $$entry[14] + 1) / $$entry[16] < $len_cutoff;
		
		# checking for undefined values #
		map{$_ = "" unless $_} @$entry[(2..3,5..7,21)];
		
		# trimming frag to just proto #
		# frag = $$entry[17];
		my $xstart = $$entry[19];
		my $xend = $$entry[20];
		my $sstart = $$entry[9];
		my $send = $$entry[10];
		
		## flipping start-end if necessary ##
		($xstart, $xend) = flip_se($xstart, $xend) if $xstart > $xend;
		($sstart, $send) = flip_se($sstart, $send) if $sstart > $send;
		
		## getting substr (protospacer) ##
		my $proto;
		if($$entry[9] <= $$entry[10]){
			$proto = substr($$entry[17], $sstart - $xstart - 1, $send - $sstart + 1);
			}
		else{
			$proto = substr($$entry[17], $sstart - $xstart, $send - $sstart + 1);
			}
		
			#print Dumper$$entry[18], $proto;		# qseq, frag(trimmed)
		
		# writing fasta #
		print join("\n",
			join("__", ">spacer" ,@$entry[(0,2,3,21,1,4..7)] ),
			$$entry[18]), "\n";
		print join("\n", 
			join("__",  ">proto", @$entry[(0,2,3,21,1,4..7)] ),
			$proto), "\n";
		}
	}


sub flip_se{
	my ($start, $end) = @_;
	return $end, $start;
	}

sub get_blast_hits{
# 3 table join #
	my ($dbh, $join_sqls, $extra_query) = @_;
	
	# basic query #
	my $query = "SELECT 
c.group_id, 
a.locus_id,
a.taxon_name,
a.taxon_id,
b.spacer_id,
c.S_Taxon_ID, 
c.S_Taxon_name,
c.S_accession,
c.sseqid,
c.sstart, 
c.send, 
c.evalue,
c.mismatch,
c.pident,
c.qstart,
c.qend,
c.qlen,
c.frag,
c.qseq,
c.xstart,
c.xend,
a.subtype
FROM Loci a, Spacers b, Blast_hits c
WHERE a.locus_id == b.locus_id
AND c.group_id == b.spacer_group
AND c.array_hit == 'no'
AND c.spacer_DR == 'Spacer'
$join_sqls";
	$query =~ s/[\n\r]+/ /g;
	
	$query = join(" ", $query, $extra_query);
	
	# status #
	print STDERR "$query\n" if $verbose;

	#print Dumper $query; exit;

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
	map{ s/"*(.+)"*/"$1"/ } @$staxon_id_r;
	map{ s/"*(.+)"*/"$1"/ } @$staxon_name_r;	
	map{ s/"*(.+)"*/"$1"/ } @$sacc_r;	
	
	if(@$staxon_id_r && @$staxon_name_r){
		return join("", " AND (c.S_taxon_id IN (", join(", ", @$staxon_id_r),
						") OR c.S_taxon_name IN (", join(", ", @$staxon_name_r),
						"))");
		}
	elsif(@$staxon_id_r){
		return join("", " AND c.S_taxon_id IN (", join(", ", @$staxon_id_r), ")");
		}
	elsif(@$staxon_name_r){
		return join("", " AND c.S_taxon_name IN (", join(", ", @$staxon_name_r), ")");	
		}
	elsif(@$sacc_r){
		return join("", " AND c.S_accession IN (", join(", ", @$sacc_r), ")");	
		}
	}

sub join_query_opts{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND a.$cat IN (", join(", ", @$vals_r), ")");
	}




__END__

=pod

=head1 NAME

CLdb_proto2fasta.pl -- Make fasta of protospacers

=head1 SYNOPSIS

CLdb_proto2fasta.pl [flags] > protospacers.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 options

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -staxon_id  <char>

Refine query to specific a subject taxon_id(s) (>1 argument allowed).

=item -staxon_name  <char>

Refine query to specific a subject taxon_name(s) (>1 argument allowed).

=item -saccession  <char>

Refine query to specific a accession(s) (>1 argument allowed).

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

perldoc CLdb_proto2fasta.pl

=head1 DESCRIPTION

Make a fasta of all protospacers.

=head2 Sequence naming (separated by "__"):

=over

=item group_id

=item taxon_name

=item taxon_id

=item subtype

=item locus_id

=item spacer_id

=item subject_taxon_id

=item subject_taxon_name

=item subject_accession

=back

=head2 WARNING!

Spacer blasting and DR filtering must be done prior!

=head1 EXAMPLES

=head2 Fasta of all protospacers

CLdb_proto2fasta.pl -d CLdb.sqlite  > all_proto_mismatch.txt

=head2 Fasta of all protospacers in subtype I-A

CLdb_proto2fasta.pl -d CLdb.sqlite -sub I-A > I-A_proto_mismatch.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

