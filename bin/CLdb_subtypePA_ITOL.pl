#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_subtypePA_ITOL.pl -- make CRISPR subtype pres-abs ITOL table

=head1 SYNOPSIS

CLdb_subtypePA_ITOL.pl [flags] > subtype_PA.meta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -colors  <char>

For providing user-defined hexidecimal colors (>= 1 argument).

=item -abundance  <bool>

Provide counts of subtypes per taxon instead of binary pres-abs? [FALSE]

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -query  <char>

Extra sql to refine which sequences are returned.

=item -group  <bool>

Get array elements de-replicated by group (ie. all uniqe sequences). [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_subtypePA_ITOL.pl

=head1 DESCRIPTION

Make a CRISPR subtype presence-absence table for plotting
with a tree in ITOL.

Subtype coloring is done in alpha-numeric order unless
the subtypes are provided via '-subtype'

=head1 EXAMPLES

=head2 Pres-abs table of all subtypes & taxa in CLdb

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite 

=head2 Subtype count table of all subtypes & taxa in CLdb (not-binary)

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite -a

=head2 Subtype count table of loci containing CAS operons

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite -q "AND operon_status!='absent'"

=head2 Subtype count table of loci containing intact CAS operons (no broken or shuffled)

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite -q "AND operon_status='intact'"

=head2 User provided colors

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite -c "#FF0000 #FF6600 #FFFF00"

=head2 Pres-abs table of specific subtypes & taxa_names

CLdb_subtypePA_ITOL.pl -d CLdb.sqlite -sub I-A I-B I-C -taxon_name Ecoli Senterica Awoodii 

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
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file, $binary_output);
my (@subtype, @taxon_id, @taxon_name);
my @colors;
my $extra_query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$extra_query,
	   "colors=s{,}" => \@colors,				# colors for plotting 
	   "abundance" => \$binary_output, 			# binary output? [TRUE]
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# querying subtypes & counting per taxon #
my ($taxa_r, $subtypes_r) = count_subtypes($dbh, $extra_query, $join_sql);

# getting colors for subtypes #
get_subtype_colors($subtypes_r, \@colors, \@subtype);

# writing ITOL metadata table #
write_ITOL_metadata($taxa_r, $subtypes_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub write_ITOL_metadata{
# writing metadata table for ITOL #
	my ($taxa_r, $subtypes_r) = @_;
	
	# conserving order #
	my @subtypes = sort keys %$subtypes_r;
	
	# header #
	## labels ##
	print join("\t", "LABELS", @subtypes), "\n";
	## colors ##
	my @colors;
	map{push @colors, $subtypes_r->{$_}} @subtypes;
	print join("\t", "COLORS", @colors), "\n";
	
	# body #
	foreach my $taxon (keys %$taxa_r){
		my @line;
		foreach my $subtype (@subtypes){
			if(exists $taxa_r->{$taxon}{$subtype}){
				if(! $binary_output){ push @line, 1;}				# by default: binary
				else{ push @line, $taxa_r->{$taxon}{$subtype}; }
				}
			else{
				push @line, 0;
				}
			}
		print join("\t", $taxon, @line), "\n";
		}
	
	}

sub get_subtype_colors{
	my ($subtypes_r, $colors_r, $subtype_r) = @_;

	# getting rainbow #
	unless(@$colors_r){		# if the user did not provide colors
		@$colors_r = qw/FF0000 FF6600 FFFF00 33CC00 00CCFF 0033FF 9900FF FF00FF 660000 CC3300 336600 006666 000066/;
		map{$_ =~ s/^/#/} @$colors_r;
		}

	# subtype order alpha-numeric or user-defined #
	my @subtype_order;
	if(@$subtype_r){
		 map{$_ =~ s/"//g} @$subtype_r; 
		 @subtype_order = @$subtype_r;  
		 }
	else{ 
		@subtype_order = sort keys %$subtypes_r; 
		}

	# connecting colors to subtypes #
	my $cnt = 0;
	foreach my $subtype (@subtype_order){
		unless($$colors_r[$cnt]){
			die " ERROR: not enough colors for all subtypes\n Number of subtypes: ", scalar keys %$subtypes_r, "\n";
			#die;
			}
		if($subtype eq "NA"){
			$subtypes_r->{$subtype} = "#000000";			# NA = black
			}
		else{
			$subtypes_r->{$subtype} = $$colors_r[$cnt];
			}
		$cnt++;
		}

		#print Dumper %$subtypes_r; exit;
	}

sub count_subtypes{
	my ($dbh, $extra_query, $join_sql) = @_;
	
	# make query #
	my $query = "SELECT taxon_name, subtype, count(subtype) FROM loci";
	
	if($join_sql || $extra_query){
		$query = join(" ", $query, "WHERE locus_id = locus_id");
		}
	$query = join(" ", $query, $join_sql) if $join_sql;
	$query = join(" ", $query, $extra_query) if $extra_query;
	$query = join(" ", $query, "group by taxon_name, subtype");
	
	# status #
	print STDERR "$query\n" if $verbose;

	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	die " ERROR: no matching entries!\n"
		unless $$ret[0];
	
	# hash by taxon=>subtype #
	my %taxa;
	my %subtypes;
	foreach my $row (@$ret){
		$$row[1] = "NA" unless $$row[1];
		$taxa{$$row[0]}{$$row[1]} = $$row[2];	# taxon_name=>subtype=>count
		$subtypes{$$row[1]} = 1;				# unique subtypes
		}
		
		#print Dumper %taxa; exit;
		#print Dumper %subtypes; exit;
	return \%taxa, \%subtypes;
	}
	
sub join_query_opts_OLD{
# joining query options for selecting loci #
	my ($vals_r, $cat) = @_;

	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND $cat IN (", join(", ", @$vals_r), ")");
	}



