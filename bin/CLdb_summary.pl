#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;

# number of arrays, operons
	# by status (+ all)
		# broken down by:
			# loci
			# subtype
			# taxon_name/id
	# total
# number of spacers, DR, genes, leaderseqs, consensus seqs
	# by cutoff (when possible)
		# broken down by:
			# loci
			# subtype
			# taxon_name/id
	# total
# spacer blast hits
	# by (array_hit/not)
		# broken down by:
			# loci
			# subtype
			# taxon_name/id
	# total

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $database_file);
my ($subtype, $taxon_id, $taxon_name, $by_total);
my $extra_query = "";
my @cutoff;
GetOptions(
	   "database=s" => \$database_file,
	   "subtype" => \$subtype,
	   "id" => \$taxon_id,
	   "name" => \$taxon_name,
	   "cutoff=f{,}" => \@cutoff, 				 #clustering cutoff(s)
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
@cutoff = (1) unless @cutoff;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# checking for tables of interest #
my $table_list_r = list_tables($dbh);
$table_list_r = entry_count($dbh, $table_list_r);

# make 'group by' sql #
my $group_by_r = make_group_by_sql($subtype, $taxon_id, $taxon_name);

# summary #
## loci ##
if ( exists $table_list_r->{"loci"} ){
	sum_loci($dbh, $group_by_r);
	}	
else{ print STDERR "...no loci table found, skipping loci summary!\n"; }


## spacers ##
if ( exists $table_list_r->{"loci"} && exists $table_list_r->{"spacers"} ){
	if( exists $table_list_r->{"spacer_hclust"} && $table_list_r->{"spacer_hclust"} > 0){
		sum_spacers_DR_hclust($dbh, "spacer", $group_by_r, \@cutoff);
		}
	else{
		sum_spacers_DR($dbh, "spacers", $group_by_r);
		}
	}	
else{ print STDERR "...no loci and/or spacer tables found, skipping spacer summary!\n"; }


## DR ##
if ( exists $table_list_r->{"loci"} && exists $table_list_r->{"directrepeats"} ){ 
	if( exists $table_list_r->{"directrepeat_hclust"} && $table_list_r->{"directrepeat_hclust"} > 0){
		sum_spacers_DR_hclust($dbh, "DirectRepeat", $group_by_r, \@cutoff);
		}
	else{
		sum_spacers_DR($dbh, "DirectRepeat", $group_by_r);
		}
	}	
else{ print STDERR "...no loci and/or DirectRepeat tables found, skipping DirectRepeat summary!\n"; }


## genes ##
if ( exists $table_list_r->{"loci"} && exists $table_list_r->{"genes"} ) {
	sum_genes($dbh, $group_by_r);
	}	
else{ print STDERR "...no loci and/or genes tables found, skipping gene summary!\n"; }


## leaderseqs ##
if ( exists $table_list_r->{"loci"} && exists $table_list_r->{"leaderseqs"} ) {
	sum_leaderseqs($dbh, $group_by_r);
	}	
else{ print STDERR "...no loci and/or leaderseqs tables found, skipping leader sequence summary!\n"; }



# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub sum_leaderseqs{
	my ($dbh, $group_by_r) = @_;

	# by operon/array status #
	my @select;
	push @select, @$group_by_r if @$group_by_r;
	push @select, qw/b.Leader_group/;
	my $select = join(",", @select);

	my $q = "SELECT $select,'NA',count(*) FROM Loci a, Leaderseqs b WHERE a.locus_id=b.locus_id GROUP BY $select";
	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "leaderseqs", @$_), "\n";
		}	

	# total #
	if(@$group_by_r){
		my $group_by = join(",", @$group_by_r);
		$q = "SELECT $group_by,count(*) FROM Loci a, Leaderseqs b WHERE a.locus_id=b.locus_id GROUP BY $group_by";
		}
	else{ $q = "SELECT count(*) FROM Loci a, Genes b WHERE a.locus_id=b.locus_id"; }

	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "leaderseqs", @$_[0..($#$_-1)], qw/Total NA/, $$_[$#$_]), "\n";
		}	
	}

sub sum_genes{
	my ($dbh, $group_by_r) = @_;

	# by operon/array status #
	my @select;
	push @select, @$group_by_r if @$group_by_r;
	push @select, qw/b.In_operon/;
	my $select = join(",", @select);

	my $q = "SELECT $select,'NA',count(*) FROM Loci a, Genes b WHERE a.locus_id=b.locus_id GROUP BY $select";
	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "genes", @$_), "\n";
		}	

	# total #
	if(@$group_by_r){
		my $group_by = join(",", @$group_by_r);
		$q = "SELECT $group_by,count(*) FROM Loci a, Genes b WHERE a.locus_id=b.locus_id GROUP BY $group_by";
		}
	else{ $q = "SELECT count(*) FROM Loci a, Genes b WHERE a.locus_id=b.locus_id"; }

	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "genes", @$_[0..($#$_-1)], qw/Total NA/, $$_[$#$_]), "\n";
		}	
	}

sub sum_spacers_DR_hclust{
	my ($dbh, $cat, $group_by_r, $cutoff_r) = @_;
	
	# group_by  #
	my @select = qw/b.cutoff/;
	my $select = join(",", @$group_by_r, @select);
	my $cutoff = join(" ", "IN (", join(",", @$cutoff_r), ")");
	my $table = $cat . "s";

	my $q = "SELECT $select,'NA',count(distinct(b.cluster_id)) 
	FROM Loci a, $cat\_hclust b 
	WHERE a.locus_id=b.locus_id 
	AND b.cutoff $cutoff
	GROUP BY $select";	
	$q =~ s/\t|\n/ /g;
	
	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", $cat . "s", @$_), "\n";
		}
	}

sub sum_spacers_DR{
	my ($dbh, $cat, $group_by_r) = @_;
	
	# group_by #
	my $select = join(",", @$group_by_r);
	my $table = $cat . "s";

	my $q = "SELECT $select,'NA','NA',count(*) FROM Loci a, $table b WHERE a.locus_id=b.locus_id GROUP BY $select";	
	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "$cat", @$_), "\n";
		}
	}

sub sum_loci{
	my ($dbh, $group_by_r) = @_;

	# by operon/array status #
	my @select;
	push @select, @$group_by_r if @$group_by_r;
	push @select, qw/a.operon_status a.CRISPR_array_status/;
	my $select = join(",", @select);

	my $q = "SELECT $select,count(*) FROM Loci a, Loci b WHERE a.locus_id=b.locus_id GROUP BY $select";
		#print Dumper $q; exit;
	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "loci", @$_), "\n";
		}	

	# total #
	if(@$group_by_r){
		my $group_by = join(",", @$group_by_r);
		$q = "SELECT $group_by,count(*) FROM Loci a, Loci b WHERE a.locus_id=b.locus_id GROUP BY $group_by";
		}
	else{ $q = "SELECT count(*) FROM Loci a, Loci b WHERE a.locus_id=b.locus_id"; }

	foreach (@{$dbh->selectall_arrayref($q)}){
		map{$_ = "NULL" unless $_} @$_;
		print join("\t", "loci", @$_[0..($#$_-1)], qw/Total NA/, $$_[$#$_]), "\n";
		}	
	}

sub make_group_by_sql{
	my ($subtype, $taxon_id, $taxon_name) = @_;
	my @group_by;
	push @group_by, "a.subtype" if $subtype;
	push @group_by, "a.taxon_id" if $taxon_id;
	push @group_by,"a.taxon_name" if $taxon_name;
	#if(@group_by){ return join(",", @group_by); }
	#else{ return ""; }
	return \@group_by;
	}
	
sub entry_count{
	my ($dbh, $table_list_r) = @_;
	my %table;
	foreach my $table (@$table_list_r){
		my $q = "SELECT count(*) FROM $table";
		my $res = $dbh->selectall_arrayref($q);
		$table =~ tr/A-Z/a-z/;
		$table{$table} = $$res[0][0];
		}
		#print Dumper %table; exit;
	return \%table;
	}

sub list_tables{
	my $dbh = shift;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}


__END__

=pod

=head1 NAME

CLdb_summary.pl -- summary stats on CLdb

=head1 SYNOPSIS

CLdb_summary.pl [flags] > summary.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -subtype

Group summary by subtype?

=item -id

Group summary by taxon_id?

=item -name

Group summary by taxon_name?

=item -cutoff

Which Spacer/DirectRepeat clustering cutoffs to summarize (>= 1 argument)? [1]

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_summary.pl

=head1 DESCRIPTION

Get summary stats on CLdb. 

The output is a tab-delimited file for easy parsing.

=head1 EXAMPLES

=head2 Total summary

CLdb_summary.pl -d CLdb.sqlite 

=head2 Summarize by subtype

CLdb_summary.pl -d CLdb.sqlite -s

=head2 Summarize by subtype & taxon_name

CLdb_summary.pl -d CLdb.sqlite -s -n

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

