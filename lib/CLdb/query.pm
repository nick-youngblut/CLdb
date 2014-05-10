package CLdb::query;

# module use #
## core ##
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use DBI;

## CLdb ##
use CLdb::seq qw/
		  revcomp
		/;


# export #
use base 'Exporter';
our @EXPORT_OK = qw/
		     table_exists
		     n_entries
		     join_query_opts
		     get_array_seq
		     get_leader_pos
		   /;

	
=head1 NAME

CLdb::query

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for querying CLdb

=head1 EXPORT_OK

=cut


=head2 queryLociTable

general query for loci table.
returning all entries & columns in loci table
that match sql 


=head3 IN

dbh =>  dbh object
refine =>  refinement query
columns => array_ref of columns (fields) to return 

=head3 OUT

@@ of matching entries.
Column names are 1st @

=cut

push @EXPORT_OK, 'queryLociTable';

sub queryLociTable{
  my %opt = @_;
  my $dbh = defined $opt{dbh} ? $opt{dbh} : 
    croak "ERROR: provide a dbh object\n";
  my $columns = defined $opt{columns} ?
    join(",", @{$opt{columns}}) : '*';


  my $sql = <<HERE;
SELECT $columns
FROM loci
WHERE locus_id IS NOT NULL
$opt{refine}
HERE
  
  my $sth = $dbh->prepare($sql);
  my $rv = $sth->execute;
  my $ret = $sth->fetchall_arrayref;

  die "ERROR: no matching entries for query:\n'$sql'\n"
    unless scalar @$ret;

  # adding column names
  unshift @$ret, $sth->{NAME};

  # making undef columns = ''
  foreach my $r (@$ret){
    map{ $_ = '' unless defined $_ } @$r;
  }

  #print Dumper @$ret; exit;
  return $ret;
}



=head2 list_columns

listing all columns in specified table (table must exist) 
lists column names

=head3 IN

# $dbh = DBI object
# $tbl = sql table of interest
# $silent_ret = no verbose & exit; return column names

=head3 OUT

=cut

push @EXPORT_OK, 'list_columns';

sub list_columns{
  my ($dbh, $tbl, $silent_ret) = @_;
  my $all = $dbh->selectall_arrayref("pragma table_info($tbl)");
  
  my %tmp;
  foreach (@$all){ 
    $$_[1] =~ tr/A-Z/a-z/;		# lower case for matching
    $tmp{$$_[1]} = 1; 
  }
  
  if(defined $silent_ret){ return \%tmp; }
  else{  print "Columns:\n", join(",\n", keys %tmp), "\n\n";  exit; }
}



=head2 orient_byArray

 orienting (rev-comp if needed) by array start-end
 
 if array_start > array_end, revcomp element
 else: revcomp

=head3 columns:
    
  # locus_ID
    # element_ID
    # eleement_seq
    # element_start
    # element_end
    # array_start
    # array_end

=cut

sub orient_byArray{    
  my ($row) = @_;
  
  if ($$row[5] <= $$row[6]) { # array start > array end
    return [@$row[0..2]];
    
  }
  else {
    $$row[2] = revcomp($$row[2]);
    return [@$row[0..2]];
    
  }
  
}


=head2 orient_byleader

  # orienting (rev-comp if needed) by leader
  # if leader comes prior to element, no rev-comp
  # else: revcomp

=head3 columns:

  # locus_ID
  # element_ID
  # eleement_seq
  # element_start
  # element_end
  # leader_stat
  # leader_end

=cut


sub orient_byleader{  
  my ($row) = @_;
  
  if ($$row[3] >= $$row[6]) { #array_start >= leader_end position
    return [@$row[0..2]];
  }
  else {
    $$row[2] = revcomp($$row[2]);
    return [@$row[0..2]];
       
  }
}



=head2 get_array_seq

=head3 description

Getting spacer or DR cluster repesentative sequence from either table

=head3 IN
 
 $dbh = DBI database object
 $opts_r = hash of options 

=head3 Options
 
 spacer_DR_b = spacer or DR [spacer]
 refine_sql = refinement sql (part of WHERE; must start with 'AND')
 cutoff = cluster sequenceID cutoff [1]

=cut

sub get_array_seq{	
  my ($dbh, $opts_r) = @_;
  
  # checking for opts #
  map{ confess "ERROR: cannot find option: '$_'" 
	 unless exists $opts_r->{$_} } qw/spacer_DR_b/;
  $opts_r->{refine_sql} = '' unless exists $opts_r->{refine_sql};

  # getting table info #
  my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
  ($tbl_oi, $tbl_prefix) = ("DRs","DR") if $opts_r->{"spacer_DR_b"};		# DR instead of spacer
  
  # columns:
  #	locus_ID
  #	spacer|DR
  #	spacer/DR_ID
  # 	clusterID
  # 	sequence

  my $query; 
  if(defined $opts_r->{by_cluster}){ # selecting by cluster
    
    # by group default options #
    $opts_r->{"cutoff"} = 1 unless exists $opts_r->{'cutoff'}; 	
    
    $query = "
SELECT 
'NA',
'$tbl_prefix',
'NA',
$tbl_prefix\_clusters.cluster_ID,
$tbl_prefix\_clusters.Rep_sequence
FROM $tbl_prefix\_clusters, loci 
WHERE loci.locus_id = $tbl_prefix\_clusters.locus_id 
AND $tbl_prefix\_clusters.cutoff = 1";
  }
  else{ # selecting all spacers
    $query = "
SELECT
$tbl_oi.Locus_ID,
'$tbl_prefix',
$tbl_oi.$tbl_prefix\_ID,
'NA',
$tbl_prefix\_clusters.Rep_sequence
FROM $tbl_oi, $tbl_prefix\_clusters, loci
WHERE loci.locus_id = $tbl_oi.locus_id
AND $tbl_oi.locus_id = $tbl_prefix\_clusters.locus_id
AND $tbl_oi.$tbl_prefix\_ID = $tbl_prefix\_clusters.$tbl_prefix\_ID 
AND $tbl_prefix\_clusters.cutoff = 1";
  }
  
  #$query =~ s/[\n\t]+/ /g;
  
  # adding group by  & order if clustering
  $query .= " GROUP BY $tbl_prefix\_clusters.cluster_ID ORDER BY $tbl_prefix\_clusters.cluster_ID"
    if defined $opts_r->{by_cluster}; # selecting by cluster
  
  #print Dumper $query; exit;

  # query db #
  my $ret = $dbh->selectall_arrayref($query);
  confess "ERROR: no matching $tbl_prefix entries!\n"
    unless $$ret[0];
  
  #print Dumper @$ret; exit;
  return $ret;
}


=head2 get_array_seq_preCluster

Selecting array sequences (spacer|DR) if clustering tables
aren't present

=head3 IN

$dbh :  dbh object
$opts :  hashref of options

=head3 Options
 
 spacer_DR_b = spacer or DR [spacer]
 refine_sql = refinement sql (part of WHERE; must start with 'AND')

=head3 OUT

=cut

push @EXPORT_OK, 'get_array_seq_preCluster';

sub get_array_seq_preCluster{
  my ($dbh, $opts_r) = @_;

  # checking for opts #
  map{ die "ERROR: cannot find option: '$_'"
         unless exists $opts_r->{$_} } qw/spacer_DR_b/;
  $opts_r->{refine_sql} = '' unless exists $opts_r->{refine_sql};

  # getting table info (spacer|DR)#
  my ($tbl_oi, $tbl_prefix);
  if( $opts_r->{spacer_DR_b} ){ # DR
    ($tbl_oi, $tbl_prefix) = ("DRs","DR");
  }
  else{  # spacer
    ($tbl_oi, $tbl_prefix) = ("spacers","spacer");
  }

  my $spacer_DR = $dbh->quote($tbl_prefix);

  my $query = "SELECT
b.Locus_ID,
$spacer_DR,
b.$tbl_prefix\_ID,
b.$tbl_prefix\_sequence,
a.array_sense_strand
FROM loci a,$tbl_oi b
WHERE a.locus_id = b.locus_id
$opts_r->{'refine_sql'}";

  # query db #
  my $sth = $dbh->prepare($query) or croak $dbh->err;
  $sth->execute;
  my $ret = $sth->fetchall_arrayref() or croak $dbh->err;
  die " ERROR: no matching entries!\n"
    unless $$ret[0];

  # making fasta #
  my %fasta;
  foreach my $row (@$ret){
    # revcomp sequence if array_sense_strand == -1
    croak "ERROR: sense strand must be 1 or -1\n"
      unless $row->[4] == 1 or $row->[4] == -1;
    $row->[3] = revcomp($row->[3]) if $row->[4] == -1;
    
    my $seq_name = join("|", @{$row}[0..2] );
    $fasta{ $seq_name } = $row->[3];
  }
  
  #print Dumper %fasta; exit;
  return \%fasta;
}


sub join_query_opts{
# making sql refinement statment
## joining query options for with 'AND' 

	my ($vals_r, $cat) = @_;
	
	return "" unless @$vals_r;	
	
	map{ s/"*(.+)"*/"$1"/ } @$vals_r;
	return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
	}

sub table_exists {
# checking for the existence of a table #
	my ($dbh, $table) = @_;
	confess "ERROR: Provide a DBI database object!\n" if ! defined $dbh;
	confess "ERROR: Provide a CLdb table name!\n" if ! defined $table;
	
	my $res = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", "tbl_name"); 
	
	confess "ERROR: '$table' table not found in CLdb!\n" 
		unless grep(/^$table$/i, keys %$res);
	}

sub n_entries {
# getting number of entries in a table #
	my ($dbh, $table) = @_;
	confess "ERROR: Provide a DBI database object!\n" if ! defined $dbh;
	confess "ERROR: Provide a CLdb table name!\n" if ! defined $table;
	
	my $q = "SELECT count(*) FROM $table";
	my $res = $dbh->selectall_arrayref($q);

	return $$res[0][0];
	}

sub get_leader_pos{
# getting spacer or DR sequence from either table #
#-- input --#
# $dbh = DBI database object
# $opts_r = hash of options 
#-- options --#
# extra_query = extra sql
# join_sql = "AND" statements 
	
	my ($dbh, $opts_r) = @_;
	
	my $query = "SELECT leaders.locus_ID, leaders.leader_start, leaders.leader_end
		FROM leaders, loci WHERE loci.locus_id = leaders.locus_id $opts_r->{'join_sql'}";
	
	$query =~ s/[\n\t]+/ /g;
	$query = join(" ", $query, $opts_r->{"extra_query"});
	
	# query db #
	my $ret = $dbh->selectall_arrayref($query);
	confess " ERROR: no matching entries!\n"
		unless $$ret[0];

	# making hash of sequences #
	my %leaders;
	foreach my $row (@$ret){
		$leaders{$$row[0]}{"start"} = $$row[1];
		$leaders{$$row[0]}{"end"} = $$row[2];
		}

		#print Dumper %leaders; exit;
	return \%leaders;	
	}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::query


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=CRISPR_db>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/CRISPR_db>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/CRISPR_db>

=item * Search CPAN

L<http://search.cpan.org/dist/CRISPR_db/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2013 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.

=cut

1; 
