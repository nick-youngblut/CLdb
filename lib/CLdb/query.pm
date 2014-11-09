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
		     n_entries
		     join_query_opts
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
  my $rv = $sth->execute() or croak $dbh->err;
  my $ret = $sth->fetchall_arrayref();

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


=head2 getLociSpacerInfo

Querying CLdb with join on loci & spacers

=head3 IN

hash of args:
dbh :  $; dbh objeckt
extra_sql :  $; refine sql query
columns  :  \@; columns to return

=cut

push @EXPORT_OK, 'getLociSpacerInfo';

sub getLociSpacerInfo{
  my %h = @_;
  my $dbh = exists $h{dbh} ? $h{dbh} : confess "Provide a dbh object";
  my $extra_sql = exists $h{extra_sql} ?
    $h{extra_sql} : "";
  my $columns_r = exists $h{columns} ? 
    $h{columns} : confess "Provide columns to select";
  
  # manditory columns
  unshift @$columns_r, 'loci.locus_id';
  unshift @$columns_r, 'spacers.spacer_id';

  my $sql = join(" ", 
		 "SELECT", 
		 join(",",  @$columns_r),
		 "FROM loci, spacers",
		 "WHERE loci.locus_id = spacers.locus_id",
		 $extra_sql
		 );
  
  # querying
  $dbh->{FetchHashKeyName} = 'NAME_lc';
  my $sth = $dbh->prepare($sql) or confess $dbh->err;
  $sth->execute or confess $dbh->err;
  my $ret_r = $sth->fetchall_hashref(['locus_id', 'spacer_id']);

  # print Dumper $ret_r; exit;
  return $ret_r;
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

etting spacer or DR cluster repesentative sequence from either table

=head3 IN
 
 $dbh = DBI database object
 $opts_r = hash of options 

=head3 Options
 
 spacer_DR_b = spacer or DR [spacer]
 refine_sql = refinement sql (part of WHERE; must start with 'AND')
 cutoff = cluster sequenceID cutoff [1]

=head3 OUT

$%{name} => seq

name = 'locus_id|element|element_id|cluster_id|cluster_cutoff'

=cut

push @EXPORT_OK, 'get_array_seq';

sub get_array_seq{	
  my $dbh = shift or confess "ERROR: provide a dbh object\n";
  my $opts_r = shift or confess "ERROR: provide a hashref of options\n";
  my $cutoff = exists $opts_r->{cutoff} ? $opts_r->{cutoff} : 1;
  my $refine_sql = exists $opts_r->{refine_sql} ? 
    $opts_r->{refine_sql} : '';
  exists $opts_r->{spacer_DR_b} or confess "ERROR: provide 'spacer_DR_b' as an arg";

  # getting table info #
  my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
  ($tbl_oi, $tbl_prefix) = ("DRs","DR") if $opts_r->{"spacer_DR_b"};		# DR instead of spacer
  
  # columns:
  #	locus_ID
  #	spacer|DR
  #	spacer/DR_ID
  # 	clusterID
  # 	sequence
  #     sense_strand

  my $query; 
  if( $opts_r->{by_cluster} ){ # selecting by cluster        
    $query = "
SELECT 
'NA',
'$tbl_prefix',
'NA',
$tbl_prefix\_clusters.cluster_ID,
$tbl_prefix\_clusters.Rep_sequence,
loci.array_sense_strand
FROM $tbl_prefix\_clusters, loci 
WHERE loci.locus_id = $tbl_prefix\_clusters.locus_id 
AND $tbl_prefix\_clusters.cutoff = $cutoff
$refine_sql
GROUP BY $tbl_prefix\_clusters.cluster_ID 
ORDER BY $tbl_prefix\_clusters.cluster_ID
";
  }
  else{  # selecting all spacers or DRs
    $query = "
SELECT
$tbl_oi.Locus_ID,
'$tbl_prefix',
$tbl_oi.$tbl_prefix\_ID,
'NA',
$tbl_oi.$tbl_prefix\_sequence,
loci.array_sense_strand
FROM $tbl_oi, loci
WHERE loci.locus_id = $tbl_oi.locus_id
$refine_sql";
  }
    
  # query db #
#  my $ret = $dbh->selectall_arrayref($query);
  my $sth = $dbh->prepare($query);
  $sth->execute() or confess $dbh->err;
  my $ret = $sth->fetchall_arrayref();
  confess "ERROR: no matching $tbl_prefix entries!\n"
    unless defined $$ret[0];
 
  # parsing output
  my %fasta;
  foreach my $row (@$ret){
    # revcomp sequence if array_sense_strand == -1
    ## ACTUALLY: don't need to revcomp because cluster rep sequences already rev-comped
    croak "ERROR: sense strand must be 1 or -1\n"
      unless $row->[5] == 1 or $row->[5] == -1;

    my $seq_name = join("|", @{$row}[0..3], 
		       $opts_r->{by_cluster} ? $cutoff : 'NA');

    $fasta{ $seq_name } = $row->[4];
  }

  return \%fasta;
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

$%{name}=>seq

name = 'locus_id|element|element_id|cluster_id|cluster_cutoff'

=cut

push @EXPORT_OK, 'get_array_seq_preCluster';

sub get_array_seq_preCluster{
  my $dbh = shift or confess "ERROR: provide a dbh object\n";
  my $opts_r = shift or confess "ERROR: provide a hashref of options\n";

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
loci.Locus_ID,
$spacer_DR,
$tbl_oi.$tbl_prefix\_ID,
'NA',
$tbl_oi.$tbl_prefix\_sequence,
loci.array_sense_strand
FROM loci,$tbl_oi
WHERE loci.locus_id = $tbl_oi.locus_id
$opts_r->{'refine_sql'}";

  # query db #
  my $sth = $dbh->prepare($query) or croak $dbh->err;
  $sth->execute() or confess $dbh->err;
  my $ret = $sth->fetchall_arrayref() or croak $dbh->err;
  die " ERROR: no matching entries!\n"
    unless defined $$ret[0];

  # making fasta #
  my %fasta;
  foreach my $row (@$ret){
    croak "\nERROR: array_sense_strand not set for all entries! Set with 'CLdb_setSenseStrand.pl'\n\n"
      unless defined $row->[5];  # array_sense_strand must exist

    # revcomp sequence if array_sense_strand == -1
    croak "ERROR: sense strand must be 1 or -1\n"
      unless $row->[5] == 1 or $row->[5] == -1;
    $row->[4] = revcomp($row->[4]) if $row->[5] == -1;
    
    my $seq_name = join("|", @{$row}[0..3], 'NA' );
    $fasta{ $seq_name } = $row->[4];
  }
  
  #print Dumper %fasta; exit;
  return \%fasta;
}


=head2 get_array_elem_pos

=head3 description

Getting the position infomation on array elemnts.
Will also get cluster info.

=head3 IN
 
 $dbh = DBI database object
 $opts_r = hash of options 

=head3 Options
 
 refine_sql = refinement sql (part of WHERE; must start with 'AND')
 cutoff = cluster sequenceID cutoff [1]
 spacer_DR_b = booleen on spacers/DRs (0 = spacers)

=head3 OUT

locus_id : spacer_id : {hashref of columns}

=cut

push @EXPORT_OK, 'get_array_elem_pos';

sub get_array_elem_pos{	
  my $dbh = shift or confess "ERROR: provide a dbh object\n";
  my $opts_r = shift or confess "ERROR: provide a hashref of options\n";

  # opts
  my $cutoff = exists $opts_r->{cutoff} ? $opts_r->{cutoff} : 1;
  my $refine_sql = exists $opts_r->{refine_sql} ? 
    $opts_r->{refine_sql} : '';

  # getting table info #
  my ($tbl_oi, $tbl_prefix) = ("spacers","spacer");	
  ($tbl_oi, $tbl_prefix) = ("DRs","DR") if $opts_r->{"spacer_DR_b"};		# DR instead of spacer
      
  # query
  my $query = <<HERE;
SELECT
loci.Locus_ID,
$tbl_oi.$tbl_prefix\_ID,
loci.array_sense_strand,
$tbl_oi.$tbl_prefix\_start,
$tbl_oi.$tbl_prefix\_end,
$tbl_prefix\_clusters.cluster_ID
FROM loci, $tbl_oi, $tbl_prefix\_clusters
WHERE loci.locus_id = $tbl_oi.locus_id
AND loci.locus_id = $tbl_prefix\_clusters.locus_id
AND $tbl_oi.$tbl_prefix\_ID = $tbl_prefix\_clusters.$tbl_prefix\_ID
AND $tbl_prefix\_clusters.cutoff = $cutoff
$refine_sql
HERE
 
  # execute
  my $sth = $dbh->prepare($query);
  $sth->execute() or confess $dbh->err;
  my $colnames_r = $sth->{NAME_lc_hash};
  my $ret = $sth->fetchall_hashref( [qw/Locus_ID Spacer_ID/] );
  ## error handling
  confess "ERROR: no matching entries for query:\n\t'$ret'\n"
    unless scalar keys %$ret > 0;
 
  return $ret;
}


=head2 get_array_start_end

Getting array_start-end for select loci

=head3 IN

$dbh -- DBI object
$loci_r -- arrayref -- locus_ids

=head3 OUT

locus_id : columeName : value

=cut

push @EXPORT_OK, 'get_array_start_end';

sub get_array_start_end{
  my $dbh = shift or confess "Provide dbh object\n";
  my $loci_r = shift or confess "Provide arrayref of locus_ids\n";

  my $sql = <<HERE;
SELECT array_start, array_end
FROM loci
WHERE locus_id = ?
HERE

  my $sth = $dbh->prepare($sql);
 
  my %tbl;
  foreach my $locus_id (@$loci_r){
    $sth->bind_param(1, $locus_id);
    $sth->execute() or confess $dbh->err;
    
    while(my $row = $sth->fetchrow_arrayref()){
      $tbl{$locus_id}{array_start} = $row->[0];
      $tbl{$locus_id}{array_end} = $row->[1];
    }
  }
 
  # check
  confess "No entries for sql:\n\t'$sql'\n"
    unless scalar keys %tbl > 0;

  # return
  return \%tbl;
}


=head2 get_leader_info

Getting all info from leaders table

=head3 Out

locus_id : columnName : value

=cut

push @EXPORT_OK, 'get_leader_info';

sub get_leader_info{
  my $dbh = shift or confess "Provide DBI object\n";
  my $locus_ids = shift or confess "Provide arrayref of locus_ids\n";
  
  my $sql = <<HERE;
SELECT *
FROM leaders
WHERE locus_id = ?
HERE

  my $sth = $dbh->prepare($sql);


  my %tbl;
  foreach my $locus_id (@$locus_ids){
    $sth->bind_param(1, $locus_id);

    $sth->execute() or confess $dbh->err;
    my $ret = $sth->fetchall_hashref('Locus_ID');

    map{ $tbl{$_} = $ret->{$_} } keys %$ret;
  }

 # print Dumper %tbl; exit;
  return \%tbl;
}


=head2 join_query_opts

Making sql refinement statement.
Joining query options for with 'AND'.

=cut

sub join_query_opts{
  my ($vals_r, $cat) = @_;
  
  return "" unless @$vals_r;	
  
  map{ s/"*(.+)"*/"$1"/ } @$vals_r;
  return join("", " AND loci.$cat IN (", join(", ", @$vals_r), ")");
}


=head2 table_exists

Determine whether table exists in CLdb

=head3 IN

$dbh :  dbh object
$table :  table name

=head3 OUT

exists ? 1 : 0

=cut

push @EXPORT_OK, 'table_exists';

sub table_exists {
  my $dbh = shift or confess "ERROR: provie a dbh object\n";
  my $table = shift or confess "ERROR: provide a table name\n";

  my $res = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", "tbl_name"); 

  grep(/^$table$/i, keys %$res) ? 1 : 0;  # 1 if exists
}

sub n_entries {
# getting number of entries in a table #
  my ($dbh, $table) = @_;
  confess "ERROR: Provide a DBI database object!\n" if ! defined $dbh;
  confess "ERROR: Provide a CLdb table name!\n" if ! defined $table;
  
  my $q = "SELECT count(*) FROM $table";
  my $res = $dbh->selectrow_arrayref($q);
  $res->[0] = 0 unless defined $res->[0];
  return $res->[0];
}


=head2 get_leader_pos

Getting spacer or DR sequence from either table

=head3 IN

dbh -- DBI database object
opts -- hash of options 

=head3 Opts

extra_query = extra sql
join_sql = "AND" statements 

=cut

sub get_leader_pos{	
  my ($dbh, $opts_r) = @_;
  
  my $query = "SELECT leaders.locus_ID, leaders.leader_start, leaders.leader_end
		FROM leaders, loci WHERE loci.locus_id = leaders.locus_id $opts_r->{'join_sql'}";
  
  $query =~ s/[\n\t]+/ /g;
  $query = join(" ", $query, $opts_r->{"extra_query"});
  
  # query db #
  my $sth = $dbh->prepare($query) or confess $dbh->err;
  $sth->execute() or confess $dbh->err;
  my $ret = $sth->fetchall_arrayref($query);
  confess " ERROR: no matching entries!\n"
    unless defined $$ret[0];
  
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
