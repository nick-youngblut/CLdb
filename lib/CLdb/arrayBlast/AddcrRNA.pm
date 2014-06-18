package CLdb::arrayBlast::AddcrRNA;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use DBI;
use IPC::Cmd qw/can_run run/;
use CLdb::seq qw/
		  read_fasta
		  revcomp
		/;

# export #
use base 'Exporter';
our @EXPORT_OK = '';

	

=head1 NAME

CLdb::arrayBlast::AddcrRNA

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

These subroutines are for extracting the crDNA (crRNA)
region from the genome fasta.

=head1 EXPORT_OK

=cut



=head2 get_query_IDs

Getting the IDs of each query sequence. 
Extracting blast_db: 'BlastOutput_db'.

=head3 IN

spacer blast hit hash ref

=head3 OUT

$%{query_ID} => blast_db

=cut 

push @EXPORT_OK, 'get_query_IDs';

sub get_query_IDs{
  my $spacer_r = shift or croak "Provide a spacer blast hit hash ref";
  my %opt = @_;

  # status
  print STDERR "Extracting all queries with blast hits from blast hit file...\n"; 

  # iterating through each hit for each blast run
  ## just selecting queries with blast hits
  ## blast hits can be to arrays 
  my %query_IDs;
  foreach my $run (keys %$spacer_r){
    next unless exists $spacer_r->{$run}{'BlastOutput_iterations'};
    # getting blast db
    my $blast_db = exists $spacer_r->{$run}{'BlastOutput_db'} ? 
      $spacer_r->{$run}{'BlastOutput_db'} :
	confess "Could not find 'BlastOutput_db' for run $run\n";

    # each iteration
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){

      next unless exists $iter->{'Iteration_hits'} and 
	$iter->{'Iteration_hits'} !~ /^\s*$/;
      # query
      my $query_id = exists $iter->{'Iteration_query-def'} ?
	$iter->{'Iteration_query-def'} :
	  die "ERROR: no query id in blast run $run\n";
      
      $query_IDs{ $query_id }++;
    }
  }


  # status
  printf STDERR "...Number of total unique query IDs:\t%i\n",
    scalar keys %query_IDs;

  #print Dumper %query_IDs; exit;
  return [keys %query_IDs];
}


=head2 detectClusteredSpacers

Detecting which spacers are clusterrs & which are not (single spacers).

Using spacer ID to determine

If single: 'CHAR-INT|spacer|INT|NA|NA'

If cluster: 'NA|spacer|NA|INT|FLOAT'

=head3 IN

{blast_db}=>{query_id}=>1 

=head3 OUT

hash of arrays: spacer IDs grouped by single|cluster

=cut

push @EXPORT_OK, 'detectClusteredSpacers';

sub detectClusteredSpacers{
  my $spacerIDs_r = shift || croak "Provide an array_ref of spacer IDs";

  my %tmp;
  foreach my $spacerID ( @{$spacerIDs_r} ){
    if( $spacerID =~ /^[^|]+\|spacer\|\d+\|NA\|[0-9-]+$/ ){  # single
      push @{$tmp{single}}, $spacerID;
    }
    elsif( $spacerID =~ /^NA\|spacer\|NA\|\d+\|[0-9-]+$/){  # cluster
      push @{$tmp{cluster}}, $spacerID;
    }
    else{
      confess "ERROR: do not recognize spacer ID: '$spacerID'";      
    }
  }
  
#  print Dumper %tmp; exit;
  return \%tmp;
}


=head2 queryBySpacer

Query CLdb by ID of spacer

=head3 IN

$dbh :  dbi connection object
$spacerIDs_r :  array_ref  of spacerIDs
$CLdb_info :  $%

=head3 OUT

#loading CLdb_info: {fasta_file}=>{scaffold}=>{locus_id}=>{spacer_id}=>{field}=>value

=cut

push @EXPORT_OK, 'queryBySpacer';

sub queryBySpacer{
  my $dbh = shift || confess "Provide a dbh object\n";
  my $spacerIDs_r = shift ||  confess "Provide an array_ref of spacerIDs\n";
  my $CLdb_info_r = shift || confess "Provide a CLdb_info arg\n";

  my $query = <<HERE;
SELECT
l.fasta_file,
l.scaffold,
l.locus_id,
l.array_sense_strand,
l.subtype,
l.cas_status,
l.array_status,
s.spacer_id,
s.spacer_start,
s.spacer_end,
s.spacer_sequence
FROM
loci l, spacers s
WHERE l.locus_id=s.locus_id
AND s.locus_id = ?
AND s.spacer_id = ?
HERE
 
  # prepare & get fields of query
  my $sth = $dbh->prepare($query) or confess $dbh->err;
  my $fields = $sth->{NAME_lc_hash} or confess $dbh->err;

  # querying CLdb for each blast query ID
  my %spacer_info;
  foreach my $ID (@$spacerIDs_r){
    my @l = split /\|/, $ID;
    confess "ERROR: cannot find clusterID and/or cluster cutoff in $ID\n"
      unless defined $l[3] and defined $l[4];

    $sth->bind_param(1, $l[0]);  # locus_id
    $sth->bind_param(2, $l[2]);  # spacer_id
    
    $sth->execute or confess $dbh->err;
    my $ret = $sth->fetchall_arrayref( ) or confess $dbh->err;
    
    # loading CLdb_info: {fasta_file}=>{scaffold}=>{locus_id}=>{spacer_id}=>{field}=>value
    foreach my $r (@$ret){
      foreach my $field (keys %$fields){
	next if $field eq 'fasta_file' or $field eq 'locus_id' 
	  or $field eq 'spacer_id' or $field eq 'scaffold';
	$CLdb_info_r->{ $r->[$fields->{fasta_file}] }
	  { $r->[$fields->{scaffold}] }{ $r->[$fields->{locus_id}] }
	  { $r->[$fields->{spacer_id}] }{ $field } = $r->[ $fields->{$field} ];
      }
      # adding query_ID to add info back to decoded blast srl
      $CLdb_info_r->{ $r->[$fields->{fasta_file}] }
	{ $r->[$fields->{scaffold}] }{ $r->[$fields->{locus_id}] }
	  { $r->[$fields->{spacer_id}] }{ query_id } = $ID;
            
    }
  }
 
#  print Dumper %$CLdb_info_r; exit;
}


=head2 queryBySpacerCluster

Query CLdb by ID of spacer cluster rep

=head3 IN

$dbh :  dbi connection object
$spacerIDs_r :  array_ref of spacerIDs
$CLdb_info :  $%

=head3 OUT

CLdb_info: {fasta_file}=>{scaffold}=>{locus_id}=>{spacer_id}=>{field}=>value

=cut

push @EXPORT_OK, 'queryBySpacerCluster';

sub queryBySpacerCluster{
  my $dbh = shift || confess "Provide a dbh object\n";
  my $spacerIDs_r = shift ||  confess "Provide an array_ref of spacerIDs\n";
  my $CLdb_info_r = shift || confess "Provide a CLdb_info arg\n";
  my %h = @_;

  my $query = <<HERE;
SELECT
l.fasta_file,
l.scaffold,
l.locus_id,
l.array_sense_strand,
l.subtype,
l.cas_status,
l.array_status,
s.spacer_id,
s.spacer_start,
s.spacer_end,
s.spacer_sequence,
sc.cluster_id
FROM
loci l, spacers s, spacer_clusters sc
WHERE l.locus_id=s.locus_id
AND l.locus_id=sc.locus_id
AND s.spacer_id=sc.spacer_id
AND sc.cluster_id = ?
AND sc.cutoff = ?
HERE
 
  # prepare & get fields of query
  my $sth = $dbh->prepare($query) or confess $dbh->err;
  my $fields = $sth->{NAME_lc_hash} or confess $dbh->err;
  #my $levels = ['Fasta_File', 'Scaffold', 'Locus_ID', 'Spacer_ID'];

  # querying CLdb for each blast query ID
  my %spacer_info;
  foreach my $ID (@$spacerIDs_r){
    my @l = split /\|/, $ID;
    confess "ERROR: cannot find clusterID and/or cluster cutoff in $ID\n"
      unless defined $l[3] and defined $l[4];

    $sth->bind_param(1, $l[3]);  # clusterID
    $sth->bind_param(2, $l[4]);  # cluster cutoff
    
    $sth->execute or confess $dbh->err;
    my $ret = $sth->fetchall_arrayref(  ) or confess $dbh->err;
 

    #loading CLdb_info: {fasta_file}=>{scaffold}=>{locus_id}=>{spacer_id}=>{field}=>value
    foreach my $r (@$ret){
      foreach my $field (keys %$fields){
    	next if $field eq 'fasta_file' or $field eq 'locus_id' 
    	  or $field eq 'spacer_id' or $field eq 'scaffold';
    	$CLdb_info_r->{ $r->[$fields->{fasta_file}] }
    	  { $r->[$fields->{scaffold}] }{ $r->[$fields->{locus_id}] }
    	  { $r->[$fields->{spacer_id}] }{ $field } = $r->[ $fields->{$field} ];
      }
    #  adding query_ID to add info back to decoded blast srl
	$CLdb_info_r->{ $r->[$fields->{fasta_file}] }
	{ $r->[$fields->{scaffold}] }{ $r->[$fields->{locus_id}] }
	  { $r->[$fields->{spacer_id}] }{ query_id } = $ID;
    }
  }
 
#  print Dumper $CLdb_info_r; exit;
}
    

=head2 getSpacerRegion

Getting spacer region (crRNA region) for each
spacer query (dereplicated spacers if blasted
with the representative sequence).

=head3 IN

$CLdb_info_r :  {fasta_file}=>{scaffold}=>{locus_id}=>{spacer_id}=>{field}=>value
$CLdb_HOME :  home directory of CLdb.sqlite
$extension :  bp to extend spacer reion on either end 

=head3 OUT

=cut

push @EXPORT_OK, 'getSpacerRegion';

sub getSpacerRegion{
  my %h = @_;
  my $CLdb_info_r = exists $h{CLdb} ? $h{CLdb} :
    confess "Provide CLdb info as $%%%";
  my $CLdb_HOME = exists $h{CLdb_HOME} ? $h{CLdb_HOME} : 
    confess "Provide \$CLdb_HOME\n";
  my $ext = exists $h{extension} ? $h{extension} : 
    confess "Provide spacer region extension";
  

  # checking for directory of genome fasta files
  my $fasta_dir = "$CLdb_HOME/fasta";
  confess "ERROR: cannot find '$fasta_dir'" unless -d $fasta_dir;


  # spacers by genome 
  my %byQuery;  # {queryID}=>[ {field} => {value} ]
  foreach my $fasta_file (keys %$CLdb_info_r){
    confess "ERROR: cannot find '$fasta_file'\n"
      unless -e "$fasta_dir/$fasta_file";

    # status
    print STDERR " Processing genome: $fasta_file...\n";

    # loading genome fasta as hashref
    my $genome = read_fasta(-file => "$fasta_dir/$fasta_file");
    

    # substr each spacer region sequence
    foreach my $scaffold (keys %{$CLdb_info_r->{$fasta_file}}){
      # checking for existence of scaffold in genome
      confess "Cannot find scaffold '$scaffold' in $fasta_file"
	unless exists $genome->{$scaffold};
      my $scaf_len = length $genome->{$scaffold};

      foreach my $locus_id (keys %{$CLdb_info_r->{$fasta_file}{$scaffold}}){
	foreach my $spacer_id (keys %{$CLdb_info_r->{$fasta_file}{$scaffold}{$locus_id}}){
	  my $info_r = $CLdb_info_r->{$fasta_file}{$scaffold}{$locus_id}{$spacer_id};
	  
	  # checking for required fields 
	  map{confess "Cannot find $_ for $locus_id->$spacer_id" unless 
		exists $info_r->{$_} } qw/spacer_start spacer_end array_sense_strand/;	  

	  #print STDERR join(",", $info_r->{spacer_start}, $info_r->{spacer_end}), "\n";

	  # floor & ceiling (1-indexing)
	  my $region_start = $info_r->{spacer_start} - $ext;
	  $region_start = 1 if $region_start < 1;

	  my $region_end = $info_r->{spacer_end} + $ext;
	  $region_end = $scaf_len if $region_end > $scaf_len; 

	  # sanity check
	  if($region_start > $region_end){
	    print STDERR "WARNING for $fasta_file, locus '$locus_id', spacer '$spacer_id': " . 
	      "region_start > region_end ($region_start > $region_end).\n\t" . 
		" Are locus_start, locus_end, array_start, & array_end set correctly?\n";
	    ($region_start, $region_end) = ($region_end, $region_start);
	  }

	  # spacer crDNA
	  my $crDNA = substr(
			     $genome->{$scaffold},
			     $region_start - 1,   # 0-indexing for substr
			     $region_end - $region_start + 1
			    );

	  # revcomp crDNA & spacer sequence if array_sense_strand == -1
	  if($info_r->{array_sense_strand} == -1){
	    $crDNA = revcomp($crDNA);
	    $info_r->{spacer_sequence} = exists $info_r->{spacer_sequence} ?
	      revcomp($info_r->{spacer_sequence}) : 
		confess "Cannot find 'spacer_sequence' for $genome->$scaffold->$locus_id\n";
	  }

	  # loading values
	  $info_r->{region_start} = $region_start;
	  $info_r->{region_end} = $region_end;
	  $info_r->{crDNA} = $crDNA;
	  $info_r->{scaffold} = $scaffold;
	  $info_r->{genome_fasta} = $fasta_file;
	  $info_r->{spacer_id} = $spacer_id;
	  $info_r->{locus_id} = $locus_id;
	  my $uID = join("|", $info_r->{locus_id}, $info_r->{spacer_id}); # unique ID for spacer query (locus_id-spacer_ID)
	  $byQuery{ $info_r->{query_id} }{$uID} = $info_r;
	}
      }
    }
  }

  #print Dumper %byQuery; exit;
  return \%byQuery;
}


=head2 addcrDNAtoBlast

Adding info on crRNA (DNA) back to decoded blast srl.

=head3 IN

$spacer_r :  decoded blast srl
$byQuery_r :  {query_id}=>[ {field}=>value ]

=head3 OUT

=cut

push @EXPORT_OK, 'addcrDNAtoBlast';

sub addcrDNAtoBlast{
  my $spacer_r = shift or confess "Provide a decoded blast srl.";
  my $byQuery_r = shift or confess "Provide the crRNA info grouped by query";

  foreach my $run (keys %$spacer_r){
    next unless exists $spacer_r->{$run}{'BlastOutput_iterations'};

    # getting blastdbfile
    my $blastdbfile = exists $spacer_r->{$run}{'BlastOutput_db'} ? 
      $spacer_r->{$run}{'BlastOutput_db'} :
	confess "Cannot find BlastOutput_db in run $run";

    # each iteration
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){

      # skipping iterations without hits
      next unless exists $iter->{'Iteration_hits'} and 
	$iter->{'Iteration_hits'} !~ /^\s*$/;
      
      # making sure query_ID found in $byQuery_r
      my $query_id = $iter->{'Iteration_query-def'};

      # adding crRNA info
      $iter->{crRNA_info} = exists $byQuery_r->{ $query_id } ?
	$byQuery_r->{ $query_id } :
	confess "Cannot find queryID '$query_id' in decoded blast srl";
      
    }
  }
 
#  print Dumper $spacer_r; exit;
}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::AddcrRNA


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
