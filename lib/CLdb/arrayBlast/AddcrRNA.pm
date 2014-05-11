package CLdb::arrayBlast::AddcrRNA;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use DBI;
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

=head3 IN

spacer blast hit hash ref

=head3 OUT

array of query IDs

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
  printf STDERR "...Number of unique query IDs extracted: %i\n\n",
    scalar keys %query_IDs;

 # print Dumper %query_IDs; exit;
  return [keys %query_IDs];
}


=head2 detect_clustered_spacers

Detecting which spacers are clusterrs & which are not (single spacers).

Using spacer ID to determine

If single: 'CHAR-INT|spacer|INT|NA|NA'

If cluster: 'NA|spacer|NA|INT|FLOAT'

=head3 IN

array of spacer IDs

=head3 OUT

hash of arrays: spacer IDs grouped by single|cluster

=cut

push @EXPORT_OK, 'detect_clustered_spacers';

sub detect_clustered_spacers{
  my $spacer_r = shift || croak "Provide an array_ref of spacer IDs";

  my %tmp;
  foreach my $spacerID ( @$spacer_r ){
    if( $spacerID =~ /^[^|]+\|spacer\|\d+\|NA\|NA$/ ){  # single
      push @{$tmp{single}}, $spacerID;
    }
    elsif( $spacerID =~ /^NA\|spacer\|NA\|\d+\|[0-9.]+$/){  # cluster
      push @{$tmp{cluster}}, $spacerID;
    }
    else{
      confess "ERROR: do not recognize spacer ID: '$spacerID'";      
    }
  }
  
  #print Dumper %tmp; exit;
  return \%tmp;
}

=head2 queryBySpacer

Query CLdb by ID of spacer

=head3 IN

$dbh :  dbi connection object
$spacerIDs_r :  array_ref of spacerIDs

=head3 OUT

@@ of CLdb entries

=cut

push @EXPORT_OK, 'queryBySpacer';

sub queryBySpacer{
  my $dbh = shift || croak "Provide a dbh object\n";
  my $spacerIDs_r = shift ||  croak "Provide an array_ref of spacerIDs\n";

  my $query = <<END;
SELECT
l.locus_id,
l.fasta_file,
l.scaffold,
l.array_sense_strand,
s.spacer_id,
s.spacer_start,
s.spacer_end
FROM
loci l, spacers s
WHERE l.locus_id=s.locus_id
AND l.locus_id=sc.locus_id
AND s.locus_id = ?
AND s.spacer_id = ?
END

  my $sth = $dbh->prepare($query) or confess $dbh->err;

  # query for each locusID->spacerID
  ## have to query each locusID->spacerID seperately
  my @ret;
  foreach my $ID (@$spacerIDs_r){
    $ID = "F|spacer|1|NA";

    # getting IDs
    my @l = split /\|/, $ID;
    die "ERROR: cannot find locus_id in '$ID'\n"
      unless defined $l[0];
    die "ERROR: cannot find spacer_id in '$ID'\n"
      unless defined $l[2];

    # adding to sql
    $sth->bind_param(1, $l[0]);  # locus_id
    $sth->bind_param(2, $l[2]);  # spacer_id
    $sth->execute() or confess $dbh->err;

    # getting results
    while(my $r = $sth->fetchrow_arrayref()){
      push @ret, $r;
    }
  }
  
  die "ERROR: no entries match query! $!"
    unless scalar @ret; 

  return \@ret;
}


=head2 queryBySpacerCluster

Query CLdb by ID of spacer cluster rep

=head3 IN

$dbh :  dbi connection object
$spacerIDs_r :  array_ref of spacerIDs

=head3 OUT

@@ of CLdb entries; 1 more column than queryBySapcer

=cut

push @EXPORT_OK, 'queryBySpacerCluster';

sub queryBySpacerCluster{
  my $dbh = shift || croak "Provide a dbh object\n";
  my $spacerIDs_r = shift ||  croak "Provide an array_ref of spacerIDs\n";

  my $query = <<END;
SELECT
l.locus_id,
l.scaffold,
l.array_sense_strand,
s.spacer_id,
s.spacer_start,
s.spacer_end,
s.spacer_sequence
FROM
loci l, spacers s, spacer_clusters sc
WHERE l.locus_id=s.locus_id
AND l.locus_id=sc.locus_id
AND s.spacer_id=sc.spacer_id
AND sc.cluster_id = ?
AND sc.cutoff = ?
END
 
  # prepare & get fields of query
  my $sth = $dbh->prepare($query) or confess $dbh->err;
  my $fields = $sth->{NAME_lc} or confess $dbh->err;

  # querying CLdb for each blast query ID
  my %spacer_info;
  foreach my $ID (@$spacerIDs_r){
    my @l = split /\|/, $ID;
    confess "ERROR: cannot find clusterID and/or cluster cutoff in $ID\n"
      unless defined $l[3] and defined $l[4];

    $sth->bind_param(1, $l[3]);  # clusterID
    $sth->bind_param(2, $l[4]);  # cluster cutoff
    
    $sth->execute or confess $dbh->err;
    my $ret = $sth->fetchall_arrayref() or confess $dbh->err;
    
    confess "ERROR: no entries match query for $ID! $!"
      unless scalar @$ret; 
    
    # loading hash of spacer info
    foreach my $r ( @$ret ){  # each spacer
      my %tmp;
      foreach my $i (0..$#$fields){
	confess "ERROR: cannot find field $i for entry from $ID\n"
	  unless defined $r->[$i];
	$tmp{ $fields->[$i] } = $r->[$i];  # field => value
      }
     push @{$spacer_info{$ID}}, \%tmp; 
    }
    
  }
 
  print Dumper %spacer_info; exit;
  ### clustering is wrong!
  ### TODO: FIX clustering 
  #return $ret;
}
    

=head2 groupByFastaFile

grouping spacers by which fasta file they belong to

=head3 IN

%{single|cluster} => @@(CLdb entries)

=head3 OUT

%{single|cluster} => genome => locus_id => spacer_id => {cat} => value

=cut

push @EXPORT_OK, 'groupByFastaFile';

sub groupByFastaFile{
  my $ret = shift || die "ERROR: provide a %@@ of CLdb spacer entries\n";

  my %byGenome;
  foreach my $cat (keys %$ret){   # cluster or single
    my $fields_r = shift @{$ret->{$cat}} or last;   # fields of entries
    #my $fasta_file = shift @$fields_r or confess $!;

    # converting each CLdb entry returned by query
    foreach my $r ( @{$ret->{$cat}} ){  
      my %entry;
      for my $i (1..$#$fields_r){  # skipping fasta file
	confess "ERROR: cannot find field %s in entry\n", $fields_r->[$i]
	  unless defined $r->[$i];
	$entry{ $fields_r->[$i] } = $r->[$i];
      }
      my $fasta_file = $r->[0];
      push @{$byGenome{$fasta_file}}, \%entry;  # {blast_db}=>[field=>value]
    }
  }

#  print Dumper %byGenome; exit;
  return \%byGenome;
}


=head2 getSpacerRegion

# adding to *srl data structure
## {BlastOutput_iterations}=>{Iteration}=>
##  [ 'Iteration_crRNA'=>{fasta_file}=>
##    {scaffold|strand|spacer_start|spacer_end|spacer_seq|
##     region_start|region_end|region_seq|ext} ]

=head3 IN

$byGenome_r :  $%{fasta_file}=>[ field => value ]
# each in array is a hit
# Datastructure produced by groubByFastaFile()
$CLdb_HOME :  home directory of CLdb.sqlite
$extension :  bp to extend spacer reion on either end 

=head3 OUT

=cut

push @EXPORT_OK, 'getSpacerRegion';

sub getSpacerRegion{
  my %h = @_;
#  my $spacers_r = exists $h{blast} ? $h{blast} :
#    confess "ERROR: provide the decoded blast .srl variable\n";
  my $byGenome_r = exists $h{byGenome} ? $h{byGenome} :
    confess "ERROR: provide a %%@% of hits by genome\n";
  my $CLdb_HOME = exists $h{CLdb_HOME} ? $h{CLdb_HOME} : 
    confess "ERROR: provide \$CLdb_HOME\n";
  my $ext = exists $h{extension} ? $h{extension} : 
    confess "ERROR: provide spacer region extension";
  

  # checking for directory of genome fasta files
  my $fasta_dir = "$CLdb_HOME/fasta";
  confess "ERROR: cannot find '$fasta_dir'" unless -d $fasta_dir;

  # parsing out spacers for each hit
  foreach my $fasta_file (keys %$byGenome_r){
    confess "ERROR: cannot find '$fasta_file'\n"
      unless -e "$fasta_dir/$fasta_file";

    # status
    print STDERR "...parsing from genome: $fasta_file\n";


    # loading genome fasta as hashref
    my $genome_fasta = read_fasta(file => "$fasta_dir/$fasta_file");
    
    # subsetting each spacer sequence
    foreach my $hit (@{$byGenome_r->{$fasta_file}}){
      my $scaf = $hit->{scaffold};
      # sanity checks
      confess "ERROR: cannot find scaffold '$scaf'\n"
	unless exists $genome_fasta->{$scaf};  # scaffold must exist in genome fasta
      confess "ERROR: spacer_start > spacer_End\n"
	unless $hit->{spacer_start} <= $hit->{spacer_end};  # assuming stary-end on + strand
      
      # floor & ceiling for extension; 1-indexed
      my $region_start = $hit->{spacer_start} - $ext;
      $region_start = 1 if $region_start < 1;
      my $region_end = $hit->{spacer_end} + $ext;
      my $scaf_len = length($genome_fasta->{$scaf});
      $region_end = $scaf_len if $region_end > $scaf_len; 

      # spacer crDNA
      my $crDNA = substr($genome_fasta->{$scaf},
			 $region_start - 1,   # 0-indexing for substr
			 $region_end - $region_start  # 0-indexing
			 );

      # revcomp if array_sense_strand == -1
      $crDNA = revcomp($crDNA) if $hit->{array_sense_strand} == -1;

      # revcomp spacer sequence if array_sense_strand == -1
      $hit->{spacer_sequence} = revcomp($hit->{spacer_sequence})
	if $hit->{array_sense_strand} == -1;

      # adding info to $byGenome
      $hit->{region_start} = $region_start;
      $hit->{region_end} = $region_end;
      $hit->{region_seq} = $crDNA;
      $hit->{ext} = $ext;

    }
  }
 
#  print Dumper %$byGenome_r; exit;
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
