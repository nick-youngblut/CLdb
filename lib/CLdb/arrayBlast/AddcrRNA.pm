package CLdb::arrayBlast::AddcrRNA;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use DBI;

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

  my @query_IDs;

  # iterating through each hit
  foreach my $run (keys %$spacer_r){
    
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      # query
      my $query_id = exists $iter->{'Iteration_query-def'} ?
	$iter->{'Iteration_query-def'} :
	  die "ERROR: no query id in blast run $run\n";
      
      push @query_IDs, $query_id;
    }
  }

  # status
  printf STDERR "%i query IDs extracted from the blast hit file\n",
    scalar @query_IDs;

  return \@query_IDs;
}


=head2 detect_clustered_spacers

Detecting which spacers are clusterrs & which are not (single spacers).

Using spacer ID to determine

If single: 'CHAR-INT|spacer|INT|NA'

If cluster: 'NA|spacer|NA|FLOAT_INT'

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
    if( $spacerID =~ /^[^|]+\|spacer\|\d+\|NA$/ ){  # single
      push @{$tmp{single}}, $spacerID;
    }
    elsif( $spacerID =~ /^NA\|spacer\|NA\|[0-9.]+_\d+$/){  # cluster
      push @{$tmp{cluster}}, $spacerID;
    }
    else{
      croak "ERROR: do not recognize spacer ID: '$spacerID'";      
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
s.spacer_id,
s.spacer_start,
s.spacer_end
FROM
loci l, spacers s
WHERE l.locus_id=s.locus_id
AND s.locus_id = ?
AND s.spacer_id = ?


END

  my $sth = $dbh->prepare($query);

  # query for each locusID->spacerID
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
    $sth->execute();

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

  # getting cluster IDs
  my @clusterIDs;
  foreach my $ID (@$spacerIDs_r){
    my @l = split /\|/, $ID;
    die "ERROR: no cluster ID in '$ID'\n"
      unless defined $l[3] and $l[3] ne 'NA';
    push @clusterIDs, $l[3];
  }
  my $clusterIDs = join(",", map{ join('', '"', $_, '"') } @clusterIDs);

  my $query = <<END;
SELECT
l.locus_id,
l.fasta_file,
l.scaffold,
s.spacer_id,
s.spacer_start,
s.spacer_end,
sc.cluster_id
FROM
loci l, spacers s, spacer_clusters sc
WHERE l.locus_id=s.locus_id
AND s.spacer_id=sc.spacer_id
AND sc.cluster_id IN ($clusterIDs)

END
  
  my $ret = $dbh->selectall_arrayref($query);
  die "ERROR: no entries match query! $!"
    unless scalar @$ret; 

  return $ret;
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

  my %group;
  foreach my $cat (keys %$ret){
    foreach my $entry ( @{$ret->{$cat}} ){
      # index: $entry[0] = locus_id, [1] = fasta_file, [2] = scaffold
      $group{$cat}{$entry->[1]}{$entry->[0]}{$entry->[3]}{scaffold} = $entry->[2];
      $group{$cat}{$entry->[1]}{$entry->[0]}{$entry->[3]}{start} = $entry->[4];
      $group{$cat}{$entry->[1]}{$entry->[0]}{$entry->[3]}{end} = $entry->[5];
      $group{$cat}{$entry->[1]}{$entry->[0]}{$entry->[3]}{cluster_id} = defined $entry->[6] ?
	$entry->[6] : undef;
    }
  }

  #print Dumper %group; exit;
  return \%group;
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
