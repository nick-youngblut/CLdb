package CLdb::arrayBlast::AddcrRNA;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

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


=head2 queryBySpacerCluster

Query CLdb by ID of spacer cluster rep

=head3 IN

$dbh :  dbi connection object
$spacerIDs_r :  array_ref of spacerIDs

=head3 OUT

=cut

push @EXPORT_OK, 'queryBySpacerCluster';

sub queryBySpacerCluster{
  my $dbh = shift || croak "Provide a dbh object\n";
  my $spacerIDs_r = shift ||  croak "Provide an array_ref of spacerIDs\n";
  
  foreach my $spacerID ( @$spacerIDs_r ){
    
  }

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
