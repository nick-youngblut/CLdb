package CLdb::arrayBlast::PAM_SEED;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;

## CLdb
use CLdb::seq qw/revcomp/;

# export #
use base 'Exporter';
our @EXPORT_OK = ();
	
=head1 NAME

CLdb::arrayBlast::PAM_SEED

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for PAMa

=cut


=head2 make_pam_index

Making start-end index for PAM region
based on user input.

=head3 IN

$PAM :  \@; start-end

=head3 OUT

\@ :  [start, end, up|downstream]; start-end are + values

=cut

push @EXPORT_OK, 'make_pam_index';

sub make_pam_index{
  my $PAM = shift || confess "provide PAM start-end";
  
  # up/down-stream
  my $stream;
  if( $PAM->[0] > 0 and $PAM->[1] > 0){
    $stream = 'down';
  }
  elsif( $PAM->[0] < 0 and $PAM->[1] < 0){
    $stream = 'up';
    map{ abs $_ } @$PAM->[0..1];  # flipping to + values
  }
  else{ croak "ERROR: -PAM values should not equal 0\n"; }

  # smaller value should be 1st (now that values are all >0)
   ($PAM->[0], $PAM->[1]) = ($PAM->[1], $PAM->[0]) 
     if $PAM->[0] > $PAM->[1];

  return [$PAM->[0], $PAM->[1], $stream].
}


=head 



=head2 get_PAM_SEED

Getting crDNA-protospacer alignments
from the blast srl DS.

See parse_outfmt for fields that may need
to be parsed from blast srl

=head3 IN

hash of args:
blast :  blast srl Ds
outfmt :  outfmt values
queries :  \%{locus_id}{spacer_id}{field} = value; used for selecting specific hits
array :  just alignments of hit to spacers in arrays? [FALSE]
crRNA_ori :  orient by crRNA? [FALSE]
verbose :  verbose [TRUE}

=head3 OUT

Writing fasta to STDOUT

=cut




=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::PAM_SEED


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
