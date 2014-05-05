package CLdb::arrayBlast::AddFullQuery;

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

CLdb::arrayBlast::AddFullQuery

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

=head1 EXPORT_OK

=cut



=head2 addFullQuery

Adding full length query sequence to blast (qseqfull).
Also adding full query start-end (qfullstart, qfullend).

=head3 IN


=head3 OUT


=cut 

push @EXPORT_OK, 'addFullQuery';

sub addFullQuery{
  my $spacer_r = shift or croak "Provide a spacer hash ref";
  my $fasta_r = shift or croak "Provide a fasta hash ref";
  my %opt = @_;

  # status
  my %status;

  # iterating through each hit
  foreach my $run (keys %$spacer_r){
    # database 
    my $db =  $spacer_r->{$run}{'BlastOutput_db'};
    $db = (File::Spec->splitpath($db))[2];
    
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      next unless exists $iter->{'Iteration_hits'} and 
	ref $iter->{'Iteration_hits'} eq 'HASH';
           
      # iterating through hits
      my $hits_ref = $iter->{'Iteration_hits'}{'Hit'}; # hits in iteration (array_ref or ref to hash)
      foreach my $hit ( @$hits_ref){
	# subjectID
	my $sseqid = $hit->{Hit_id};
	my $subj = join("__", $db, $sseqid);
	
	# iterating through hsp
	my $hsp_ref = $hit->{Hit_hsps}{Hsp};	
	foreach my $hsp (@{$hsp_ref}){
	  
	  # hsp values
	  my $sstart = $hsp->{'Hsp_hit-from'};  # hit location on subject
	  my $send = $hsp->{'Hsp_hit-to'};
	  my $evalue = $hsp->{'Hsp_evalue'};
	  my $aln_len = $hsp->{'Hsp_align-len'};
	  
	  # strand
	  my $strand = $sstart <= $send ? 1 : -1;	 
	  ($sstart, $send) = ($send, $sstart) if 
	    $sstart > $send; # start must be <= end

	    
	}
      }
    }
  }

  # status

}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::AddFullQuery


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
