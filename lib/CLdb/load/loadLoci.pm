package CLdb::load::loadLoci;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

## export #
use base 'Exporter';
our @EXPORT_OK = '';

## CLdb
use CLdb::utilities qw/
			lineBreaks2unix
			/;



=head1 NAME

CLdb::load::loadLoci

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing the tab-delimited loci 
file and loading the entries into CLdb

=head1 EXPORT_OK

=cut

=head1 SUBROUTINES


=head2 unix_line_breaks

=head3 IN

$locis_r :  hash_ref of loci table
$db_path :  database path

=head3 OUT

=cut 

push @EXPORT_OK, 'unix_line_breaks';

sub unix_line_breaks{
  my ($loci_r, $db_path) = @_;

  print STDERR "### checking line breaks for all external files (converting to unix) ###\n";

  #my @file_columns = qw/genbank_file fasta_file array_file/;
  my %file_cols = (
                   "genbank_file" => "genbank",
                   "fasta_file" => "fasta",
                   "array_file" => "array");

  foreach my $locus_id (keys %$loci_r){
    foreach my $file_col (keys %file_cols){
      if(exists $loci_r->{$locus_id}{$file_col}){
        next unless $loci_r->{$locus_id}{$file_col};                            # if no file; nothing to check
        my $file_path = join("/", $db_path, $file_cols{$file_col},
                             $loci_r->{$locus_id}{$file_col});

        print STDERR " processing: $file_path\n";
        lineBreaks2unix($file_path, 1);
      }
    }
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

    perldoc CLdb::arrayBlast::loadLoci


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
