package CLdb::genbank;

use 5.006;
use strict;
use warnings FATAL => 'all';
use Data::Dumper;
use Carp qw/carp confess/;

use Bio::SeqIO;


use base 'Exporter';
our @EXPORT_OK = '';



=head1 NAME

CLdb - The great new CLdb!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT

=head1 SUBROUTINES/METHODS


=head2 get_loci_scaf

Getting scaffold info for select CRISPR loci.

=head3 IN

$dbh -- DBI connection object
$dir -- directory to sequence files associated with CLdb (eg., fasta or genbank)
$loci_r -- arrayref of locus_ids
%opts -- options

=head3 OUT

hashref -- locus_id : category : info

=cut 

push @EXPORT_OK, 'get_loci_scaf';

sub get_loci_scaf{
  my $dbh = shift or confess "Provide dbh\n";
  my $fileDir = shift or confess "Provide seq file directory\n";
  my $loci_r = shift or confess "Provide loci\n";
  my %opts = @_;

  # input check
  $opts{-seqType} ||= 'fasta';
  map{ $opts{$_} =~ tr/A-Z/a-z/ } keys %opts;  # values to lower case

  # getting fasta/genbank files associated with each locus
  my $sql = <<HERE;
SELECT locus_id, genbank_file, fasta_file, scaffold
FROM loci
WHERE locus_id = ?
HERE
  
  my $sth = $dbh->prepare($sql);
  my $colnames = $sth->{NAME_lc_hash};

  ## query
  my @tbl;
  foreach my $locus_id (@$loci_r){
    $sth->bind_param(1, $locus_id);
    $sth->execute();
    confess "$DBI::err\n" if $DBI::err;
    my $ret = $sth->fetchall_arrayref();
    map{ push @tbl, $_ } @$ret;
  }
  

  # unpack
  map{ exists $colnames->{$_} or confess "KeyError: '$_'\n" }
    qw/fasta_file scaffold locus_id/;
 
  my $file_idx;
  if($opts{-seqType} eq 'fasta'){
    $file_idx = $colnames->{fasta_file};
  }
  elsif($opts{-seqType} eq 'genbank'){
    $file_idx = $colnames->{genbank_file};
  }
  else{ confess "seqType ".$opts{-seqType}." not supported\n" }
  my $scaf_idx = $colnames->{scaffold}; 
  my $locus_id_idx = $colnames->{locus_id};
  
  # grouping by fasta file and scaffold
  ## fasta_file : scaffold : locus_id : 1
  my %ffscaf;
  foreach my $row (@tbl){
    $ffscaf{$row->[$file_idx]}{$row->[$scaf_idx]}{$row->[$locus_id_idx]} = 1;
  }

  # loading fastas and getting info
  my %retInfo;
  foreach my $seq_file (keys %ffscaf){
    printf STDERR "Processing: '%s'\n", $seq_file;
    my $filepath = File::Spec->catfile($fileDir, $seq_file);
    my $infh = Bio::SeqIO->new(
			       -file => $filepath,
			       -format => $opts{-seqType}
			       );
    while(my $seq = $infh->next_seq){
      if( exists $ffscaf{$seq_file}{ $seq->id } ){
	# scaffold found, returning info for each locus_id
	foreach my $locus_id (keys %{$ffscaf{$seq_file}{$seq->id}}){
	  $retInfo{$locus_id}{-len} = $seq->length if $opts{-len};
	  $retInfo{$locus_id}{-seq} = $seq->seq if $opts{-seq};
	}
      }
    }
  }

  #print Dumper %retInfo; exit;
  return \%retInfo;
}




=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.




=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CRISPR_db


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

1; # End of CRISPR_db
