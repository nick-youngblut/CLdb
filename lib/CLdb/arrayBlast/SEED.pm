package CLdb::arrayBlast::SEED;

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

CLdb::arrayBlast::SEED

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for SEED sequences

=cut


=head2 read_proto_aln

Reading a crDNA-protospacer alignment.
Assuming that the aligned protospacer and 
crRNA sequences are adjacent in the fasta.

=head3 IN

fh :  filehandle of fasta
file :  file name (instead of file handle)

=head3 OUT

\%{pair}{seqNaem} => seq

=cut

push @EXPORT_OK, 'read_proto_aln';

sub read_proto_aln{
  my %h = @_;
  confess "Provide either 'fh' or 'file'"
    unless exists $h{fh} or exists $h{file};

  my $fh;
  exists $h{fh} ? $fh = $h{fh} :
    open $fh, $h{file} or confess $!;
    
  my (%fasta, $tmpkey);
  my $pair = -1;
  while(<$fh>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;

    if(/>.+/){
      $pair++;
      s/>//;
      $fasta{int($pair / 2)}{$_} = "";
      $tmpkey = $_;     # changing key
    }
    else{$fasta{int($pair / 2)}{$tmpkey} .= $_; }
  }

  #print Dumper %fasta; exit;
  return \%fasta;
}


=head2 parseProtoBySEED

Parsing the protospacer sequence
into 'total', 'seed', and 'non-seed'

=head3 IN

$fasta_r :  \%{pair}{seqname} => protspacer
$seed_index_r :  

=head3 OUT

=cut

push @EXPORT_OK, 'parseProtoBySEED';

sub parseProtoBySEED{
  my $fasta_r = shift || confess "Provide fasta\n";
  my $seed_index = shift || confess "Provide seed index";

  my %mismatch_sum;
  foreach my $pair (keys %$fasta_r){    
    my $mismatch_r = makeMismatchIndex($fasta_r->{$pair}, $seed_index);
    my $SE_r = getProtoCoords($fasta_r->{$pair});
  }
  

}

sub getProtoCoords{
  my $fasta_r = shift or confess "Provide a fasta";
  
  my @start_end;
  foreach my $seqname (keys %$fasta_r){
    next if $seqname =~ /^>*cr[DR]NA\|/;

    my $seq = $fasta_r->{$seqname};
    $seq =~ /^([a-z.-]+)([A-Z][A-Z.-]+[A-Z])/;
    @start_end = ( length($1) -1 , length($1) + length($2) -1 );
  }

  #print Dumper @start_end;
  \@start_end;
}

sub makeMismatchIndex{
  my $fasta_r = shift or confess "Provide a fasta";
  my $seed_index = shift or confess "Provde a seed index";

  # check
  my @seqs = values %$fasta_r;
  confess "ERROR: don't have 2 sequences to compare!\n$!"
    unless scalar @seqs == 2;
  my $seq_len = length $seqs[0];
  confess "Sequences are not the same length!"
    unless length $seqs[1] == $seq_len;

  # exploding sequences
  for my $i (0..$#seqs){    
    $seqs[$i] = [split //, $seqs[$i]];
  }

  # determining location of mismatchs
  my @mismatches;
  for my $i (0..($seq_len-1)){
    $mismatches[$i] = ($seqs[0][$i] eq $seqs[1][$i] and
      $seqs[0][$i] !~ /[.-]/ and $seqs[1][$i] !~ /[.-]/) ? 0 : 1;
  }
  

  #print Dumper @mismatches;
  return \@mismatches;
}




=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::SEED


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
