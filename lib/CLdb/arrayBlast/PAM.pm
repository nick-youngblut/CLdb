package CLdb::arrayBlast::PAM;

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

CLdb::arrayBlast::PAM

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
  if( $PAM->[0] > 0 and $PAM->[1] > 0 ){
    $stream = 'down';
  }
  elsif( $PAM->[0] < 0 and $PAM->[1] < 0){
    $stream = 'up';
    $PAM->[0] = abs $PAM->[0];  # all values positive
    $PAM->[1] = abs $PAM->[1];    
  }
  else{ croak "ERROR: -PAM values should not equal 0\n"; }

  # smaller value should be 1st (now that values are all >0)
   ($PAM->[0], $PAM->[1]) = ($PAM->[1], $PAM->[0]) 
     if $PAM->[0] > $PAM->[1];

  return [$PAM->[0], $PAM->[1], $stream];
}


=head2 read_protoTable

Reading in table of protospacer sequences
produced by CLdb_arrayBlastGetProto.pl

=head3 IN

hash of args
fh :  file handle
file :  file name (instead of file handle)

=head3 OUT

\%{seq_name} => seq

=cut

push @EXPORT_OK, 'read_protoTable';

sub read_protoTable{
  my %h = @_;
  confess "Either provide 'fh' or 'file'\n"
    unless exists $h{fh} or exists $h{file};

  my $fh;
  exists $h{fh} ? $fh = $h{fh} :
    open $fh, $h{file} or confess $!;

  my %proto;  # fasta hash format
  while(<$fh>){
    chomp;
    next if /^\s*$/;
    next if /^protoFull/;   # skipping header

    my @l = split /\t/, $_, 2;    
    $l[1] = "Seq$." unless defined $l[1];   # default
    confess "ERROR: '$l[0]' does not appear to be a valid sequence\n"
      unless $l[0] =~ /^[A-Za-z.-]+$/;
    
    $proto{ $l[1] } = $l[0];
  }
  close $fh or confess $!;  

  return \%proto;
}


=head2 getPAM

Getting crDNA-protospacer alignments
from the blast srl DS.

See parse_outfmt for fields that may need
to be parsed from blast srl

=head3 IN

$fasta_r :  fasta of protospacer (& possibly crRNA)
$pam_index_r :  [pam_start, pam_end, stream]
$fasta_in :  fasta file name
$table_in :  table file name

=head3 OUT

\%{seq_name} => seq

=cut

push @EXPORT_OK, 'getPAM';

sub getPAM{
  my ($fasta_r, $pam_index_r, $revcomp_b) = @_;

  # IO check
  scalar @$pam_index_r == 3 or confess "PAM index length != 3";

  # each protospacer
  foreach my $seqName (keys %$fasta_r){
    # screening out crRNA
    if($seqName =~ /^>*cr[DR]NA\|/){
      delete $fasta_r->{$seqName};
      next;
    }

    # revcomp if needed
    $fasta_r->{$seqName} = revcomp $fasta_r->{$seqName}
      if $revcomp_b;

    # removing gaps
    $fasta_r->{$seqName} =~ s/[.-]+//g;    

    # splitting by CAPS letters
    my @ext = split /[A-Z]+/, $fasta_r->{$seqName}, 2;  # [0] = up, [1] = down

    # selecting upstream or downstream extension
    my $ext;
    if( $pam_index_r->[2] eq 'up'){  # upstream sequence, reversing before substr
      $ext = $ext[0];
      $ext = reverse $ext;
    }
    elsif( $pam_index_r->[2] eq 'down' ){
      $ext = $ext[1];
    }
    
    # ceiling & floor for index
    my $pam_start = $pam_index_r->[0];
    my $pam_end = $pam_index_r->[1];
    my $ext_len = length $ext;
    $pam_start = 1 if $pam_start < 1;   # shouldn't happen
    $pam_start = $ext_len if $pam_start > $ext_len;
    $pam_end = $ext_len if $pam_end > $ext_len; 
  
    # substr out pam
    $fasta_r->{$seqName} = substr($ext,
				  $pam_start - 1,  # 0 index
				  $pam_end - $pam_start + 1);    

    # flipping back if 'down'
    $fasta_r->{$seqName} = reverse $fasta_r->{$seqName}
      if $pam_index_r->[2] eq 'up';

  }

#  print Dumper %$fasta_r; exit;
  return $fasta_r;
}


=head2 writePAM

Writing out PAM sequence either as a fasta
or a tab-delim table, depending on how the
data was provided. 

=head3 IN

$pam_r :  \%{seqName} => pam_sequence
$fasta_in :  defined ? out as fasta
$table_in :  defined ? out as table

=head3 OUT

writing to STDOUT

=cut

push @EXPORT_OK, 'writePAM';

sub writePAM{
  my ($pam_r, $fasta_in, $table_in) = @_;
  defined $pam_r || confess "Provide PAM sequence as {seqName}{sequence}\n";

  foreach my $seqName (keys %$pam_r){
    if( defined $fasta_in){
      print join("\n", ">$seqName", $pam_r->{$seqName}), "\n"
    }
    elsif( defined $table_in){
      print join("\t", $pam_r->{$seqName}, $seqName), "\n";
    }
    else{ confess "LOGIC ERROR $!\n"; }
      
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

    perldoc CLdb::arrayBlast::PAM


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
