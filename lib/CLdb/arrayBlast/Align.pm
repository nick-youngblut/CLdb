package CLdb::arrayBlast::Align;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use List::MoreUtils qw/all/;
use Bio::Tools::Run::Alignment::Clustalw;
use IO::String;
use Data::Uniqid qw/ uniqid /;
use Clone  qw/ clone /;

# CLdb
use CLdb::seq qw/ read_fasta /;

# export #
use base 'Exporter';
our @EXPORT_OK = ();
	
=head1 NAME

CLdb::arrayBlast::Align

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for aligning protospacer and crDNA

=cut


=head2 alignProto

Aligning crDNA and protospacer sequences
using bioperl implementation of clustalw.

=head3 IN

hash of args:
blast =>  decoded blast srl
verbose => verbose [false]

=head3 OUT

edited blast srl structure

=cut

push @EXPORT_OK, 'alignProto';

sub alignProto{
  my %h = @_;
  my $spacer_r = exists $h{blast} ? 
    $h{blast} : confess "Provide 'blast'";
  my $verbose = $h{verbose};

    # getting blast db file name
  foreach my $run (keys %$spacer_r){
    next unless exists $spacer_r->{$run}{'BlastOutput_iterations'};  # must have hit

    # getting blastdbfile
    my $blastdbfile = exists $spacer_r->{$run}{'BlastOutput_db'} ?
      $spacer_r->{$run}{'BlastOutput_db'} :
        confess "Cannot find BlastOutput_db in run $run";

    # status
    my @parts = File::Spec->splitpath($blastdbfile);
    print STDERR "Aligning protospacers from blast db: '$parts[2]'\n"
      unless $verbose;

    # each iteration
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){

      # skipping iterations without hits or crRNA_info
      next unless exists $iter->{Iteration_hits} and
        $iter->{Iteration_hits} !~ /^\s*$/;
      next unless exists $iter->{Iteration_hits}{Hit};
      my $crRNA_info = exists $iter->{crRNA_info} ?
	$iter->{crRNA_info} : next;

      # iterating through hits
      foreach my $hit ( @{$iter->{Iteration_hits}{Hit}} ){
        next unless exists $hit->{Hit_hsps}{Hsp};

        # iterating through each hsp
        ## adding protospacer & info to $hsp
        while( my($hspUID, $hsp) = each %{$hit->{Hit_hsps}{Hsp}} ){
	  # checking for needed info
	  next unless exists $hsp->{protoFullXSeq};

	  # aligning for each crDNA
	  foreach my $q (keys %$crRNA_info){
	    confess "Cannot find 'crDNA'" unless exists $crRNA_info->{$q}{crDNA};
	    
	    # call clustalw
	    my $fasta_aln = alignSeqs( 'crDNA', $crRNA_info->{$q}{crDNA}, 
		       'proto', $hsp->{protoFullXSeq});
	    
	    # set protospacer extensions to lower case 
	    setProtoXLower($fasta_aln, $hsp);

	    # set DR portions of crDNA to lower case
	    setcrDNALower($fasta_aln, $crRNA_info->{$q});

	    # adding alignment to hsp
	    $hsp->{$q} = $fasta_aln; 
	  }
	}
      }
    }
  }
}


=head2 setcrDNALower

Setting DR portions of crDNA to lower case

=head3 IN

$asta_aln :  \% of alignment: $fasta_aln->{crDNA} must exist
$q  :  \%; CLdb_info entry

=head3 OUT

=cut

sub setcrDNALower{
  my ($fasta_aln, $q) = @_;

  # IO check
  map{ confess "Cannot find '$_'" unless exists $q->{$_} }
    qw/spacer_start spacer_end region_start region_end 
       array_sense_strand crDNA/;
  exists $fasta_aln->{crDNA} || confess "Cannot find 'crRNA' in aligned fasta";

  # determining lenght of extension up & downstream of proto (pos strand)
  my $upX = $q->{spacer_start} - $q->{region_start};
  my $downX = $q->{region_end} - $q->{spacer_end};

  # flipping up & down if neg strand (all to + strand orientation)
  ($upX, $downX) = ($downX, $upX) if $q->{array_sense_strand} == -1;

  # to lower case 
  ## upstream
  $fasta_aln->{crDNA} =~ s/([A-Z])/\L$1/ for(1..$upX);
  ## downstream 
  my $rev = reverse $fasta_aln->{crDNA};
  $rev =~ s/([A-Z])/\L$1/ for(1..$downX);
  $fasta_aln->{crDNA} = reverse $rev;

#  print Dumper $fasta_aln; 
}

=head2 setProtoXLower

Setting protospacer extensions to lower case.
Determining extension by proto start-end info.

=head3 IN

$fasta_aln :  \% of alignment; $fasta_aln->{proto} must exists
$hsp  :   \%hsp data structure

=head3 OUT

=cut

sub setProtoXLower{
  my ($fasta_aln, $hsp) = @_;

  # IO check
  map{ confess "Cannot find '$_'" unless  exists $hsp->{$_} } 
    qw/protoFullStart protoFullEnd protoFullXStart protoFullXEnd
       subjectStrand/;
  exists $fasta_aln->{proto} || confess "Cannot find 'proto' in aligned fasta";
    
  # determining lenght of extension up & downstream of proto (pos strand)
  my $upX = $hsp->{protoFullStart} - $hsp->{protoFullXStart};
  my $downX = $hsp->{protoFullXEnd} - $hsp->{protoFullEnd};

  # flipping up & down if neg strand (all to + strand orientation)
  ($upX, $downX) = ($downX, $upX) if $hsp->{subjectStrand} == -1;

  # to lower case 
  ## upstream
  $fasta_aln->{proto} =~ s/([A-Z])/\L$1/ for(1..$upX);
  ## downstream 
  my $rev = reverse $fasta_aln->{proto};
  $rev =~ s/([A-Z])/\L$1/ for(1..$downX);
  $fasta_aln->{proto} = reverse $rev;
  
}


=head2 alignSeqs

Aligning sequences using bioperl implementation of clustalw

=head3 IN

$name1 :  Sequence1 name
$seq1 :   Sequence1
$name2 :  Sequence2 name
$seq2 :   Sequence2

=head3 OUT

\% :  fasta of alignment

=cut 

sub alignSeqs{
  my ($name1, $seq1, $name2, $seq2) = @_;

  map{ s/^>// } ($name1, $name2);

  # making sequence string
  my $seq_str = join("\n", ">$name1", $seq1,
		     ">$name2", $seq2);
  # making seqIO object
  my $strfh = IO::String->new($seq_str);
  my @seqIOs;
  my $seqo = Bio::SeqIO->new( -fh => $strfh, -format => 'fasta' );
  while(my $seq = $seqo->next_seq()) { push @seqIOs, $seq; }

  # alignment
  my $factory = Bio::Tools::Run::Alignment::Clustalw->new(quiet => 1);
  my $aln = $factory->align(\@seqIOs);  
  
  ## changing back sequence names
  my @seq_names = ($name1, $name2);
  my %fasta;
  foreach my $seq ( $aln->each_seq() ){
    $fasta{ $seq->display_id} = $seq->seq();
  }

#  print Dumper %fasta; exit;
  return \%fasta;
}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::Align


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
