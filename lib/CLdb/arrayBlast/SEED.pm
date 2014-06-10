package CLdb::arrayBlast::SEED;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use Set::IntervalTree;
use List::Util qw/sum/;

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
    else{ $fasta{int($pair / 2)}{$tmpkey} .= $_; }
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
# 1-indexing everything!
  my $fasta_r = shift || confess "Provide fasta\n";
  my $seed_index = shift || confess "Provide seed index\n";
  my $no_gap = shift; 
  
  my %mismatch_sum;  # pair => field => value
  my %mismatchByPos;
  foreach my $pair (keys %$fasta_r){    
    # getting mismatches by position
    my ($itree, $seq_len) = makeMismatch_itree($fasta_r->{$pair}, 
					       $seed_index, $no_gap);

    # getting position of the protospacer (& extension)
    my $protoSE_r = getProtoCoords($fasta_r->{$pair});  # proto start-end

    # setting seed coords relative to this protospacer (& non-seed region)
    my ($seedSE_r, $NseedSE_r) = setSEEDcoords($protoSE_r, $seed_index);

    # summing mismatches by position (sep for total, seed, non-seed, and extensions)
    $mismatch_sum{$pair} = sumMismatch($itree, $protoSE_r, $seedSE_r, $NseedSE_r);
    $mismatch_sum{$pair}{seqs} = [sort keys %{$fasta_r->{$pair}}];

    # getting mismatch by position
    mismatchByPos(\%mismatchByPos, $itree, $protoSE_r, $seedSE_r, $NseedSE_r);

  }

#  print Dumper %mismatch_sum; exit;
#  print Dumper %mismatchByPos; exit;
  return \%mismatch_sum, \%mismatchByPos;
}


=head2 mismatchByPos

Summing mismatches with postion set relative to the start of the seed sequence

=cut

sub mismatchByPos{
  my $mismatchByPos_r = shift or confess "Provide mismatch \%\n";
  my $itree = shift or confess "Provide interval tree of mismatches\n";
  my $protoSE_r = shift or confess "Provide the protospacer [start, end]\n";
  my $seedSE_r = shift or confess "Provide the seed [start,end]\n";
  my $NseedSE_r = shift or confess "Provide the non-seed [start,end]\n";

  
  # relative to seed start
  my $seed_start = $seedSE_r->[0];

  
  # proto
  for my $i ($protoSE_r->[0]..$protoSE_r->[1]){
    my $ret_r = $itree->fetch($i-1, $i+1);
    $mismatchByPos_r->{protospacer}{$i - $seed_start}{mismatch} += $ret_r->[0];
    $mismatchByPos_r->{protospacer}{$i - $seed_start}{count}++;
    $mismatchByPos_r->{protospacer}{$i - $seed_start}{proto_pos} = $i;  # 'absolute' position
  }
  
  # seed
  for my $i ($seedSE_r->[0]..$seedSE_r->[1]){
    my $ret_r = $itree->fetch($i-1, $i+1);
    $mismatchByPos_r->{SEED}{$i - $seed_start}{mismatch} += $ret_r->[0];
    $mismatchByPos_r->{SEED}{$i - $seed_start}{count}++;
    $mismatchByPos_r->{SEED}{$i - $seed_start}{proto_pos} = $i;  # 'absolute' position
  }

  # non-seed
  if($NseedSE_r->[0] <= $NseedSE_r->[1]){
    for my $i ($NseedSE_r->[0]..$NseedSE_r->[1]){
      my $ret_r = $itree->fetch($i-1, $i+1);
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{mismatch} += $ret_r->[0];
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{count}++;
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{proto_pos} = $i;  # 'absolute' position
    }
  }
  if($NseedSE_r->[2] <= $NseedSE_r->[3]){
    for my $i ($NseedSE_r->[2]..$NseedSE_r->[3]){
      my $ret_r = $itree->fetch($i-1, $i+1);
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{mismatch} += $ret_r->[0];
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{count}++;
      $mismatchByPos_r->{nonSEED}{$i - $seed_start}{proto_pos} = $i;  # 'absolute' position
    }
  }

    
  #print Dumper %$mismatchByPos_r; exit;
}


=head2 sumMismatch

summing mismatches by postion.
Grouping by seed, non-seed, extensions, & total

=cut

sub sumMismatch{
  my $itree = shift or confess "Provide interval tree of mismatches\n";
  my $protoSE_r = shift or confess "Provide the protospacer [start, end]\n";
  my $seedSE_r = shift or confess "Provide the seed [start,end]\n";
  my $NseedSE_r = shift or confess "Provide the non-seed [start,end]\n";

  my %mismatch;

  # summing across all protospacer
  my $ret_r = $itree->fetch( $protoSE_r->[0] - 1,  # fetch is no inclusive of min-max values
			     $protoSE_r->[1] + 1);
  $mismatch{protospacer}{mismatch} = sum @$ret_r;
  $mismatch{protospacer}{count} = scalar @$ret_r;
  

  # summing across just seed
  $ret_r = $itree->fetch($seedSE_r->[0] - 1,
			 $seedSE_r->[1] + 1);
  $mismatch{SEED}{mismatch} = sum @$ret_r;
  $mismatch{SEED}{count} = scalar @$ret_r;

  
  # summing across non-seed
  $mismatch{nonSEED} = {mismatch =>  0,
			count => 0};
  
  if($NseedSE_r->[0] <= $NseedSE_r->[1]){
    $ret_r = $itree->fetch($NseedSE_r->[0] - 1,
			   $NseedSE_r->[1] + 1);
    $mismatch{nonSEED}{mismatch} += sum @$ret_r;
    $mismatch{nonSEED}{count} += scalar @$ret_r;
  }
  if($NseedSE_r->[2] <= $NseedSE_r->[3]){
    $ret_r = $itree->fetch($NseedSE_r->[2] - 1,
			   $NseedSE_r->[3] + 1);
    $mismatch{nonSEED}{mismatch} += sum @$ret_r;
    $mismatch{nonSEED}{count} += scalar @$ret_r;
  }

#  print Dumper carp, %mismatch; exit;
  return \%mismatch;
}

=head2 setSEEDcoords

setting SEED coords relative to the protospacer SE

=cut

sub setSEEDcoords{
  my $protoSE_r = shift or confess "provide protospacer [start,end]\n";
  my $seed_index = shift or confess "provide seed_index [first,last,stream]\n";

  # if up stream: 'left' side of proto
  my ($seed_start, $seed_end);
  my ($Nseed_up_start, $Nseed_up_end);
  my ($Nseed_down_start, $Nseed_down_end);
  if($seed_index->[2] eq 'up'){   
    $seed_start = $protoSE_r->[0] + $seed_index->[0] - 1;
    $seed_end = $protoSE_r->[0] + $seed_index->[1] - 1;

    $Nseed_up_start = $protoSE_r->[0];
    $Nseed_up_end = $seed_start -1; 
    $Nseed_down_start = $seed_end + 1;
    $Nseed_down_end = $protoSE_r->[1];
  }
  # else: 'right' side of proto
  elsif($seed_index->[2] eq 'down'){
    $seed_start = $protoSE_r->[1] - $seed_index->[1] + 1;
    $seed_end = $protoSE_r->[1] - $seed_index->[0] + 1;    

    $Nseed_up_start = $protoSE_r->[0];
    $Nseed_up_end = $seed_start -1; 
    $Nseed_down_start = $seed_end + 1;
    $Nseed_down_end = $protoSE_r->[1];
  }
  else{ die "Logic error: $!\n"; }


  return [$seed_start, $seed_end], 
    [$Nseed_up_start, $Nseed_up_end,
    $Nseed_down_start, $Nseed_down_end]; # Nseed start can be < end; those should be skipped
}


sub getProtoCoords{
  my $fasta_r = shift or confess "Provide a fasta";
  
  my @start_end;
  foreach my $seqname (keys %$fasta_r){
    next if $seqname =~ /^>*cr[DR]NA\|/;

    my $seq = $fasta_r->{$seqname};
    $seq =~ /^([a-z.-]+)([A-Z][A-Z.-]+[A-Z])/;
    @start_end = ( length($1) , length($1) + length($2) );  
  }

  #print Dumper @start_end;
  \@start_end;
}

=head2 makeMismatch_itree

Loading mixmatches into an interval tree.

=cut

sub makeMismatch_itree{
  my $fasta_r = shift or confess "Provide a fasta";
  my $seed_index = shift or confess "Provde a seed index";
  my $no_gap = shift;

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

  # initialize itree
  my $itree = Set::IntervalTree->new();

  # determining location of mismatchs
  for my $i (0..($seq_len-1)){
    my $mismatch;
    if( $no_gap and ($seqs[0][$i] =~ /[.-]/ or
		     $seqs[1][$i] =~ /[.-]/) ){  # skipping gaps
      $mismatch = 0;
    }
    else{
      $mismatch = ($seqs[0][$i] eq $seqs[1][$i]) ? 0 : 1;  # match : mismatch
    }
    
    $itree->insert($mismatch, $i, $i);
  }
  
  return $itree, $seq_len;
}


=head2 write_sum_table

Writing summary mismatch table

=cut

push @EXPORT_OK, 'write_sum_table';

sub write_sum_table{
  my $prefix = shift or confess "Provide file prefix\n";
  my $MMSum_r = shift or confess "Provide mismatchSum as \%\n";

  my $outfile = "$prefix\_SEEDstats-summary.txt";
  open OUT, ">$outfile" or confess $!;

  # header
  print OUT join("\t", qw/alignment crDNA protospacer region count mismatches mismatch_norm/), "\n";

  # body
  foreach my $pair (keys %$MMSum_r){
    foreach my $region (keys %{$MMSum_r->{$pair}}){
      next if $region eq 'seqs';

      print OUT join("\t", $pair, 
		 @{$MMSum_r->{$pair}{seqs}},
		 $region, 
		 $MMSum_r->{$pair}{$region}{count},
		 $MMSum_r->{$pair}{$region}{mismatch},
		 $MMSum_r->{$pair}{$region}{mismatch} /
		 $MMSum_r->{$pair}{$region}{count}
		 ), "\n";	    
    }	       
  }
  close OUT;

  print STDERR "SEED summary file written: '$outfile'\n";
}


=head2 write_byPos_table

Writing mismatches by position

=cut

push @EXPORT_OK, 'write_byPos_table';

sub write_byPos_table{
  my $prefix = shift or confess "Provide file prefix\n";
  my $MMByPos_r = shift or confess "Provide mismatch ByPos as \%\n";


  my $outfile = "$prefix\_SEEDstats-byPos.txt";
  open OUT, ">$outfile" or confess $!;

  # header
  print OUT join("\t", qw/region pos_rel_SEED pos_rel_align
			  count mismatches mismatch_norm/), "\n";

  # body
  foreach my $region (sort keys %$MMByPos_r){
    foreach my $pos (sort{$a<=>$b} keys %{$MMByPos_r->{$region}}){

      print OUT join("\t", $region, 
		     $pos + 1,   # 1-indexed; relative to the start of the SEED
		     $MMByPos_r->{$region}{$pos}{proto_pos},
		     $MMByPos_r->{$region}{$pos}{count},
		     $MMByPos_r->{$region}{$pos}{mismatch},
		     $MMByPos_r->{$region}{$pos}{mismatch} /
		     $MMByPos_r->{$region}{$pos}{count}
		    ), "\n";	    
    }	       
  }
  close OUT;
 
  print STDERR "SEED by-position table written: '$outfile'\n";
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
