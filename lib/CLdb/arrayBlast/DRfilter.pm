package CLdb::arrayBlast::DRfilter;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use Sereal qw/ decode_sereal /;
use Set::IntervalTree;

# export #
use base 'Exporter';
our @EXPORT_OK = '';

	

=head1 NAME

CLdb::arrayBlast::DRfilter

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

=head1 EXPORT_OK

=cut


=head2 decode_file

Sereal file decoding

=head3 IN

Sereal file name

=head3 OUT

decoded data structure

=cut

push @EXPORT_OK, 'decode_file';

sub decode_file{
  my $file = shift;
  open IN, $file or die $!;
  my $str = '';
  $str .= $_ while <IN>;
  close IN;
  return decode_sereal($str);
}

=head2 make_DR_itree

Making iterval tree of all DR hits (start-end)

For each hsp getting: subject (db,sseqid), start, end, strand

=head3 IN

DR hashref

hashref of options

=head3 OUT

=cut

push @EXPORT_OK, 'make_DR_itree';

sub make_DR_itree{
  my $DR_r = $_[0];
  my $range = $_[1]->{range} || croak $!;
  my $len_cut = $_[1]->{len_cut} || croak $!;
  my $evalue_cut = $_[1]->{evalue_cut} || croak $!;

  # status;
  my %filter;
  my %itrees;

  # iterating through each hit
  foreach my $run (keys %$DR_r){
    # database 
    my $db =  $DR_r->{$run}{'BlastOutput_db'};
    $db = (File::Spec->splitpath($db))[2];

    foreach my $iter ( @{$DR_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
            
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
	  my $sstart = $hsp->{'Hsp_hit-from'};
	  my $send = $hsp->{'Hsp_hit-to'};
	  my $evalue = $hsp->{'Hsp_evalue'};
	  my $aln_len = $hsp->{'Hsp_align-len'};

	  # strand
	  my $strand = $sstart <= $send ? 1 : -1;	 
	  ($sstart, $send) = ($send, $sstart) if 
	    $sstart > $send; # start must be <= end

	  # filtering DR: not loading DR if poor hit
	  $filter{'total'}++;
	  if (defined $evalue_cut and  $evalue > $evalue_cut){
	    $filter{'evalue'}++;
	    next;
	  }
	  if (defined $len_cut and $aln_len < $len_cut){
	    $filter{'length'}++;
	    next;
	  }
	  
	  # itree
	  ## initializing itree
	  unless( exists $itrees{$subj}{$strand} ){
	    $itrees{$subj}{$strand} = Set::IntervalTree->new();
	  }
	  
	  ## load itree
	  $sstart = $sstart - $range -1 < 1 ? 1 : $sstart - $range;
	  $send = $send + $range +2;
	  $itrees{$subj}{$strand}->insert( 1, $sstart, $send );
	  $filter{'used'}++;
	}
      }
    }
  }
  
  # status
  map{ $filter{$_} = 0 unless exists $filter{$_}} qw/total evalue length used/;
  print STDERR "### filter of 'bad' DR hits ###\n";
  print STDERR "Total DR hits: $filter{'total'}\n";
  print STDERR "DR hits with evalues < $evalue_cut: $filter{'evalue'}\n" if $evalue_cut;
  print STDERR "DR hits with hit lengths (fraction of total) < $len_cut: $filter{'length'}\n";
  print STDERR "DR hits used for filtering: $filter{'used'}\n";

  return \%itrees;
}


=head2 DR_filter_blast

filtering spacer hits based on adjacent DR hits 

=head3 IN

spacer_r = blast xml hash ref
itees_r = itrees of DR hits

=head3 OUT

just editing spacer_r 

=cut 

push @EXPORT_OK, 'DR_filter_blast';

sub DR_filter_blast{
  my $spacer_r = shift or croak $!;
  my $itrees_r = shift or croak $!;
  my $DR_cnt = defined $_[0]->{DR_cnt} 
    ? $_[0]->{DR_cnt} : croak "Provide DR count cutoff";

#  print Dumper $spacer_r; exit;

  # status
  my %filter;

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
	  $filter{'total'}++;

	  # hsp values
	  my $sstart = $hsp->{'Hsp_hit-from'};  # hit location on subject
	  my $send = $hsp->{'Hsp_hit-to'};
	  my $evalue = $hsp->{'Hsp_evalue'};
	  my $aln_len = $hsp->{'Hsp_align-len'};

	  # strand
	  my $strand = $sstart <= $send ? 1 : -1;	 
	  ($sstart, $send) = ($send, $sstart) if 
	    $sstart > $send; # start must be <= end

	  # filtering by adjacencies to DR hits (DR itree)
	  if( exists $itrees_r->{$subj}{$strand} ){ # hits to same subject for spacer & DR
	    # number of adjacent DRs
	    ## upstream hits
	    my $up_adj = $itrees_r->{$subj}{$strand}->fetch( $sstart, $sstart + ($send - $sstart)); 
	    ## downstream hits
	    my $down_adj = $itrees_r->{$subj}{$strand}->fetch($sstart + ($send - $sstart) +1, $send);
	    my $DR_adj = 0;
	    $DR_adj += scalar @$up_adj > 0 ? 1 : 0;
	    $DR_adj += scalar @$down_adj > 0 ? 1 : 0;
	   
	    if( $DR_adj >= $DR_cnt){ # hits 'array'
	      $filter{'array'}++;
	      $hsp->{'CLdb_array-hit'} = 1;
	    }
	    else{  # not hitting an 'array'; keeping
	      $filter{'proto'}++;
	      $hsp->{'CLdb_array-hit'} = 0;
	    }
	    
	  }
	  else{ # no hits; keeping 
	    $filter{'proto'}++;
	    $hsp->{'CLdb_array-hit'} = 0;
	  }	  
	}
      }
    }
  }

  # status
  map{ $filter{$_} = 0 unless exists $filter{$_}} qw/total array proto/;
  print STDERR "### DR filtering of spacer blast hits ###\n";
  print STDERR "Total spacer blast hits: $filter{'total'}\n";
  print STDERR "Spacer blast hits hitting an array (filtered out): $filter{'array'}\n";
  print STDERR "Spacer blast hits hitting a protospacer: $filter{'proto'}\n";

}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::DRfilter


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
