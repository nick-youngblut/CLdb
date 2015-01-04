package CLdb::arrayBlast::DRfilter;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
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
    
   # print Dumper $DR_r->{$run} unless exists $DR_r->{$run}{'BlastOutput_iterations'}{'Iteration'};

    next unless exists $DR_r->{$run}{'BlastOutput_iterations'}{'Iteration'}; 
    foreach my $iter ( @{$DR_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      next unless exists $iter->{'Iteration_hits'} &&
	$iter->{'Iteration_hits'} !~ /^\s*$/;  # next unless hits to subject
           
      # iterating through hits
      ## hits in iteration (array_ref or ref to hash)
      my $hits_ref = $iter->{'Iteration_hits'}{'Hit'}; 
      foreach my $hit ( @$hits_ref){
	# subjectID
	my $sseqid = $hit->{Hit_id};
	my $subj = join("__", $db, $sseqid);

	# iterating through hsp
	my $hsp_ref = $hit->{Hit_hsps}{Hsp};	
	while( my($hspUID, $hsp) = each %$hsp_ref){
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
  print STDERR "### Summary for filtering of 'bad' DR hits ###\n";
  print STDERR "Total DR hits: $filter{'total'}\n";
  print STDERR "DR hits with E-values < $evalue_cut: $filter{'evalue'}\n" if $evalue_cut;
  print STDERR "DR hits with hit lengths (fraction of total) < $len_cut: $filter{'length'}\n";
  print STDERR "DR hits used for filtering: $filter{'used'}\n";

  return \%itrees;
}


=head2 make_CRISPR_itrees

Making iterval trees of all CRISPRs

=head3 IN

hashref -- {fasta_file : scaffold : locus_id : featID : featVal}

=head3 OUT

hashref -- {fasta_file : scaffold : itree}

=cut

push @EXPORT_OK, 'make_CRISPR_itrees';

sub make_CRISPR_itrees{
  my $CRISPR_se_r = shift or confess "Provide CRISPR start-end hashref\n";

  # status;
  my %itrees;

  foreach my $fasta_file (keys %$CRISPR_se_r){
    foreach my $scaf_ID (keys %{$CRISPR_se_r->{$fasta_file}}){
      # initilize itree
      $itrees{$fasta_file}{$scaf_ID} = Set::IntervalTree->new()
	unless exists $itrees{$fasta_file}{$scaf_ID};
    
      # iterating by locusID
      while (my ($locus_ID,$feats_r) = each %{$CRISPR_se_r->{$fasta_file}{$scaf_ID}}){

	# assertions
	my $start = exists $feats_r->{Array_Start} ?
	  $feats_r->{Array_Start} :
	    confess "KeyError: Array_Start\n";
	my $end = exists $feats_r->{Array_End} ?
	  $feats_r->{Array_End} : 
	    confess "KeyError: Array_End\n";
        next if $start !~ /^\d+$/ or $end !~ /^\d+$/;

	# making sure start-end is for + strand and start != end
	if ($start > $end){
	  ($start,$end) = ($end,$start);
	}
	elsif ($start == $end){
	  next;
	}
	elsif ($start > $end){ 
	  continue;
	}
		
	# setting itree ranges
	$itrees{$fasta_file}{$scaf_ID}->insert($locus_ID, $start, $end);	          
      }
    }
  }

  #print Dumper %itrees; exit;
  return \%itrees;
}


=head2 DR_filter_blast

filtering spacer hits based on adjacent DR hits (and CRISPR array locations)

=head3 IN

spacer_r -- blast xml hash ref
itees_r -- itrees of DR hits
keep -- bool; 
CRISPR -- intervalTree object on CLdb CRISPR array locations

=head3 OUT

edited spacer_r object

=cut 

push @EXPORT_OK, 'DR_filter_blast';

sub DR_filter_blast{
  my $spacer_r = shift or confess $!;
  my $itrees_r = shift or confess $!;
  my %h = @_;

  # assertions & optional args
  my $DR_cnt = exists $h{DR_cnt} ?
    $h{DR_cnt} : croak "Provide DR count cutoff\n";
  my $keep = exists $h{keep} ?
    $h{keep} : croak "Provide keep\n";
  my $CRISPR_itrees = (exists $h{CRISPR} and defined $h{CRISPR}) ?
    $h{CRISPR} : undef;

  # status
  my %filter = (CRISPR => 0,
		proto => 0,
		array => 0);


  # iterating through each hit
  foreach my $run (keys %$spacer_r){
    # database 
    my $db =  $spacer_r->{$run}{'BlastOutput_db'};
    $db = (File::Spec->splitpath($db))[2];

    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      next unless exists $iter->{'Iteration_hits'} and 
	ref $iter->{'Iteration_hits'} eq 'HASH';
           
      # iterating through spacer hits
      # hits in iteration (array_ref or ref to hash)
      my $hits_ref = $iter->{'Iteration_hits'}{'Hit'}; 
      foreach my $hit ( @$hits_ref){
	# subjectID
	my $sseqid = $hit->{Hit_id};
	my $subj = join("__", $db, $sseqid);

	# iterating through hsp
	my $hsp_ref = $hit->{Hit_hsps}{Hsp};	
	while( my ($hspUID, $hsp) = each %$hsp_ref ){
	  $filter{'total'}++;

	  # hsp values of spacer hits
	  my $sstart = $hsp->{'Hsp_hit-from'};  # hit location on subject
	  my $send = $hsp->{'Hsp_hit-to'};
	  my $evalue = $hsp->{'Hsp_evalue'};
	  my $aln_len = $hsp->{'Hsp_align-len'};

	  # strand
	  my $strand = $sstart <= $send ? 1 : -1;	 
	  ($sstart, $send) = ($send, $sstart) if 
	    $sstart > $send; # start must be <= end

	  # filtering any hits to known CRISPR arrays (ID by CLdb)
	  ## $db must == genome fasta file name for array
	  ## using start-stop for array
	  if( defined $CRISPR_itrees and exists $CRISPR_itrees->{$db}{$sseqid}){
	    my $ret_r = $CRISPR_itrees->{$db}{$sseqid}->fetch($sstart, $send);
	    if( @$ret_r ){
	      $filter{CRISPR}++;
	      $hsp->{'CLdb_array-hit'} = 1;
	      delete $hsp_ref->{$hspUID};
	      next;
	    }
	  }

	  # filtering by adjacencies to DR hits (DR itree)
	  if( exists $itrees_r->{$subj}{$strand} ){ # hits to same subject for spacer & DR
	    # number of adjacent DRs
	    ## upstream hits
#	    my $up_adj = $itrees_r->{$subj}{$strand}->fetch( $sstart - ,
#							     $sstart + ($send - $sstart)); 
	    ## downstream hits
#	    my $down_adj = $itrees_r->{$subj}{$strand}->fetch($sstart + ($send - $sstart) +1, 
#							      $send);
	    #my $DR_adj = 0;
	    #$DR_adj += scalar @$up_adj > 0 ? 1 : 0;
	    #$DR_adj += scalar @$down_adj > 0 ? 1 : 0;


	    # hits that overlap with DR-hits (DR-hits extened by '-length')
	    my $DR_adj = $itrees_r->{$subj}{$strand}->fetch($sstart,$send);
	    $DR_adj = scalar @$DR_adj;
	   
	    if( $DR_adj >= $DR_cnt){   # hits CRISPR array
	      $filter{array}++;  # number filtered
	      if($keep){        #keeping hit, just marking as hit to array
		$hsp->{'CLdb_array-hit'} = 1;
	      }
	      else{  # deleting hsp
		delete $hsp_ref->{$hspUID};
		next;
	      }
	    }
	    else{    # not hitting an 'array'; keeping
	      $filter{proto}++;
	      $hsp->{'CLdb_array-hit'} = 0;
	    }
	    
	  }
	  else{    # no hits; keeping 
	    $filter{proto}++;
	    $hsp->{'CLdb_array-hit'} = 0;
	  }
	}
      }
    }
  }

  # status
  map{ $filter{$_} = 0 unless exists $filter{$_}} qw/total array proto/;
  print STDERR "### DR filtering of spacer blast hits ###\n";
  printf STDERR "Total spacer blast hits: %i\n", $filter{total};
  printf STDERR "Spacer blast hits hitting CLdb-defined array: %i\n",
    $filter{CRISPR} if $CRISPR_itrees;
  printf STDERR "Spacer blast hits hitting an array defined just by DR-blast hits: %i\n",
    $filter{array};
  $keep ? print STDERR " NOTE: Keeping the blast hits to arrays and writing to the srl file\n" :
    print STDERR " NOTE: Deleting the blast hits to arrays\n";
  printf STDERR "Spacer blast hits hitting a protospacer: %i\n", $filter{proto};
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
