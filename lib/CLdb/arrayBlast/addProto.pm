package CLdb::arrayBlast::addProto;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use List::Util qw/max/;

# export #
use base 'Exporter';
our @EXPORT_OK = ();
	
=head1 NAME

CLdb::arrayBlast::addProto

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for editing sequence data

=head get_proto_seq

getting protospacer/extension start-stop & strand

=cut

push @EXPORT_OK, 'get_proto_seq';

sub get_proto_seq{
  my ($blast_r, $fields_r, $extend, 
      $len_cut, $pam_r, $verbose) = @_;

  my %filter_sum;
  foreach my $query (sort keys %$blast_r){
    print STDERR "Processing query: '$query'\n" unless $verbose;
    foreach my $db (keys %{$blast_r->{$query}}){
      foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
        # checking existence of start-stop, qlen, slen #
        check_values($row, $fields_r);

        # filtering #
        $filter_sum{"total"}++;
        ## deleting  hits with gaps ##
        if ($$row[$fields_r->{"BTOP"}] =~ /-/){         # gaps in blast alignment
          $filter_sum{"gaps"}++;
          next;
        }
        ## deleting if hit_len/query_len < $len_cut ##
        if( ($$row[$fields_r->{"q. end"}] - $$row[$fields_r->{"q. start"}]) /
            $$row[$fields_r->{"query length"}] < $len_cut ){
          $filter_sum{"length"}++;
          next;
        }

        # getting subject start-stop, strand #
        my $Sstart = $$row[$fields_r->{"s. start"}];
        my $Send = $$row[$fields_r->{"s. end"}];
        my $strand;
        if($Sstart > $Send){ $strand = 'minus'; }
        else{ $strand = 'plus'; }

        # flipping Sstart-Send to + if - strand #
        ($Sstart, $Send) = ($Send, $Sstart) if $strand eq 'minus';

        # calling blastdbcmd #
        ## sequence will be - strand if $strand eq 'minus' #
        my $proto_seq = call_blastdbcmd($$row[$fields_r->{"subject id"}],
                                        $Sstart, $Send, $strand, $db, 1);  # inverted so it is protospacer match

        # extending hit to full length & getting sequence #
        ## sstart-send, qstart-qend full length (as much as extension would allow) ##
        my ($Sstart_full, $Send_full, $Qstart_full, $Qend_full) =
          extend_full_len($Sstart, $Send, $strand,
                          $$row[$fields_r->{"q. start"}],
                          $$row[$fields_r->{"q. end"}],
                          $$row[$fields_r->{"query length"}],
                          $$row[$fields_r->{"subject length"}]);

        #print Dumper $Sstart_full, $Send_full, $Qstart_full, $Qend_full; exit;

        my $proto_seq_full = call_blastdbcmd($$row[$fields_r->{"subject id"}],
                                             $Sstart_full, $Send_full, $strand, $db, 1); # not inverted; used for alignment

        # extending hit beyond full length (-x) #
        my ($Sstart_fullx, $Send_fullx) = extend_beyond_full($Sstart_full, $Send_full, $strand,
                                                             $$row[$fields_r->{"subject length"}], $extend);

        # calling blastdbcmd #
        ## sequence will be - strand if $strand eq 'minus' ##
        my $proto_seq_full_x = call_blastdbcmd($$row[$fields_r->{"subject id"}],
                                               $Sstart_fullx, $Send_fullx, $strand, $db, 1);

        # making extensions lower case #
        $proto_seq_full_x = ext_lower_case($proto_seq_full_x,
                                           $Sstart_fullx, $Send_fullx,
                                           $Sstart_full, $Send_full);


        # extracting the pam regions #
        my ($fiveP_pam, $threeP_pam) = pam_extract($proto_seq_full_x, $pam_r,
                                                   $Sstart_fullx, $Send_fullx,
                                                   $Sstart_full, $Send_full);

        # adding sequences to blast table #
        ## proto sequence extension start-end flip if hits to minus strand #
        ($Sstart_full, $Send_full) = ($Send_full, $Sstart_full) if $strand eq "minus";
        ($Sstart_fullx, $Send_fullx) = ($Send_fullx, $Sstart_fullx) if $strand eq "minus";

        # adding to table #
        push @$row, $proto_seq, $Qstart_full, $Qend_full,
          $Sstart_full, $Send_full, $proto_seq_full,
            $Sstart_fullx, $Send_fullx, $proto_seq_full_x,
              $extend, $fiveP_pam, $threeP_pam;
      }
    }
  }

  # subject_seq = aln to spacer; subject_full_seq = aln to spacer; subject_full_x_seq = revcomp (for PAMs)
  #print Dumper %$blast_r; exit;
  my @new_fields = ("proto_seq",
                    "q_start_full", "q_end_full", "s_start_full", "s_end_full", "proto_seq_full",
                    "s_start_fullx", "s_end_fullx", "proto_seq_fullx", "x_length", "5p_pam_region", "3p_pam_region");

  # status #
  map{$filter_sum{$_} = 0 unless exists $filter_sum{$_}} qw/total gaps length/;
  print STDERR "\n### filtering summary ###\n";
  print STDERR "Total number of hits:\t\t\t\t$filter_sum{'total'}\n";
  print STDERR "Number removed due to gaps in alignment: \t$filter_sum{'gaps'}\n";
  print STDERR "Number removed due to short hit length:\t\t$filter_sum{'length'}\n";
  print STDERR "Number remaining:\t\t\t\t", $filter_sum{'total'} -
    ($filter_sum{'length'} + $filter_sum{'gaps'}), "\n\n";

  return \@new_fields;
}

=head1 check_values

Checking values in blast hit data structure to make sure they exist.

=cut

sub check_values{
  my ($row, $fields_r) = @_;

  my @chk = ("query id", "subject id", "s. start", "s. end",
             "q. start", "q. end",
             "evalue", "query length", "subject length", "BTOP");
  map{die " ERROR: '$_' not found!\n" unless $$row[$fields_r->{$_}]} @chk;
}

=head1 ext_lower_cas

Changing extension sequences to lower case

=cut

sub ext_lower_case{
  # making extensions to protospacer lower case #
  ## substr extensions, and full length; lower case extensions; stitching together ##
  my ($proto_seq_full_x, $Sstart_fullx, $Send_fullx,
      $Sstart_full, $Send_full) = @_;

  # getting extension lengths #
  my $start_x_len = $Sstart_full - $Sstart_fullx;
  my $end_x_len = $Send_fullx - $Send_full;
  my $proto_len = $Send_full - $Sstart_full + 1;                # inclusive

  # substringing extensions and full #
  my $startx_seq = substr($proto_seq_full_x, 0, $start_x_len);
  my $proto_seq = substr($proto_seq_full_x, $start_x_len, $proto_len);
  my $endx_seq = substr($proto_seq_full_x, $start_x_len + $proto_len, $end_x_len);

  # lower case extensions #
  $startx_seq =~ tr/[A-Z]/[a-z]/;
  $endx_seq =~ tr/[A-Z]/[a-z]/;


  #print Dumper $startx_seq, $proto_seq, $endx_seq;
  return join("", $startx_seq, $proto_seq, $endx_seq);
}

=head extend_beyond_full

extending proto? beyond the full length

=cut 

sub extend_beyond_full{
  # extending beyond full length #
  my ($Sstart, $Send, $strand, $Slen, $extend) = @_;

  $Sstart -= $extend;
  $Sstart = 1 if $Sstart < 1;
  $Send += $extend;
  $Send = $Slen if $Send > $Slen;

  return $Sstart, $Send;
}


=head1 extend_full_len

extending spacer hit to full length if needed

=cut

sub extend_full_len{
  my ($Sstart, $Send, $strand, $Qstart, $Qend, $Qlen, $Slen) = @_;

  # missing lengths in spacer blast (if partial hit) #
  my $Q5px_len = $Qstart - 1;                   # extension off 5' end of spacer (+ strand); number of bp added
  my $Q3px_len = $Qlen - $Qend;         # extension off 3' end of spacer (+ strand); number of bp added
  ($Q5px_len, $Q3px_len) = ($Q3px_len, $Q5px_len) if $strand eq "minus";        # applyin the extension for other strand

  # changing subject start-end to full length (as much as allowed) #
  my ($Sstart_full, $Send_full);
  $Sstart_full = $Sstart - $Q5px_len;
  $Sstart_full = 1 if $Sstart_full < 1;                 # border issues
  $Send_full = $Send + $Q3px_len;
  $Send_full = $Slen if $Send_full > $Slen;             # border issues

  # getting actual subject extension (acounting for border issues) #
  my $S5px_len = $Sstart - $Sstart_full;                # number of bp added
  my $S3px_len = $Send_full - $Send;                    # number of bp added

  # making gaps based on difference between intended extend and actual extend #
  #my $5px_gap = "-" x ($Q5px_len - $S5px_len);
  #my $3px_gap = "-" x ($Q3px_len - $S3px_len);
  #map{$_ = "" unless $_} ($5px_gap, $3px_gap);
  #($5px_gap, $3px_gap) = ($3px_gap, $5px_gap) if $strand eq "minus";           # gap applied to the other end

  # getting qstart/qend full (as must as subject extension would allow) #
  ($S5px_len, $S3px_len) = ($S3px_len, $S5px_len) if $strand eq "minus";        # applying the extension for other strand
  my $Qstart_full = $Qstart - $S5px_len;
  my $Qend_full = $Qend + $S3px_len;

  # sanity check #
  die " LOGIC ERROR $!" unless $Qstart_full > 0 && $Qend_full > 0;
  die " LOGIC ERROR $!" unless $Sstart_full > 0 && $Send_full > 0;

  return $Sstart_full, $Send_full, $Qstart_full, $Qend_full; #$5px_gap, $3px_gap;
}


=head1 pam_extract

Extracting pam sequence

=cut

sub pam_extract{
  my ($proto_seq_full_x, $pam_r, $Sstart_fullx, $Send_fullx,
      $Sstart_full, $Send_full) = @_;

  # getting extension lengths #
  my $start_x_len = $Sstart_full - $Sstart_fullx;
  my $end_x_len = $Send_fullx - $Send_full;
  my $proto_len = $Send_full - $Sstart_full + 1;                # inclusive

  ## 5' end ##
  my $FiveP_pam_len = abs($$pam_r[0]) - abs($$pam_r[1]) + 1;
  die " LOGIC ERROR: 5' pam length is < 1!\n" unless $FiveP_pam_len > 0;

  my $FiveP_start = $start_x_len + $$pam_r[0];
  my $FiveP_pam;
  if($FiveP_start < 0){         # not enough extension, returning ""
    $FiveP_pam = "";
  }
  else{
    $FiveP_pam = substr $proto_seq_full_x, $FiveP_start, $FiveP_pam_len;
  }


  ## 3' end #
  my $ThreeP_pam_len = abs($$pam_r[3]) - abs($$pam_r[2]) + 1;
  die " LOGIC ERROR: 3' pam length is < 1!\n" unless $ThreeP_pam_len > 0;

  my $ThreeP_start = $proto_len + $start_x_len + $$pam_r[2] -1;
  my $ThreeP_pam;
  if($ThreeP_start < 0){                # not enough extension, returning ""
    $ThreeP_pam = "";
  }
  else{
    $ThreeP_pam = substr $proto_seq_full_x, $ThreeP_start, $ThreeP_pam_len;
  }
  #print Dumper $proto_len, $ThreeP_start, $ThreeP_pam_len, $ThreeP_pam;

  #print Dumper $FiveP_pam, $ThreeP_pam; exit;
  return $FiveP_pam, $ThreeP_pam;
}


=head1 add_new_fields

Adding new fields to blast output

=cut

push @EXPORT_OK, 'add_new_fields';

sub add_new_fields{
  my ($fields_r, $new_fields_r) = @_;

  my $max_i = max(values %$fields_r);   # starting after last fields index

  for my $i (0..$#$new_fields_r){
    $fields_r->{$$new_fields_r[$i]} = $i+$max_i+1;
  }
  #print Dumper %$fields_r; exit;
}


=head1 call_blastdbcmd

Calling blastdbcmd to extract subject sequence (full protospacer & extensions)

=cut

sub call_blastdbcmd{
  my ($subject_id, $Sstart, $Send, $strand, $db) = @_;

  #($Sstart, $Send) = ($Send, $Sstart) if $strand eq 'minus';

  # calling #
  my $cmd = "blastdbcmd -db $db -entry '$subject_id' -range '$Sstart-$Send' -strand $strand |";
  #print Dumper $cmd; exit;
  open PIPE, $cmd or die $!;
  my $seq;
  while(<PIPE>){
    chomp;
    die " ERROR: >1 sequence returned by blastdbcmd!\n"
      if /^>/ && $seq;
    next if /^>/;
    $seq .= $_;
  }

  #print Dumper $seq; exit;
  return $seq;
}


=head1 write_editted_blast

Writing edited blast table. Making sure to update fields with appended ones.

=cut

push @EXPORT_OK, 'write_editted_blast';

sub write_editted_blast{
  my ($blast_r, $fields_r) = @_;

  foreach my $query (sort keys %$blast_r){
    foreach my $db (sort keys %{$blast_r->{$query}}){
      
      # getting nubmer of hits passing filtering  # # skipping hit if no adding info (did not pass filtering) #
      my $hit_cnt = 0;
      foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
        $hit_cnt++ if $$row[$fields_r->{"proto_seq_fullx"}];
      }
      next if $hit_cnt == 0;

      # modifying & writing header #
      my $header_r;
      foreach my $row (@{$blast_r->{$query}{$db}{"header"}}){
        $row = new_header($row, $fields_r) if $row =~ /^# Field/;
        print $row, "\n";
      }

      # writing hits #
      foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
        next unless $$row[$fields_r->{"proto_seq_fullx"}];
        print join("\t", @$row), "\n";
      }
    }
  }
}

=head1 new_header

New header of output blast hits

=cut

sub new_header{
  my ($row, $fields_r) = @_;
  my $fields = join(", ", sort{$fields_r->{$a}<=>$fields_r->{$b}} keys %$fields_r);
  return join(" ", "# Fields:", $fields);
}



=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::addProto


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
