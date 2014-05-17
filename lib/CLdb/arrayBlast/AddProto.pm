package CLdb::arrayBlast::AddProto;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use List::Util qw/max/;
use IPC::Cmd qw/can_run run/;
use Parallel::ForkManager;

# CLdb
use CLdb::seq qw/ read_fasta /;

# export #
use base 'Exporter';
our @EXPORT_OK = ();
	
=head1 NAME

CLdb::arrayBlast::AddProto

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for editing sequence data

=cut


=head2 queryBlastDBs

Querying blast DBs for spacer hits.
Using blastdbcmd to query and retrieve proteospacer
sequence.

=head3 IN

hash of args:
blast =>  decoded blast srl
extension => extension beyond protospacer

=head3 OUT

=cut

push @EXPORT_OK, 'queryBlastDBs';

sub queryBlastDBs{
  my %h = @_;
  my $spacer_r = exists $h{blast} ? $h{blast} : 
    confess "Provide a decoded blast srl as $%";
  my $ext = exists $h{extension} ? $h{extension} :
    confess "Provide the protospacer extension (bp)";
  my $fork = exists $h{fork} ? $h{fork} : confess "Provide 'fork'";
  my $verbose = $h{verbose};

  # checking for blastdbcmd in path
  can_run('blastdbcmd') or confess "Cannot find 'blastdbcmd' in \$PATH";

  # parallel queries by hit
  my $pm = new Parallel::ForkManager($fork);

  # getting blast db file name
  foreach my $run (keys %$spacer_r){
    next unless exists $spacer_r->{$run}{'BlastOutput_iterations'};  # must have hit

    # getting blastdbfile
    my $blastdbfile = exists $spacer_r->{$run}{'BlastOutput_db'} ?
      $spacer_r->{$run}{'BlastOutput_db'} :
        confess "Cannot find BlastOutput_db in run $run";

    # status
    my @parts = File::Spec->splitpath($blastdbfile);
    print STDERR "Retrieving sequences from blast db: '$parts[2]'\n"
      unless $verbose;

    # each iteration
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      
      # skipping iterations without hits
      next unless exists $iter->{Iteration_hits} and
        $iter->{Iteration_hits} !~ /^\s*$/;
      next unless exists $iter->{Iteration_hits}{Hit};

      # getting crRNA info if provided
      # my $crRNA_info = $iter->{crRNA_info} if exists $iter->{crRNA_info};

      # getting query length
      my $query_len = exists $iter->{'Iteration_query-len'} ?
	$iter->{'Iteration_query-len'} :
	  confess "Cannot find Iteration_query-len";

      # iterating through hits
      foreach my $hit ( @{$iter->{Iteration_hits}{Hit}} ){
	next unless exists $hit->{Hit_hsps}{Hsp};

	# collecting updated Hsps 
	my @Hsps;
	$pm->run_on_finish(
	   sub{ 
	     my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $hsp) = @_;
	     push @Hsps, $hsp;			     
	   });


	# getting subject scaffold
	my $sub_scaf = exists $hit->{Hit_id} ? 
	  $hit->{Hit_id} : confess "Cannot find 'Hit_id'";
	my $sub_scaf_len = exists $hit->{Hit_len} ? 
	  $hit->{Hit_len} : confess "Cannot find 'Hit_len'";
	
	# iterating through each hsp
	## adding protospacer & info to $hsp
	foreach my $hsp ( @{$hit->{Hit_hsps}{Hsp}} ){
	  $pm->start and next;

	  addProtoCoords( hsp => $hsp, 
			  query_len => $query_len, 
			  extension => $ext,
			  subject_scaffold => $sub_scaf,
			  scaffold_len => $sub_scaf_len);
	  getProtoFromDB( hsp => $hsp,
			  blastdb => $blastdbfile );	 
	  #print Dumper $hsp; exit;
	  $pm->finish(0, $hsp);
	}	
	$pm->wait_all_children;

	# replacing old Hsps with udpates
	$hit->{Hit_hsps}{Hsp} = \@Hsps;
      }      
    }
    # time testing: exit after 1st run
    exit; 
  }
}


=head2 addProtoCoords

Adding full-length protospacer start-end
and strand.

To get full length protospacer, extending
back from spacer query.

Hsp_qseq = query hit
Hsp_hit-from = subject hit start
Hsp_hit-to = subject hit end


=head3 IN

hash of args:
hsp :  hsp $%
query_len :  query length
extension :  extension length
subject_scaffold :  subject scaffold name

=head3 OUT

=cut

sub addProtoCoords{
  my %h = @_;
  my $hsp = exists $h{hsp} ? $h{hsp} : confess "Provide an hsp";
  my $query_len = exists $h{query_len} ?
    $h{query_len} : confess "Provide the query length";
  my $ext = exists $h{extension} ? 
    $h{extension} : confess "Provide proteo spacer extension";
  my $sub_scaf = exists $h{subject_scaffold} ?
    $h{subject_scaffold} : confess "Provide subject scaffold ID";
  my $sub_scaf_len = exists $h{scaffold_len} ? 
    $h{scaffold_len} : confess "Provide 'scaffold_len'";

  # determining strand of subject sequence
  my $sub_start = exists $hsp->{'Hsp_hit-from'} ?
    $hsp->{'Hsp_hit-from'} : confess $!;
  my $sub_end = exists $hsp->{'Hsp_hit-to'} ?
    $hsp->{'Hsp_hit-to'} : confess $!;

  # strand of subject hit; moving subject start-end to + strand
  my $sub_strand = $sub_start <= $sub_end ? 1 : -1;
  ($sub_start, $sub_end) = ($sub_end, $sub_start)
    if $sub_strand == -1;

  # determing length not covered by query hit
  ## perspective of query
  my $q_start = exists $hsp->{'Hsp_query-from'} ?
    $hsp->{'Hsp_query-from'} : confess $!;
  my $q_end = exists $hsp->{'Hsp_query-to'} ?
    $hsp->{'Hsp_query-to'} : confess $!;
 
  my $query_up_miss = $q_start - 1;  # upstream from query
  my $query_down_miss = $query_len - $q_end;  # downstream from query

  # flipping up & down miss if strand hit is '-'
  ## still will use subject start-end from + strand, but the extension length matters
  ($query_up_miss, $query_down_miss) = ($query_down_miss, $query_up_miss)
    if $sub_strand == -1;

  # extending subject hit start-end to full length
  my $protoFullStart = $sub_start - $query_up_miss;
  my $protoFullEnd = $sub_end + $query_down_miss;

  # extending subject hit start-end to '-extension' length
  my $protoFullXStart = $protoFullStart - $ext;
  my $protoFullXEnd = $protoFullEnd + $ext;

  # floor to start (1-indexed)
  map{ $_ = 1 if $_ < 1 } ($protoFullStart, $protoFullXStart);

  # ceiling to end (1-indexed); ceiling = scaffold length
  map{ $_ = $sub_scaf_len if $_ > $sub_scaf_len } ($protoFullEnd, $protoFullXEnd);

  # adding info to hsp
  ## will use location info to query blast DB for sequence info    
  $hsp->{protoFullStart} = $protoFullStart;
  $hsp->{protoFullEnd} = $protoFullEnd;
  $hsp->{protoFullXStart} = $protoFullXStart;
  $hsp->{protoFullXEnd} = $protoFullXEnd;
  $hsp->{subjectStrand} = $sub_strand;
  $hsp->{subjectScaffold} = $sub_scaf;
  $hsp->{protoX} = $ext;

#  print Dumper $query_len;
#  print Dumper $hsp; exit;
}


=head2 getProtoFromDB

Getting protospacer sequence (full length)
and sequence of proto+extensions. 

Using blastdbcmd to query blast db

=head3 IN

hash of args:
hsp :  hash $%
blastdb :  blast database file

=head3 OUT

=cut

sub getProtoFromDB{
  my %h = @_;
  my $hsp = exists $h{hsp} ? 
    $h{hsp} : confess "Provide 'hsp'";
  my $blastdbfile = exists $h{blastdb} ?
    $h{blastdb} : confess "Provide 'blastdb'";

  # 'entry' = subject scaffold name
  my $entry = exists $hsp->{subjectScaffold} ?
    $hsp->{subjectScaffold} : confess "Cannot find subjectScaffold";

  # strand
  my $strand = exists $hsp->{subjectStrand} ?
    $hsp->{subjectStrand} : confess "Cannot find subjectStrand";
  if($strand == 1){ $strand = 'plus'; }
  elsif($strand == -1){ $strand = 'minus'; }
  else{ confess "strand '$strand' should be 1 or -1"; }

  # range 
  ## proto
  my $proto_range = join("-", $hsp->{protoFullStart}, $hsp->{protoFullEnd});
  ## proto & extension
  my $protoX_range = join("-", $hsp->{protoFullXStart}, $hsp->{protoFullXEnd});

  ## sub for querying for full proto & extension
  sub call_blastdbcmd{
    my ($blastdbfile, $entry, $proto_range, $strand) = @_;
    my $cmd = "blastdbcmd -db $blastdbfile -entry '$entry' -range '$proto_range' -strand $strand";

    ## running cmd
    open PIPE, "$cmd |" or confess $!;
    my $fasta_r = read_fasta(fh => \*PIPE);
    close PIPE;
    confess "No sequence retrieved from command: '$cmd'" unless %$fasta_r;

    ## returning just 1 sequence
    foreach my $name (keys %$fasta_r){
      return $fasta_r->{$name};
      last;
    }
  }

  # adding sequence to hsp
  $hsp->{protoFullSeq} = call_blastdbcmd($blastdbfile, $entry, $proto_range, $strand);
  $hsp->{protoFullXSeq} = call_blastdbcmd($blastdbfile, $entry, $protoX_range, $strand);
  
  #print Dumper $hsp; exit;
}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::AddProto


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
