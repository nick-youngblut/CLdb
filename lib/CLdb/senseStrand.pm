package CLdb::senseStrand;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

# export #
use base 'Exporter';
our @EXPORT_OK = '';

	

=head1 NAME

CLdb::senseStrand

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for determining sense strand of array

Sense strand is either determined by array start-end
or loci

=head1 EXPORT_OK

=cut


=head2 setSenseByArraySE

Setting sense_strand by array_start & array_end

=head3 IN

$loci_tbl_r :  @@_ref of loci entries
$index :  %_Ref of field indices

=head3 OUT

%sense : {locus_id} => [start, end, sense_strand]

=cut

push @EXPORT_OK, 'setSenseByArraySE';

sub setSenseByArraySE{
  my $loci_tbl_r = shift || croak "ERROR: provide loci table";
  my $index_r = shift;

  unless(defined $index_r){
    $index_r = {locus_id => 0, array_start => 1, array_end => 2};
  }

  my %sense;
  for my $i (1..$#$loci_tbl_r){  # skipping header
    my $locus_id = $loci_tbl_r->[$i][$index_r->{locus_id}];
    croak "ERROR: cannot find locus_id\n" unless defined $locus_id;
    my $start = $loci_tbl_r->[$i][$index_r->{array_start}];
    my $end = $loci_tbl_r->[$i][$index_r->{array_end}];
    croak "ERROR: cannot find start or end\n"
      unless defined $start and defined $end;
    
    my $sense = $start <= $end ? 1 : -1;
    $sense{ $locus_id } = [$start, $end, $sense];
  }
  
  return \%sense;
}


=head2 getLeaderCoords

Getting start-end for leaders of interest

=head3 IN

dbh :  dbh object
loci_tbl :  @@_ref of loci entries
locus_id_column :  which column has locus_id

=head3 OUT

%leaderSE : {locus_id} => [leader_start, leader_end]

=cut

push @EXPORT_OK, 'getLeaderCoords';

sub getLeaderCoords{
  my %opts = @_;
  my $dbh = $opts{dbh} or croak "ERROR: provide dbh object";
  my $loci_tbl = $opts{loci_tbl} or croak "ERROR: provide a loci table";
  my $id_column = exists $opts{locus_id_column} ?
    $opts{locus_id_column} : 0;

  my $sql = <<HERE; 
SELECT locus_id, leader_start, leader_end
FROM leaders
WHERE locus_id = ?
HERE

  my %leaderCoords;
  my $sth = $dbh->prepare($sql);
  for my $i (1..$#$loci_tbl){  # moving through loci table, skipping header
    my $locus_id = defined $loci_tbl->[$i][$id_column] ?
      $loci_tbl->[$i][$id_column] : croak $!;
    
    $sth->bind_param(1, $locus_id);
    my $rv = $sth->execute;
    my $ret = $sth->fetchall_arrayref( );
    next unless @$ret;

    $leaderCoords{ $ret->[0][0] } = [ @{$ret->[0]}[1..2] ];  # locus_id => [leader_start, leader_end]

  }

  return \%leaderCoords;
}


=head2 senseByLeaderLoc

Setting sense strand based on lociatio of leader sequences

=head3 IN

$sense_r : hashref of locus_id => [start, end, sense]
$leaders_r : hashref of locus_id => [start, end]

=head3 OUT

sense_r edited 

=cut

push @EXPORT_OK, 'senseByLeaderLoc';

sub senseByLeaderLoc{
  my ($sense_r, $leaders_r) = @_;

  foreach my $locus_id (keys %$sense_r){
    # skipping locus unless leader present for locuz
    next unless exists $leaders_r->{$locus_id};

    # making sure start <= end (+ strand)
    ## array
    my $array_start = $sense_r->{$locus_id}->[0];
    my $array_end = $sense_r->{$locus_id}->[1];
    ($array_start, $array_end) = ($array_end, $array_start)
      if $array_start > $array_end;
    ## leader
    my $l_start = $leaders_r->{$locus_id}->[0];
    my $l_end = $leaders_r->{$locus_id}->[1];
    ($l_start, $l_end) = ($l_end, $l_start)
      if $l_start > $l_end;

    # setting sense strand by leader location
    ## if leader upstream (on + strand): sense is + strand
    ## if leader downstream (on - strand): sense is - strand
    if($l_end <= $array_start){
      $sense_r->{$locus_id}->[2] = 1;   # + strand
    }
    elsif($l_start >= $array_end){
      $sense_r->{$locus_id}->[2] = -1;   # - strand
    }
    else{
      warn "WARNING: it appears that the leader & CRISPR array overlap for locus: '$locus_id'\n";
      next;
    }    
  }

  return $sense_r;
}


=head2 parse_sense_file

Parsing sense_file

=head3 IN

$sense_file :  string with file name

=head3 OUT

\%sense : {locus_id} => {strand}

=cut 

push @EXPORT_OK, 'parse_sense_file';

sub parse_sense_file{
  my $sense_file = shift || croak "ERROR: provide a sense file\n";
  
  open IN, $sense_file or die $!;
  my %sense;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    my @l = split /\t/;
    croak "ERROR: sense file has < 2 columns\n"
      unless scalar @l >= 2;

    $sense{$l[0]} = $l[1];
  }
  close IN;
  
  return \%sense;
}


=head2 updateSenseStrand

updating the sense_strand field in loci table

=head3 IN

$dbh :  dbh object
$sense_r :  \%{locus_id} => [start,end,sense]

=head3 OUT

=cut

push @EXPORT_OK, 'updateSenseStrand';

sub updateSenseStrand{
  my $dbh = shift or croak "ERROR: provide a dbh object\n$!\n";
  my $sense_r = shift or croak "ERROR: provide a sense_r ds\n$!\n";

  my $sth = $dbh->prepare("UPDATE loci set array_sense_strand=? WHERE locus_id=?")
    or croak $dbh->err;

  foreach my $locus_id (keys %$sense_r){
    #print Dumper $sense_r->{$locus_id}->[0]; exit;
    $sth->bind_param(1, $sense_r->{$locus_id}->[2]) or croak $dbh->err;		     
    $sth->bind_param(2, $locus_id) or croak $dbh->err;
    $sth->execute or croak $dbh->err;
  }
  
  $dbh->commit() or croak $dbh->err;
}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::senseStrand


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
