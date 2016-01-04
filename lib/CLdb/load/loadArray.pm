package CLdb::load::loadArray;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;

## export #
use base 'Exporter';
our @EXPORT_OK = '';



=head1 NAME

CLdb::load::loadArray

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing array files (CRISPRFinder format)

=head1 EXPORT_OK

=cut

=head1 SUBROUTINES


=head2 make_headers

indexing headers for CLdb entry loading

=head3 IN

=head3 OUT

hash_hash_refs of indexes

=cut 

push @EXPORT_OK, 'make_headers';

sub make_headers{
# headers for CLdb entry loading #
  my %dr_header = (
                   locus_id => 1,
                   dr_id => 2,
                   dr_start => 3,
                   dr_end => 4,
                   dr_sequence => 5
                  );
  my %sp_header = (
                   locus_id => 1,
                   spacer_id => 2,
                   spacer_start => 3,
                   spacer_end => 4,
                   spacer_sequence => 5
                  );

  #print Dumper %dr_header; exit;
  return \%dr_header, \%sp_header;
}


=head2 parse_array_file

Parsing a array file (CRISPRFinder format)

=head3 IN

path to array (string)
array file (string)
locus id (string)
arrays (hash ref)

=head3 OUT

=cut

push @EXPORT_OK, 'parse_array_file';

sub parse_array_file{
  my ($array_path, $array_file, $locus_id, $arrays_r) = @_;

  die " ERROR: '$array_path/$array_file' not found!\n"
    unless -e "$array_path/$array_file";
  open IN, "$array_path/$array_file" or die $!;

  my $seq;
  my $repeat_cnt = 0;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    s/^\s+//;
    my @line = split /\t/;
    map{$_ =~ s/\s+//g} @line;

    # sanity check #
    die " ERROR: table not formatted correctly (>=4 columns required)\n"
      if scalar(@line) < 4;
    die " ERROR: $array_file not formatted correctely!\n"
      unless $line[0] =~ /^\d+$/;
    die " ERROR: $array_file not formatted correctely!\n"
      unless $line[1] =~ /^[A-Z-]+$/;

    # repeat count #
    $repeat_cnt++;

    # getting positions from array file #
    my $pos_r = determine_positions(\@line, 1, $seq);

    # loading hash #
    my $uID = join("_", $locus_id, "dri.$repeat_cnt");
    $arrays_r->{"DR"}{$uID}{"locus_id"} = $locus_id;
    $arrays_r->{"DR"}{$uID}{"dr_id"} = "$repeat_cnt";
    $arrays_r->{"DR"}{$uID}{"dr_start"} = $$pos_r[0];
    $arrays_r->{"DR"}{$uID}{"dr_end"} = $$pos_r[1];
    $arrays_r->{"DR"}{$uID}{"dr_sequence"} = $line[1];

    # loading spacer #
    if($$pos_r[2] && $$pos_r[3]){
      my $uID = join("_", $locus_id, "$repeat_cnt");
      $arrays_r->{"spacer"}{$uID}{"locus_id"} = $locus_id;
      $arrays_r->{"spacer"}{$uID}{"spacer_id"} = "$repeat_cnt";
      $arrays_r->{"spacer"}{$uID}{"spacer_start"} = $$pos_r[2];
      $arrays_r->{"spacer"}{$uID}{"spacer_end"} = $$pos_r[3];
      $arrays_r->{"spacer"}{$uID}{"spacer_sequence"} = $line[2];
    }

    # concat array sequence #
    $seq .= join("", @line[1..2]);       # total array sequence
    $seq =~ s/ +//g;
  }

  close IN;

  #print Dumper %$arrays_r; exit;
}


=head2 determine_positions

determining the start-end for each spacer & DR from a row in an array file

=head3 IN

$line_r :  array_ref for array_file row
$location : 
$seq :  sequence string

=head3 OUT

array_ref of position info

=cut

sub determine_positions{
# determing start-end for each spacer & DR from array file (1 row)
  my ($line_r, $location, $seq) = @_;

  my @pos;
  my $seq_length = length($seq);
  $seq_length = 0 if ! $seq_length;

  $$line_r[0] = $seq_length + 1 unless $location;  # start = sequence length

  $pos[0] = $$line_r[0] if $$line_r[1];                          
  
  # repeat start
  $pos[1] = $$line_r[0] + length($$line_r[1]) if $$line_r[1];      
  $pos[2] = $$line_r[0] + length($$line_r[1]) if $$line_r[2];   
  $pos[3] = $$line_r[0] + length($$line_r[1]) + length($$line_r[2]) 
    if $$line_r[2];             # spacer end

  die "HERE" unless $$line_r[0];

  return \@pos;
}





=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::loadArray

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
