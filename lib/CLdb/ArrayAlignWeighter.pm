package ArrayAlignWeighter;

use strict;
use warnings;
use Data::Dumper;

use File::Basename;
use lib dirname (__FILE__);

use base 'ArrayAlign';

#use Array::Align;


=head2 weighter

Setting weight for aligner.
See https://github.com/jkahn/Array-Align/blob/master/lib/Array/Align.pm
for more info

=cut

sub weighter {
  my ($self, $left, $right) = @_;
  return 1 if not defined $left;
  return 1 if not defined $right;
  return 1.5 if $left ne $right;
  return 0;
}


=head1 Example

my @left = qw/111 4453 324 12 1 1234 999 99 9 222 22 2/;
my @right = qw/324 12 1 1234 99 9 222 22 2/;


my $aligner = Aligner->new(left => \@left, right => \@right);

for my $pair ($aligner->pairwise) {
  map{ $pair->[$_] = '' unless defined $pair->[$_]} 0..1;
  printf "%s - %s\n", @$pair;
}

=cut 


1;

