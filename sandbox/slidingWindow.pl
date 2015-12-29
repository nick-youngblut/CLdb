#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

slidingWindow.pl -- Sliding window summary of a table

=head1 VERSION

This is version 0.0.1

=head1 USAGE

slidingWindow.pl [options]

=head1 REQUIRED ARGUMENTS

=over

=item <file.txt>

tab-delimited file ('-' if provided via STDIN)

=item -s[eries] <s> 

Column containing numeric values; used to slide the window across.

Default: s.default

=for Euclid:
s.type: int > 0
s.default: 2

=item -v[alue] <v>

Column containing numeric values; used for calculations.

Default: v.default

=for Euclid:
v.type: int > 0
v.default: 2

=back

=head1 OPTIONS

=over

=item -w[indow] <w>

Window size.

Default: w.default

=for Euclid:
w.type: num > 0
w.default: 100

=item -j[ump] <j>

Jump size. 

Default: j.default

=for Euclid:
j.type: num > 0
j.default: 10

=item -g[roup] <g>...

Grouping columns(s) for faceting the sliding windows. (>1 argument allowed).

=for Euclid:
g.type: int > 0

=item -f[ork] <f>

Number of groups to process in parallel (0 = no forking).

Default: f.default

=for Euclid;
f.type: int >= 0
f.default: 0

=item --debug [<log_level>]

Set the log level. Default is log_level.default but if you provide --debug,
then it is log_level.opt_default.

=for Euclid:
    log_level.type:        int
    log_level.default:     0
    log_level.opt_default: 1

=item --quiet

=item --version

=item --usage

=item --help

=item --man

Print the usual program information

=back

=head1 DESCRIPTION

Perform a sliding window analysis on a column
of numeric values in a table.

Each group (if designated) can be processed in parallel
using multiple cores.

=head2 Output columns:

=over

=item group(s)

=item window

=item count

=item sum

=item mean

=item stdev

=back

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut (ndy2@cornell.edu)

=head1 BUGS

There are undoubtedly serious bugs lurking somewhere in this code.
Bug reports and other feedback are most welcome.

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

use Data::Dumper;
use Getopt::Euclid;
use Hash::MultiKey; 
use Statistics::Descriptive;
use Parallel::ForkManager;
#use List::Util qw/min max/;
use List::MoreUtils qw/minmax/;


#--- I/O error ---#
my $fh;
$ARGV{'<file.txt>'} eq '-' ? $fh = \*STDIN :
  open $fh, $ARGV{'<file.txt>'}; 
$ARGV{'-s'}--;
$ARGV{'-v'}--;
map{ $_ -= 1 } @{$ARGV{'g'}} if exists $ARGV{'g'};

#--- MAIN ---#

# loading table
my %tbl;
tie %tbl, 'Hash::MultiKey';
while(<$fh>){
  chomp;
  next if /^\s*$/;
  my @l = split /\t/;

  # I/O check 
  die "ERROR: table must have >= 2 columns (<2 at line $.)\n"
    unless @l >= 2;
  die "ERROR: cannot find series column value (line $.)\n"
    unless defined $l[$ARGV{'-s'}];
  die "ERROR: cannot find value column value (line $.)\n"
    unless defined $l[$ARGV{'-v'}];

  # loading by group
  my $group_r = exists $ARGV{'-g'} ? [@l[@{$ARGV{'-g'}}]] : ['all'];
  push @{$tbl{$group_r}{$l[$ARGV{'-s'}]}},  $l[$ARGV{'-v'}]; 
}

# sliding window on each group
my $pm = Parallel::ForkManager->new($ARGV{'-f'});

## fork finish
$pm->run_on_finish(
  sub{
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $ret) = @_;
    # writing out summary by window
    foreach my $window (sort{$a<=>$b} keys %$ret){
      #print Dumper @{$ret->{$window}}{qw/count sum/}; exit;
      map{ $_ = 'NA' unless defined $_ } @{$ret->{$window}}{qw/count sum mean stdev/};

      print join("\t", @{$ret->{$window}{group}}, $window,
		 @{$ret->{$window}}{qw/count sum mean stdev/}), "\n"
    }
  }
);

## getting min-max of all values
#print Dumper %tbl; exit;
my @pos;
map{ push @pos, keys $_ } @tbl{keys %tbl}; 
my ($min, $max) = minmax(@pos);

## each group
foreach my $group (keys %tbl){
  $pm->start and next;

  my %stats;
  for(my $i=$min; $i<=$max; $i+=$ARGV{'-j'}){ # each window
    # loading values for each window
    my $stat = Statistics::Descriptive::Full->new();
    for my $x (@{$tbl{$group}}{$i..($i+$ARGV{'-w'})} ){
      next unless defined $x;
      map{ $stat->add_data($_) } @$x;
    }
    
    # loading stats hash
    $stats{$i}{count} = $stat->count();
    $stats{$i}{sum} = $stat->sum();
    $stats{$i}{mean} = $stat->mean();
    $stats{$i}{stdev} = $stat->standard_deviation();
    $stats{$i}{group} = $group;
  }

  $pm->finish(0, \%stats);
}
$pm->wait_all_children;
