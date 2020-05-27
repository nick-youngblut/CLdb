#!/usr/bin/env perl
use strict;
use warnings;

=pod

=head1 NAME

arrayBlast -- blast of CRISPR array sequences

=head1 SYNOPSIS

arrayBlast [options] subcommand -- [subcommand_options]

=head2 Options

=over

=item --list

List all subcommands.

=item --perldoc

Get perldoc of subcommand.

=item -v        Verbose output

=item -h        This help message

=back

=head2 For more information:

CLdb --perldoc -- arrayBlast

=head1 DESCRIPTION

This a the top-level commands for running all other commands
associated with array blast-ing.
See --list for all subcommands.

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut



use Data::Dumper;
use Getopt::Long;
use Pod::Usage;
use FindBin;
use lib "$FindBin::RealBin/../lib";
use IPC::Cmd qw/run can_run/;
use File::Spec;


#--- option parsing ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));
sub Pod2Usage{ pod2usage({@_, -exitval => 0}) }
my ($verbose, $listSubcmds, $getPerlDoc, $database);
GetOptions (
	    "--list" => \$listSubcmds,
	    "--perldoc" => \$getPerlDoc,
	    "--verbose" => \$verbose,
	    "database=s" => \$database,  # 
	    "-help|?" => \&Pod2Usage
	   );


#--- I/O error ---#
my $bindir = File::Spec->catdir($FindBin::RealBin, 'bin_arrayBlast');

# check that script exists 
my $scripts_r = list_scripts($bindir, $listSubcmds);
(my $subcmd = $ARGV[0]) =~ s/\.pl$//i;
if (! grep(/^$subcmd(\.pl)*$/, @$scripts_r)){
  die "ERROR: '" . $ARGV[0] . "' is not a a valid subcommand\n";
}


#--- MAIN ---#
call_subcommand($bindir, \@ARGV, $database, $getPerlDoc);


#--- subroutines ---#
sub call_subcommand{
  my $bindir = shift or die $!;
  my $argv_r = shift or die $!;
  my $database = shift;
  my $getPerlDoc = shift;
  
  # check exists
  (my $subcmd = $argv_r->[0]) =~ s/(\.pl)*$/.pl/i;
  $subcmd = File::Spec->join($bindir, $subcmd);  

  unless (can_run($subcmd)){
    die "ERROR: '$subcmd' is not executable.\n"
  }

  # adding quotes to any args w/ spaces
  @$argv_r = map{/ / ? join('', '"', $_, '"') : $_} @$argv_r;


  # calling subcommand
  if ($getPerlDoc){
    print `perldoc -T $subcmd`;
  }
  else{
    $subcmd .= " -database $database" if $database;
    $subcmd = join(" ", $subcmd, @$argv_r[1..$#$argv_r]);
    print `$subcmd`;
  }
}

sub list_scripts{
  my $bindir = shift or die $!;
  my $listSubcmds = shift;

  # getting scripts
  opendir INDIR, $bindir or die $!;
  my @scripts = grep(/\.pl$/i, readdir INDIR);
  closedir INDIR or die $!;

  # listing or returning
  if ($listSubcmds){
    map(s/\.pl//i, @scripts);
    print join("\n", sort  @scripts), "\n";
    exit;
  }
  else{
    return \@scripts;
  }
}
