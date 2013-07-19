#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;

### args/flags
#pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $replace);
GetOptions(
	   "replace" => \$replace,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$ARGV[0] = "CLdb.sqlite" unless $ARGV[0];

### MAIN


### Subroutines


__END__

=pod

=head1 NAME

CLdb_makeDB.pl -- Initial DB construction

=head1 SYNOPSIS

CLdb_makeDB.pl [options] [DATABASE_name]

=head2 options

=over

=item -r 	Replace existing database.

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_makeDB.pl

=head1 DESCRIPTION

Make all of the CRISPR_db tables.

=head1 EXAMPLES

=head2 Naming database 'CRISPR_db1'

CLdb_makeDB.pl CRISPR_db1

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

