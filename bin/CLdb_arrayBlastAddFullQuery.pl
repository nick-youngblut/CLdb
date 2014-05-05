#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddFullQuery -- adding the full query sequence to the blast table

=head1 SYNOPSIS

CLdb_arrayBlastAddFullQuery.pl [flags] < blast_results.srl > blast_results_edit.srl

=head2 Required flags

=over

=item -fasta  <char>

fasta file of spacer query sequences used for BLASTn

=back

=head2 Optional flags

=over

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddFullQuery.pl

=head1 DESCRIPTION

Adding the full query sequence to the spacer blast hits table.
Just using the query sequences provided during the blastn run(s).
New field with full sequence: 'qseqfull'

=head1 EXAMPLES

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

#--- modules ---#
# core #
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Sereal qw/ encode_sereal /;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::seq qw/read_fasta/;
use CLdb::utilities qw/file_exists/;
use CLdb::arrayBlast::sereal qw/decode_file/;
use CLdb::arrayBlast::AddFullQuery qw/
				       addFullQuery
				     /;

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $fasta_in);
GetOptions(
	   "fasta=s" => \$fasta_in,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
die "ERROR: provide a fasta file of the query sequences ('-fasta')\n"
  unless defined $fasta_in;
file_exists($fasta_in, "fasta");


#--- MAIN ---#
# loading files #
## fasta of query sequences ##
my $fasta_r = read_fasta($fasta_in);
## blast ##
my $blast_r = decode_file(fh => \*STDIN);

# adding full query to blast table
addFullQuery($blast_r, $fasta_r);

print Dumper $blast_r; exit;

# encoding blast
print encode_sereal( $blast_r );


#--- Subroutines ---#


