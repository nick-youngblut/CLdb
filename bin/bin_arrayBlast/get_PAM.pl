#!/usr/bin/env perl

=pod

=head1 NAME

get_PAM -- getting PAM for each protospacer

=head1 SYNOPSIS

get_PAM [flags] > PAMs

=head2 Required flags

=over

=head3 NOTE: Use either -fasta OR -table

=item -fasta  <char>

Fasta file of the protospacers (& possibly crDNA). 
'-' if from STDIN.

=item -table  <char>

Tab-delim table containing the protospacer sequence 
(1st column of table). '-' if from STDIN.

=back

=head2 Optional flags

=over

=item -PAM  <int>

start-stop of the PAM region (2 values required). 
See DESCRIPTION for details. [-3 -1]

=item -revcomp  <bool>

Reverse complement protospacer sequence
before extracting PAM (if protospacer
is not oriented on the correct strand). [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

perldoc get_PAM

=head1 DESCRIPTION

Pull out the PAM regions for each hit to a protospacer. 
You can either provide the output fasta from
arrayBlastGetAlign.pl (-fasta) or the table
from arrayBlastGetProto.pl (-table). 
The output format will match the input.

=head2 -PAM

The flag designates the region (relative to the
protospacer) containing the PAM. Negative 
values = upstream from the proto, while
positive values = downstream. 
So, the default [-3 -1] will select the 3bp
immediately upstream of the protospacer.
'-PAM 1 2' would select the 2bp immediately
upstream of the protospacer.


=head1 EXAMPLES

=head2 Basic Usage:


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use FindBin;
use lib "$FindBin::RealBin/../../lib/";

### CLdb
use CLdb::arrayBlast::PAM qw/ make_pam_index
			      read_protoTable
			      getPAM
			      writePAM/;
use CLdb::seq qw/read_fasta/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $revcomp_b);
my @PAM;
my ($fasta_in, $table_in); 
GetOptions(
	   "PAM=i{2,2}" => \@PAM,	   
	   "fasta=s" => \$fasta_in,
	   "table=s" => \$table_in,
	   "revcomp" => \$revcomp_b,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# pam index
@PAM = (-3, -1) unless @PAM;
my $pam_index_r = make_pam_index(\@PAM);
# fasta or table from STDIN?
die "Provide either -fasta or -table\n"
  unless defined $fasta_in or defined $table_in;
my $fh;
foreach my $file ( $fasta_in, $table_in ){
  if(defined $file){
    $file eq '-' ? $fh = \*STDIN : 
      open $fh, $file or die $!;
  }
}


#--- MAIN ---#
# getting fasta if fasta
my $fasta_r;
if(defined $fasta_in){ 
  $fasta_r = read_fasta(-fh => $fh); 
}
elsif(defined $table_in){ 
  $fasta_r = read_protoTable(-fh => $fh); 
}
else{ die "Logic error $!\n"; }

# getting PAM region from proto
my $pam_r = getPAM($fasta_r, $pam_index_r, $revcomp_b);

# writing out pam
writePAM($pam_r, $fasta_in, $table_in);

