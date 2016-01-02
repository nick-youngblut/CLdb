#!/usr/bin/env perl

=pod

=head1 NAME

xml2srl -- converting blast output (xml format) to binary serialization format

=head1 SYNOPSIS

xml2srl [flags] < blast_output.xml > blast_output.srl

=head2 Required flags

NONE

=head2 Optional flags

=over

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

CLdb_perldoc xml2srl

=head1 DESCRIPTION

Simple script for converting blast output in xml format
(-outfmt 5) to a perl binary serialization format.

Multiple concatenated blast xml output files can be 
provided via STDIN. 

Each blast hit->hsp will get a unique identifier.

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

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
use XML::LibXML::Simple qw/XMLin/;
use Sereal qw/ encode_sereal /;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../../lib";
use CLdb::arrayBlast::sereal qw/
				 blast_all_array
			       /;

#--- parsing args ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database);
GetOptions(
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#


#--- MAIN ---#
my %blast;
my $run_cnt = 0;
while(<>){
  if( exists $blast{$run_cnt} and (/<\?xml version="1.0"\?>/ or eof)){
    $blast{$run_cnt} = XMLin( $blast{$run_cnt} );
    $run_cnt++;
  }
  $blast{$run_cnt} .= $_;
}

# making each hit & hsp an array (even if only 1)
foreach my $run (keys %blast){
  unless (ref $blast{$run} eq 'HASH'){
    delete $blast{$run};
    next;
  }
  blast_all_array($blast{$run});
}

# serializing 
my $encoder = Sereal::Encoder->new();
print $encoder->encode( \%blast );

