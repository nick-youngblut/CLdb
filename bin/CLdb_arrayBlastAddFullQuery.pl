#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddFullQuery -- adding the full query sequence to the blast table

=head1 SYNOPSIS

CLdb_arrayBlastAddFullQuery.pl [flags] < blast_results.txt > blast_results_info.txt

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
The blast table should be formatted as '-outfmt 7'


=head2 Fields added to blast table:

=over

=item query_seq_full

Full length query sequence

=back

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
use DBI;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::utilities qw/
	file_exists/;
use CLdb::seq qw/
	read_fasta/;
use CLdb::blast qw/
	read_blast_file
	write_blast_file/;


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
my ($lines_r) = read_blast_file();

## adding full length query ##
add_full_query_seq($lines_r, $fasta_r);

# writing edited fasta #
write_blast_file($lines_r);




#--- Subroutines ---#
sub add_full_query_seq{
# adding full length query sequence to blast table #
  my ($lines_r, $fasta_r) = @_;
  
  foreach my $query ( keys %$lines_r ){
	  
    # checking for existence of query in fasta #
    (my $name = $query) =~ s/# Query: //i;
    unless (exists $fasta_r->{$name}){
      warn "WARNING: '$name' not found in fasta. Skipping!\n";
      next;
    }
    
    # adding sequence to blast table #
    foreach my $db ( keys %{$lines_r->{$query}} ){
      
      foreach my $blast ( keys %{$lines_r->{$query}{$db}} ){
	next unless exists $lines_r->{$query}{$db}{$blast}{'fields'};
	
	$lines_r->{$query}{$db}{$blast}{'fields'} .= ", query_seq_full";
	$lines_r->{$query}{$db}{$blast}{'fields_sep'}{"query_seq_full"} = 
	  scalar keys %{$lines_r->{$query}{$db}{$blast}{'fields_sep'}};
		    
	foreach my $hit ( @{$lines_r->{$query}{$db}{$blast}{'hits'}} ){
	  $hit .= "\t$fasta_r->{$name}";
	}
      }
    }
  }
  #print Dumper $lines_r; exit
}








