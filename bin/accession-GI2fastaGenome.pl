#!/usr/bin/env perl

=pod

=head1 NAME

accession-GI2fastaGenome.pl -- get genome fasta (or genbank) from NCBI accession/GI numbers

=head1 SYNOPSIS

accession-GI2fastaGenome.pl [flags] < accessions_GIs.txt

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -format

Output file format (see BioPerl SeqIO). [fasta]

=item -trial

Number of trials to attempt for downloading a genome. [10]

=item -prefix

Output file prefix. [.]

=item -forks

Number of genomes to download in parallel (0 = no forking). [0]

=item -h	This help message

=back

=head2 For more information:

CLdb_perldoc accession-GI2fastaGenome.pl

=head1 DESCRIPTION

A wrapper for Bio::DB::GenBank designed
for downloading genomes or other large
sequence files. 

Each accession or GI provided will be written
to a different file.

Provide a list of accession numbers.
Format:

=over

=item accession or GI    

=item name (optional)

=back

Accession numbers, versioned accession numbers,
or GI numbers can be used.

Downloading is not always successful, so
sequences are downloaded in batches. If
not all sequences in the batch are downloaded
successfully, the batch is re-tried.

=head1 EXAMPLES

=head2 GI of ecoli strain

printf "170079663\nEcoli" | accession-GI2fastaGenome.pl > ecoli.fasta

=head2 genbank file output

printf "170079663\nEcoli" | accession-GI2fastaGenome.pl -f genbank > ecoli.gbk

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


#--- modules ---#
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::SeqIO;
use Bio::DB::GenBank;
use Parallel::ForkManager;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my $trial_limit = 10;
my $format = 'fasta';
my $prefix;
my $cpus = 0;
GetOptions(
	   "trial=i" => \$trial_limit,
	   "format=s" => \$format,
	   "prefix=s" => \$prefix,
	   "forks=i" => \$cpus,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
# output files
$prefix = File::Spec->rel2abs($prefix);
die "ERROR: cannot find $prefix\n"
  unless -d $prefix;
print STDERR "Writing files to '$prefix'\n";

my $ext = get_ext($format);

#--- MAIN ---#
# getting accessions|GIs 
my $acc_r = load_accessions();
# initialize
my $gb = new Bio::DB::GenBank;

# fork setup
my $pm = Parallel::ForkManager->new($cpus);

foreach my $acc_GI ( keys %$acc_r ){
  # forking
  $pm->start and next;

  # status
  printf STDERR "Attempting to stream: %s (accession/GI = %s)\n", 
    $acc_r->{$acc_GI}, $acc_GI;

  my $trials = 0;
  while(1){
    $trials++;
    print STDERR "Starting trial: $trials\n"
      unless $trials == 1;
    
    # streaming in sequences
    my $seqio = $gb->get_Stream_by_acc( [$acc_GI] );
    
    my $seq_out = Bio::SeqIO->new( -file => join("", ">", "$prefix/", $acc_r->{$acc_GI}, $ext), 
				   -format => $format);
	
    # checking to make sure seqio is defined
    if(defined $seqio){ # all sequences appear to have been written out; writing from tempFile to STDOUT
      # writing out sequences
      my $seq_cnt = 0;
      while( my $seq = $seqio->next_seq ){
	$seq_out->write_seq($seq);
	$seq_cnt++;
      }
      if($seq_cnt > 0){ 
	last;
	}
      elsif( $trials >= $trial_limit ){
	print STDERR "WARNING: exceded number of trials for $acc_GI. Skipping\n";
	last;
      }
      else{
	printf STDERR "WARNING: nothing streamed for $acc_GI.Retrying\n";      
	next;
      }
    }
    elsif( $trials >= $trial_limit ){
      print STDERR "WARNING: exceded number of trials for $acc_GI. Skipping\n";
      last;
    }
    else{
      printf STDERR "WARNING: nothing streamed for $acc_GI.Retrying\n";
      next;
    }
  }
  $pm->finish;
}
$pm->wait_all_children;

#--- Subroutines ---#
sub load_accessions{
  my %acc;
  while(<>){
    chomp;
    next if /^\s*$/;

    my @line = split /\t/;
    
    # name associated with GI|accession
    my $name = defined $line[1] ? $line[1] : $line[0];
    $name =~ s/[ .*&:!@_()#+[\]]+/_/g;
    $name =~ s/_$//;
    $acc{$line[0]} = $name;
  }

  #print Dumper %acc; exit;
  return \%acc;
}


sub get_ext{
# making a file extension based on output file format
  my ($format) = @_;

  my $ext;
  if($format eq 'fasta'){ $ext = '.fasta'; }
  elsif($format eq 'genbank'){ $ext = '.gbk'; }
  else{ die "ERROR: only 'fasta' or 'genbank' supported for output file format\n"; }

  return $ext;
}

