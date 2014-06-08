#!/usr/bin/env perl

=pod

=head1 NAME

crt_parse.pl -- batch run of the CRISPR Recognition Tool (CRT); convert to CRISPRFinder format; works on draft genomes.

=head1 SYNOPSIS

crt_parse.pl [flags] *.fna 

=head2 Required flags

=over

NONE

=back

=head2 Optional flags

=over

=item -crt  <str>

CRT executable name. [CRT1.2-CLI]

=item -opts  <str>

Any opts to provide to crt (eg., '-minNR 4 -minRL 25')

=item -coord  <bool>

Array start-end coords added to each array file name? [FALSE]

=item -forks  <int>

Max number of genomes to process at the same time (0 = no forking). [0]

=item -verbose  <bool>

Verbose output

=item -h	This help message

=back

=head2 For more information:

perldoc crt_parse.pl

=head1 DESCRIPTION

This script make CRT a little easier to use for
many genomes, especially if some/all are draft genomes.

Provide >=1 genome. For each scaffold/chromosome in
each genome, crt is run. If >=1 CRISPR array is found
it is written to an individual file in a directory
labeled as the input genome fasta file.

The genome and scaffold names are comment lines in the
written array files.

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

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
use File::Temp;
use File::Path qw/rmtree/;
use IPC::Cmd qw/run can_run/;
use Parallel::ForkManager;

#--- args/flags ---#
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $coords);
my $crt_x = 'CRT1.2-CLI';
my $opts = "";
my $forks = 0;
GetOptions(
	   "crt=s" => \$crt_x,
	   "coords" => \$coords,
	   "opts=s" => \$opts,
	   "forks=i" => \$forks,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	  );

#--- I/O error ---#
can_run($crt_x) || die "ERROR: $crt_x is not in your \%PATH\n";


#--- MAIN ---#
my $pm = Parallel::ForkManager->new($forks);

foreach my $infile (@ARGV){
  die "ERROR: canot find $infile\n" unless -e $infile;
  
  # making output direcotry
  my @parts = File::Spec->splitpath($infile);
  (my $outdir = $parts[2]) =~ s/\.[^.]+$//;
  rmtree $outdir if -d $outdir;
  mkdir $outdir or die $!;

  # load genome as bio::seqio object
  my $in = Bio::SeqIO->new( -file => $infile,
			    -format => 'fasta' );

  # foreach scaffold/chromosome:
  while(my $seqo = $in->next_seq){
    $pm->start and next;
    
    # creating temp file
    my $tmpdir = File::Temp->newdir();
    my $dirname = $tmpdir->dirname;
    my $tmpfile = File::Spec->catfile($dirname, $seqo->id);

    open OUT, ">", $tmpfile or die $!;
    
    # writing sequence to temp file
    print OUT join("\n", ">" . $seqo->id, $seqo->seq), "\n";
    close OUT;

    # calling crt
    my $arrays_r = call_crt( -exe => $crt_x, 
			     -opts => $opts,	      
			     -input => $tmpfile,
			     -verbose => $verbose);
    
    # status
    printf STDERR  "Genome:\t\t%s\nScaf/chromo:\t%s %s\n", $infile,
		      $seqo->id, $seqo->desc unless $verbose;

    printf STDERR "Arrays found:\t%i\n", scalar keys $arrays_r
      unless $verbose;
    
    # writing out arrays
    write_arrays( 
		 -genome => $parts[2],
		 -arrays => $arrays_r,
		 -outdir => $outdir,
		 -coords => $coords,
		 -scafid => join("_", $seqo->id, $seqo->desc),
		 -verbose => $verbose
		  );
    
    # end fork
    $pm->finish;
  }
  
}
$pm->wait_all_children;


#--- Subroutines ---#
sub write_arrays{
  my %h = @_;

  my $arrays_r = $h{-arrays};
  my $scafid = $h{-scafid};
  $scafid =~ s/^_+//;
  $scafid =~ s/_+$//;


  foreach my $arrayNum (keys %$arrays_r){
    # outfile
    ## coords in name?
    my $start_end = $h{-coords} ? 
      join("-", '', $arrays_r->{$arrayNum}{start},
	   $arrays_r->{$arrayNum}{end}, '' ) : "";	   
      

    my $outfile = join("", $scafid, $start_end, $arrayNum, ".txt");
    $outfile = File::Spec->catfile($h{-outdir}, $outfile);

    open OUT, ">", $outfile or die $!;
    printf OUT "# %s\n", $h{-genome};
    print OUT "# $scafid\n";

    foreach my $line (@{$arrays_r->{$arrayNum}{'array'}}){
      print OUT join("\t", @$line), "\n";    
    }

    close OUT;

    # status
    print STDERR "Array file written:\t$outfile\n"
      unless $h{-verbose};
  }

}


sub call_crt{
  my %h = @_;

  my $cmd = join(" ", $h{-exe}, 
		 "-screen 1", 
		 $h{-opts}, 
		 $h{-input}, "|");


  open my $pipe, $cmd or die $!;
  my $arrays_r = parse_crt($pipe, %h);
  close $pipe;

  return $arrays_r;
}


sub parse_crt{
  # parsing CRT output

  my $pipe = shift or die "Provide a pipe fh\n";
  my %h = @_;

  my %res;
  while(<$pipe>){
    chomp;

    # status of output
    if(/^invalid argument/i){
      print;
      while(<$pipe>){ print; }
      exit;
    }

    # if(/^CRISPR \d+/){
    #   # array number

    #   while(<$PIPE>){
    # 	chomp;
    # 	next if /^-/;
    # 	last if /^Repeats:/;


    #   }
      
    # }
    
    # loading arrays
    if(/^CRISPR \d+/){
      my @head = split / +/;    # 1=num; 3=start; 5=end
      die "ERROR: line $. is not formatted correctly!\n"
	unless scalar @head == 6;

      $res{$head[1]}{'start'} = $head[3];
      $res{$head[1]}{'end'} = $head[5];
      die "ERROR: start > end at line $.\n"
	unless $head[3] <= $head[5];

      # each array
      while(<$pipe>){
	chomp;
	last if /^Repeats:/;
	next if /^(POSITION|^-)/;

	my @line = split /\t+/;
	die "ERROR: line $. does not have 2 or 4 tab-delimited columns!\n"
          unless scalar @line == 2 || scalar @line == 4;

        # calculating end #
        my ($DR_len, $spacer_len) = (0,0);
                                $DR_len = length $line[1] if defined $line[1];
        $spacer_len = length $line[2] if defined $line[2];
        $line[3] = $line[0] + $DR_len + $spacer_len - 1;

        $line[2] = "\t" unless defined $line[2];
        push @{$res{$head[1]}{'array'}}, \@line;
      }
    }
  }

  return \%res;
}
  
