#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlast2Align.pl -- align sequence fields in an arrayBlast file

=head1 SYNOPSIS

CLdb_arrayBlast2Align.pl [flags] < spacer_blast_hits.txt

=head2 Optional flags

=over

=item -seq_fields

Fields containing the sequences to align (2 field names required). 
['query_seq_full', 'proto_seq_fullx']

=item -names

Fields to use for sequence names in the alignment. ['query id', 'subject id']

=item -v  <bool>

Verbose output. [TRUE]

=item -h  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlast2Align.pl

=head1 DESCRIPTION

Make an alignment of fields in arrayBlast output (blastn -outfmt 7).
The fields must contain nucleotide sequence and just sequence.

clustalw is used to align the sequences.

=head2 Output

'-name' fields are separated by '__'.

The aligments are writen next to the sequence names (name\tsequence).
A blank line separates each alignment.

=head1 EXAMPLES

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
use Bio::SeqIO;
use IO::String;
use Bio::Tools::Run::Alignment::Clustalw;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::utilities qw/
	file_exists 
	get_file_path/;
use CLdb::blast qw/
	read_blast_file/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose);
my @seq_fields = ('query_seq_full', 'proto_seq_fullx');
my @name_fields = ('query id', 'subject id');
GetOptions(
	   "seq_fields=s{2,2}" => \@seq_fields,
	   "name=s{1,}" => \@name_fields, 		# output by query
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
#--- MAIN ---#
my $blast_r = read_blast_file();

#print Dumper $blast_r; exit;

blast_2_align($blast_r, \@seq_fields, \@name_fields);



#--- Subroutines ---#
sub blast_2_align{
  my ($blast_r, $seq_fields_r, $name_fields_r) = @_;

  foreach my $query (keys %$blast_r){
    foreach my $db (keys %{$blast_r->{$query}}){
      foreach my $blast (keys %{$blast_r->{$query}{$db}}){
	# getting fields
	my $field_i = $blast_r->{$query}{$db}{$blast}{'fields_sep'};
	map{ die "ERROR: cannot find field: '$_' in $query -> $db\n" 
	       unless exists $field_i->{$_} } (@$seq_fields_r, @$name_fields_r);

	foreach my $hit (@{$blast_r->{$query}{$db}{$blast}{hits}}){
	  my @hit = split /\t/, $hit;

	  # getting sequence fields 
	  ## check fields exist in blast hit	  
	  map{ die "ERROR: cannot find %i column in blast table\n",
		 $_+1 unless defined $hit[$_] }  @{$field_i}{ @{$name_fields_r} };
	  (my $tmp = $db) =~ s/.+\///;
	  my $name_str = join("__", 
			      @hit[ @{$field_i}{ @{$name_fields_r} } ],
			      $tmp
			     );
	  
	  # pulling out sequences
	  ## check fields exist in blast hit
	  map{ die "ERROR: cannot find %i column in blast table\n",
		 $_+1 unless defined $hit[$_] }  @{$field_i}{ @{$seq_fields_r} };

	  ## combining sequence strings
	  my @seq_names = (join("__", ">$seq_fields_r->[0]", $name_str),
			   join("__", ">$seq_fields_r->[1]", $name_str));
	  my $seq_str = join("\n", $seq_names[0],
			     $hit[$field_i->{$seq_fields_r->[0]}],
			     $seq_names[1],
			     $hit[$field_i->{$seq_fields_r->[1]}]);
	 
	  ## making seqIO object
	  my $strfh = IO::String->new($seq_str);
	  my @seqIOs;
	  my $seqo = Bio::SeqIO-> new(
                             -fh     => $strfh,
                             -format => 'fasta',
                             );
	  while(my $seq = $seqo->next_seq()){ push @seqIOs, $seq; }

	  ## alignment
	  my $factory = Bio::Tools::Run::Alignment::Clustalw->new(quiet => 1);
	  my $aln = $factory->align(\@seqIOs);

	  ## changing back sequence names
	  my %seqs; 
	  foreach my $seq ( $aln->each_seq() ){
	    #$seq->display_id( shift @seq_names );
	    $seqs{shift @seq_names} = $seq->seq();
	  }

	  # writing out alignment	  
	  foreach my $seq (sort{$b cmp $a} keys %seqs){
	    $seqs{$seq} =~ s/\./-/g;
	    print join("\t", $seq, $seqs{$seq}), "\n";
	  }
	  print "\n";
	}
      }
    }
  }
}

