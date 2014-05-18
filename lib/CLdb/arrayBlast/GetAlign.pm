package CLdb::arrayBlast::GetAlign;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use List::MoreUtils qw/any/;

## CLdb
use CLdb::seq qw/revcomp/;

# export #
use base 'Exporter';
our @EXPORT_OK = ();
	
=head1 NAME

CLdb::arrayBlast::GetAlign

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for getting alignments of protospacer and crDNA
form decoded blast srl.

=cut


=head2 parse_outfmt

Parsing the comma-delimited arg 
provided via the '-outfmt' flag.

Binning columns by where they can be found.

=head3 IN

$outfmt :  $; comma-delimited string

=head3 OUT

=cut

push @EXPORT_OK, 'parse_outfmt';

sub parse_outfmt{
  my $outfmt = shift or return undef;

  $outfmt =~ tr/A-Z/a-z/;
  my @outfmt = split /\s*,\s*/, $outfmt;
  my %outfmt_order;
  for my $i (0..$#outfmt){ $outfmt_order{$outfmt[$i]} = $i; }

# supported columns:
## from CLdb:
  # subtype
  # taxon_name
  # taxon_id
  # scaffold
## from blast srl
### from run
  # blastdb  (BlastOutput_db)
### from CLdb_info
  # crDNA_start
  # crDNA_end
  # array_sense_strand
### from hsp
  # subject_scaffold
  # proto_start
  # proto_end
  # protoX_start
  # protoX_end
  # proto_strand
  # identity
  # evalue
  # bitscore

  my %index;       # user_input => name_in_db
  $index{CLdb} = {   # only loci.* fields supported
		  'subtype' => 'loci.subtype',
		  'taxon_name' => 'loci.taxon_name',
		  'taxon_id' => 'loci.taxon_id',
		  'scaffold' => 'loci.scaffold' 
		 };
  $index{run} = { blastdb => 'BlastOutput_db' };
  $index{crRNA_info} = {
			crdna_start => 'region_start',
			crdna_end => 'region_end',
			array_sense_strand => 'array_sense_strand'
		       };
  $index{hsp} = { 
		 subject_scaffold => 'subjectScaffold',
		 proto_start => 'protoFullStart',
		 proto_end => 'protoFullEnd',
		 protox_start => 'protoFullXStart',
		 protox_end => 'protoFullXEnd',
		 proto_strand => 'subjectStrand',
		 identity => 'Hsp_identity',
		 evalue => 'Hsp_evalue',
		 bitscore => 'Hsp_bit-score'
		};


  # checking for existence of columns
  foreach my $field (keys %outfmt_order){
    my $cnt = 0;
    foreach my $cat (keys %index){
      $cnt++ if exists $index{$cat}{$field};
    }
    croak "ERROR: -outfmt field '$field' not supported\n"
      unless $cnt == 1;
  }

  # loading outfmt index
  my %outfmt;
  sub set_index{  # setting index for which fields needed and where they are
    my $index_r = shift;
    my $outfmt_r = shift;
    my $outfmt_order_r = shift;
    my $cat = shift;
    map{ $outfmt_r->{ $cat }{ $index_r->{$cat}{$_} } = $outfmt_order_r->{$_} 
       if exists $index_r->{$cat}{$_} }  keys %$outfmt_order_r;
  }

  set_index(\%index, \%outfmt, \%outfmt_order,  'CLdb', );
  set_index(\%index, \%outfmt,  \%outfmt_order, 'run');
  set_index(\%index, \%outfmt, \%outfmt_order,  'crRNA_info');
  set_index(\%index, \%outfmt, \%outfmt_order, 'hsp');

 
  #print Dumper %outfmt; exit;
  return \%outfmt;
}


=head2 get_alignProto

Getting crDNA-protospacer alignments
from the blast srl DS.

See parse_outfmt for fields that may need
to be parsed from blast srl

=head3 IN

hash of args:
blast :  blast srl Ds
outfmt :  outfmt values
queries :  \%{locus_id}{spacer_id}{field} = value; used for selecting specific hits
array :  just alignments of hit to spacers in arrays? [FALSE]
crRNA_ori :  orient by crRNA? [FALSE]
verbose :  verbose [TRUE}

=head3 OUT

Writing fasta to STDOUT

=cut

push @EXPORT_OK, 'get_alignProto';

sub get_alignProto{
  my %h = @_;
  my $spacer_r = exists $h{blast} ? 
    $h{blast} : confess "Provide 'blast'";
  my $outfmt_r = exists $h{outfmt} ?
    $h{outfmt} : confess "Provide 'outfmt'";
  my $queries_r = exists $h{queries} ? 
    $h{queries} : confess "Provide 'queries'";
  my $array = exists $h{array} ? 
    $h{array} : confess "Provide 'array'";
  my $crDNA_ori = exists $h{crDNA_ori} ?
    $h{crDNA_ori} : confess "Provide 'crDNA_ori'";
  my $verbose = $h{verbose};
  
  foreach my $run (keys %$spacer_r){
    next unless exists $spacer_r->{$run}{'BlastOutput_iterations'};  # must have hit

    # getting blastdbfile 
    my $blastdbfile = exists $spacer_r->{$run}{'BlastOutput_db'} ?
	$spacer_r->{$run}{'BlastOutput_db'} :
	  confess "Cannot find BlastOutput_db in run $run";
    my @parts = File::Spec->splitpath($blastdbfile);
    $blastdbfile = $parts[2];

    # status
    print STDERR "Retrieving sequences from blast db: '$blastdbfile'\n"
      if $verbose;

    # each iteration
    foreach my $iter ( @{$spacer_r->{$run}{'BlastOutput_iterations'}{'Iteration'}} ){
      # skipping iterations without hits
      next unless exists $iter->{Iteration_hits} and
        $iter->{Iteration_hits} !~ /^\s*$/;
      next unless exists $iter->{Iteration_hits}{Hit};

      # getting crRNA info if needed
      ## come back to this for each hsp (using locus_id & spacer_id)
	
      # iterating through hits
      foreach my $hit ( @{$iter->{Iteration_hits}{Hit}} ){
        next unless exists $hit->{Hit_hsps}{Hsp};

	foreach my $hspUID (keys %{$hit->{Hit_hsps}{Hsp}}){
	  # filtering 	  
	  ## next unless alignment
	  my $aln_r = exists $hit->{Hit_hsps}{Hsp}{$hspUID}{aln} ?
	    $hit->{Hit_hsps}{Hsp}{$hspUID}{aln} : next;	    
	  ## if CLdb_array-hit: skip if 
	  if( exists $hit->{Hit_hsps}{Hsp}{$hspUID}{'CLdb_array-hit'} ){
	    # if -array; just keeping hits to array
	    next if $array and $hit->{Hit_hsps}{Hsp}{$hspUID}{'CLdb_array-hit'} == -1;
	    # if ! -array; just keeping hits to protospacer
	    next if ! $array and $hit->{Hit_hsps}{Hsp}{$hspUID}{'CLdb_array-hit'} == 1;
	  }	  	  	  

	  # foreach alignment, writing sequence
	  foreach my $locus_spacer (keys %$aln_r){
	    my ($locus_id, $spacer_id) = split /\|/, $locus_spacer, 2;
	    # check to make sure crRNA_info exists
	    exists $iter->{crRNA_info}{$locus_spacer} ||
	      confess "Cannot find crRNA_info for $locus_spacer\n";


	    # collecting columns of iterest for sequence names
	    my %columns;
	    ## run fields
	    $columns{$blastdbfile} = $outfmt_r->{run}{BlastOutput_db}
	      if exists $outfmt_r->{run}{BlastOutput_db};


	    ## CLdb
	    if( %$queries_r){
	      # skip if selecting subset and this locus-spacer is not 1 of them
	      next unless exists $queries_r->{$locus_id}{$spacer_id};
	      foreach my $field (keys %{$queries_r->{$locus_id}{$spacer_id}}){
		next if $field eq 'locus_id' or $field eq 'spacer_id';
		my $val = $queries_r->{$locus_id}{$spacer_id}{$field};
		$field = join(".", 'loci', $field);
		$columns{ $val } = $outfmt_r->{CLdb}{$field};  # field_value => order in seq name		
	      }
	    }

	    ## crRNA_info fields
	    if(exists $outfmt_r->{crRNA_info}){
	      foreach my $field (keys %{$outfmt_r->{crRNA_info}} ){    
		my $val = exists $iter->{crRNA_info}{$locus_spacer}{$field} ?
		  $iter->{crRNA_info}{$locus_spacer}{$field} :
		    confess "ERROR: cannot find field '$field' for $locus_spacer\n";
		$columns{ $val } = $outfmt_r->{crRNA_info}{$field};   # field_value => order
	      }
	    }

	    ## hsp fields
	    if(exists $outfmt_r->{hsp}){
	      foreach my $field (keys %{$outfmt_r->{hsp}} ){    
		my $val = exists $hit->{Hit_hsps}{Hsp}{$hspUID}{$field} ?
		  $hit->{Hit_hsps}{Hsp}{$hspUID}{$field} :
		    confess "ERROR: cannot find field '$field' for $locus_spacer\n";
		$columns{ $val } = $outfmt_r->{hsp}{$field};   # field_value => order
	      }
	    }

	    # writing out the alignment
	    foreach my $ele (sort keys %{$aln_r->{$locus_spacer}}){  # $ele= crDNA|protospacer
	      my $seqName = join("", ">", 
				 join("|", $ele,
				      $locus_id, 
				      $spacer_id,
				      $hspUID,
				      sort{$columns{$a} cmp $columns{$b}}
					     keys %columns),
				);

	      # by default, sequences oriented by blast hit (crDNA)
	      ## orienting by protospacer unless $crRNA_ori 
	      $aln_r->{$locus_spacer}{$ele} = revcomp($aln_r->{$locus_spacer}{$ele})
		unless $crDNA_ori;

	      print join("\n", $seqName, 
			 $aln_r->{$locus_spacer}{$ele}), "\n";
	    }
	  }
	}
      }
    }
  }
}



=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::GetAlign


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=CRISPR_db>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/CRISPR_db>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/CRISPR_db>

=item * Search CPAN

L<http://search.cpan.org/dist/CRISPR_db/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2013 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See L<http://dev.perl.org/licenses/> for more information.

=cut

1; 
