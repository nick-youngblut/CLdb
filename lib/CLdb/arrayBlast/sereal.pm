package CLdb::arrayBlast::sereal;

use 5.006;
use strict;
use warnings FATAL => 'all';
use Data::Dumper;
use Carp qw/carp croak confess/;
use Clone qw/clone/;
use Data::Uniqid qw/uniqid/;

use base 'Exporter';
our @EXPORT_OK = '';

use Sereal qw/decode_sereal/;

=head1 NAME

sereal - The great new sereal!

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

subroutines for editing blast xml converted to json format

=head1 EXPORT


=head1 SUBROUTINES/METHODS


=head2 decode_file

Sereal file decoding

=head3 IN

hash with either: file => file_name  or fh => filehandle_ref

=head3 OUT

decoded data structure

=cut

push @EXPORT_OK, 'decode_file';

sub decode_file{
  my %opt = @_;

  my $fh;
  if(exists $opt{file}){
    open $fh, $opt{file} or die $!;
  }
  elsif(exists $opt{fh}){
    $fh = $opt{fh};
  }
  else{ die "ERROR: provide a file or filehandle\n"; }
  
  my $str = '';
  $str .= $_ while <$fh>;
  close $fh;

  return decode_sereal($str);
}


=head2 blast_all_array

Making sure that 'Hit' and 'Hsp' are arrays, 
even if just 1 present

=head3 IN

blast table hash_ref. Parsed xml blast output

=head3 Output

blast table hash_ref altered

=cut

push @EXPORT_OK, 'blast_all_array';

sub blast_all_array{
  my ($blast_r) = @_;

  # iterating through each hit
  foreach my $iter ( @{$blast_r->{'BlastOutput_iterations'}{'Iteration'}} ){
    next unless exists $iter->{'Iteration_hits'};
    next if ! ref $iter->{'Iteration_hits'} or 
      ! exists $iter->{'Iteration_hits'}{'Hit'};
    my $hits_ref = $iter->{'Iteration_hits'}{'Hit'}; # hits in iteration (array_ref or ref to hash)

    # multiple hits or just 1?
    if(ref $hits_ref eq 'ARRAY'){ # multiple
      foreach my $hit ( @$hits_ref){
        make_hsp_hash($blast_r, $iter, $hit);
      }
    }
    else{ # making array
      make_hsp_hash($blast_r, $iter, $hits_ref);
      $iter->{'Iteration_hits'}{'Hit'} = [$hits_ref]
    }

    # making a unique ID for each hsp (each 'blast hit' has a unique ID)
    sub make_hsp_hash{
      my ($blast_r, $iter, $hit) = @_;

      # multiple hsp or just 1?
      if( ref $hit->{Hit_hsps}{Hsp} eq 'ARRAY' ){
	my %hsp_byUID;
	foreach my $hsp (@{$hit->{Hit_hsps}{Hsp}}){
	  my $UID = uniqid;
	  $hsp_byUID{$UID} = clone( $hsp );	  
	}
	$hit->{Hit_hsps}{Hsp} = \%hsp_byUID;
      }
      else{ # just hash ref
	my $UID = uniqid;
	$hit->{Hit_hsps}{Hsp} = {$UID => $hit->{Hit_hsps}{Hsp}};
      }
    }
  }
# print Dumper $blast_r; exit;
}


=head2 parse_outfmt

Parsing '-outfmt' argument.
Formatted like blast -outfmt.
Only 6 or 7 supported

=head3 Input

space-delim string 

=head3 Output
{
 {comment: [1|0]}
 {fields: [fields]}
}

=cut

push @EXPORT_OK, 'parse_outfmt';

sub parse_outfmt{
  my ($outfmt) = @_;

  my @l = split / +/, $outfmt;
  croak "-outfmt argument must start with '6' or '7'"
    unless $l[0] == 6 || $l[0] == 7;


  my %fields;
 
  # comments 
  $fields{comments} = $l[0] == 7 ? 1 : 0;

  # all lower case
#  map{ tr/A-Z/a-z/ } @l[1..$#l];

  # fields
  $fields{unclassified_fields} = [@l[1..$#l]];

 
#  print Dumper %fields; exit;
  return \%fields;
}


=head2 classify_fields

Classifying blast++ fields based on location in blast xml
file structure

=head3 Input

parse fields hash_ref

=head3 Output

edited fields hash_ref

%: classification => [fields]


=cut

push @EXPORT_OK, 'classify_fields';

sub classify_fields {
  my ($fields_r) = @_;

  # classification table
  my %class = (
	       # run
	       blastdb => ['run', 'BlastOutput_db'],
	       # iteration
	       qseqid => ['iter', 'Iteration_query-def'],
	       qlen => ['iter', 'Iteration_query-len'], 
	       # hit
	       sseqid => ['hit', 'Hit_id'],
	       sacc => ['hit', 'Hit_accession'],
	       slen => ['hit', 'Hit_len'],
	       # hsp
	       qstart => ['hsp', 'Hsp_query-from'],
	       qend => ['hsp', 'Hsp_query-to'],
	       sstart => ['hsp', 'Hsp_hit-from'],
	       send => ['hsp', 'Hsp_hit-to'],
	       qseq => ['hsp', 'Hsp_qseq'],
	       sseq => ['hsp', 'Hsp_hseq'],
	       evalue => ['hsp', 'Hsp_evalue'],
	       bitscore => ['hsp', 'Hsp_bit-score'],
	       score => ['hsp', 'Hsp_score'],
	       length => ['hsp', 'Hsp_align-len'],
	       pident => ['hsp', 'Hsp_identity'],
	       mismatch => ['hsp', undef],
	       gaps => ['hsp', 'Hsp_gaps'],
	       frames => ['hsp', 'Hsp_query-frame', 'Hsp_hit-frame'],
	       qframe => ['hsp', 'Hsp_query-frame'],
	       sframe => ['hsp', 'Hsp_hit-frame'],
	       ## CLdb hsp
	       arrayHit => ['hsp', 'CLdb_array-hit'],
	       subjectScaffold => ['hsp', 'subjectScaffold'],
	       protoFullStart => ['hsp', 'protoFullStart'],
	       protoFullXSeq => ['hsp', 'protoFullXSeq'],
	       subjectStrand => ['hsp', 'subjectStrand'],
	       protoFullXStart => ['hsp', 'protoFullXStart'],
	       protoFullSeq => ['hsp', 'protoFullSeq'],
	       protoX => ['hsp', 'protoX'],
	       protoFullEnd => ['hsp', 'protoFullEnd'],
	       protoFullXEnd => ['hsp', 'protoFullXEnd'],
	       
	       # param
	       gapopen => ['param', 'Parameters_gap-open'],
	       # crRNA
	       spacer_start => ['crRNA_info', 'spacer_start'],
	       spacer_end => ['crRNA_info', 'spacer_end'],
	       region_start => ['crRNA_info', 'region_start'],
	       region_end => ['crRNA_info', 'region_end'],
	       spacer_scaffold => ['crRNA_info', 'scaffold'],
	       cluster_id => ['crRNA_info', 'cluster_id'],
	       genome_fasta => ['crRNA_info', 'genome_fasta'],
	       locus_id => ['crRNA_info', 'locus_id'],
	       spacer_id => ['crRNA_info', 'spacer_id'],	       
	       query_id => ['crRNA_info', 'query_id'],
	       array_sense_strand => ['crRNA_info', 'array_sense_strand'],
	       spacer_seq =>['crRNA_info', 'spacer_sequence'],
	       crDNA_seq => ['crRNA_info', 'crDNA']
	      );


  # classification
  foreach my $field ( @{$fields_r->{unclassified_fields}} ){
#    (my $tmp = $field) =~ tr/A-Z/a-z/;

    # adding classification for field
    if( exists $class{$field} ){   # field can be classified
      if( $class{$field}->[0] eq 'crRNA_info' ){
	$fields_r->{crRNA_index}{ $field } = $class{$field};
	push @{$fields_r->{crRNA_fields}}, $field;
      }
      else{
	$fields_r->{index}{ $field } = $class{$field};
	push @{$fields_r->{classified_fields}}, $field;
      }
    }
    else{
      print STDERR "ERROR: field '$field' not supported\n";
      print STDERR "\nSUPPORTED FIELDS:\n" . 
	join(",\n", sort keys %class) . "\n"; 
      exit(1);
    }    
  }

#  print Dumper $fields_r; exit;
}


=head2 blast_xml2txt

convert blast from parsed xml format to a blast table

=head3 Input (hashref)

'blast': hash_ref; blast xml parsed to hash 

'fields': fields: fields structure from parse_outfmt & classify_fields

=head3 Output

tab-delimited blast table written to STDOUT

=cut

push @EXPORT_OK, 'blast_xml2txt';

sub blast_xml2txt {
  my %opts = @_; 
  my $blast_r = $opts{blast} || croak "No 'blast' arg provided $!\n";
  my $fields_r = $opts{fields} || croak "No 'fields' arg provided $!\n";


  # iterating through each hit
  foreach my $iter ( @{$blast_r->{'BlastOutput_iterations'}{'Iteration'}} ){
    # skipping iterations without hits
    next unless exists $iter->{Iteration_hits} and
      $iter->{Iteration_hits} !~ /^\s*$/;
    next unless exists $iter->{Iteration_hits}{Hit};
    
    # comments for the query (if needed)
    if( $fields_r->{comments} == 1){
      print join(" ", '# blast:', $blast_r->{'BlastOutput_version'}), "\n";
      print join(" ", '# Query:', $blast_r->{'BlastOutput_query-def'}), "\n";
      print join(" ", '# Database:', $blast_r->{'BlastOutput_db'}), "\n";
      my $crRNA_fields = exists $fields_r->{crRNA_fields} ? 
	join(", ",  @{$fields_r->{crRNA_fields}} ) : "";
      print join(" ", '# Fields:',
		 join(", ", 
		      @{$fields_r->{classified_fields}}, 
		      $crRNA_fields ),
		), "\n";
    }

    # each hit
    foreach my $hit ( @{$iter->{Iteration_hits}{Hit}} ){
      next unless exists $hit->{Hit_hsps}{Hsp};

      # each hsp
      while( my($hspUID, $hsp) = each %{$hit->{Hit_hsps}{Hsp}} ){
	
	my @row; # output row	
	foreach my $field ( @{$fields_r->{classified_fields}} ){
	  # index to variable: [level, key]
	  push @row, getFieldValue( 
				   $field,
				   $fields_r->{index}{$field},
				   $blast_r,
				   $iter,
				   $hit,
				   $hsp
				  );
	}
       
       	## crRNA_info: if provided writing row for each crRNA entry
	if( exists $fields_r->{crRNA_fields}  ){
	  foreach my $crRNA_entry ( keys %{$iter->{crRNA_info}} ){  # each entry
	    my @crRNA_row;	    
	    # getting info for the entry
	    foreach my $field ( @{$fields_r->{crRNA_fields}} ){
	      push @crRNA_row, getcrRNAFieldValue(
						  $field,
						  $fields_r->{crRNA_index}{$field},
						  $iter->{crRNA_info}{$crRNA_entry}
						 );
	    }
	    print join("\t", @row, @crRNA_row), "\n";
	  }
	}
	## else: just 1 row
	else{	
	  print join("\t", @row), "\n";
	}
      }
    }
  }  
}


sub getcrRNAFieldValue{
  # input
  my $field = shift or confess "Provide field\n";
  my $index = shift or confess "Provide index\n";
  my $entry_r = shift or confess "Provide entry_r\n";


  if( exists $entry_r->{$index->[1]} ){
    return $entry_r->{$index->[1]};
  }
  else{
    warn "WARNING: Cannot find find field '$field'\n";
    return 'undef';
  }  
}


sub getFieldValue{
  # input
  my $field = shift or confess "Provide field\n";
  my $index = shift or confess "Provide index\n";
  my $blast_r = shift or confess "Provide blast_r\n";
  my $iter = shift or confess "Provide iter\n";
  my $hit = shift or confess "Provide hit\n";
  my $hsp = shift or confess "Porvide hsp\n";
	  
	  

  # 'undef' for variables not actually in blast data
  unless(defined $index->[1]){
    return 'undef';  
  }
  
  
  # adding field values to @row
  ## run
  if( $index->[0] eq 'run'){
    if( exists $blast_r->{$index->[1]} ){
      return $blast_r->{$index->[1]};
    }
    else{
      warn "WARNING: Cannot find field '$field'\n";
      return 'undef';
    }	    
  }
  
  ## iteration
  if( $index->[0] eq 'iter' ){
    if( exists $iter->{$index->[1]} ){
      return $iter->{$index->[1]};
    }
    else{
      warn "WARNING: Cannot find field '$field'\n";
      return'undef';
    }	    
  }
  
  ## hit
  if( $index->[0] eq 'hit' ){
    if( exists $hit->{$index->[1]} ){
      return $hit->{$index->[1]};
    }
    else{
      warn "WARNING: Cannot find field '$field'\n";
      return 'undef';
    }	    
  }
  
  ## hsp
  if( $index->[0] eq 'hsp' ){
    if( exists $hsp->{$index->[1]} ){
      # hsp -> frames field
      if (defined $index->[2]){
	return join("/", $hsp->{$index->[1]},
			$hsp->{$index->[2]});
      }
      # hsp_identity
      elsif( $index->[1] eq 'Hsp_identity' ){
	return $hsp->{Hsp_identity} / 
	  $hsp->{'Hsp_align-len'} * 100;
      }
      # other fields
      else{ 
	return $hsp->{$index->[1]};
      }	      
    }
    else{
      warn "WARNING: Cannot find field '$field'\n";
      return 'undef';
    }	    
  }
  
  ## param
  if( $index->[0] eq 'param' ){
    if( exists $blast_r->{BlastOutput_param}{Parameters}{$index->[1]} ){
      return $blast_r->{'BlastOutput_param'}{'Parameters'}{$index->[1]};
    }
    else{
      warn "WARNING: Cannot find field '$field'\n";
      return 'undef';
    }
  }

  # if made it this far, die
  confess "Internal error\n";
}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CRISPR_db


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

1; # End of CL_db
