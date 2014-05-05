package CLdb::arrayBlast::sereal;

use 5.006;
use strict;
use warnings FATAL => 'all';
use Data::Dumper;
use Carp qw/carp croak/;

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
        make_hsp_array($blast_r, $iter, $hit);
      }
    }
    else{ # making array
      make_hsp_array($blast_r, $iter, $hits_ref);
      $iter->{'Iteration_hits'}{'Hit'} = [$hits_ref]
    }

    # make blast output row of field values
    sub make_hsp_array{
      my ($blast_r, $iter, $hit) = @_;

      # multiple hsp or just 1?
      if( ref $hit->{Hit_hsps}{Hsp} eq 'ARRAY' ){
	# nothing needed
      }
      else{ # making array
	$hit->{Hit_hsps}{Hsp} = [$hit->{Hit_hsps}{Hsp}];
      }
    }
  }
 #print Dumper $blast_r; exit;
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

  # fields
  $fields{fields} = [@l[1..$#l]];
 
 # print Dumper %fields; exit;
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

=head3 Index

# header
blast: 'BlastOutput_version'
Query: 'BlastOutput_query-def'
Database: 'BlastOutput_db'
Fields: user defined
Hits: 'BlastOutput_iterations' => 'Iteration' => scalar []

# body
qseqid: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_query-def' 
qgi: NA
qacc: NA
qaccver: NA
qlen: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_query-len' 
sseqid: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_id'
sallseqid: NA
sgi: NA
sallgi: NA
sacc: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_accession'
accver: NA
sallacc: NA
slen: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_len'
qstart: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_query-from'
qend: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_query-to'
sstart: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_hit-from'
send: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_hit-to'
qseq: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_qseq'
sseq: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_hseq'
evalue: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_evalue'
bitscore: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_bit-score'
score: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_score'
length: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_align-len'
pident: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_identity'
nident: NA
mismatch: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => ?
positive: NA
gapopen: 'BlastOutput_iterations' => 'Parameters' => 'Parameters_gap-open'
gaps: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_gaps'
ppos: NA
frames: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_query-frame' /
'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_hit-frame'
qframe: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_query-frame'
sframe: 'BlastOutput_iterations' => 'Iteration' => [] => 'Iteration_hits' => 'Hit' => [] => 'Hit_hsp' => 'Hsp' => 'Hsp_hit-frame'
btop: NA
staxids: NA
sscinames: NA
scomnames: NA
sblastnames: NA
sskingdoms: NA
stitle: NA
salltitles: NA
sstrand: NA
qcovs: NA
qcovhsp: NA

=cut

push @EXPORT_OK, 'classify_fields';

sub classify_fields {
  my ($fields_r) = @_;

  # classification table
  my %class = (
	       qseqid => ['run', 'Iteration_query-def'],
	       qlen => ['run', 'Iteration_query-len'], 
	       sseqid => ['hit', 'Hit_id'],
	       sacc => ['hit', 'Hit_accession'],
	       slen => ['hit', 'Hit_len'],
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
	       gapopen => ['param', 'Parameters_gap-open']
	      );


  # classification
  foreach my $field ( @{$fields_r->{fields}} ){
    (my $tmp = $field) =~ tr/A-Z/a-z/;

    if( exists $class{$tmp} ){
      $fields_r->{index}{ $tmp } = $class{$tmp};
    }
    else{
      die "ERROR: field '$tmp' not supported\n";
    }    
  }

 # print Dumper $fields_r; exit;
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
    my $hits_ref = $iter->{'Iteration_hits'}{'Hit'}; # hits in iteration (array_ref or ref to hash)

    # comments for the query (if needed)
    if( $fields_r->{comments} == 1){
      print join(" ", '# blast:', $blast_r->{'BlastOutput_version'}), "\n";
      print join(" ", '# Query:', $blast_r->{'BlastOutput_query-def'}), "\n";
      print join(" ", '# Database:', $blast_r->{'BlastOutput_db'}), "\n";
      print join(" ", '# Fields:', 
		 join(", ", @{$fields_r->{fields}})
		), "\n";
    }
    
    # multiple hits or just 1?
    if(ref $hits_ref eq 'ARRAY'){ # multiple
      foreach my $hit ( @$hits_ref){
	make_blast_row($blast_r, $iter, $hit, $fields_r);
      }
    }
    else{ 
      make_blast_row($blast_r, $iter, $hits_ref, $fields_r);
    }
    
    # make blast output row of field values
    sub make_blast_row{
      my ($blast_r, $iter, $hit, $fields_r) = @_;

      # parsing for each hsp
      if( ref $hit->{Hit_hsps}{Hsp} eq 'ARRAY' ){
	foreach my $hsp (@{$hit->{Hit_hsps}{Hsp}}){
	  row_by_hsp($blast_r, $iter, $hit, 
				$hsp, $fields_r);
	}
      }
      else{
	row_by_hsp($blast_r, $iter, $hit, 
			      $hit->{Hit_hsps}{Hsp}, $fields_r);
      }

      sub row_by_hsp{
	my ($blast_r, $iter, $hit, $hsp, $fields_r) = @_;
	
	my @row; # output row
	foreach my $field ( @{$fields_r->{fields}} ){
	  my $index = $fields_r->{index}{$field};
	  
	  # if variable undefined
	  unless(defined $index->[1]){
	    push @row, 0;  # really, is 'undef'
	    next;
	  }
	  
	  # determing path to variable 
	  ## run
	  if( $index->[0] eq 'run' ){
	    push @row, $iter->{$index->[1]};
	    
	  }
	  ## hit
	  if( $index->[0] eq 'hit' ){
	    push @row, $hit->{$index->[1]};
	  }
	  ## hsp
	  if( $index->[0] eq 'hsp' ){
	    # hsp -> frames field
	    if (defined $index->[2]){
	      push @row, join("/", $hsp->{$index->[1]},
			  $hsp->{$index->[2]});
	    }
	    elsif( $index->[1] eq 'Hsp_identity' ){
	      push @row, $hsp->{Hsp_identity} / 
		$hsp->{'Hsp_align-len'} * 100;
	    }
	    else{ # other fields
	      push @row, $hsp->{$index->[1]};
	    }
	  }
	  ## param
	  if( $index->[0] eq 'param' ){
	    push @row, $blast_r->{'BlastOutput_param'}{'Parameters'}{$index->[1]};
	  }
	}
	print join("\t", @row), "\n";
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
