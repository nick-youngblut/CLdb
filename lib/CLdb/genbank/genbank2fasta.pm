package CLdb::genbank::genbank2fasta;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use List::Util qw/max/;

# export #
use base 'Exporter';
our @EXPORT_OK = '';

	

=head1 NAME

CLdb::genbank::genbank2fasta

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

=head1 EXPORT_OK

=cut

=head1 SUBROUTINES


=head2 genbank2fasta

getting fasta from genbank unless genome fasta file
already exists

=head3 IN 

$loci_r :  hash_ref of loci table
$db_path :  CLdb.sqlite path
$header_r :  loci header index

=head3 OUT


=cut 

push @EXPORT_OK, 'genbank2fasta';

sub genbank2fasta{
# getting fasta from genbank unless -e fasta #
  my ($loci_r, $db_path, $header_r) = @_;

  print STDERR "### Checking for existence of genome fasta ###\n";

  foreach my $locus_id (keys %$loci_r){
    print STDERR "### Processing locus: \"$locus_id\" ###\n";

    if(! exists $loci_r->{$locus_id}{'fasta_file'} || $loci_r->{$locus_id}{'fasta_file'} =~ /^\s*$/ ){
      print STDERR "  No genome fasta for locus_id: \"$locus_id\"! Trying to extract sequence from genbank...\n";
      $loci_r->{$locus_id}{'fasta_file'} =
        genbank2fasta_extract($loci_r->{$locus_id}{'genbank_file'}, 
			      $db_path, "$db_path/fasta/", $locus_id);
    }
    elsif( ! -e "$db_path/fasta/$loci_r->{$locus_id}{'fasta_file'}"){
      my $fasta_name = $loci_r->{$locus_id}{'fasta_file'};
      print STDERR "  WARNING: Cannot find $fasta_name! Trying to extract sequence from genbank...\n";
      $loci_r->{$locus_id}{'fasta_file'} =
        genbank2fasta_extract($loci_r->{$locus_id}{'genbank_file'}, 
			      $db_path, "$db_path/fasta/", $locus_id);
    }
    else{ 
      printf STDERR " The fasta_file field %s was provided for locus %s\n",
	    $loci_r->{$locus_id}{'fasta_file'}, $locus_id; 
    }
  }

  $header_r->{'fasta_file'} = max(values %$header_r) + 1
    unless exists $header_r->{'fasta_file'};
}


=head2 genbank2fasta_extract

Extracting genome fasta sequence from genbank.
Used by genbank2fasta subroutine

=head3 IN

$genbank_file :  genbank file name

=head3 OUT

=cut

sub genbank2fasta_extract{
  my ($genbank_file, $db_path, $fasta_dir, $locus_id) = @_;

  # genbank file must exist in order to extract fasta from it
  unless( defined $genbank_file ){
    print STDERR " WARNING: no genbank file name provided for locus '$locus_id'.";
    print STDERR " Cannot extract the genome sequence to make a fasta. Skipping\n";
    return undef;
  }            # cannot do w/out genbank


  # making fasta dir if not present #
  mkdir $fasta_dir unless -d $fasta_dir;
  
  # checking for existence of fasta #
  my @parts = File::Spec->splitpath($genbank_file);
  $parts[2] =~ s/\.[^.]+$|$/.fasta/;
  my $fasta_out = "$fasta_dir/$parts[2]";
  if(-e $fasta_out){
    print STDERR "  Success! Found fasta in '$fasta_out'\n";
    print STDERR "  Adding fasta to loci table\n";
    return $parts[2];
  }
  
  # sanity check #
  croak " ERROR: cannot find $genbank_file!\n"
    unless -e "$db_path/genbank/$genbank_file";
  
  # I/O #
  my $seqio = Bio::SeqIO->new(-file => "$db_path/genbank/$genbank_file", -format => "genbank");
  open OUT, ">$fasta_out" or die $!;
  
  # writing fasta #
  my $seq_cnt = 0;
  while(my $seqo = $seqio->next_seq){
    $seq_cnt++;
    
    # seqID #
    my $scafID = $seqo->display_id;
    print OUT join("\n", ">$scafID", $seqo->seq), "\n";
  }
  close OUT;
  
  # if genome seq found and fasta written, return fasta #
  if($seq_cnt == 0){
    print STDERR "\tWARNING: no genome sequnece found in Genbank file: $genbank_file!\nSkipping BLAST!\n";
    unlink $fasta_out;
    return 0;
  }
  else{
    print STDERR "   Success! fasta file extracted from $genbank_file.\n   File written to: $fasta_out\n";
    return $parts[2];                   # just fasta name
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

    perldoc CLdb::genbank::genbank2fasta


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
