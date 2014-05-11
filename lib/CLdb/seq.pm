package CLdb::seq;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use DBI;

# export #
use base 'Exporter';
our @EXPORT_OK = qw/
revcomp
read_fasta
seq_from_genome_fasta
/;

	
=head1 NAME

CLdb::seq

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';

=head1 SYNOPSIS

Subroutines for editing sequence data

=head1 EXPORT_OK

revcomp

=cut

=head2 write_fasta

writing out fasta to file or file handle

=head3 IN

fasta => hashref %{name}=>seq
file => file name  ('>file_name' or other write operator)
fh => file handle 

=head3 OUT

=cut

push @EXPORT_OK, 'write_fasta';

sub write_fasta{
  my %h = @_;
  
  # I/O
  croak "ERROR: provide a fasta as hash_ref"
    unless exists $h{fasta};
  my $outfh;
  if(exists $h{file}){
    open $outfh, $h{file} or croak 
      "ERROR: could not write to $h{'file'}";
  }
  elsif(exists $h{fh}){ $outfh = $h{fh}; }
  else{ croak "ERROR: provide a file name or file handle"; }

  # writing
  foreach my $name (keys %{$h{fasta}}){
    $name =~ s/^>//;
    print $outfh join("\n", ">$name", $h{fasta}{$name} ), "\n";
  }
  close $outfh;
}



=head2 seq_from_genome_fasta

=cut

sub seq_from_genome_fasta{
#-- Description --#
# extracting leader sequence from genome fasta
#-- Input --#
# $fasta_r = %$; {seq_name}=>seq
# $coords_r = @$; [scaf, start, end]
# $db_path = file path

	my ($fasta_r, $coords_r) = @_;
	
	# IO check #
	confess "ERROR: coords must be [scaffold,start,end]\n"
		unless scalar @$coords_r == 3;
	
	# getting fasta #
	#my $fasta_r = load_fasta("$db_path/fasta/$locus_r->{'fasta_file'}");
	
	# checking for existence of scaffold #
	my $scaf = $$coords_r[0];
	unless(exists $fasta_r->{$scaf}){
		#print STDERR "WARNING: '$scaf' not found . Not loading leader sequence!\n";
		return "";			# no sequence
		}
		
	# start-stop #
	my $start = $$coords_r[1];
	my $end = $$coords_r[2];
	
	my $leader_seq;
	if($start < $end){		# neg strand
		$leader_seq = substr($fasta_r->{$scaf}, $start -1, $end-$start+1);
		}
	else{
		$leader_seq = substr($fasta_r->{$scaf}, $end -1, $start-$end+1);
		$leader_seq = revcomp($leader_seq);
		}
		
		#print Dumper $leader_seq; exit;
	return $leader_seq;
	}


=head2 read_fasta

Reading in fasta file

=head3 IN

hash of args:
file :  file name
fh  : file handle

=head3 OUT

$%{name}=>seq

=cut

sub read_fasta{
# loading fasta file as a hash #
  my %h = @_;
  confess "ERROR: a file name or file handle must be provided\n"
    unless exists $h{file} or $h{fh};

  # file or file handle
  my $fh;
  exists $h{fh} ? $fh = $h{fh} : open $fh, $h{file} or confess $!;

  my (%fasta, $tmpkey);
  while(<$fh>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;	

    if(/>.+/){
      s/>//;
      $fasta{$_} = "";
      $tmpkey = $_;	# changing key
    }
    else{$fasta{$tmpkey} .= $_; }
  }
  close $fh;
  return \%fasta;
} 


sub revcomp{
# reverse complement of a sequence
# caps invariant 
	my $seq = shift;
	$seq = reverse($seq);
	#$seq =~ tr/[a-z]/[A-Z]/;
	$seq =~ tr/ACGTNBVDHKMRYSWacgtnbvdhkmrysw\.-/TGCANVBHDMKYRSWtgcanvbhdmkyrsw\.-/;
	return $seq;
	}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::seq


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
