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


sub seq_from_genome_fasta{
#-- Description --#
# extracting leader sequence from genome fasta
#-- Input --#
# $fasta_r = %$; {seq_name}=>seq
# $coords_r = @$; [scaf, start, end]
# $db_path = file path

	my ($fasta_r, $coords_r) = @_;
	
	# IO check #
	die "ERROR: coords must be [scaffold,start,end]\n"
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

sub read_fasta{
# loading fasta file as a hash #
	my $fasta_in = shift;
	die " ERROR: cannot find $fasta_in!" unless -e $fasta_in || -l $fasta_in;

	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
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
	close IN;
		#print Dumper %fasta; exit;
	return \%fasta;
	} 

sub revcomp{
# reverse complement of a sequence
# caps invariant 
	my $seq = shift;
	$seq = reverse($seq);
	$seq =~ tr/[a-z]/[A-Z]/;
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
