#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::DB::GenBank;


### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $subject_file, $database_file);
GetOptions(
	   "database=s" => \$database_file,
	   "subject=s" => \$subject_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: cannot find $subject_file!\n"
	if $subject_file && ! -e $subject_file;


### MAIN
my $blast_r = load_blast();
query_genbank($blast_r);

### Subroutines
sub load_blast{
	my @blast;
	while(<>){
		chomp;
		next if /^\s*$/;
		push @blast, [split /\t/];
		}
	return \@blast;
	}

sub query_genbank{
	my ($blast_r) = @_;
	my $gb = new Bio::DB::GenBank;

	foreach my $line (@$blast_r){

		#print Dumper $$line[$#$line]; exit;
		
		my $seqio = $gb->get_Stream_by_acc( [$$line[$#$line]] );

		while(my $seq = $seqio->next_seq){
			# seqfeature #
			my @line = (
				$seq->primary_seq->desc,
				$seq->primary_seq->seq,
				);
			print join("\t", @line), "\n"; exit;
			#print Dumper $seq; exit;
			}
		}
	}




__END__

=pod

=head1 NAME

CLdb_loadSpacerBlast.pl -- Parse blast output (xml format)

=head1 SYNOPSIS

CLdb_loadSpacerBlast.pl [options] < file.blast.xml > file.blast.txt

=head2 options

=over

=item -taxa_summary 		

Print a summary of taxa hit

=item -evalue 			

e-value cutoff (examples: '10' or '1e-5'). [10]

=item -hsp_length 		

Hit length (gaps included)

=item -query_length 		

Query length (gaps included)

=item -percent_identity 	

Percent identity (values from 0-100)

=item -bit_score 		

Bit score

=item -header 			

Print header? [TRUE]

=item -sequence

Subject sequence column? [FALSE]

=item -type_strain

Just type strain hit (no 'uncultured')? [FALSE

=item -help			

This help message

=back

=head2 For more information:

perldoc CLdb_loadSpacerBlast.pl

=head1 DESCRIPTION

Parse the xml output from blast+ (and probably blast).

A summary of taxa hit can also be produced. The summary can be grouped by different categories.

=head1 EXAMPLES

=head2 Parse blast output

CLdb_loadSpacerBlast.pl < file.blast.xml > file.blast.txt

=head2 Make a taxa-hit summary count (grouping by hitID & queryID)

CLdb_loadSpacerBlast.pl -taxa 1 < file.blast.xml > file.taxa-hit.txt

=head2 Make a taxa-hit summary stats (grouping by hitID)

CLdb_loadSpacerBlast.pl -taxa 2 < file.blast.xml > file.taxa-hit.txt

=head2 Make a taxa-hit summary stats (grouping by hitID & queryID)

CLdb_loadSpacerBlast.pl -taxa 3 < file.blast.xml > file.taxa-hit.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/NY_misc_perl/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

