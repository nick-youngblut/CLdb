package CLdb::blast;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

# export #
use base 'Exporter';
our @EXPORT_OK = qw/
parse_blast_hits
read_blast_file
write_blast_file
/;

	
=head1 NAME

CLdb::blast

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

=head1 EXPORT_OK

file_exists
connect2db


=cut

sub write_blast_file{
# writing a blast file read in format from read_blast_file subroutine #
	my ($lines_r) = @_;
	
	foreach my $query ( sort keys %$lines_r ){
		foreach my $db ( sort keys %{$lines_r->{$query}} ){
			foreach my $blast ( keys %{$lines_r->{$query}{$db}} ){
				print $blast, "\n";
				print $query, "\n";
				print $db, "\n";
				print $lines_r->{$query}{$db}{$blast}{'fields'}, "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'fields'};
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'comments'}}), "\n";
				# printing all hits #
				print join("\n", @{$lines_r->{$query}{$db}{$blast}{'hits'}}), "\n"
					if exists $lines_r->{$query}{$db}{$blast}{'hits'};
				}
			}
		}
	}

sub read_blast_file{
# reading in each blast entry (-outfmt 7) & extracting names and line numbers #
	my %lines;
	my $blast;
	my $query;
	my $db;
	while(<>){
		chomp;
		
		if(/^# BLAST/i){
			$blast = $_;
			}
		elsif(/^# Query/i){
			$query = $_;
			}
		elsif(/^# Database/i){
			$db = $_;
			}
		elsif(/# Fields/i){
			$lines{$query}{$db}{$blast}{'fields'} = $_;
			
			my @l = split /[:,] /;
			shift @l;
			for my $i (0..$#l){
				$lines{$query}{$db}{$blast}{'fields_sep'}{$l[$i]} = $i;
				}
			}
		elsif(/^# /){
			push @{$lines{$query}{$db}{$blast}{'comments'}}, $_;
			}
		else{
			push @{$lines{$query}{$db}{$blast}{'hits'}}, $_;
			}
			
		}		
		#print Dumper %lines; exit;	
	return \%lines;
	}

sub parse_blast_hits{
# REQUIRED: outfmt = 7;
	my ($blast_in) = @_;
	
	my %blast;
	my %fields;
	my @header;
	my ($db, $query);
	while(<>){
		chomp;		
		next if /^\s*$/;		# skipping blank lines
		if(/^#/){ 
			@header = () if /^# BLASTN/;		# start of next set of comment lines 
			push @header, $_; 
			}
		else{
			if(/^DR_/){
				@header = () if @header;
				next;
				}
			
			if(@header){ 				# need to parse header #
				# parsing header #
				($db, $query, %fields) = parse_header(\@header, \%blast, \%fields);
				my @tmp = @header;
				$blast{$query}{$db}{"header"} = \@tmp;
				
				# reseting header #
				@header = ();
				}
			
			my @line = split /\t/;	
			
			push @{$blast{$query}{$db}{"hits"}}, [split /\t/];
			}

		}
		
		confess " ERROR: not blast hits found in input file!\n" unless %blast;
	return \%blast, \%fields;
	}

sub parse_header{
# parsing each header of blast output #
	my ($header_r, $blast_r, $fields_r) = @_;
	
	my ($db, $query);
	foreach (@$header_r){
		if(/^# Database:/){
			($db = $_) =~ s/^# Database: //; 
			
			}
		elsif(/^# Fields:/){
			$fields_r = fields2header($_) unless %$fields_r;
			check_fields($fields_r);
			}
		elsif(/^# Query:/){
			($query = $_) =~ s/^# Query: //; 
			}
		}
	return $db, $query, %$fields_r;
	}
	
sub fields2header{
# making header hash from fields #
	my ($line) = @_;
	
	$line =~ s/# Fields: //;
	my @l = split /\s*,\s*/, $line;		
	
	my %header;
	for my $i (0..$#l){
		$header{$l[$i]} = $i;
		}
		
		#print Dumper %header; exit;
	return \%header;
	}
	
sub check_fields{
# checking to make sure that the correct fields are available #
	my ($header_r) = @_;
	
	my @needed = ("query id", "subject id", "s. start", "s. end",
				"q. start", "q. end", 
				"evalue", "query length", "subject length", "BTOP");
	
	map{confess " ERROR: can't find '$_'!" unless exists $header_r->{$_}} @needed;
	}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::blast


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
