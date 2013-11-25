package CLdb::load;

# module use #
## core ##
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use DBI;

## CLdb ##


# export #
use base 'Exporter';
our @EXPORT_OK = qw/
load_db_table
/;

	
=head1 NAME

CLdb::load

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for querying CLdb

=head1 EXPORT_OK

=cut

sub load_db_table{
#-- description --#
# adding table entries; assuming redundant entries will replace or be ignored #
#-- input --#
# dbh = DBI object
# table = table to update
# $header_r = %$; all columns to add values to; ${column_name} => column_index
# $vals_r = %@$;	all entries

	my ($dbh, $table, $header_r, $vals_r, $vb) = @_;

	# sanity check #
	croak "ERROR: no categories provided!\n"
		unless defined $header_r || scalar keys %$header_r > 0;

	# header order #
	my @header = keys %$header_r;
	
	# making query #
	my $l = join(",", @header);
	my $Qmrk = join(",", ("?") x scalar @header);
	

	my $q = "INSERT INTO $table($l) values ($Qmrk)";	
	my $sql = $dbh->prepare($q);
	
	my $cnt = 0;
	foreach my $entry (keys %$vals_r){
		$sql->execute(@{$vals_r->{$entry}}[@{$header_r}{@header}]);
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr for $entry\n";
			}
		else{ $cnt++; }
		}
	$dbh->commit;

	print STDERR "...Number of entries added/updated to '$table' table: $cnt\n"
		unless $vb;
	}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::load


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
