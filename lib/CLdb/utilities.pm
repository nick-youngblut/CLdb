package CLdb::utilities;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

# export #
use base 'Exporter';
our @EXPORT_OK = qw/
file_exists
connect2db
lineBreaks2unix
get_file_path
/;

	
=head1 NAME

CLdb::utilities

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Utility subroutines

=head1 EXPORT_OK

file_exists
connect2db


=cut

sub get_file_path{
#-- Description --#
# converting line breaks to unix #
#-- input --#
# $file = file name 
	my ($file) = @_;
	die "ERROR: no file provided!\n" unless defined $file;
	my @parts = File::Spec->splitpath(File::Spec->rel2abs($file));
	return $parts[1];
	}

sub lineBreaks2unix{
#-- Description --#
# converting line breaks to unix #
#-- input --#
# $fh = filehandle or file name
# $ext = external file? 
	my ($fh, $ext) = @_;
	
	if($ext){			# external file 
		die "ERROR: cannot find $fh!\n" unless -e $fh || -l $fh;
		# mac to unix #
		my $cmd = "perl -pi -e 's/\\r/\\n/g' $fh";
		`$cmd`;
		# win to unix #
		$cmd = "perl -pi -e 's/\\r\$//g' $fh";
		`$cmd`;					
		}
	else{				# internal; filehandle
		my @table = <$fh>;
		my @tmp;
		map{ s/\r$//; s/\r/\n/g; s/\n+/\n/g; push @tmp, split /\n/;  } @table;
		@table = @tmp;
		return \@table;
		}
	}

sub connect2db{
# connecting to CLdb 
# $db_file = database file
	my $db_file = shift;

	my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
	my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
		or die " Can't connect to $db_file!\n";
		
	return $dbh;
	}

sub file_exists{
# checking to see if file name was provided and if it exists #
# $file = file name
# $cat = type of file (eg. "database")
	my ($file, $cat) = @_;
	
	$cat = "" unless defined $cat;
	
	die "ERROR: provide a $cat file name!\n"
		unless defined $file;
	die "ERROR: cannot find $cat file: '$file'\n"
		unless -e $file || -l $file;

	}

=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::utilities


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
