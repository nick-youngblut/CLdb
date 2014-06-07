package CLdb::utilities;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use IPC::Cmd qw/run/;

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


=cut

=head2 get_file_path

Getting just path of file. Using File::Spec

=head3 IN

$file :  file name

=head3 OUT

$ :  file path

=cut

sub get_file_path{
  my ($file) = @_;
  confess "ERROR: no file provided!\n" unless defined $file;
  my @parts = File::Spec->splitpath(File::Spec->rel2abs($file));
  return $parts[1];
}


=head2 lineBreaks2unix

Converting line breaks to unix

=head3 IN

$fh :  file handle
$ext :  external file

=head3 OUT

=cut

sub lineBreaks2unix{
  my ($fh, $ext) = @_;
  
  if($ext){			# external file 
    confess "ERROR: cannot find $fh!\n" unless -e $fh || -l $fh;
    # mac to unix #
    my $cmd = "perl -pi -e 's/\\r/\\n/g' $fh";
    my $ret_r = run( command => $cmd, verbose => 0 );
    `$cmd`;
    if ($? == -1) {  confess "ERROR: failed to execute: $!\n"; }
    # win to unix #
    $cmd = "perl -pi -e 's/\\r\$//g' $fh";
    `$cmd`;					
    if ($? == -1) {  confess "ERROR: failed to execute: $!\n"; }
  }
  else{				# internal; filehandle
    my @table = <$fh>;
    my @tmp;
    map{ s/\r$//; s/\r/\n/g; s/\n+/\n/g; push @tmp, split /\n/;  } @table;
    @table = @tmp;
    return \@table;
  }
}


=head2 connect2db

Connecting to sqlite3 database

=head3 IN

$db_file :  sqlite3 database file name

=head3 OUT

$ :  dbh object

=cut

sub connect2db{
  my $db_file = shift or confess "Provide a database file name";
  
  my %attr = (RaiseError => 0, PrintError=>1, AutoCommit=>0 );
  my $dbh = DBI->connect("dbi:SQLite:dbname=$db_file", '','', \%attr) 
    or confess " Can't connect to $db_file!\n";
  
  return $dbh;
}


=head2 file_exists

Checking to see if file name was provided and if it exists.
Die if it does not.

=head3 IN

$file :  file name 
$cat  :  Alias for file name (eg., 'database'). Used for error reporting

=head3 OUT

=cut

sub file_exists{
  my ($file, $cat) = @_;
  
  $cat = "" unless defined $cat;
  
  confess "ERROR: provide a $cat file name!\n"
    unless defined $file;
  confess "ERROR: cannot find $cat file: '$file'\n"
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
