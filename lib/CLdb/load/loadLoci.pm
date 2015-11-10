package CLdb::load::loadLoci;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use Parallel::ForkManager;
use IPC::Cmd qw/run/;

## export #
use base 'Exporter';
our @EXPORT_OK = '';

## CLdb
use CLdb::utilities qw/
			lineBreaks2unix
		      /;
use CLdb::seq qw/read_fasta
		 seq_from_genome_fasta/;


=head1 NAME

CLdb::load::loadLoci

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing the tab-delimited loci 
file and loading the entries into CLdb

=head1 EXPORT_OK

=cut

=head1 SUBROUTINES


=head2 unix_line_breaks

=head3 IN

$locis_r :  hash_ref of loci table
$db_path :  database path

=head3 OUT

=cut 

push @EXPORT_OK, 'unix_line_breaks';

sub unix_line_breaks{
  my $loci_r = shift || confess "Provide loci hashref\n";
  my $db_path = shift || confess "Provide database path\n";
  my $forks = shift;
  $forks = 0 unless defined $forks;

  print STDERR "### checking line breaks for all external files (converting to unix) ###\n";

  #my @file_columns = qw/genbank_file fasta_file array_file/;
  my %file_cols = (
                   "genbank_file" => "genbank",
                   "fasta_file" => "fasta",
                   "array_file" => "array");

  # fork object
  my $pm = Parallel::ForkManager->new($forks);

  # processing each file
  foreach my $locus_id (keys %$loci_r){
    foreach my $file_col (keys %file_cols){
      if(exists $loci_r->{$locus_id}{$file_col}){
        next unless $loci_r->{$locus_id}{$file_col};                            # if no file; nothing to check

	$pm->start and next;

	my $dirpath = File::Spec->catdir($db_path, $file_cols{$file_col});
	my $filepath = File::Spec->catdir($dirpath, $loci_r->{$locus_id}{$file_col});

        print STDERR " processing: $filepath\n";
#        lineBreaks2unix($filepath, 1);

	# system call for command line perl (hopefully, will prevent file truncations upon crashing)
	my $cmd = "perl -pe 's/\r/\n/g; s/\r\$//g' $filepath";
	my( $success, $error_message, $full_buf, 
	    $stdout_buf, $stderr_buf ) =
	      run( command => $cmd, verbose => 0 );
	die "System call failed: '$cmd'" unless $success;

	$pm->finish;
      }
    }
  }
  $pm->wait_all_children;
}


=head2 just_table_columns

Loading just the columns (fields) needed for CLdb

=cut

push @EXPORT_OK, 'just_table_columns';

sub just_table_columns{
#-- Description --#
# loading just the columns of interest #
  my ($header_r, $table) = @_;
  $table =~ tr/A-Z/a-z/;

  my %col_lists;
  @{$col_lists{'loci'}} = qw/Locus_ID Taxon_ID Taxon_Name
                             Subtype Scaffold
                             Locus_start Locus_end
                             CAS_Start CAS_End
                             Array_Start Array_End
                             array_sense_strand
                             CAS_Status Array_Status
                             Genbank_File Fasta_File Array_File
                             Scaffold_count
                             File_Creation_Date Author/;

  @{$col_lists{'leader'}} = qw/Locus_ID Leader_start Leader_end Leader_end Leader_sequence/;
  @{$col_lists{'pam'}} = qw/Locus_ID PAM_seq PAM_start PAM_end PAM_sequence/;

  die "ERROR: '$table' not found in column lists!\n"
    unless exists $col_lists{$table};

  map{ $_ =~ tr/A-Z/a-z/ } @{$col_lists{$table}};               # table of interest

  my %header_parse;
  map{$header_parse{$_} = $header_r->{$_} if exists $header_r->{$_}} @{$col_lists{$table}};

  #print Dumper %header_parse; exit;
  return \%header_parse;
}


=head2 make_external_file_dirs

Moving files to directory in $CLdb_home

=cut

push @EXPORT_OK, 'make_external_file_dirs';

sub make_external_file_dirs{
  my ($loci_r, $header_r, $dir, $quiet) = @_;

  print STDERR "### Checking/copying external files specfied in loci table ###\n";

  my %fileTypes = (genbank_file => 'genbank',
		   fasta_file => 'fasta',
		   array_file => 'array');


  foreach my $locus_id (keys %$loci_r){
    print STDERR "Processing locus '$locus_id'...\n";
    while (my ($fileType, $ext_dir) = each %fileTypes){
      # next unless file exists in loci table
      my $fileName = exists $loci_r->{$locus_id}{$fileType} ?
	$loci_r->{$locus_id}{$fileType} : next;
      next if $fileName eq ''; 
      
      # making directory to store files unless exists
      $ext_dir = make_CLdb_dir($dir, $ext_dir);
      
      # copying external files to ext_dir if not in ext_dir
      $loci_r->{$locus_id}{$fileType} = copy_external_files($fileName, $ext_dir, $fileType);
    }
  }
}


=head2 make_CLdb_dir

making directory in $CLdb_home

Args:
dir -- name of CLdb_home dir
name -- file of directory to store external files

=cut

sub make_CLdb_dir{
  my ($dir, $name) = @_;

  my $newdir = File::Spec->catdir($dir, $name);
  mkdir $newdir unless -d $newdir;
  
  return $newdir;
}


=head2 copy_ext_files

identifying location of files (if found) and copying if needed

Args:
fileName -- name of file as specified in locus table
ext_dir -- name of exteranl file directory
fileType -- file type checking on

Return:
str -- basename of existing file || undef

=cut

sub copy_external_files{
  my $fileName = shift or confess $!;
  my $ext_dir = shift or confess $!;
  my $fileType = shift or confess $!;

  # file basenamne
  my $baseName = (File::Spec->splitpath($fileName))[2];
  my $fileInExt = File::Spec->catfile($ext_dir, $baseName);


  # check for file in external directory  
  print STDERR "  Checking: '$fileInExt'\n";
  if (-e $fileInExt){
    print STDERR "\tFile found: '$fileInExt'\n";
    return $baseName; #$fileInExt;
  }
  else{
    print STDERR "\tFile NOT found: '$fileInExt'\n";
  }
  
  # check $fileInExt
  print STDERR "  Checking: '$fileName'\n";
  if( -e $fileName){   # checking for file in loci_table-specified location
    print STDERR "\t'$baseName' not in $ext_dir of \$CLdb_home. Copying file from '$fileName'.\n";
    copy($fileName, $fileInExt) or die $!;    
    return $baseName; #$fileInExt;
  }
  else{
    print STDERR "\tFile NOT found: '$fileName'\n";
  }

  # checking basename location
  print STDERR "  Checking: '$baseName'\n";
  if( -e $baseName){
    print STDERR "\t'$baseName' not in $ext_dir of \$CLdb_home. Copying file from '$baseName'.\n";
    copy($baseName, $fileInExt) or die $!;    
    return $baseName; #$fileInExt;
  }
  else{
    print STDERR "\tFile NOT found: '$baseName'\n";    
  }

  # file cannot be found; return undef
  if($fileType eq 'genbank_file'){  # genbank required
    die "ERROR: cannot find required genbank file: '$baseName'\n";
  }
  else{
    print STDERR  "\tWARNING: Cannot find '$baseName' in any location. Skipping.\n";
  }
  return undef;  
}


=head2 get_leader_seq

Loading leader. Pulling out sequence from genome if available

=cut

push @EXPORT_OK, 'get_leader_seq';
sub get_leader_seq{
# loading leader; pulling out sequence from genome if available #
  my ($loci_r, $db_path) = @_;

  # status #
  print STDERR "### Leader start-end provided. Loading values into leader table ###\n";

#  print Dumper $loci_r; exit;

  # getting leader sequence if possible #
  my %cp;
  foreach my $locus_id (keys %$loci_r){
    # next unless leader_start && leader_end #
    next unless exists $loci_r->{$locus_id}{'leader_start'}
      && exists $loci_r->{$locus_id}{'leader_end'}
        && exists $loci_r->{$locus_id}{'scaffold'}
          && exists $loci_r->{$locus_id}{'fasta_file'};

    # checking for existence of genome fasta, if yes, extract leader sequence #
    if(exists $loci_r->{$locus_id}{'fasta_file'}){
      my $catdir = File::Spec->catdir( $db_path, "fasta" );
      my $catfile = File::Spec->catfile( $catdir, $loci_r->{$locus_id}{fasta_file} );
      my $fasta_r = read_fasta(-file => $catfile);

      $loci_r->{$locus_id}{'leader_sequence'} =
        seq_from_genome_fasta( $fasta_r,
                               [$loci_r->{$locus_id}{'scaffold'},
                                $loci_r->{$locus_id}{'leader_start'},
                                $loci_r->{$locus_id}{'leader_end'}]
                             ) unless exists $loci_r->{$locus_id}{'leader_sequence'};

      if($loci_r->{$locus_id}{'leader_sequence'} eq ""){
        print STDERR "WARNING: '", $loci_r->{$locus_id}{'scaffold'},
          "' not found for ", $loci_r->{$locus_id}{'fasta_file'},
            ". Not loading leader sequence!\n";
      }
      else{             # just entries w/ leader sequence #
        $cp{$locus_id} = $loci_r->{$locus_id};
      }
    }
  }

  #print Dumper %cp; exit;
  return \%cp;
}


=head2 get_pam_seq

Getting pam sequence if location provided.

=cut

push @EXPORT_OK, 'get_pam_seq';
sub get_pam_seq{
  my ($loci_r, $db_path, $pam_header_r) = @_;

  if(exists $pam_header_r->{"pam_start"} && exists $pam_header_r->{"pam_end"}
     && exists $pam_header_r->{"pam_sequence"}){
    print STDERR "### PAM sequence & start-end columns provided. Getting PAM sequence if needed. Loading table ###\n"
  }
  elsif(exists $pam_header_r->{"pam_start"} && exists $pam_header_r->{"pam_end"}){
    print STDERR "### PAM start-end columns provided. Getting PAM sequence if needed. Loading table ###\n"
  }
  elsif(exists $pam_header_r->{"pam_sequence"}){
    print STDERR "### PAM sequence column provided. Loading values into PAM table ###\n"
  }
  else{
    print STDERR "### no PAM info provided. Skipping PAM loading ###\n"
  }

  # getting pam sequence if possible #
  my %cp;
  foreach my $locus_id (keys %$loci_r){
    if(exists $loci_r->{$locus_id}{'pam_sequence'}){
      $cp{$locus_id} = $loci_r->{$locus_id};
    }
    elsif( exists $loci_r->{$locus_id}{'pam_start'}
           && exists $loci_r->{$locus_id}{'pam_start'}
           && exists $loci_r->{$locus_id}{'scaffold'}
           && exists $loci_r->{$locus_id}{'fasta_file'}){               # geting sequence
      my $catdir = File::Spec->catdir( $db_path, "fasta" );
      my $catfile = File::Spec->catfile( $catdir, $loci_r->{$locus_id}{fasta_file} );
      my $fasta_r = read_fasta(-file => $catfile);

      $loci_r->{$locus_id}{'pam_sequence'} =
        seq_from_genome_fasta( $fasta_r,
                               [$loci_r->{$locus_id}{'scaffold'},
                                $loci_r->{$locus_id}{'pam_start'},
                                $loci_r->{$locus_id}{'pam_end'}]
                             );
    }
    else{ next; }

    # checking for existence of genome fasta, if yes, extract leader sequence #
    if(exists $loci_r->{$locus_id}{'fasta_file'}){
      my $catdir = File::Spec->catdir( $db_path, "fasta" );
      my $catfile = File::Spec->catfile( $catdir, $loci_r->{$locus_id}{fasta_file} );
      my $fasta_r = read_fasta(-file => $catfile);

      $loci_r->{$locus_id}{'leader_sequence'} =
        seq_from_genome_fasta( $fasta_r,
                               [$loci_r->{$locus_id}{'scaffold'},
                                $loci_r->{$locus_id}{'leader_start'},
                                $loci_r->{$locus_id}{'leader_end'}]
                             ) unless exists $loci_r->{$locus_id}{'leader_sequence'};

      if($loci_r->{$locus_id}{'pam_sequence'} eq ""){
        print STDERR "WARNING: '", $loci_r->{$locus_id}{'scaffold'},
          "' not found for ", $loci_r->{$locus_id}{'fasta_file'},
            ". Not loading pam sequence!\n";
      }
      else{             # just entries w/ leader sequence #
        $cp{$locus_id} = $loci_r->{$locus_id};
      }
    }
  }

  #print Dumper %cp; exit;
  return \%cp;
}


=head1 AUTHOR

Nick Youngblut, C<< <nyoungb2 at illinois.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-crispr_db at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=CRISPR_db>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.


=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::arrayBlast::loadLoci


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
