package CLdb::load::loadLoci;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;

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
  my ($loci_r, $db_path) = @_;

  print STDERR "### checking line breaks for all external files (converting to unix) ###\n";

  #my @file_columns = qw/genbank_file fasta_file array_file/;
  my %file_cols = (
                   "genbank_file" => "genbank",
                   "fasta_file" => "fasta",
                   "array_file" => "array");

  foreach my $locus_id (keys %$loci_r){
    foreach my $file_col (keys %file_cols){
      if(exists $loci_r->{$locus_id}{$file_col}){
        next unless $loci_r->{$locus_id}{$file_col};                            # if no file; nothing to check

	my $dirpath = File::Spec->catdir($db_path, $file_cols{$file_col});
	my $filepath = File::Spec->catdir($dirpath, $loci_r->{$locus_id}{$file_col});

        print STDERR " processing: $filepath\n";
        lineBreaks2unix($filepath, 1);
      }
    }
  }
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

Moving files to directory in CLdb_home

=cut

push @EXPORT_OK, 'make_external_file_dirs';

sub make_external_file_dirs{

  my ($loci_r, $header_r, $dir) = @_;

  foreach my $locus_id (keys %$loci_r){
    $loci_r->{$locus_id}{"genbank_file"} =
      make_CLdb_dir($dir, 'genbank', $loci_r->{$locus_id}{"genbank_file"})
        if exists $loci_r->{$locus_id}{"genbank_file"};
    $loci_r->{$locus_id}{"array_file"} =
      make_CLdb_dir($dir, 'array', $loci_r->{$locus_id}{"array_file"})
        if exists $loci_r->{$locus_id}{"array_file"};
    $loci_r->{$locus_id}{"fasta_file"} =
      make_CLdb_dir($dir, 'fasta', $loci_r->{$locus_id}{"fasta_file"})
        if exists $loci_r->{$locus_id}{"fasta_file"};
  }

  #print Dumper %$loci_r; exit;
  return $dir;
}

=head2 make_CLdb_dir

making directory in $CLdb_home

=cut

sub make_CLdb_dir{
  my ($dir, $name, $infile) = @_;
  my @parts = File::Spec->splitpath( $infile );

  # making dir; copying files #
  if(File::Spec->rel2abs($parts[1]) ne "$dir/$name"){
    my $newdir = File::Spec->catdir($dir, $name);
    mkdir $newdir unless -d $newdir;

    my $focal_file = File::Spec->catfile($newdir, $parts[2]);
    unless(-e $focal_file){          # can't find file in needed directory, is it in specified directory?
      croak " ERROR: cannot find $infile.\n Checked:\n\t'$focal_file'\n\t'$infile'\n\n"
        unless -e $infile;                      # specified file can't be found anywhere

      print STDERR "'$infile' not in $name of \$CLdb_home. Moving to \n";
      copy($infile, $focal_file) or die $!;
      print STDERR "\tCopied '$infile' to '$focal_file'\n";
    }
  }

  return $parts[2];
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
