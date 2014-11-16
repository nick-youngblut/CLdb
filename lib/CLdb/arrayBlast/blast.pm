package CLdb::arrayBlast::blast;

# module use #
use strict;
use warnings FATAL => 'all';
use Carp  qw( carp confess croak );
use Data::Dumper;
use File::Spec;
use Parallel::ForkManager;
use File::Path qw/rmtree/;
use IPC::Cmd qw/run can_run/;

# export #
use base 'Exporter';
our @EXPORT_OK = '';

	

=head1 NAME

CLdb::arrayBlast::blast

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

Subroutines for parsing & editing spacer/DR blast files

=head1 EXPORT_OK

=cut


=head2 make_blast_dir

Making blast database directory at $CLdb_HOME/arrayBlast/blast_db/

Any existing directory will be deleted.

=head3 IN

=head3 OUT

=cut

push @EXPORT_OK, 'make_blast_dir';

sub make_blast_dir{
  my $db_path = shift;

  my $dir = join("/", $db_path, "arrayBlast");
  mkdir $dir unless -d $dir;

  $dir = "$dir/blast_db/";
  rmtree $dir if -d $dir;
  mkdir $dir;

  return $dir;
}


=head2 call_blastn_short

Calling blastn -task 'blastn-short'

=cut

push @EXPORT_OK, 'call_blastn_short';

sub call_blastn_short{
  my ($blast_dbs_r, $query_file, $blast_params, $fork) = @_;

  # check for blastn
  can_run('blastn') or die "ERROR: 'blastn' executable not found!\n";

  # blast output format
  my $outfmt = '5';

  # sanity check #
  map{die " ERROR: cannot find $_!\n" unless -e $_} @$blast_dbs_r;

  # forking
  my $pm = Parallel::ForkManager->new($fork);

  $pm->run_on_finish(
    sub{
      my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $ret_r) = @_;
      # status #
      foreach my $out (sort keys %$ret_r){
	my @parts = File::Spec->splitpath($out);   # just db file name needed
        print STDERR "Number of blast hits to $parts[2]:\t",
          scalar @{$ret_r->{$out}}, "\n";
      }

      # writing output #
      foreach my $out (sort keys %$ret_r){
        print join("", @{$ret_r->{$out}});
      }
    }
   );

  my %output;
  foreach my $blast_db (@$blast_dbs_r){
    $pm->start and next;

    my %output;
    my $cmd = join(" ", "blastn -task 'blastn-short'",
		   "-query $query_file",
		   "-db $blast_db",
		   "-outfmt '$outfmt'",
		   "$blast_params");
    open PIPE, "$cmd |" or die $!;
    while(<PIPE>){
      push @{$output{$blast_db}}, $_;
    }
    close PIPE;
    $pm->finish(0, \%output);
  }
  $pm->wait_all_children;
}


=head2 make_blast_db

Making a blastn database for array blast

=head3 Output

array of blast files

=cut

push @EXPORT_OK, 'make_blast_db';

sub make_blast_db{
# foreach fasta, making a blast DB #
  my ($subject_loci_r, $db_path, 
      $blast_dir, $verbose) = @_;

  # checking for makeblastdb in path
  can_run('makeblastdb') or 
    croak "ERROR: 'makeblastdb' not in \$PATH." .
      " Add the blast+ toolkit to your path\n";

  # dir
  my $fasta_dir = "$db_path/fasta";

  # sanity check #
  die " ERROR: cannot find $fasta_dir!\n" unless -d $fasta_dir;
  die " ERROR: cannot find $blast_dir!\n" unless -d $blast_dir;

  # status #
  print STDERR "...making blast databases in $blast_dir\n";

  # making blast dbs #
  my @blastdbs;
  foreach my $row (@$subject_loci_r){
    unless($$row[2]){
      map{$_ = "" unless $_} @$row;
      print STDERR " Skipping: '", join(",", @$row), "'! 'fasta_file' value empty!\n";
      next;             # next unless fasta file present
    }

    # sanity check #
    die " ERROR: cannnot find $fasta_dir/$$row[2]!"
      unless -e "$fasta_dir/$$row[2]";

    # making symlink in blast directory #
    unless(-e "$blast_dir/$$row[2]" || -l "$blast_dir/$$row[2]"){
      symlink("$fasta_dir/$$row[2]", "$blast_dir/$$row[2]") or die $!;
    }

    # making blast db #
    my $cmd = "makeblastdb -dbtype nucl -parse_seqids -in $blast_dir/$$row[2]";
    print STDERR "$cmd\n" unless $verbose;
    `$cmd`;

    push @blastdbs, "$blast_dir/$$row[2]";
  }

  return \@blastdbs;
}



=head2 write_blast_file

writing out blast file

=cut

push @EXPORT_OK, 'write_blast_file';

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


=head2 read_blast_file

Reading blast file (format: '-outfmt 7') from STDIN

=cut

push @EXPORT_OK, 'read_blast_file';

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
      confess "ERROR: no '# Query' comment found for an entry in the blast table!\n"
	unless defined $query;
      confess "ERROR: no '# Database' comment found for an entry in the blast table!\n"
	unless defined $db;
      confess "ERROR: no '# BLAST' comment found for an entry in the blast table!\n"
	unless defined $blast;
      
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


=head2 parse_blast_hits

Parsing blast hit file (file in '-outfmt 7' format)

Input: $(blast hit file), [$(OID field)]

Output: %{query} => {db} => {category} => [hits]

=cut

push @EXPORT_OK, 'parse_blast_hits';

sub parse_blast_hits{
  my ($blast_in, $OID_field) = @_;

  $OID_field-- if defined $OID_field;

  my $fh;
  $blast_in == 0 ?  $fh = *STDIN : 
    open $fh, $blast_in or croak $!;
	
  my %blast;
  my %fields;
  my @header;
  my ($db, $query);
  while(<$fh>){
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
	#@header = ();
      }
      # parsing hit line
      my @line = split /\t/;	      
      # setting OID field
      my $subjectID_i = exists $fields{'subject id'} ?
	$fields{'subject id'} : croak $!;
      $line[$subjectID_i] = (split /\|/, $line[$subjectID_i])[$OID_field];      

      # loading hash with hit
      push @{$blast{$query}{$db}{"hits"}}, \@line;
    }
    
  }
  close $fh;

  confess " ERROR: not blast hits found in input file!\n" unless %blast;
  return \%blast, \%fields;
}


=head2 parse_header

Parsing header from comments of blast in outfmt 7 format

=cut

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

    perldoc CLdb::arrayBlast::blast


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
