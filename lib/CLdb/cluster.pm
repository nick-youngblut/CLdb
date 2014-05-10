package CLdb::cluster;

use 5.006;
use strict;
use warnings;
use IPC::Cmd qw/can_run/;

=head1 NAME

CLdb::cluster - subroutines for clustering sequences

=head1 VERSION

Version 0.01

=cut

our $VERSION = '0.01';


=head1 SYNOPSIS

=head1 EXPORT


=head1 SUBROUTINES/METHODS

=cut

use base 'Exporter';
our @EXPORT_OK = '';
use Carp qw/ confess carp croak/;
use Data::Dumper;


=head2 delete_cluster_entries

Deleting all existing entries in cluster table.
Resetting tables before redoing clustering.

=head3 IN

$dbh :  dbh object
$tbl :  name of table to delete

=head3 OUT

=cut

push @EXPORT_OK, 'clear_cluster_table';

sub clear_cluster_table{
  my $dbh = shift or croak "ERROR: provide a dbh object";
  my $tbl = shift or croak "ERROR: provide a table to clear";

  my $cmd = "DELETE FROM $tbl";
  my $sql = $dbh->prepare($cmd) or croak $dbh->err;
  $sql->execute();
  $dbh->commit();

  # status
  print STDERR "Deleted all entries in '$tbl' table\n";
}


=head2 insert_cluster

insert entries into cluster table

=head3 IN

hash of variables:

dbh => dbh object
cluster => cluster DS
element => spacer|DR

=head3 OUT

=cut

push @EXPORT_OK, 'insert_cluster';

sub insert_cluster{
  my %h = @_;
  my $dbh = exists $h{dbh} ? $h{dbh} :
    croak "ERROR: provide a dbh object\n";
  exists $h{cluster} or croak "ERROR: provide a cluster hashref\n";
  my $element = exists $h{element} ? $h{element} : 
    croak "ERROR: element = spacer or DR\n";

  # sql #
  my $cmd = "INSERT INTO 
$element\_clusters(locus_id, $element\_id, cutoff, cluster_id, Rep_sequence) 
values(?,?,?,?,?)";
  my $sql = $dbh->prepare($cmd) or croak $dbh->err;

  # insert #
  my $entry_cnt = 0;
  foreach my $cutoff ( keys %{$h{cluster}} ){
    foreach my $seq_name ( keys %{$h{cluster}{$cutoff}} ){
      $sql->bind_param(1, $h{cluster}{$cutoff}{$seq_name}{locus_id});
      $sql->bind_param(2, $h{cluster}{$cutoff}{$seq_name}{element_id});
      $sql->bind_param(3, $cutoff);
      $sql->bind_param(4, $h{cluster}{$cutoff}{$seq_name}{cluster_id});
      $sql->bind_param(5, $h{cluster}{$cutoff}{$seq_name}{rep_seq});
     
      $sql->execute() or croak $dbh->err;
      $entry_cnt++;
    }
  }
  $dbh->commit or croak $dbh->err;

  print STDERR "...$entry_cnt $element cluster entries added/updated\n" 
    unless $h{verbose};
}


=head2 cluster_cutoffs 

Example subroutine

=head3 IN

hash of variables: 

cluster_cutoff => array_ref of clustering cutoffs
fasta_file => fasta file name
element => spacer or DR
threads => number of threads

=head3 OUT

$%cluster{cluster_cutoff}=>{seq_name}=>{cat}=>value

=cut

push @EXPORT_OK, 'cluster_cutoffs';

sub cluster_cutoffs{
  my %h = @_;
  exists $h{cluster_cutoff} or croak "ERROR: provide a cluster cutoff array ref [start,end,jump]";
  exists $h{fasta_file} or croak "ERROR: provide a fasta file name";
  exists $h{element} or croak "ERROR: defined element [spacer|DR]";
  $h{threads} = 1 unless exists $h{threads};
  
  # checking for existing cd hit
  can_run("cd-hit-est") or croak "ERROR: cannot find cd-hit-est\n";

  # workflow for each clustering cutoff #  
  my %cluster;
  for (my $i=$h{cluster_cutoff}->[0];
       $i<=$h{cluster_cutoff}->[1];                 # rounding issues
       $i+=$h{cluster_cutoff}->[2]){
    $i = sprintf('%.2f', $i);
    printf STDERR "...Clustering %s at cutoff: $i\n", $h{element} unless $h{verbose};

    ## call cd-hit-est ##
    my $cdhit_out = call_cdhit( \%h, $i );

    ## parsing cdhit clustering file
    my $cluster_index = ID2clusterID("$cdhit_out.clstr", $h{verbose});
    
    ## parse cd-hit-est cluster representative sequence output ##
    $cluster{$i} = parse_cdhit_seq($cdhit_out, $cluster_index);    # saving for each cluster cutoff
  }
  
  #print Dumper %cluster; exit;
  return \%cluster;
  
}


=head2 call_cdhit

Calling ed-hit-est on spacers or DRs

Strand-specific clustering.

Sequences must be the same length to cluster together.

=head3 IN

$fasta :  fasta file
$cluster :  cluster cutoff
$threads : number of theads used  

=head3 OUT

=cut

sub call_cdhit{
# calling cd-hit-est on spacers/DRs #
## sequences must be same length to cluster together
  my ($h, $cluster) = @_;
  
  (my $cdhit_out = $h->{fasta_file}) =~ s/\.fna/.cdhit/;
  
  my $cmd = join(" ", "cd-hit-est -i",  $h->{fasta_file}, 
		 "-o $cdhit_out -c $cluster -n 8 -s 1 -r 0 -d 0 -T",
		 $h->{threads}, "-s", $h->{length} ) ;
  
  print STDERR "$cmd\n" if $h->{verbose};
  `$cmd`;
  if ($? == -1) {
    die "ERROR: failed to execute: $!\n";
  }
  elsif ($? & 127) {
    die "ERROR: child died with signal %d, %s coredump\n",
      ($? & 127),  ($? & 128) ? 'with' : 'without';
  }
  else {
    #die "ERROR: child exited with value %d\n", $? >> 8;
  }
  
  #exit;
  return $cdhit_out;
}


=head2 ID2clusterID

Parsing cd-hit cluster report output to make rep-seqid => clusterID index

=head3 IN

$cdhit_out_clst :  cluster report output file name

=head3 OUT

=cut

sub ID2clusterID{
  my $cdhit_out_clst = shift or croak "ERROR: provide a cdhit cluster output file\n";
  my $verbose = shift;

  open IN, $cdhit_out_clst or croak $!;
  my %N_clusters;
  my %cluster_index;
  my $cluster_id;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    if(/>Cluster (\d+)/){
      $cluster_id = $1;
      $N_clusters{$cluster_id} = 1;
    }
    elsif( /nt,.+>(.+)\.\.\. [a*]/ ){
      croak "ERROR: cd-hit cluster file not formatted corrrectly\n"
	unless defined $cluster_id;
      $cluster_index{ $1} = $cluster_id;
    }
    else{ next; }
  }
  close IN;

  # status
  printf STDERR "\tNumber of clusters produced: %i\n",
    scalar keys %N_clusters;

  #print Dumper %cluster_index; exit;
  return \%cluster_index;
}


=head2 parse_cdhit_seq

Parsing the output from cd-hit-est

Assuming sequence name in format: 'locus_id|element|element_id'

=head3 IN

$cdhit_out :  output cluster file for cd-hit-est run
$cluster_index :  seq_name => cluster_ID  # no '>' 

=head3 OUT

seq_name => {locus_id|element|element_id|cluster_ID|rep_seq} => value 

=cut

sub parse_cdhit_seq{
  my ($cdhit_out, $cluster_index) = @_;
 
  open IN, $cdhit_out or die $!;
  my $seq_name;
  my %seq_cluster;  # seq_name => {locus_id|element|element_id|cluster_ID|rep_seq} => value 
  while(<IN>){
    chomp;
    next if /^\s*$/;
    
    if(/^>(.+)/){
      $seq_name = $1;
      my @seq_name = split /\|/, $seq_name;
      croak "ERROR: sequence name '$seq_name' is not formated correctly\n"
	unless scalar @seq_name == 3;
      $seq_cluster{$seq_name}{locus_id} = $seq_name[0];
      $seq_cluster{$seq_name}{element} = $seq_name[1];
      $seq_cluster{$seq_name}{element_id} = $seq_name[2];
      $seq_cluster{$seq_name}{cluster_id} = exists $cluster_index->{$seq_name} ?
	$cluster_index->{$seq_name} : croak "ERROR: cannot find cluster ID for $seq_name\n";
    }
    else{
      $seq_cluster{ $seq_name }{rep_seq} .= $_;
    }
  }
  close IN;
  
  #print Dumper %seq_cluster; exit;
  return \%seq_cluster;
}


=head1 AUTHOR

Nick Youngblut, C<< <ndy2 at cornell.edu> >>

=head1 BUGS

Please report any bugs or feature requests to C<bug-cluster at rt.cpan.org>, or through
the web interface at L<http://rt.cpan.org/NoAuth/ReportBug.html?Queue=cluster>.  I will be notified, and then you'll
automatically be notified of progress on your bug as I make changes.

=head1 SUPPORT

You can find documentation for this module with the perldoc command.

    perldoc CLdb::cluster


You can also look for information at:

=over 4

=item * RT: CPAN's request tracker (report bugs here)

L<http://rt.cpan.org/NoAuth/Bugs.html?Dist=cluster>

=item * AnnoCPAN: Annotated CPAN documentation

L<http://annocpan.org/dist/cluster>

=item * CPAN Ratings

L<http://cpanratings.perl.org/d/cluster>

=item * Search CPAN

L<http://search.cpan.org/dist/cluster/>

=back


=head1 ACKNOWLEDGEMENTS


=head1 LICENSE AND COPYRIGHT

Copyright 2014 Nick Youngblut.

This program is free software; you can redistribute it and/or modify it
under the terms of either: the GNU General Public License as published
by the Free Software Foundation; or the Artistic License.

See http://dev.perl.org/licenses/ for more information.


=cut

1; # End of CLdb::cluster
