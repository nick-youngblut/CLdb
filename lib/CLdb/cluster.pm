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
use Carp qw/ confess carp confess/;
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
  my $dbh = shift or confess "ERROR: provide a dbh object";
  my $tbl = shift or confess "ERROR: provide a table to clear";

  my $cmd = "DELETE FROM $tbl";
  my $sql = $dbh->prepare($cmd) or confess $dbh->err;
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
    confess "ERROR: provide a dbh object\n";
  exists $h{cluster} or confess "ERROR: provide a cluster hashref\n";
  my $element = exists $h{element} ? $h{element} : 
    confess "ERROR: element = spacer or DR\n";

  # sql #
  my $cmd = "INSERT INTO 
$element\_clusters(locus_id, $element\_id, cutoff, cluster_id, Rep_sequence) 
values(?,?,?,?,?)";
  my $sql = $dbh->prepare($cmd) or confess $dbh->err;

  # insert #
  my %debug;
  my $entry_cnt = 0;
  foreach my $cutoff ( keys %{$h{cluster}} ){
    foreach my $seq_name (sort  keys %{$h{cluster}{$cutoff}} ){
      $sql->bind_param(1, $h{cluster}{$cutoff}{$seq_name}{locus_id});
      $sql->bind_param(2, $h{cluster}{$cutoff}{$seq_name}{element_id});
      $sql->bind_param(3, $cutoff);
      $sql->bind_param(4, $h{cluster}{$cutoff}{$seq_name}{cluster_id});
      $sql->bind_param(5, $h{cluster}{$cutoff}{$seq_name}{rep_seq});
     
      $sql->execute() or confess $dbh->err;
      $entry_cnt++;

      $debug{$cutoff}{ $h{cluster}{$cutoff}{$seq_name}{cluster_id} }++;
    }
  }
  $dbh->commit or confess $dbh->err;

  #print Dumper %debug; exit;
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
cat = {locus_id|element|element_id|cluster_ID|rep_seq}
 

=cut

push @EXPORT_OK, 'cluster_cutoffs';

sub cluster_cutoffs{
  my %h = @_;
  exists $h{cluster_cutoff} or 
    confess "ERROR: provide a cluster cutoff array ref [start,end,jump]";
  exists $h{fasta_file} or confess "ERROR: provide a fasta file name";
  exists $h{element} or confess "ERROR: defined element [spacer|DR]";
  $h{threads} = 1 unless exists $h{threads};
  
  # checking for existing cd hit
  can_run("cd-hit-est") or confess "ERROR: cannot find cd-hit-est\n";

  # workflow for each clustering cutoff #  
  my %cluster;
  for (my $i=$h{cluster_cutoff}->[0];
       $i<=$h{cluster_cutoff}->[1];                 # rounding issues
       $i+=$h{cluster_cutoff}->[2]){
    $i = sprintf('%.2f', $i);
    printf STDERR "...Clustering %s at cutoff: $i\n", $h{element} unless $h{verbose};

    ## call cd-hit-est ##
    my $cdhit_out = call_cdhit( \%h, $i );

    ## parsing cdhit cluster report file
    my $clustReport_r = parse_cdhit_clstr("$cdhit_out.clstr",  # name of cluster file
		      $h{verbose});

    ## parse cd-hit-est cluster representative sequence output ##
    my $clustRepSeq_r = parse_cdhit_seq($cdhit_out, $clustReport_r);    

    ## merging fasta file and cluster report data
    ## saving for each cluster cutoff
    $cluster{$i} = merge_cdhit_output($clustRepSeq_r, $clustReport_r);
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
  }
  
  #exit;
  return $cdhit_out;
}


=head2 merge_cdhit_output

Merging data from fasta and clsuter report output from cd-hit.
Assuming seq_name formated as: 
 'locus_id|element|element|element_id|cluster_id|cutoff'

=head3 IN

$clustRepSeq_r :  cluster_id => [seq_names]
$clustReport_r :  cluster_id => rep_sequence

=head3 OUT

$%seq_name => {locus_id|element|element_id|cluster_ID|rep_seq} => value 

=cut

sub merge_cdhit_output{
  my $clustRepSeq_r = shift or confess "Provide a $% of rep sequences";
  my $clustReport_r = shift or confess "Provide a $% of cluster report entries";

  # making index associating rep_sequence with cluster
  my %merged;
  foreach my $clusterID (keys %{$clustReport_r->{clust_name}} ){
    confess "Cannot find $clusterID in rep sequences"
      unless exists $clustRepSeq_r->{$clusterID};

    # loading %merged
    foreach my $seq_name ( @{$clustReport_r->{clust_name}{$clusterID}} ){
      # checking name
      confess "Do not recognize $seq_name format"
	unless $seq_name =~ /[^|]+\|(spacer|DR)\|\d+\|NA\|NA/;
      my @l = split /\|/, $seq_name;  # locus_id|element|element_id
      confess "Do not recognize $seq_name format"
	unless scalar @l == 5;

      $merged{$seq_name}{locus_id} = $l[0];
      $merged{$seq_name}{element} = $l[1];
      $merged{$seq_name}{element_id} = $l[2];
      $merged{$seq_name}{cluster_id} = $clusterID;
      $merged{$seq_name}{rep_seq} = $clustRepSeq_r->{$clusterID};      
    }
  }
 
  #print Dumper %merged; exit;
  return \%merged;
}


=head2 parse_cdhit_seq

Parsing the output from cd-hit-est

Associating rep sequence ot cluster

=head3 IN

$cdhit_out :  output cluster file for cd-hit-est run
$cluster_index :  seq_name => cluster_ID  # no '>' 

=head3 OUT

$%cluster_id => rep_sequence
# sequence_id is just the name of the rep sequene for the cluster

=cut

sub parse_cdhit_seq{
  my ($cdhit_out, $clustReport_r) = @_;
 
  open IN, $cdhit_out or die $!;
  my $clust_name;
  my %seq_cluster;  
  while(<IN>){
    chomp;
    next if /^\s*$/;
    
    if(/^>(.+)/){
      my $seq_name = $1;
      confess "Cannot find $seq_name in cluster report"
	unless exists $clustReport_r->{name_clust}{$seq_name};
      $clust_name = $clustReport_r->{name_clust}{$seq_name}; # cluster_id=>rep_seq
    }
    else{
      $seq_cluster{ $clust_name } .= $_;
    }
  }
  close IN;
  
  #print Dumper %seq_cluster; exit;
  return \%seq_cluster;
}


=head2 parse_cdhit_clstr

Parsing cd-hit cluster report output to make rep-seqid => clusterID index

=head3 IN

$cdhit_out_clst :  cluster report output file name
$verbose :  verbose bool

=head3 OUT

$%{clust_name}{cluster} => [seq_names]
$%{name_clust}{seq_name} => clust_id

=cut

sub parse_cdhit_clstr{
  my $cdhit_out_clst = shift or confess "ERROR: provide a cdhit cluster output file\n";
  my $verbose = shift;

  open IN, $cdhit_out_clst or confess $!;
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
      my $seq_name = $1;
      confess "ERROR: cd-hit cluster file not formatted corrrectly\n"
	unless defined $cluster_id;
      push @{$cluster_index{clust_name}{ $cluster_id }}, $seq_name;
      $cluster_index{name_clust}{$seq_name} = $cluster_id;
    }
    else{ next; }
  }
  close IN;

  # stats
  printf STDERR "\tNumber of clusters produced: %i\n",
    scalar keys %N_clusters;

#  print Dumper %cluster_index; exit;
  return \%cluster_index;
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
