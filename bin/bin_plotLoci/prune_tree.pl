#!/usr/bin/env perl

=pod

=head1 NAME

prune_tree -- prune tree to just taxa in provided dna_segs or xlims table.

=head1 SYNOPSIS

=head2 Ordering 'dna_segs' table

prune_tree [flags] -d dna_segs.txt < tree_file > tree_file_prn

=head2 Ordering 'xlims' table

prune_tree [flags] -x xlims.txt < tree_file > tree_file_prn

=head2 Required flags

=over

=item -dna_segs  <char>

dna_segs table (NOT required if -x used).

=item -xlims  <char>

xlims table (NOT required if -d used).

=back

=head2 Optional flags

=over

=item -format  <char>

Tree file format. [newick]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb -- plotLoci --perldoc -- prune_tree

=head1 DESCRIPTION

Prune a tree to just taxa found in 
the provided dna_segs or xlims table.

=head2 Tree pruning

Tree pruning will be done with the 'nw_prune' from
the newick utilities package, if the 'nw_prune' script
is in your $PATH. If not, BioPerl will be used for tree
pruning, which may cause some warnings/errors with downstream
uses of the pruned tree.

=head2 Output

The tree is written to STDOUT in newick format (required for plotting).

=head1 EXAMPLES

=head2 Basic usage

=head3 Pruning by a dna_segs table

CLdb -- plotLoci -- prune_tree -d dna_segs.txt < tree.nwk > tree_prn.nwkk

=head3 Pruning by an xlims table

CLdb -- plotLoci -- prune_tree -x xlims.txt < tree.nwk > tree_prn.nwk

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

https://github.com/nyoungb2/CLdb

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut


### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use IPC::Cmd qw/can_run run/;
use IO::String;
use File::Temp;


# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../../lib";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $format, $database);
my ($dna_segs_in, $xlims_in);
GetOptions(
	   "format=s" => \$format,
	   "dna_segs=s" => \$dna_segs_in,
	   "xlims=s" => \$xlims_in,    # xlims input? 
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
# checking input files
$format = check_tree_format($format);
print "ERROR: provide -dna_segs or -xlims flag\n"
  unless ($dna_segs_in or $xlims_in);
# checking for nw_prune in path
my $nw_prune = can_run('nw_prune') ? 1 : 0;



#--- Main ---#
# loading tree #
my $treeio = Bio::TreeIO -> new(-fh => \*STDIN,
				-format => $format);
my $treeo = $treeio->next_tree;


# loading dna_segs|xlims table #
my $infile;
if($dna_segs_in){ $infile = $dna_segs_in; }
elsif($xlims_in){ $infile = $xlims_in; }
else{ die "LOGIC ERROR!\n"; }

my ($dna_segs_r, $dna_segs_order_r, 
    $header_r, $dna_segs_ids_r, 
    $name_index_r) = load_input_table($infile);



# prune tree by dna_segs (just taxon_names of interest) #
if($nw_prune and $format eq 'newick'){
  $treeo = nw_prune_tree($dna_segs_r, $treeo);
}
else{
  $treeo = prune_tree($dna_segs_r, $treeo);
}

	
## writing editted tree ##
tree_write($treeo);


#--- Subroutines ---#

=head2 tree_write

Writing tree object using bioperl

=cut

sub tree_write{
  my ($treeo) = @_;
  
  my $out = new Bio::TreeIO(-fh => \*STDOUT,
			    -format => "newick");
  $out->write_tree($treeo);
}


=head2 nw_prune_tree

Using nw_prune_tree for pruning.

=cut

sub nw_prune_tree{
  my ($dna_segs_r, $treeo) = @_;
  
  # getting taxa not in tree #
  my %rm;
  for my $node ($treeo->get_leaf_nodes){
    $rm{$node->id} = 1 unless exists $dna_segs_r->{$node->id};
  }

  # status
  print STDERR "#--- Pruning with nw_prune (from newick utilities) ---#\n";
  map{ print STDERR "Removing node: $_\n" } keys %rm;

  # writing tree to temp file
  my $tmpo = File::Temp->new();
  my $tmpfilename = $tmpo->filename;

  my $outo = Bio::TreeIO->new(-file => ">$tmpfilename",
			      -format => 'newick');
  $outo->write_tree($treeo);
  $outo = ();

  # calling nw_prune 
  my $cmd = join(" ", "nw_prune $tmpfilename", keys %rm, "|");

  open PIPE, $cmd or die $!;
  my $prune_tree = <PIPE>;
  close PIPE or die $!;

  
  # loading tree using bio::treeIO
  my $strIO = IO::String->new($prune_tree);
  my $treeio = Bio::TreeIO->new(-fh => $strIO,
				-format => 'newick');

  # return 1st tree (should only be 1)
  my $treeo_prn = $treeio->next_tree;
  return $treeo_prn;
}


=head2 prune_tree

Using Bio::Tree for pruning.

WARNING: appears to create odd newick files.

=cut

sub prune_tree{
  # pruning tree by dna_segs, just dna_segs taxa remaining #
  my ($dna_segs_r, $treeo) = @_;
  
  # getting taxa not in tree #
  my %rm;
  for my $node ($treeo->get_leaf_nodes){
    $rm{$node->id} = 1 unless exists $dna_segs_r->{$node->id};
  }
  
  # sanity check #	
  die " ERROR: < 2 leaves found in dna_segs table; no need to prune!\n"
    if ( (scalar $treeo->get_leaf_nodes) - (scalar keys %rm) ) <= 2;
  
  # editing internal node labels #
  for my $node ($treeo->get_nodes){
    next if $node->is_Leaf;
    $node->id("") unless $node->id;
    $node->id(join("__", "INTERNAL", $node->id));
  }
  
  # pruning tree #
  for my $node ($treeo->get_leaf_nodes){
    if (exists $rm{$node->id}){
      print STDERR "Pruning node: ", $node->id, "\n";
      #$treeo->remove_Node($node);
      $treeo->splice(-remove_id => $node->id, -preserve_lengths => 1);
    }
  }
  
  # pruning any non-labeled leaves #
  while(1){
    my $intern_cnt = 0;
    for my $node ($treeo->get_leaf_nodes){
      if($node->id =~ /^INTERNAL__/){
	#$treeo->remove_Node($node);
	$treeo->splice(-remove_id => $node->id, -preserve_lengths => 1);
	$intern_cnt++;
      }
    }
    last unless $intern_cnt;		# continue until all internal node leaves are removed
  }	

  # splicing out all singleton internal nodes
  my $root = $treeo->get_root_node;
  for my $node ($treeo->get_nodes){
    next if $node->is_Leaf;
   # next if $node->id eq $root->id;  # skipping root
    my @Desc = $node->each_Descendent;
    if(scalar @Desc  == 1){ # singleton internal node
      if($node->id eq $root->id){
	# if root = singleton: reroot on desc; removing root 
	$treeo->reroot($Desc[0]);
	$treeo->remove_Node($node);
      }
      else{  # not root
	$treeo->splice(-remove_id => $node->id, -preserve_lengths => 1);
      }
    }
  }  

  # editing internal node labels (setting names back) #
  for my $node ($treeo->get_nodes){
    next if $node->is_Leaf;
    next unless $node->id =~ /^INTERNAL__/;
    (my $tmp = $node->id) =~ s/^INTERNAL__//;
    $node->id( $tmp );
  }	
  
  return $treeo;
}

=head2 load_input_table

Loading input table: either dna_segs or xlims

=head3 IN

infile -- table file name

=cut

sub load_input_table{
  my $infile = shift or die "ERROR: provide infile\n";

  open IN, $infile or die $!;
  
  my @dna_segs_order;
  my %dna_segs;
  my %header;
  my %dna_segs_ids;
  my %name_index;
  while(<IN>){
    chomp;
    next if /^\s*$/;
    
    # header #
    if($.==1){
      # indexing header #
      my @tmp = split /\t/;
      for my $i (0..$#tmp){
	$header{$tmp[$i]} = $i;
      }
      
      # checking for existence of columns #
      my @check = qw/taxon_name locus_id subtype dna_segs_id/;
      map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
	     unless exists $header{$_} } @check;
    }
    else{
      my @line = split /\t/;
      my $taxon_name = $line[$header{"taxon_name"}];
      my $locus_id = $line[$header{"locus_id"}];
      my $subtype = $line[$header{"subtype"}];
      my $dna_seg_id = $line[$header{"dna_segs_id"}];
      
      # order of loci#
      push( @dna_segs_order, $taxon_name)
	unless exists $dna_segs{$taxon_name};
      
      # dna_segs  #
      push( @{$dna_segs{$taxon_name}{$locus_id}{"entries"}}, \@line );
      $dna_segs{$taxon_name}{$locus_id}{"subtype"} = $subtype;
      
      # dna_segs_id index #
      $dna_segs_ids{$locus_id} = $dna_seg_id;
      
      # name index #
      $name_index{$dna_seg_id} = $taxon_name;
    }
  }
  
  close IN or die $!;

  #print Dumper @dna_segs_order;
  #print Dumper %dna_segs; exit;
  #print Dumper %dna_segs_ids; exit;
  #print Dumper %name_index; exit;
  return \%dna_segs, \@dna_segs_order, \%header, \%dna_segs_ids ,\%name_index;
}

sub check_tree_format{
# determining tree format based on user input (if provided) or setting default
  my $format = shift;
  $format = "newick" if ! $format;
  $format =~ s/^new$/newick/i;
  $format =~ s/^nex$/nexus/i;
  die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
  return $format;
}


