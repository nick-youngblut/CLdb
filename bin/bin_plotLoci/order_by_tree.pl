#!/usr/bin/env perl

=pod

=head1 NAME

order_by_tree -- needed if plotting a tree with CRISPR loci

=head1 SYNOPSIS

=head2 Ordering 'dna_segs' table

order_by_tree [flags] < dna_segs.txt > dna_segs_ordered.txt

=head2 Ordering 'xlims' table

order_by_tree [flags] < xlims.txt > xlims_ordered.txt

=head2 Required flags

=over

=item -tree  <char>

Tree file name (nexus or newick).

=back

=head2 Optional flags

=over

=item -format  <char>

Tree file format. [newick]

=item -name  <char>

Output file name for pruned tree. ['-tree' + '_prn.nwk']

=item -xlims  <bool>

Ordering xlims table instead of dna_segs. An editted tree will not be written. [FALSE]

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

CLdb -- plotLoci --perldoc -- order_by_tree

=head1 DESCRIPTION

Formatting tree and 'dna_segs' or 'xlims' table for plotting. 

The dna_segs/xlims tables are ordered to match the provided
tree ordering. 

=head2 Output

If the dna_segs/xlims table has the same taxon present twice
(eg., if 2 CRISPR loci from the same taxon), then 
extra leaves (with edited taxon_names) are added to the tree.
This editted tree is written to a new file.

=head2 Warnings

Extra taxa can be in the tree, but a pruned version of the
tree is needed for plotting!

Both xlims and dna_segs tables must be ordered by the tree
if the tree is used for plotting!


=head1 EXAMPLES

=head2 Basic usage

=head3 Ordering a dna_segs table

CLdb -- plotLoci -- order_by_tree -t tree.nwk < dna_segs.txt > dna_segs_ordered.txt

=head3 Ordering the complemenary xlims table

CLdb -- plotLoci -- order_by_tree -t tree.nwk -x < xlims.txt > xlims_ordered.txt

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
use DBI;
use Bio::TreeIO;


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


my ($verbose, $tree_in, $format, $tree_name, $xlims_bool, $database);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "name=s" => \$tree_name,	
	   "xlims" => \$xlims_bool, # xlims input? 
	   "database=s" => \$database, # unused
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
# checking input files
file_exists($tree_in, "tree");
$format = check_tree_format($format);
($tree_name = $tree_in) =~ s/\.[^.]+$|$/_edit.nwk/ unless $tree_name;



#--- Main ---#
# loading tree #
my $treeio = Bio::TreeIO -> new(-file => $tree_in,
				-format => $format);
my $treeo = $treeio->next_tree;

# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r, 
	$dna_segs_ids_r, $name_index_r) = load_dna_segs();

# mutli-loci / multi-subtype problem ##
## adding leaves if multiple taxon entries; also adding subtype ##
$treeo = add_leaves($dna_segs_r, $treeo, $dna_segs_ids_r);
	
## writing editted tree ##
tree_write($treeo, $tree_name) unless $xlims_bool;

# ordering dna_segs by tree #
## getting tree order ##
my $tree_order_r = get_tree_order($treeo);

## ordering table by tree and writing ##
order_dna_segs($dna_segs_r, $tree_order_r, $header_r, $name_index_r);



#--- Subroutines ---#
sub tree_write{
  my ($treeo, $tree_name) = @_;
  
  my $out = new Bio::TreeIO(-file => ">$tree_name",
			    -format => "newick");
  $out->write_tree($treeo);

  print STDERR "Editted tree written: $tree_name\n";
}

sub order_dna_segs{
  # writing out dna segs with new tree order #
  my ($dna_segs_r, $tree_order_r, $header_r, $name_index_r) = @_;
  
  # header #
  print join("\t", sort{$header_r->{$a}<=>$header_r->{$b}} 
	     keys %$header_r), "\n";


  # body #
  foreach my $leaf (@$tree_order_r){
    my $taxon_name = $name_index_r->{$leaf};
    next unless defined $taxon_name;
    unless (exists $dna_segs_r->{$taxon_name}){
      warn "Warning: leaf->$taxon_name not found in table. Skipping\n";
      next;
    }
      
    if($leaf =~ /__[^_]+/){   	## if lociID already in name (prevents duplicates)
      (my $locus_id = $leaf) =~ s/^.*?__([^_]+).*/$1/; 	#if $leaf =~ /\|\d+$/;
      if(exists $dna_segs_r->{$taxon_name}{$locus_id}){ 
	foreach my $row (@{$dna_segs_r->{$taxon_name}{$locus_id}{"entries"}}){
	  print join("\t", @$row), "\n";
	}
      }		
      else{  # not locus_id, should then be subtype
	foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}}){
	  foreach my $row (@{$dna_segs_r->{$taxon_name}{$locus_id}{"entries"}}){
	    print join("\t", @$row), "\n";
	  }
	}
      }
    }
    else{     ## if 1 lociID per taxon_name
      foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}}){
	foreach my $row (@{$dna_segs_r->{$taxon_name}{$locus_id}{"entries"}}){
	  print join("\t", @$row), "\n";
	}
      }
    }
  }	
}


sub get_tree_order{
  my ($treeo) = @_;
  
  my @tree_order;
  for my $node ($treeo->get_leaf_nodes){
    push @tree_order, $node->id;
  }
  
  map{ s/['"]//g } @tree_order;
  
  return \@tree_order;
}


sub add_leaves{
  # adding leaves to any nodes that have multiple subtypes #
  my ($dna_segs_r, $treeo, $dna_segs_ids_r) = @_;
	
  # splitting taxa if needed #
  for my $node ($treeo->get_leaf_nodes){
    my @loci = keys %{$dna_segs_r->{$node->id}};  # locus_ids for each taxon_name 
    
    for my $i (0..$#loci){
      if ($#loci > 0){		# adding nodes if needed 	
	# sanity check #
	die " ERROR: locus_id -> $loci[$i] not found in dna_segs_ids!\n"
	  unless exists $dna_segs_ids_r->{$loci[$i]};
	
	# getting node info #
	my $anc_node = $node->ancestor;
	my $brlen = $node->branch_length();
	(my $new_node_id = $node->id) =~ s/__[^_]+$//g;
	
	my $new_node = new Bio::Tree::Node(
					   -id => $dna_segs_ids_r->{$loci[$i]},
					   -branch_length => 0);      
	$node->add_Descendent($new_node);
	
	$node->id(100) if $i == $#loci;
      }
      else{	# appending locus_id #				
	# changing node id #
	$node->id( $dna_segs_ids_r->{$loci[$i]} );
      }
    }
  }
  return $treeo;
}
	
sub check_multi{
  # checking for multiple entries per taxon #
  my ($dna_segs_r) = @_;
  
  my $multi_loci = 0;				# mutliple loci per taxon_name
  my $multi_subtype = 0;			# multiple subtypes total 
  my %subtypes; 
  foreach my $taxon_name (keys %$dna_segs_r){
    
    $multi_loci = 1 if scalar keys %{$dna_segs_r->{$taxon_name}} > 1;
    
    foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}}){
      $subtypes{$dna_segs_r->{$taxon_name}{$locus_id}{"subtype"}}++;
    }
  }
  $multi_subtype = 1 if scalar keys %subtypes > 1;
  
  #print Dumper $multi_loci, $multi_subtype; exit;
  # status #
  print STDERR join(" ", "...Found multiple loci for 1 or more taxa.",
		    "Adding leaves to the tree!", 
		    "Adding loci_ids to leaves & dna_segs table!"), "\n"
		      if $multi_loci;
  print STDERR join(" ", "...Found multiple subtypes.",
		    "Adding subtype to names in tree & dna_segs table!"), "\n"
		      if $multi_subtype;
  
  return $multi_loci, $multi_subtype;
}


=head2 nw_prune_tree

Using nw_prune_tree for pruning.

=cut

sub nw_prune_tree{
  my ($dna_segs_r, $tree_in) = @_;
  
  # getting taxa not in tree #
  my %rm;
  for my $node ($treeo->get_leaf_nodes){
    $rm{$node->id} = 1 unless exists $dna_segs_r->{$node->id};
  }

  # status
  print STDERR "#--- Pruning with nw_prune (from newick utilities) ---#\n";
  map{ print STDERR "Removing node: $_\n" } keys %rm;

  # calling nw_prune
  my $cmd = join(" ", "nw_prune $tree_in", keys %rm, "|");
  
  open PIPE, $cmd or die $!;
  my $prune_tree = <PIPE>;
  close PIPE or die $!;

  # loading tree using bio::treeIO
  my $strIO = IO::String->new($prune_tree);
  my $treeio = Bio::TreeIO->new(-fh => $strIO,
				-format => 'newick');

  # return 1st tree (should only be 1)
  my $treeo = $treeio->next_tree;
  return $treeo;
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

sub load_dna_segs{
  # loading dna_segs from stdin #
  
  my @dna_segs_order;
  my %dna_segs;
  my %header;
  my %dna_segs_ids;
  my %name_index;
  while(<>){
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


