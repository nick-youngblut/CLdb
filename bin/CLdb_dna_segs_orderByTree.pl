#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;
use Bio::TreeIO;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $tree_in, $format, $tree_name);
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "name=s" => \$tree_name,	   
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: cannot find tree file!\n" unless -e $tree_in;
$format = check_tree_format($format);
($tree_name = $tree_in) =~ s/\.[^.]+$|$/_prn.nwk/ unless $tree_name;

### Main
# loading tree #
my $treeio = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
my $treeo = $treeio->next_tree;

# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r) = load_dna_segs();

# prune tree by dna_segs #
$treeo = prune_tree($dna_segs_r, $treeo);

# mutli-loci / multi-subtype problem
## determining if multiple loci/subtypes per taxon_name ##
my ($multi_loci, $multi_subtype)  = check_multi($dna_segs_r);

## adding leaves if multiple taxon entries; also adding sub ##
$treeo = add_leaves($dna_segs_r, $treeo) if $multi_loci;
	
## adding subtype if needed ##
$treeo = add_subtype($dna_segs_r, $treeo) if $multi_subtype;
	
## writing editted tree ##
tree_write($treeo, $tree_name);

## adding loci & subtypes to DNA_segs names ##
$dna_segs_r = edit_dna_segs_taxon_name($dna_segs_r, $header_r, $multi_loci, $multi_subtype);

# getting tree order #
my $tree_order_r = get_tree_order($treeo);

# ordering table by tree and writing #
order_dna_segs($dna_segs_r, $tree_order_r, $header_r);



### Subroutines
sub tree_write{
	my ($treeo, $tree_name) = @_;
	
	my $out = new Bio::TreeIO(-file => ">$tree_name",
						-format => "newick");
	$out->write_tree($treeo);

	#exit;
	}

sub order_dna_segs{
	my ($dna_segs_r, $tree_order_r, $header_r) = @_;
	
	# header #
	print join("\t", sort{$header_r->{$a}<=>$header_r->{$b}} keys %$header_r);
	print  "\ttaxon_name_original\n";
	
	# body #
	foreach my $leaf (@$tree_order_r){
		die " ERROR: leaf -> \"$leaf\" not found in dna_segs table!\n"
			unless exists $dna_segs_r->{$leaf};
		foreach my $locus_id (keys %{$dna_segs_r->{$leaf}}){
			foreach my $row (@{$dna_segs_r->{$leaf}{$locus_id}{"entries"}}){
				print join("\t", @$row), "\n";
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

sub edit_dna_segs_taxon_name{
	my ($dna_segs_r, $header_r, $multi_loci, $multi_subtype) = @_;
	
	my %new_dna_segs;
	foreach my $taxon_name (keys %$dna_segs_r){
		foreach my $locus_id (keys %{$dna_segs_r->{$taxon_name}}){
			foreach my $row (@{$dna_segs_r->{$taxon_name}{$locus_id}{"entries"}}){				
				# appending new column of original taxon_names #
				push(@$row, $taxon_name);
				
				# editing taxon_name in row #
				my $new_taxon_name = $taxon_name;
				$new_taxon_name = join("__", $taxon_name, "cli$locus_id")
							if $multi_loci;
				$new_taxon_name = join("__", $new_taxon_name, $dna_segs_r->{$taxon_name}{$locus_id}{"subtype"})
							if $multi_subtype;
								
				# appending to new hash #
				$$row[$header_r->{"taxon_name"}] = $new_taxon_name;
				
				push(@{$new_dna_segs{$new_taxon_name}{$locus_id}{"entries"}}, $row);
				}
			}
		}
	
		#print Dumper %new_dna_segs; exit;
	return \%new_dna_segs;
	}

sub add_subtype{
# adding subtype label to leaves if needed #
	my ($dna_segs_r, $treeo, $multi_subtype) = @_;
	for my $node ($treeo->get_leaf_nodes){
		my @parts = split(/__/, $node->id);
		my $taxon_id  = join("__", @parts[0..($#parts-1)]);
		(my $locus_id = $parts[$#parts]) =~ s/^cli//;
		my $subtype = $dna_segs_r->{$taxon_id}{$locus_id}{"subtype"};
		
		$node->id( join("__", $node->id, $subtype)); 
		}
		#print Dumper $treeo; exit;
	return $treeo;
	}

sub add_leaves{
# adding leaves to any nodes that have multiple subtypes #
	my ($dna_segs_r, $treeo, $multi_loci) = @_;
	
	# splitting taxa if needed #
	for my $node ($treeo->get_leaf_nodes){
		my @loci = keys %{$dna_segs_r->{$node->id}};		# locus_ids 
		
		for my $i (0..$#loci){
			if ($i > 0){		# adding nodes if needed 	
				# getting node info #
				my $anc_node = $node->ancestor;
				my $brlen = $node->branch_length();
				(my $new_node_id = $node->id) =~ s/__[^_]+$//g;
									
				my $new_node = new Bio::Tree::Node(
							-id => $node->id,
							-branch_length => $brlen);					
				
				$new_node->id( join("__", $new_node_id, "cli$loci[$i]") );					
				$anc_node->add_Descendent($new_node);
				}
			else{	# appending locus_id #
				$node->id( join("__", $node->id, "cli$loci[$i]") );
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
	print STDERR "...Found multiple loci for 1 or more taxa. Adding leaves to the tree! Adding loci_ids to leaves & dna_segs table!\n"
		if $multi_loci;
	print STDERR "...Found multiple subtypes. Adding subtype to names in tree & dna_segs table!\n"
		if $multi_subtype;
		
	return $multi_loci, $multi_subtype;
	}

sub prune_tree{
# pruning tree by dna_segs, just dna_segs taxa remaining #
	my ($dna_segs_r, $treeo) = @_;

	# getting taxa no in tree #
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
			$treeo->remove_Node($node);
			}
		}
		
		
	# pruning any non-labeled leaves #
	while(1){
		my $intern_cnt = 0;
		for my $node ($treeo->get_leaf_nodes){
			if($node->id =~ /^INTERNAL__/){
				$treeo->remove_Node($node);
				$intern_cnt++;
				}
			}
		last unless $intern_cnt;		# continue until all internal node leaves are removed
		}	
	
	# editing internal node labels (names back) #
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
	while(<>){
		chomp;
		next if /^\s*$/;

		# header #
		if($.==1){
			my @tmp = split /\t/;
			for my $i (0..$#tmp){
				$header{$tmp[$i]} = $i;
				}
			die " ERROR: 'taxon_name' not found in dna_seg header!\n"
				unless exists $header{"taxon_name"};
			die " ERROR: 'locus_id' not found in dna_seg header!\n"
				unless exists $header{"locus_id"};				
			die " ERROR: 'locus_id' not found in dna_seg header!\n"
				unless exists $header{"subtype"};	
			}
		else{
			my @line = split /\t/;
			my $taxon_name = $line[$header{"taxon_name"}];
			my $locus_id = $line[$header{"locus_id"}];
			my $subtype = $line[$header{"subtype"}];
			
			# order of loci#
			push( @dna_segs_order, $taxon_name)
				unless exists $dna_segs{$taxon_name};
			
			# dna_segs  #
			push( @{$dna_segs{$taxon_name}{$locus_id}{"entries"}}, \@line );
			$dna_segs{$taxon_name}{$locus_id}{"subtype"} = $subtype;
			}
		}
		
		#print Dumper @dna_segs_order;
		#print Dumper %dna_segs; exit;
	return \%dna_segs, \@dna_segs_order, \%header;
	}
	
sub check_tree_format{
	my $format = shift;
	$format = "newick" if ! $format;
	$format =~ s/^new$/newick/i;
	$format =~ s/^nex$/nexus/i;
	die " Designated tree format ($format) not recognized.\n" if $format !~ /newick|nexus/;
	return $format;
	}



__END__

=pod

=head1 NAME

CLdb_dna_segs_orderByTree.pl -- needed if plotting a tree with CRISPR loci

=head1 SYNOPSIS

=head2 Ordering 'dna_segs' table

CLdb_dna_segs_orderByTree.pl [flags] < dna_segs.txt > dna_segs_ordered.txt

=head2 Ordering 'xlims' table

CLdb_dna_segs_orderByTree.pl [flags] < xlims.txt > xlims_ordered.txt


=head2 Required flags

=over

=item -tree

Tree file name (nexus or newick).

=back

=head2 Optional flags

=over

=item -format

Tree file format. [newick]

=item -name

Output file name for pruned tree. ['-tree' + '_prn.nwk']

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_dna_segs_orderByTree.pl

=head1 DESCRIPTION

Formatting tree and 'dna_segs' or 'xlims' table for plotting. 

The tree file is pruned to just the taxa found in
the dna_segs/xlims table.

The dna_segs/xlims tables are ordered to match the pruned
tree ordering.

Both xlims and dna_segs tables must be ordered by the tree
if the tree is used for plotting!

=head1 EXAMPLES

=head2 Basis usage

=head3 Ordering a dna_segs table

CLdb_dna_segs_orderByTree.pl -t tree.nwk < dna_segs.txt > dna_segs_ordered.txt

=head3 Ordering the complemenary xlims table

CLdb_dna_segs_orderByTree.pl -t tree.nwk < xlims.txt > xlims.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

