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
tree_write($treeo, $tree_name);

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

sub order_dna_segs{
	my ($dna_segs_r, $tree_order_r, $header_r) = @_;
	
	# header #
	print join("\t", sort{$header_r->{$a}<=>$header_r->{$b}} keys %$header_r), "\n";
	
	# body #
	foreach my $leaf (@$tree_order_r){
		die " ERROR: leaf -> \"$leaf\" not found in dna_segs table!\n"
			unless exists $dna_segs_r->{$leaf};
		foreach my $row (@{$dna_segs_r->{$leaf}}){
			print join("\t", @$row), "\n";
			}
		}

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
			}
		else{
			my @line = split /\t/;
			my $taxon_name = $line[$header{"taxon_name"}];
			
			# order #
			push( @dna_segs_order, $taxon_name)
				unless exists $dna_segs{$taxon_name};
			
			# dna_segs to %@ #
			push( @{$dna_segs{$taxon_name}}, \@line );
			}
		}
		
		#print Dumper @dna_segs_order;
		#print Dumper %dna_segs; exit;
	return \%dna_segs, \@dna_segs_order, \%header;
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

CLdb_dna_segs_formatColor.pl -- needed if plotting a tree with CRISPR loci

=head1 SYNOPSIS

CLdb_dna_segs_formatColor.pl [flags] > dna_segs_ordered.txt

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

perldoc CLdb_dna_segs_formatColor.pl

=head1 DESCRIPTION

Formatting tree and dna_segs table for plotting. 

The tree file is pruned to just the taxa found in
the dna_segs table.

The dna_segs file is ordered to match the pruned
tree ordering.

=head1 EXAMPLES

=head2 Basic usage

CLdb_dna_segs_formatColor.pl -t tree.nwk < dna_segs.txt > dna_segs_ordered.txt

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

