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
use Color::Mix;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $tree_in, $format, $adjacency_opt, $write_brlen_stats);
my $gene_color_opt = 0;
my $brlen_cutoff = 0;
my @default_colors = ("#666666", "#FF9933");
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "adjacency" => \$adjacency_opt, 
	   "gene=i" => \$gene_color_opt,
	   "branch=f" => \$brlen_cutoff,			# branch-length cutoff (>= max branch length)	   
	   "stats" => \$write_brlen_stats,
	   "default=s{1,2}" => \@default_colors,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$format = check_tree_format($format) if $tree_in;

### Main
# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r) = load_dna_segs();

# tree-based color formatting #
## using brlen distance to remove similar colors ##
if($tree_in){
	my $treeio = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeo = $treeio->next_tree;
	
	}




### Subroutines
sub load_dna_segs{
# loading dna_segs from stdin #
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
			die " ERROR: 'col' not found in dna_seg header!\n"
				unless exists $header{"col"};				
			}
		else{
			my @line = split /\t/;
			my $taxon_name = $line[$header{"taxon_name"}];
			my $col = $line[$header{"col"}];
			
			
			# dna_segs to %@ #
			push {@{$dna_segs{}} = \@line );		# taxon_name => row
			}
		}
		#print Dumper %dna_segs; exit;
	return \%dna_segs, \%header;
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

