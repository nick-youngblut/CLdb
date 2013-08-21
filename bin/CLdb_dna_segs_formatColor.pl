#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use Bio::TreeIO;
use Color::Mix;
use List::Util qw/min max sum/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $tree_in, $format, $adjacency_opt, $write_brlen_stats);
my $gene_color_opt = 0;
my $spacer_color_opt = 3;
my $brlen_cutoff = 0;
my @default_colors = ("#666666", "#666666");
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "spacer=i" => \$spacer_color_opt, 		# 0 = none, 1 = tree, 2 = adjacency, 3 = both
	   "gene=i" => \$gene_color_opt,			# 0 = none, 1 = tree, 2 = adjacency, 3 = both
	   "branch=f" => \$brlen_cutoff,			# branch-length cutoff (>= max branch length)	   
	   "stats" => \$write_brlen_stats,
	   "default=s{1,2}" => \@default_colors,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$format = check_tree_format($format) if $tree_in;

# default colors #
my %default_colors;
$default_colors{"gene"} = $default_colors[0];
$default_colors{"spacer"} = $default_colors[1];


### Main
my %same_color; 		# colors thats will be all the same color

# loading dna_segs table #
my ($dna_segs_r, $dna_segs_order_r, $header_r) = load_dna_segs();

# tree-based color formatting #
## using brlen distance to remove similar colors ##
if($tree_in){
	my $treeio = Bio::TreeIO -> new(-file => $tree_in,
								-format => $format);
	my $treeo = $treeio->next_tree;
	
	# getting brlens #
	my $brlen_mat_r = make_brlen_matrix($treeo, $write_brlen_stats);
	add_upper($brlen_mat_r);
	
	# finding colors where max brlen between taxa are < cutoff #
	apply_brlen_cutoff($dna_segs_r, $brlen_mat_r, \%same_color, "spacer")
		if $spacer_color_opt == 1 || $spacer_color_opt == 3;
		
	apply_brlen_cutoff($dna_segs_r, $brlen_mat_r, \%same_color, "gene")
		if $gene_color_opt == 1 || $gene_color_opt == 3;	
	}

# adjacency-based color formatting #
check_adjacency($dna_segs_r, $dna_segs_order_r, \%same_color, "spacer")
	if $spacer_color_opt == 2 || $spacer_color_opt == 3;
check_adjacency($dna_segs_r, $dna_segs_order_r, \%same_color, "gene")
	if $gene_color_opt == 2 || $gene_color_opt == 3;

# editting colors #
## applying default 'same' color ##
($dna_segs_r, my $new_color_r) = apply_same_color($dna_segs_r, \%same_color);

## adding 'good' colors to rest ##
$dna_segs_r = apply_rainbow($dna_segs_r, $new_color_r);

# writing editted table #
write_dna_segs($dna_segs_r, $header_r, $dna_segs_order_r);



### Subroutines
sub write_dna_segs{
# writing dna_segs table #
	my ($dna_segs_r, $header_r, $dna_segs_order_r) = @_;
	
	# header #
	print join("\t", sort{$header_r->{$a}<=>$header_r->{$b}} keys %$header_r), "\n";
	
	# body #
	foreach my $dna_segs_id (@$dna_segs_order_r){
		foreach my $feat (keys %$dna_segs_r){
			foreach my $col (keys %{$dna_segs_r->{$feat}}){
				foreach my $row (@{$dna_segs_r->{$feat}{$col}{$dna_segs_id}}){
					$$row[$header_r->{"col"}] = $col;
					print join("\t", @$row), "\n";
					}
				}
			}
		}
	}

sub apply_rainbow{
	my ($dna_segs_r, $need_color_r) = @_;
	
	my %new_dna_segs;
	foreach my $feat (keys %$dna_segs_r){
		next if $feat eq "directrepeat";
		
		# getting hexa colors for color wheel #
		my $color = Color::Mix->new;
		my @cols = keys %{$need_color_r->{$feat}};
		my $ncol = scalar @cols;
		my @hexa = $color->analogous('0000ff', $ncol, $ncol);
		
		# formatting hexa #
		map{ $_ = "#$_"} @hexa;
		
		# color index #
		my %col_index;
		for my $i (0..$#cols){ $col_index{$cols[$i]} = $hexa[$i]; }
		
		# applying colors #
		foreach my $col (keys %{$dna_segs_r->{$feat}}){
			my $col2 = $col;
			$col2 = $col_index{$col} if exists $col_index{$col};
			$new_dna_segs{$feat}{$col2} = $dna_segs_r->{$feat}{$col};
			}
		}
		
		#print Dumper %new_dna_segs; exit;
	return \%new_dna_segs;
	}

sub apply_same_color{
# making colors all the same if needed #
	my ($dna_segs_r, $same_color_r) = @_;
		#print Dumper $same_color_r; exit;
		# same_color: feat=>color=>hexadecimal
	my %new_dna_segs;
	my %new_color;
	foreach my $feat (keys %$dna_segs_r){
		next if $feat eq "directrepeat";
		foreach my $col (keys %{$dna_segs_r->{$feat}}){
			my $col2 = $col;
			if(exists $same_color_r->{$feat}{$col}){		# applying non-descriminant color
				$col2 = $same_color_r->{$feat}{$col};
				}
			else{ $new_color{$feat}{$col} = 1; } 
				
			$new_dna_segs{$feat}{$col2} = $dna_segs_r->{$feat}{$col};	
			}	
		}
		#print Dumper %new_dna_segs; exit;
		#print Dumper %new_color; exit
	return \%new_dna_segs;
	}

sub check_adjacency{
# making adjacency list for taxa by color #
	my ($dna_segs_r, $dna_segs_order_r, $same_color_r, $feat) = @_;
	
	# making and checking adjacency matrix #
	foreach my $col (keys %{$dna_segs_r->{$feat}}){
		next if exists $same_color_r->{$feat}{$col};
		
		my @adjmtx;
		for my $i (0..$dna_segs_order_r){
			last if $i == $#$dna_segs_order_r;
			if(exists $dna_segs_r->{$feat}{$col}{$$dna_segs_order_r[$i]}
			   && exists $dna_segs_r->{$feat}{$col}{$$dna_segs_order_r[$i+1]}){
				$adjmtx[$i][$i+1] = 1;
				}
			else{ $adjmtx[$i][$i+1] = 0; }
			}
		# checking adjacencies #
		for my $i (0..$dna_segs_order_r){
			last if $i == $#$dna_segs_order_r;
			if($adjmtx[$i][$i+1]) {
				my $cnt = transverse_diag(\@adjmtx, $i);
				$same_color_r->{$feat}{$col} = $default_colors{$feat}
					if $cnt == scalar @adjmtx - $i;
				}
			}
		}
	
		#print Dumper $same_color_r; exit;
	}
	
sub transverse_diag{
# transversing diag of matrix; end when hit a zero/undef #
	my ($mtx_r, $i) = @_;
	
	my $cnt = 0;
	for my $ii ($i..$#$mtx_r){
		if($$mtx_r[$ii][$ii+1]){
			$cnt++;
			}
		else{
			last;
			}
		}
	return $cnt;
	}

sub apply_brlen_cutoff{
	my ($dna_segs_r, $brlen_mat_r, $same_color_r, $feat) = @_;
	
	#print Dumper %$brlen_mat_r; exit;
	
	foreach my $col (keys %{$dna_segs_r->{$feat}}){
		# finding max brlen distance #
		my @dna_segs_ids = keys %{$dna_segs_r->{$feat}{$col}};
		
		# getting brlen distances among all taxa in a color #
		my @dists;
		for my $i (0..$#dna_segs_ids){
			for my $ii (0..$#dna_segs_ids){
				next if $ii >= $i;
				die " ERROR: $dna_segs_ids[$i], $dna_segs_ids[$ii] not found in branch length matrix!\n"
					unless exists $brlen_mat_r->{$dna_segs_ids[$i]}{$dna_segs_ids[$ii]}; 
				push(@dists, $brlen_mat_r->{$dna_segs_ids[$i]}{$dna_segs_ids[$ii]});
				}
			}
		push(@dists, 0) unless $dists[0];
		
		if(max (@dists) < $brlen_cutoff){
			$same_color_r->{$feat}{$col} = $default_colors{$feat};
			}
		}
		#print Dumper %$same_color_r; exit;
	}

sub make_brlen_matrix{
# making a pairwise distance matrix of all taxa in tree #
	my ($treeo, $write_stats) = @_;
	
	my @nodes = $treeo->get_leaf_nodes;
	
	my %brlens;
	my @vals;
	for my $i (0..$#nodes){
		for my $ii (0..$#nodes){
			next if $ii >= $i;
			my $dist = $treeo->distance(@nodes[($i, $ii)]);
			$brlens{$nodes[$i]->id}{$nodes[$ii]->id} = $dist;
			push @vals, $dist;
				
			}
		}
		
	# stats #
	if($write_stats){
		print "Min branch length:\t", min(@vals), "\n";
		print "Mean branch length:\t", sum(@vals) / scalar @vals, "\n";
		print "Max branch length:\t", max(@vals), "\n";	
		print "Stdev branch length:\t", stdev(\@vals), "\n";
		exit;
		}
		
		#print Dumper %brlens; exit;
	return \%brlens;
	}

sub add_upper{
# adding upper half of matrix #
	my ($mat_r) = @_;
	foreach my $i (keys %$mat_r){
		foreach my $ii (keys %{$mat_r->{$i}}){
			$mat_r->{$ii}{$i} = $mat_r->{$i}{$ii};
			}
		}
		#print Dumper $mat_r; exit;
	}

sub load_dna_segs{
# loading dna_segs from stdin #
	my %dna_segs_order;
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
			
			my @check = qw/dna_segs_id col feat/;
			map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
				unless exists $header{$_} } @check;
			}
		else{
			my @line = split /\t/;
			my $dna_segs_id = $line[$header{"dna_segs_id"}];
			my $col = $line[$header{"col"}];
			my $feat = $line[$header{"feat"}];
	
			# order of loci#
			push( @dna_segs_order, $dna_segs_id)
				unless exists $dna_segs_order{$dna_segs_id};	
			$dna_segs_order{$dna_segs_id} = 1;
	
			
			# dna_segs to %@ #
			push ( @{$dna_segs{$feat}{$col}{$dna_segs_id}}, \@line );
			}
		}
		#print Dumper %dna_segs; exit;
		#print Dumper @dna_segs_order; exit;
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

sub average{
        my($data) = @_;
        if (not @$data) {
                die("Empty array\n");
        }
        my $total = 0;
        foreach (@$data) {
                $total += $_;
        }
        my $average = $total / @$data;
        return $average;
}

sub stdev{
        my($data) = @_;
        if(@$data == 1){
                return 0;
        }
        my $average = &average($data);
        my $sqtotal = 0;
        foreach(@$data) {
                $sqtotal += ($average-$_) ** 2;
        }
        my $std = ($sqtotal / (@$data-1)) ** 0.5;
        return $std;
}

__END__

=pod

=head1 NAME

CLdb_dna_segs_formatColor.pl -- reduce the number of colors needed for plotting

=head1 SYNOPSIS

CLdb_dna_segs_formatColor.pl [flags] < dna_segs_order.txt > dna_segs_order_color.txt

=head2 Optional flags

=over

=item -tree

Tree file name (nexus or newick).

=item -format

Tree file format. [newick]

=item -spacer

Methods used to reduce colors needed for spacer clusters (see DESCRIPTION).

=item -gene

Methods used to reduce colors needed for gene clusters (see DESCRIPTION).

=item -branch

The branch length cutoff for not coloring features just in that clade.

=item -stats

Write branch length stats for the provided tree, then exit? [FALSE]

=item -default

The default color used when discrimant coloring is not needed (needs 2 values: 'spacer_color' 'gene_color').

=item -v 	Verbose output. [FALSE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_dna_segs_formatColor.pl

=head1 DESCRIPTION


=head1 EXAMPLES

=head2 Basic usage

CLdb_dna_segs_formatColor.pl 

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

