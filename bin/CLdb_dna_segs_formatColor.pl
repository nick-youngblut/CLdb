#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_dna_segs_formatColor.pl -- reduce the number of colors needed for plotting

=head1 SYNOPSIS

CLdb_dna_segs_formatColor.pl [flags] < dna_segs_order.txt > dna_segs_order_color.txt

=head2 Optional flags

=over

=item -tree  <char>

Tree file name (nexus or newick).

=item -format  <char>

Tree file format. [newick]

=item -spacer  <int>

Methods used to reduce colors needed for spacer clusters (see DESCRIPTION).

=item -gene  <int>

Methods used to reduce colors needed for gene clusters (see DESCRIPTION).

=item -branch  <float>

The branch length cutoff for not coloring features just in that clade.

=item -stats  <bool>

Write branch length stats for the provided tree, then exit? [FALSE]

=item -default  <char> 

The default colors for when coloring isn't needed (needs 3 values: 'gene_color' 'spacer_color' 'DR_color').

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_dna_segs_formatColor.pl

=head1 DESCRIPTION

Add descriminatory coloring to genes and spacers.
Genes and/or spacers that do not need coloring
because: 
=over 

=item i) the gene/spacer always adjacent
in the plot, so 'blast' connections in the 
plot will provide the needed information.

=item ii) the gene/spacer is only found in
a very similar clade, so no coloring is needed
(if a tree & branch length cutoff is provided).

=back

Discrimantory coloring can be based on:

=over 

=item Continuous adjacency in the plot: [1]

=item Max branch length: [2]

=item Both: [3]

=back

By default: genes = 1, spacers = 3

The default non-descriminatory colors are:

=over

=item gene: #666666

=item spacer: #CCCCCC

=item direct repeat: #000000

=back 

=head2 WARNING

If the dna_segs table has been ordered by a tree,
use the pruned/edited tree.

=head1 EXAMPLES

=head2 Basic usage:

CLdb_dna_segs_formatColor.pl < dna_segs_order.txt > dna_segs_order_col.txt

=head2 Branch length cutoff (<= 0.1):

CLdb_dna_segs_formatColor.pl -t tree.nwk -b 0.1 < dna_segs_order.txt > dna_segs_order_col.txt


=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

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
use Color::Mix;
use List::Util qw/min max sum/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
	table_exists
	n_entries
	join_query_opts/;
use CLdb::utilities qw/
	file_exists 
	connect2db/;
	

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));


my ($verbose, $tree_in, $format, $write_brlen_stats, $compare_in);
my $gene_color_opt = 0;
my $spacer_color_opt = 3;
my $brlen_cutoff = 0;
my @default_colors = ("#666666", "#CCCCCC", "#000000");		# gene, spacer, DR
GetOptions(
	   "tree=s" => \$tree_in,
	   "format=s" => \$format,
	   "compare=s" => \$compare_in,
	   "spacer=i" => \$spacer_color_opt, 		# 0 = none, 1 = tree, 2 = adjacency, 3 = both
	   "gene=i" => \$gene_color_opt,			# 0 = none, 1 = tree, 2 = adjacency, 3 = both
	   "branch=f" => \$brlen_cutoff,			# branch-length cutoff (>= max branch length)	   
	   "stats" => \$write_brlen_stats,
	   "default=s{3}" => \@default_colors,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
$format = check_tree_format($format) if $tree_in;

# default colors #
my %default_colors;
$default_colors{"gene"} = $default_colors[0];
$default_colors{"spacer"} = $default_colors[1];
$default_colors{"directrepeat"} = $default_colors[2];

### Main
my %color_mod; 			# how colors will be modified

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
	apply_brlen_cutoff($dna_segs_r, $brlen_mat_r, \%color_mod, "spacer")
		if $spacer_color_opt == 1 || $spacer_color_opt == 3;
		
	apply_brlen_cutoff($dna_segs_r, $brlen_mat_r, \%color_mod, "gene")
		if $gene_color_opt == 1 || $gene_color_opt == 3;	
	}

# adjacency-based color formatting #	
check_adjacency($dna_segs_r, $dna_segs_order_r, \%color_mod, "spacer", $header_r)
	if $spacer_color_opt == 2 || $spacer_color_opt == 3;
check_adjacency($dna_segs_r, $dna_segs_order_r, \%color_mod, "gene", $header_r)
	if $gene_color_opt == 2 || $gene_color_opt == 3;
#	}

# editting colors #
## finding colors that need descriminating hexadecimal coloring ##
my $descrim_cnt_r = find_descrim_color($dna_segs_r, \%color_mod);

## adding 'good' colors to rest ##
apply_rainbow(\%color_mod, $descrim_cnt_r);

## applying default color to repeats ##
apply_DR_color(\%color_mod);

# writing editted table #
write_dna_segs($dna_segs_r, $header_r, $dna_segs_order_r, \%color_mod);


### Subroutines
sub write_dna_segs{
# writing dna_segs table #
	my ($dna_segs_r, $header_r, $dna_segs_order_r, $color_mod_r) = @_;
	
	# header #
	print join("\t", sort{$header_r->{$a}<=>$header_r->{$b}} keys %$header_r), "\n";
	
	# body #
	foreach my $dna_segs_id (@$dna_segs_order_r){
		foreach my $feat (keys %$dna_segs_r){
			foreach my $col (keys %{$dna_segs_r->{$feat}}){
				foreach my $row (@{$dna_segs_r->{$feat}{$col}{$dna_segs_id}}){
					# sanity check #
					die " ERROR: $feat -> $col not found in color index! $!\n"
						unless exists  $color_mod_r->{$feat}{$$row[$header_r->{"col"}]};
					
					# applying modified colors (hexadecimal) #
					$$row[$header_r->{"col"}] = $color_mod_r->{$feat}{$$row[$header_r->{"col"}]};
					print join("\t", @$row), "\n";
					}
				}
			}
		}
	}

sub apply_DR_color{
# applying default to direct repeat #
	my $color_mod_r = shift;	
	foreach my $col ( keys %{$color_mod_r->{"directrepeat"}} ){
		$color_mod_r->{"directrepeat"}{$col} = $default_colors{"directrepeat"};
		}
	}

sub apply_rainbow{
# applying rainbow to color mod #
	my ($color_mod_r, $descrim_cnt_r) = @_;
		
		
		#print Dumper $color_mod_r; exit;
		
	my %new_dna_segs;
	foreach my $feat (keys %$color_mod_r){
		next if $feat eq "directrepeat";
		
		# getting hexa colors for color wheel #
		my $color = Color::Mix->new;
		my @hexa = $color->analogous('0000ff', 
							$descrim_cnt_r->{$feat}, 
							$descrim_cnt_r->{$feat});
		
		# formatting hexa #
		map{ $_ = "#$_"} @hexa;
		
		# applying colors #
		my $cnt = 0;
		if($feat eq 'gene'){
			foreach my $col (sort {$a cmp $b} keys %{$color_mod_r->{$feat}}){
				unless($color_mod_r->{$feat}{$col}){
					$color_mod_r->{$feat}{$col} = $hexa[$cnt];
					$cnt++;
					}
				}			
			}
		else{
			foreach my $col (sort {$a<=>$b} keys %{$color_mod_r->{$feat}}){
				unless($color_mod_r->{$feat}{$col}){
					$color_mod_r->{$feat}{$col} = $hexa[$cnt];
					$cnt++;
					}
				}
			}
		}
		#print Dumper "here", $color_mod_r; exit;		
	}

sub find_descrim_color{
	my ($dna_segs_r, $color_mod_r) = @_;
	
	my %descrim_cnt;
	foreach my $feat (keys %$dna_segs_r){
		foreach my $col (keys %{$dna_segs_r->{$feat}}){
			unless (exists $color_mod_r->{$feat}{$col}){
				$color_mod_r->{$feat}{$col} = 0;
				$descrim_cnt{$feat}++;
				}
			}
		}
		#print Dumper $color_mod_r; exit; 
	return \%descrim_cnt;
	}

sub check_adjacency{
# making adjacency list for taxa by color #
	my ($dna_segs_r, $dna_segs_order_r, $color_mod_r, $feat, $header_r) = @_;

	#print Dumper $dna_segs_r->{"spacer"}; exit;

	# making and checking adjacency list #
	foreach my $col (keys %{$dna_segs_r->{$feat}}){			# for each spacer, gene, or DR
		#next if exists $color_mod_r->{$feat}{$col};		# next if color already assigned
		
		my %adjlist;
		for my $i (0..($#$dna_segs_order_r-1)){				# order
			my $dna_segs_id1 = $$dna_segs_order_r[$i];
			my $dna_segs_id2 = $$dna_segs_order_r[$i+1];			
			
			if(exists $dna_segs_r->{$feat}{$col}{$dna_segs_id1} &&		# adjacent taxa found in same color
				exists $dna_segs_r->{$feat}{$col}{$dna_segs_id2}){
				$adjlist{$dna_segs_id1} = $dna_segs_id2;				# color must be found in adjacent 
				#print Dumper $dna_segs_id1, $dna_segs_id2
				#	if $dna_segs_id1 eq "3.F.A.2.12__cli74" || $dna_segs_id2 eq "3.F.A.2.12__cli74";
				}
			elsif(exists $dna_segs_r->{$feat}{$col}{$dna_segs_id1}){	# if 'loner'
				$adjlist{$dna_segs_id1} = 0;
				}
			}
			#print Dumper $col, %adjlist;
			
		# checking adjacencies #
		my $adjcnt=0;			# number of adjacency strings
		my $i = 0;
		my $adjlen = 0;											# total adjacency length; for checking 'loners'
		while(1){												# finding 1st adjacency; moving through dna_seg_order
			last if $i >= $#$dna_segs_order_r;
			if(exists $adjlist{$$dna_segs_order_r[$i]}){		# found adjacency or loner; following
				$adjlen++;										# length of adjacency
				
				# following adjacency #
				while(1){										# last in adjacency
					$i++;
					if(exists $adjlist{$$dna_segs_order_r[$i]}){
						$adjlen++ if $adjlist{$$dna_segs_order_r[$i]};				
						}
					else{
						last;
						}
					}
				# counting adjacency if len > 1 #
				$adjcnt++;
				}
			$i++;
			}	
		#print Dumper "adjcnt: $adjcnt; adjlen: $adjlen";
		$color_mod_r->{$feat}{$col} = $default_colors{$feat} if $adjcnt <= 1;
		}
		
		#exit;
		#print Dumper %{$color_mod_r->{"spacer"}}; exit;
		#print Dumper %debug; exit;
		#print Dumper scalar keys %{$color_mod_r->{"spacer"}}; exit;
	}

sub check_adjacency_old{
# making adjacency list for taxa by color #
	my ($dna_segs_r, $dna_segs_order_r, $compare_r, $color_mod_r, $feat, $header_r) = @_;
	
		#print Dumper $dna_segs_order_r; exit;
	
	# making and checking adjacency list #
	foreach my $col (keys %{$dna_segs_r->{$feat}}){		# for each spacer, gene, or DR
		#next if exists $color_mod_r->{$feat}{$col};		# next if color already assigned
		
		# making adjacency list #
		my %adjlist;
		for my $i (0..($#$dna_segs_order_r-1)){			# order
			my $dna_segs_id1 = $$dna_segs_order_r[$i];
			my $dna_segs_id2 = $$dna_segs_order_r[$i+1];			
			
			# checking comparisons #
			foreach my $row (@{$dna_segs_r->{$feat}{$col}{$dna_segs_id1}}){
				my $feat_id = $$row[$header_r->{"feat_id"}];
				
				if(exists $compare_r->{$feat}{$dna_segs_id1}{$dna_segs_id2}{$feat_id}){ 	# feature must exist in comparison 
					$adjlist{$dna_segs_id1} = $dna_segs_id2;								# color must be found in adjacent 
					}
				}
			}
		
			#print Dumper %adjlist;		# zero means no adjacency
		
		# checking adjacencies #
		my $adjcnt=0;			# number of adjacency strings
		my $i = 0;
		my $adjlen = 0;											# total adjacency length
		while(1){												# finding 1st adjacency; moving through dna_seg_order
			last if $i >= $#$dna_segs_order_r;
			if(exists $adjlist{$$dna_segs_order_r[$i]}){				# found adjacency or loner; following
					#print Dumper "here", $adjlist{$$dna_segs_order_r[$i]};
					#$debug{$col} = 1 if $adjlist{$$dna_segs_order_r[$i]} eq "3.F.A.2.12__cli74";
					#$alone++ unless exists $adjlist{$$dna_segs_order_r[$i+1]};		# counting 'loner' arrays; no adjacency
				$adjlen++;										# length of adjacency
				
				# following adjacency #
				while(1){										# last in adjacency
					$i++;
					if(exists $adjlist{$$dna_segs_order_r[$i]}){
						$adjlen++;								
						}
					else{
						last;
						}
					}
				# counting adjacency if len > 1 #
				$adjcnt++;
				
				}
			$i++;
				#$debug{$col} = $i if $debug{$col};
			}
			#print Dumper "adj: $adjcnt";
			#print Dumper "alone: $alone";
		
		# modifying to same color (one continued adjacency will be colored default color) #
		## 1 adjacency & adjacency length must be >1 ##
		$color_mod_r->{$feat}{$col} = $default_colors{$feat} if $adjcnt == 1 && $adjlen > 1;
		}
		
		#exit;
		#print Dumper %debug; exit;
		#print Dumper scalar keys %{$color_mod_r->{"spacer"}}; exit;
	}

sub load_compare{
	my ($compare_in) = @_;
	
	open IN, $compare_in or die $!;
	my %compare;
	my %header;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		if($.==1){ 	#header
			my @tmp = split /\t/;
			for my $i (0..$#tmp){
				$header{$tmp[$i]} = $i;
				}
			
			my @check = qw/dna_segs_id1 dna_segs_id2 feat feat_id/;
			map{ die " ERROR: \"$_\" not found in dna_seg header!\n"
				unless exists $header{$_} } @check;
			}
		else{   	#body
			my @line = split /\t/;
			my $dna_segs_id1 = $line[$header{"dna_segs_id1"}];
			my $dna_segs_id2 = $line[$header{"dna_segs_id2"}];
			my $feat = $line[$header{"feat"}];
			my $feat_id = $line[$header{"feat_id"}];
	
			# compare #
			push ( @{$compare{$feat}{$dna_segs_id1}{$dna_segs_id2}{$feat_id}}, \@line );
			}	
		}
	close IN;
	
		#print Dumper %compare; exit;
	return \%compare;
	}

sub apply_brlen_cutoff{
	my ($dna_segs_r, $brlen_mat_r, $color_mod_r, $feat) = @_;
	
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
			$color_mod_r->{$feat}{$col} = $default_colors{$feat};
			}
		}
		#print Dumper %$color_mod_r; exit;
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



