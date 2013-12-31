#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_positionByScaffold.pl -- Convert positional values from a scaffold-merged (1 scaffold) to multi-scaffold 

=head1 SYNOPSIS

CLdb_positionByScaffold.pl [Flags] < loci/Loci.txt > loci/Loci_byScaf.txt

=head2 Required flags

=over

=item -genbank  <char>

List of genbank files (3-column; tab-delim; 'taxon_name\tgenbank_unmerged\tgenbank_merged')

=back

=head2 Optional flags

=over

=item -array  <char>

Directory of the CLdb array files.

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>	

This help message

=back

=head2 For more information:

perldoc CLdb_positionByScaffold.pl

=head1 DESCRIPTION

Convert merged genbank position information (merged with EMBOSS 'union')
to scaffold-based (multiple scaffold) position info. 
This is designed to be used before loading the loci table, genbank files,
and array files into CLdb.
It is only necessary if you would like to have scaffold position information
in CLdb.

=head2 Input

=head3 Loci table

All of the position info will be changed to scaffold (as long as
the corresponding merged & unmerged genbank files are provided).
Also, the scaffold ID for the position information will be added
to the editted loci table.

=head3 genbank list

3-column; tab-delim; 'taxon_name\tgenbank_unmerged\tgenbank_merged'

Example of unmerged genbank file format: RAST genbank output for a draft genome.

=head3 array directory

This is used to determine where the array files listed in the loci table
are located. The edited array files will be written to a seperate directory. 

=head2 Matching unmerged and merged scaffolds

All scaffold start-end info is obtained form the 'source' tags
in both the unmerged and merged genbank files. 1-to-1
connections between scaffolds are determined in a number of steps:

=over 

=item 1) Comparing scaffold lengths.

=item 2) Comparing scaffold sequences. No duplicate scaffolds allowed!

=back

=head1 EXAMPLES

=head2 Normal usage (array files editted):

CLdb_positionByScaffold.pl < loci.txt -genbank genbank_file_list.txt -array > loci_byScaf.txt

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
use Bio::SeqIO;
use Set::IntervalTree;
use File::Path qw/rmtree/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::utilities qw/
	file_exists 
	connect2db/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose_b, $array_dir, $gbk_list_in);
GetOptions(
	   "array=s" => \$array_dir,
	   "genbank=s" => \$gbk_list_in,
	   "verbose_b" => \$verbose_b,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($gbk_list_in, "genbank list");

# array directory #
my $array_dir_out;
if(defined $array_dir){
	die "ERROR: cannot find directory of array files!\n"
	unless -d $array_dir;
	($array_dir_out = $array_dir) =~ s|(.+)[\/]$|$1_byScaf|;
	die "REGEX ERROR: $!\n" if $array_dir_out eq $array_dir;
	rmtree $array_dir_out or die $! if -d $array_dir_out;
	mkdir $array_dir_out;
	}

#--- MAIN ---#
# loading input lists/tables #
my $gbk_list_r = load_genbank_list($gbk_list_in);
my $loci_tbl_r = load_loci_table();

# edited loci table header #
my @header_order = sort{$loci_tbl_r->{'header'}{$a} <=> $loci_tbl_r->{'header'}{$b}} 
						keys %{$loci_tbl_r->{'header'}};
print join("\t", @header_order), "\n";
my @indices = @{$loci_tbl_r->{'header'}}{@header_order};

# processing each taxon #
foreach my $taxon_name (keys %$gbk_list_r){
	# genbank #
	## skipping unless present in loci file ##
	next unless exists $loci_tbl_r->{'body'}{$taxon_name};
	
	## getting genbank bioperl objects #
	my $gbk_r = load_genbank_io($gbk_list_r->{$taxon_name});
	
	## getting scaffold lengths & sequence #
	my $unmerged_r = parse_genbank($gbk_r->{'unmerged'}, 
									'unmerged', $gbk_list_r->{$taxon_name}{'unmerged'});
	my $merged_r = parse_genbank($gbk_r->{'merged'}, 
									'merged', $gbk_list_r->{$taxon_name}{'merged'});

	# position association #
	## finding merged-unmerged contig 1-to-1 assocations in scaffold length #
	my ($scaf_index_r, $needs_blasting_r) = find_same_length($unmerged_r, $merged_r);

	## determine if sequences are the same (either identical or rev-comp)
	find_same_sequence($scaf_index_r, $needs_blasting_r, $unmerged_r, $merged_r);

	## checking to make sure all scaffolds accounted for #
	if(%{$needs_blasting_r->{"merged"}}){
		print STDERR " ERROR: not all 1-to-1 connects made. Ambiguous contigs:\n\t";
		print STDERR join(",\n\t", keys %{$needs_blasting_r->{"unmerged"}}), "\n\n";
		exit(1);
		}
	
	## making an interval tree relating absolute location to scaffold-based location #
	my $itree = make_loc_itree($scaf_index_r, $merged_r, $unmerged_r);

	# location editing #
	## loci table ##	
	foreach my $locus ( @{$loci_tbl_r->{'body'}{$taxon_name}} ){		# locus = @$row
		# loci #
		edit_loci_tbl_loc( $loci_tbl_r->{'header'}, $itree, $locus, 
							$taxon_name, 'locus_start', 'locus_end', 'scaffold');
		# array #
		edit_loci_tbl_loc( $loci_tbl_r->{'header'}, $itree, $locus,
							$taxon_name, 'array_start', 'array_end', '');
		# CAS #
		edit_loci_tbl_loc( $loci_tbl_r->{'header'}, $itree, $locus,
							$taxon_name, 'CAS_start', 'CAS_end', '');
		# leader #
		edit_loci_tbl_loc( $loci_tbl_r->{'header'}, $itree, $locus,
							 $taxon_name, 'leader_start', 'leader_end', '');
				
		# writing out new loci file #
		print join("\t", @$locus[ @indices ] ), "\n";
						 
		# array file #
		edit_array_file($array_dir, $array_dir_out, $locus, $loci_tbl_r->{'header'}, $itree, $taxon_name)
				if defined $array_dir;
		}
	
	}


### Subroutines
sub edit_array_file{
# editing array file (changing positions to scaffold-relative positions) #
	my ($array_dir, $array_dir_out, $locus, $header_r, $itree, $taxon_name) = @_;
	
	# cannot edit array file unless file exists #
	return 0 unless exists $header_r->{'array_file'};
	my $array_i = $header_r->{'array_file'};
	
	# check for array file value #
	unless (defined $$locus[$array_i]){
		warn "WARNING: cannot find array file value for $taxon_name loci\n";
		return 0;
		}
		
	# making an array file edit #
	my @parts = File::Spec->splitpath($$locus[$array_i]);
	
	open IN, "$array_dir/$parts[2]" or die $!;
	open OUT, ">$array_dir_out/$parts[2]" or die $!;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;			# start, DR, spacer, end 
		
		# getting by-scaffold position #
		my $ret = $itree->fetch($l[0], $l[3]);		# ret = scaffold that positions span
		die "ERROR: could not find position info for: $taxon_name, ", join(", ", @l), "\n"
			unless defined $ret;
		
		# updating position info #
		## new start = scaf_start + (current_pos - scaf_gw_pos)
		## new end = scaf_end - (scaf_gw_pos - current_pos)
		$l[0] = $$ret[0]->{"start"} + ($l[0] - $$ret[0]->{"gw_start"});
		$l[3] = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $l[3]);
					
		print OUT join("\t", @l), "\n";
		}
	close IN;
	close OUT;
	
	# status #
	print STDERR "Edited array file for written: '$array_dir_out/$parts[2]'\n"
		unless $verbose_b;
	}

sub edit_loci_tbl_loc{

	my ($header, $itree, $locus, $taxon_name, $start_cat, $end_cat, $scaf_cat) = @_;
		
	# locus start-end #
	## index #
	my $start_i = $header->{$start_cat}
			if exists $header->{$start_cat};
	my $end_i = $header->{$end_cat}
			if exists $header->{$end_cat};
	my $scaf_i = $header->{$scaf_cat}
			if exists $header->{$scaf_cat};
	
	return 0 unless exists $header->{$start_cat} && exists $header->{$end_cat};
	
	## edit start-end values #
	if( defined $$locus[ $start_i ] && defined $$locus[ $end_i ] ){
		my $start = $$locus[ $start_i ];
		my $end = $$locus[ $end_i ];
	
		if( defined $scaf_i ){			# scaffold info
			($$locus[ $start_i ], 	
 			 $$locus[ $end_i ],
 			 $$locus[ $scaf_i ]) = 
 					edit_loc($itree, $taxon_name, $locus,
 							$start, $end
 							);
 				}
 		else{				# no scaffold info 
 			($$locus[ $start_i ], 
 			 $$locus[ $end_i ]) = 
 					edit_loc($itree, $taxon_name, $locus,
 							$start, $end
 							);	
 			}
 		}
		#print Dumper $loci_tbl_r; exit;
	}

sub edit_loc{
# editing start-end location # 
	my ($itree, $taxon_name, $locus, $start, $end) = @_;
	
	# checking start-end orientation; flipping if - strand for itree fetching #
	my $flip_b = 0;
	if($start > $end){
		$flip_b = 1;
		($start, $end) = ($end, $start);
		}
	
	# getting position info (in unmerged) #
	my $ret = $itree->fetch($start, $end);				# ret = scaffold that positions span
	die "ERROR: no scaffold span for taxon_name->$taxon_name, start->$start, end->$end\n"
		unless defined $ret;
	die "ERROR: >1 scaffold spanned for taxon_name->$taxon_name, start->$start, end->$end\n"
		if scalar @$ret > 1;

	# updating position info #
	## new start = scaf_start + (current_pos - scaf_gw_pos)
	## new end = scaf_end - (scaf_gw_pos - current_pos)
	$start = $$ret[0]->{"start"} + ($start - $$ret[0]->{"gw_start"});
	$end = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $end);
	
	# flipping start-end back if necessary #
	($start, $end) = ($end, $start) if $flip_b;

	# check start stop #
	die " ERROR: scaffold start or end < 0 for taxon_name->$taxon_name, start->$start, end->$end\n"
		if $start < 0 || $end < 0;
			
	# scaffold #
	my $scaffold_name = $$ret[0]->{"id"};

		#print Dumper $start, $end, $scaffold_name; exit;
	return $start, $end, $scaffold_name;
	}

sub merged_to_unmerged_pos{
# getting scaffold that positions fall into #
	my ($itree, $loci_r, $tables_oi_r) = @_;
	
	# updating per-loci #
	foreach my $locus (@$loci_r){			# $$locus[0] = locus_id
		my $scaffold_name;
		# locus table #
		for (my $i=2; $i<=6; $i+=2){
			# skipping if no values #
			next unless $$locus[$i] && $$locus[$i+1];
			
			# checking start-end orientation #
			my $flip_bool = 0;
			$flip_bool = 1 if $$locus[$i] > $$locus[$i+1];
			($$locus[$i], $$locus[$i+1]) = flip_se($$locus[$i], $$locus[$i+1]) if $flip_bool;
			
			# getting position info #
			my $ret = $itree->fetch($$locus[$i], $$locus[$i+1]);		# ret = scaffold that positions span
			check_scaffold_span($ret, $locus, $$locus[$i], $$locus[$i+1], "loci");						# positions should span 1 scaffold

			# updating position info #
			## new start = scaf_start + (current_pos - scaf_gw_pos)
			## new end = scaf_end - (scaf_gw_pos - current_pos)
			$$locus[$i] = $$ret[0]->{"start"} + ($$locus[$i] - $$ret[0]->{"gw_start"});
			$$locus[$i+1] = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $$locus[$i+1]);
			
			# flipping start-end back if necessary #
			($$locus[$i], $$locus[$i+1]) = flip_se($$locus[$i], $$locus[$i+1]) if $flip_bool;
			
			# check start stop #
			die " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->loci; start = $$locus[$i], end = ", $$locus[$i+1], "\n"
				#print STDERR " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->loci; start = $$locus[$i], end = ", $$locus[$i+1], "\n"
				if $$locus[$i] < 0 || $$locus[$i+1] < 0;
			
			# scaffold #
			$scaffold_name = $$ret[0]->{"id"} unless $scaffold_name;
			}
		$$locus[1] = $scaffold_name;
		
		# other tables #
		foreach my $table (keys %$tables_oi_r){
			foreach my $row ( @{$tables_oi_r->{$table}{$$locus[0]}} ){		 # 2 value
				# skipping if not present #
				next unless $$row[0] && $$row[1];

					#print Dumper @$row;
				
				# checking start-end orientation #
				my $flip_bool = 0;
				$flip_bool = 1 if $$row[0] > $$row[1];
				($$row[0], $$row[1]) = flip_se($$row[0], $$row[1]) if $flip_bool;
				
				# getting position info #
				my $ret = $itree->fetch($$row[0], $$row[1]);		# ret = scaffold that positions span
				check_scaffold_span($ret, $locus, $$row[0], $$row[1], $table);						# positions should span 1 scaffold

				# updating position info #
				## new start = scaf_start + (current_pos - scaf_gw_pos)
				## new end = scaf_end - (scaf_gw_pos - current_pos)
				$$row[0] = $$ret[0]->{"start"} + ($$row[0] - $$ret[0]->{"gw_start"});
				$$row[1] = $$ret[0]->{"end"} - ($$ret[0]->{"gw_end"} - $$row[1]);
			
				# flipping start-end back if necessary #
				($$row[0], $$row[1]) = flip_se($$row[0], $$row[1]) if $flip_bool;

				# check start stop #
				die " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->$table; start = $$row[0], end = $$row[1]\n"
					#print STDERR " ERROR: scaffold start or end < 0 for locus->$$locus[0], table->$table; start = $$row[0], end = $$row[1]\n"
					if $$row[0] < 0 || $$row[1] < 0;
		
					#print Dumper @$row;
				}
			}
		}
	
		#print Dumper @$loci_r;
		#print Dumper %$tables_oi_r; 
		#exit;
	}

sub make_loc_itree{
# making a location itree #
## input: genome-wide loc; output: unmerged scaf loc ##
	my ($scaf_index_r, $merged_r, $unmerged_r) = @_;
	
	my $itree = Set::IntervalTree->new();

	foreach my $merged_contig (keys %$scaf_index_r){
		foreach my $unmerged_contig (@{$scaf_index_r->{$merged_contig}}){
			# addig genome-wide position info to unmerged contigs #
			$unmerged_r->{$unmerged_contig}{"gw_start"} = $merged_r->{$merged_contig}{"start"};
			$unmerged_r->{$unmerged_contig}{"gw_end"} = $merged_r->{$merged_contig}{"end"};
			
			# deleting sequence info #
			delete $unmerged_r->{$unmerged_contig}{"seq"};
			
			# adding to interval tree #
			$itree->insert($unmerged_r->{$unmerged_contig},
						$merged_r->{$merged_contig}{"start"},
						$merged_r->{$merged_contig}{"end"});
			}
		}

	#print_itree($itree, 1, 10000);
	#my $res_r = $itree->fetch(200, 300);
	#print Dumper $res_r; 
	
	return $itree;
	}

sub print_itree{
	my ($itree, $start, $end) = @_;
	
	for my $i ($start..$end){
		my $res = $itree->fetch($i, $i);
		
		print join("\t", "pos:$i", join("::", @$res)), "\n";
		}
	}

sub find_same_sequence{
# checking contigs by same sequence; should match all contigs #
## sequences should match exactly & 1-to-1 ##
	my ($scaf_index_r, $needs_blasting_r, $unmerged_r, $merged_r) = @_;

	foreach my $merged_contig (keys %{$needs_blasting_r->{"merged"}}){
		foreach my $unmerged_contig (keys %{$needs_blasting_r->{"unmerged"}}){
			if ($merged_r->{$merged_contig}{"seq"} eq $unmerged_r->{$unmerged_contig}{"seq"}){
				push( @{$scaf_index_r->{$merged_contig}}, $unmerged_contig);
				}
			}
		}
	
	# checking problem contigs for 1 to 1 associations #
	foreach my $merged_contig (keys %{$needs_blasting_r->{"merged"}}){
		if( exists $scaf_index_r->{$merged_contig} ){
			if(scalar @{$scaf_index_r->{$merged_contig}} == 1){
				# deleting merged #
				delete $needs_blasting_r->{"merged"}{$merged_contig};
			
				# deleting unmerged #
				foreach my $unmerged_contig (@{$scaf_index_r->{$merged_contig}}){
					delete $needs_blasting_r->{"unmerged"}{$unmerged_contig};
					}
				}
			else{
				die " ERROR: Scaffold->\"$merged_contig\" has same sequence as >=2 scaffolds in unmerged genbank!\n";
				}
			}
		else{		# no sequenc hit
			die " ERROR: Scaffold->\"$merged_contig\" has no identical sequence in unmerged genbank!\n";
			}
		}
	}

sub find_same_length{
# foreach merged-unmerged, find scaffolds that match in length #
## finding 1-to-1 associations in length ##
	my ($unmerged_r, $merged_r) = @_;
	
	my %scaf_index;
	foreach my $mcontig (keys %$merged_r){				# merged => unmerged
		# comparing merged & unmerged contig lengths #
		foreach my $ucontig (keys %$unmerged_r){
			push( @{$scaf_index{$mcontig}}, $ucontig) if
				$merged_r->{$mcontig}{"length"} == $unmerged_r->{$ucontig}{"length"};
			}
		# checking that contig had at least on hit #
		unless (exists $scaf_index{$mcontig} ){
			die " ERROR: Contig$mcontig did not have a match in the merged genbank!\n";
			}
		}		

	# making list of taxa that need blasting #
	my %needs_blasting;
	foreach my $contig (keys %scaf_index){
		if(scalar @{$scaf_index{$contig}} > 1){			# if to many associations
			$needs_blasting{"merged"}{$contig} = 1;
			foreach ( @{$scaf_index{$contig}} ){
				$needs_blasting{"unmerged"}{$_} = 1;
				}
			delete $scaf_index{$contig};				# deleting from index; non 1-to-1
			}
		elsif( scalar @{$scaf_index{$contig}} < 1){
			die " LOGIC ERROR: $!\n";
			}
		}
	
		#print Dumper %needs_blasting; exit;
		#print Dumper %scaf_index; exit;	
	return \%scaf_index, \%needs_blasting;
	}

sub parse_genbank{
# pulling out start-stop & sequence for all 'source' features #
	my ($gbko, $cat, $gbk_file) = @_;
	
	my %loc;			
	print STDERR "...loading $cat genbank file: '$gbk_file'\n" unless $verbose_b;
	my $source_cnt = 0;
	while(my $seqo = $gbko->next_seq){
		my @sources = grep{ $_->primary_tag eq "source" } $seqo->get_SeqFeatures;
		foreach my $source (@sources){
			# umerged contigs: ID by count; merged: ID by ID #
			my $contig_id;
			if($cat eq 'merged'){
				$contig_id = $source_cnt; 
				}
			else{
				$contig_id = $source->seq->display_id; 
				}
			$loc{$contig_id}{"start"} = $source->location->start;
			$loc{$contig_id}{"end"} = $source->location->end;
			$loc{$contig_id}{"length"} = abs($loc{$contig_id}{"end"} -
															$loc{$contig_id}{"start"});
			$loc{$contig_id}{"seq"} = $source->seq->seq;
			$loc{$contig_id}{"seq"} =~ tr/A-Z/a-z/;
			$loc{$contig_id}{"id"} = $contig_id;
			
			$source_cnt++;
			}
		}

	
		#print Dumper %loc; exit;
	return \%loc;

	}

sub load_genbank_io{
# parsing the merged genbank into scaffold lengths & sequence #
	my ($gbk_list_r) = @_;
	
		#print Dumper $gbk_list_r; exit;
	
	my %gbk;
	foreach my $cat (keys %$gbk_list_r){		
		$gbk{$cat}	= 	Bio::SeqIO->new(-format => "genbank", 
						-file => $gbk_list_r->{$cat});
		}
		#print Dumper %gbk; exit;
	return \%gbk;
	}

sub load_array_list{
# loading list of array files #
	my ($array_in) = @_;
	
	open IN, $array_in or die $!;
	my %array_list;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @l = split /\t/;
		die "ERROR: array list file must be in 2-column format: ('taxon_name\\tarray_file')\n"
			unless scalar @l == 2;
		
		$array_list{$l[0]} = $l[1];
		}
	close IN;
	
		#print Dumper %array_list;
	return \%array_list;
	}

sub load_loci_table{
# loci loci table file used for CLdb #	
	my %loci;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @l = split /\t/;
		
		if($.==1){	# header
			map{ tr/A-Z/a-z/ } @l; 		# caps-invariant
			for my $i (0..$#l){
				$loci{'header'}{$l[$i]} = $i;
				}
			die "ERROR: cannot find taxon_name column in header of loci table!\n"
				unless exists $loci{'header'}{'taxon_name'};
			}
		else{		# body 
			die "ERROR: no taxon_name value for: $_\n"
				unless defined $l[ $loci{'header'}{'taxon_name'} ];
			my $taxon_name = $l[ $loci{'header'}{'taxon_name'} ];
			push @{$loci{'body'}{$taxon_name}}, \@l;
			}
		}

	# checking for genbank file #
	die "ERROR: no 'genbank_file' column found in loci table! Possitions in array files will not be editted!\n"
		unless exists $loci{'header'}{'genbank_file'};

	# checking for array file #
	warn "WARNING: no 'array_file' column found in loci table! Possitions in array files will not be editted!\n"
		unless exists $loci{'header'}{'array_file'};
		
	# adding scaffold column if not present #
	$loci{'header'}{'scaffold'} = scalar keys %{$loci{'header'}}
		unless exists $loci{'header'}{'scaffold'};

		#print Dumper %loci; exit;
	return \%loci;
	}

sub load_genbank_list{
# getting list of genbanks #
## 3 column: taxon_name\tunmerged\tmerged ##
	my ($gbk_list_in) = @_;
	
	open IN, $gbk_list_in or die $!;
	
	my %list;
	while(<IN>){
		chomp;
		next if /^\s*$/;
		
		my @tmp = split /\t/;
		die "ERROR: the genbank list must have 3 columns (taxon_name\\tunmerged\\tmerged)!\n"
			unless scalar @tmp == 3;
		map{ die "ERROR: cannot find file '$_'\n" unless -e $_ } @tmp[1..2];
		
		$list{$tmp[0]}{'unmerged'} = $tmp[1];
		$list{$tmp[0]}{'merged'} = $tmp[2];
		}
	close IN;
	
		#print Dumper %list; exit;
	return \%list;
	}


