#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_arrayBlastAddProto.pl -- getting the protospacers for spacer blast hits (& adjacent regions)

=head1 SYNOPSIS

CLdb_arrayBlastAddProto.pl [flags] < spacerBlast.txt > spacerBlast_proto.txt

=head2 flags

=over

=item -x  <int>

Extension beyond spacer blast to check for PAMs (bp). [10]

=item -length  <float>

Length cutoff for blast hit (>=; fraction of spacer length). [0.66]

=item -pam  <int>

Columns containing just he supposed PAM region. This is designated with 4 values to determine 
5' & 3' PAM region adjacent to the protospacer. Example: -pam -3 -1 1 3 will get the 3 bp 
adjacent up and downstream of the protospacer. [-3 -1 1 3] 

=item -revcomp  <bool>

Protospacer as opposite strand of spacer blast hit? [FALSE] (default: protospacer oriented to subject + strand)

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>	

This help message

=back

=head2 For more information:

perldoc CLdb_arrayBlastAddProto.pl

=head1 DESCRIPTION

Get sequences of protospacers (found via 
blast hits) plus adjacent nucleotides
to look for PAMs.

The input is a blast table with comment lines ('-outfmt 7').
The output is the same blast table with appended fields & 
spacers not meeting the cutoffs removed.

Adjacent regions to the protospacer ('proto_seq_fullx' field)
are in lower case, while the protospacer is upper case.

By default, the protospacer is oriented to subject + strand.

=head2 Extra fields appended to the spacer blast table:

=over

=item * "proto_seq_rel2query" = protospacer sequence (just blast match); oriented to query + strand

=item * "proto_seq_rel2subject" = protospacer sequence (just blast match); oriented to subject + strand 

=item * "q_start_full" = query start (full length spacer)

=item * "q_end_full" = query end (full length spacer)

=item * "s_start_full" = protospacer start (full length spacer hit)

=item * "s_end_full" = protospacer end (full length spacer hit)

=item * "proto_seq_full" = protospacer sequence (full length spacer hit) 

=item * "s_start_fullx" = protospacer sequence start (full length spacer hit) & extension 

=item * "s_end_fullx" = protospacer sequence end (full length spacer hit) & extension 

=item * "proto_seq_fullx" = protospacer sequence (full length spacer hit) & extension 

=item * "x_length" = length of extension beyond full protospacer sequence (each side), ('-x' flag)

=item * "5p_pam_region" = 5' pam sequence (by subject + strand). ('-pam' flag)

=item * "3p_pam_region" = 3' pam sequence (by subject + strand). ('-pam' flag)

=back 

If using '-revcomp', the protospacer sequences will be 
the reverse-complement of spacer blast hits.
Their start-end values should reflect this.

If a spacer blast hits the end of a scaffold/chromosome, 
extending the spacer hit to the full spacer length ("proto_seq_full")
or protospacer+extension length ("proto_seq_fullx") will be 
limited to the end of the genomic sequence.

=head2 WARNING

The spacer blast table must have comment lines!

The blast databases must be in the location specified
in the "# Database:" comment lines! 

=head1 EXAMPLES

=head2 Protospacers for all hits (by spacer group)

CLdb_arrayBlastAddProto.pl < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

=head2 Protospacers for just full length spacer blast hits

CLdb_arrayBlastAddProto.pl -l 1 < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

=head2 Protospacers with 5 bp extensions 

CLdb_arrayBlastAddProto.pl -x 5 < spacerBlast_DR-filtered.txt > spacerBlast_proto.txt

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
use DBI;
use Bio::SeqIO;
use List::Util qw/max/;

# CLdb libs #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::blast qw/
	parse_blast_hits/;
use CLdb::seq qw/
	revcomp/;



### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $full_len, $revcomp_b, @pam);
my $extend = 10;							# length to extend off of each side
my $len_cut = 0.66;							# cutoff for length of blast hit relative to query length
GetOptions(
	   "x=i" => \$extend,
	   "length=f" => \$len_cut,				# min length of spacer hit (fraction of total)
	   "revcomp" => \$revcomp_b,			# reverse complement protospacer relative to blast hit? [FALSE]
	   "pam=i{4,4}" => \@pam, 				# pam region to write out
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
if(@pam){
	map{die " ERROR: the pam start-end should be negative values!\n" unless $_ < 0} @pam[0..1];
	map{die " ERROR: the pam start-end should be positive values!\n" unless $_ > 0} @pam[2..3];
	die " ERROR: the 5' pam start should be <= pam end (default: -3 -1)\n" unless $pam[0] <= $pam[1];
	die " ERROR: the 3' pam start should be <= pam end (default: 1 3)\n" unless $pam[2] <= $pam[3];
	}
else{ @pam = (-3,-1,1,3); }

### main
# parsing blast table #
my ($blast_r, $fields_r) = parse_blast_hits();

# getting protospacer & extensions #
## making list of protospacer/extension start-stop, & strand ##
my $new_fields_r = get_proto_seq($blast_r, $fields_r, $extend, $len_cut);

# writing editted table #
## adding new fields ##
add_new_fields($fields_r, $new_fields_r);
## writing table ##
write_editted_blast($blast_r, $fields_r);


#--- Subroutines ---#
sub write_editted_blast{
# writing editted blast table #
## making sure to update fields with appended ones ##
	my ($blast_r, $fields_r) = @_;
	
	foreach my $query (sort keys %$blast_r){
		foreach my $db (sort keys %{$blast_r->{$query}}){
			# getting nubmer of hits passing filtering  # # skipping hit if no adding info (did not pass filtering) #
			my $hit_cnt = 0;
			foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
				$hit_cnt++ if $$row[$fields_r->{"proto_seq_fullx"}];
				}
			next if $hit_cnt == 0;

			# modifying & writing header #
			my $header_r;
			foreach my $row (@{$blast_r->{$query}{$db}{"header"}}){
				$row = new_header($row, $fields_r) if $row =~ /^# Field/;
				print $row, "\n";
				}			

			# writing hits #
			foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
				next unless $$row[$fields_r->{"proto_seq_fullx"}];
				print join("\t", @$row), "\n";
				}			
			}
		}	
	}
	
sub new_header{
	my ($row, $fields_r) = @_;
	my $fields = join(", ", sort{$fields_r->{$a}<=>$fields_r->{$b}} keys %$fields_r);
	return join(" ", "# Fields:", $fields);
	}
	
sub add_new_fields{
# adding new fields for output #
	my ($fields_r, $new_fields_r) = @_;
	
	my $max_i = max(values %$fields_r);	# starting after last fields index
	
	for my $i (0..$#$new_fields_r){
		$fields_r->{$$new_fields_r[$i]} = $i+$max_i+1;
		}
		#print Dumper %$fields_r; exit;
	}

sub get_proto_seq{
# getting protospacer/extension start-stop & strand #
	($blast_r, $fields_r, $extend, $len_cut) = @_;
	
	my %filter_sum;
	foreach my $query (sort keys %$blast_r){
		print STDERR "Processing query: '$query'\n" unless $verbose;
		foreach my $db (keys %{$blast_r->{$query}}){
			foreach my $row (@{$blast_r->{$query}{$db}{"hits"}}){
				# checking existence of start-stop, qlen, slen #
				check_values($row, $fields_r);	
				
				# filtering #
				$filter_sum{"total"}++;
				## deleting  hits with gaps ##
				if ($$row[$fields_r->{"BTOP"}] =~ /-/){ 	# gaps in blast alignment
					$filter_sum{"gaps"}++;
					next;
					}
				## deleling if hit_len/query_len < $len_cut ##
				if( ($$row[$fields_r->{"q. end"}] - $$row[$fields_r->{"q. start"}]) / 
							$$row[$fields_r->{"query length"}] < $len_cut ){
					$filter_sum{"length"}++;
					next;
					}
				
				# getting subject start-stop, strand #
				my $Sstart = $$row[$fields_r->{"s. start"}];
				my $Send = $$row[$fields_r->{"s. end"}];
				my $strand = "plus";

				# flipping if needed #
				if($Sstart > $Send){
					$strand = "minus";
					($Sstart,$Send) = ($Send, $Sstart);		# start-stop to + strand
					}
				
				#$strand = "plus" unless $revcomp_b; 		# by default; protospacer oriented to subject + strand
					
				# calling blastdbcmd #
				my $proto_seq = call_blastdbcmd($$row[$fields_r->{"subject id"}],
								$Sstart, $Send, $strand, $db, 1);			 # inverted so it is protospacer match
				
				# extending hit to full length & getting sequence #
				## sstart-send, qstart-qend full length (as much as extension would allow) ##
				my ($Sstart_full, $Send_full, $Qstart_full, $Qend_full) = 
						extend_full_len($Sstart, $Send, $strand,
									$$row[$fields_r->{"q. start"}],
									$$row[$fields_r->{"q. end"}],
									$$row[$fields_r->{"query length"}],
									$$row[$fields_r->{"subject length"}]);

				my $proto_seq_full = call_blastdbcmd($$row[$fields_r->{"subject id"}],
								$Sstart_full, $Send_full, $strand, $db, 1);		# not inverted; used for alignment
				
				# extending hit beyond full length (-x) #	
				my ($Sstart_fullx, $Send_fullx) = extend_beyond_full($Sstart_full, $Send_full, $strand,
													$$row[$fields_r->{"subject length"}], $extend);
				
				# calling blastdbcmd #
				my $proto_seq_full_x = call_blastdbcmd($$row[$fields_r->{"subject id"}],
								$Sstart_fullx, $Send_fullx, $strand, $db, 1);				
				
				# making extensions lower case #
				$proto_seq_full_x = ext_lower_case($proto_seq_full_x,
										$Sstart_fullx, $Send_fullx,
										$Sstart_full, $Send_full);					
				
				
				# extracting the pam regions #
				my ($fiveP_pam, $threeP_pam) = pam_extract($proto_seq_full_x, \@pam, 
										$Sstart_fullx, $Send_fullx,
										$Sstart_full, $Send_full);
				
				# adding sequences to blast table #
				## proto sequence extension start-end flip if hits to minus strand #
				($Sstart_full, $Send_full) = ($Send_full, $Sstart_full) if $strand eq "minus";		
				($Sstart_fullx, $Send_fullx) = ($Send_fullx, $Sstart_fullx) if $strand eq "minus";	
				
				
				# bring subject (protospacer) to + strand unless revcomp_b #
				my $proto_seq_orig = $proto_seq;				# proto sequence relative to query (not potentially revcomp)
				if($strand eq "minus" && ! $revcomp_b){
					$proto_seq = revcomp($proto_seq);
					($Sstart_full, $Send_full) = ($Send_full, $Sstart_full);
					$proto_seq_full = revcomp($proto_seq_full);
					($Sstart_fullx, $Send_fullx) = ($Send_fullx, $Sstart_fullx);
					$proto_seq_full_x = revcomp($proto_seq_full_x);
					$fiveP_pam = revcomp($fiveP_pam);
					$threeP_pam = revcomp($threeP_pam);
					($fiveP_pam,$threeP_pam) = ($threeP_pam, $fiveP_pam);
					}
				
				# adding to table #
				push @$row, $proto_seq_orig, $proto_seq, $Qstart_full, $Qend_full,
						$Sstart_full, $Send_full, $proto_seq_full,
						$Sstart_fullx, $Send_fullx, $proto_seq_full_x,
						$extend, $fiveP_pam, $threeP_pam; 
				}
			}
		}
	
		# subject_seq = aln to spacer; subject_full_seq = aln to spacer; subject_full_x_seq = revcomp (for PAMs)
		#print Dumper %$blast_r; exit;
	my @new_fields = ("proto_seq_rel2query", "proto_seq_rel2subject", 
						"q_start_full", "q_end_full", "s_start_full", "s_end_full", "proto_seq_full", 
						"s_start_fullx", "s_end_fullx", "proto_seq_fullx", "x_length", "5p_pam_region", "3p_pam_region");
	
	# status #
	map{$filter_sum{$_} = 0 unless exists $filter_sum{$_}} qw/total gaps length/;
	print STDERR "\n### filtering summary ###\n";
	print STDERR "Total number of hits:\t\t\t\t$filter_sum{'total'}\n";
	print STDERR "Number removed due to gaps in alignment: \t$filter_sum{'gaps'}\n";
	print STDERR "Number removed due to short hit length:\t\t$filter_sum{'length'}\n";
	print STDERR "Number remaining:\t\t\t\t", $filter_sum{'total'} - 
						($filter_sum{'length'} + $filter_sum{'gaps'}), "\n\n";
	
	return \@new_fields;
	}
	
sub pam_extract{
	my ($proto_seq_full_x, $pam_r, $Sstart_fullx, $Send_fullx, 
		$Sstart_full, $Send_full) = @_;
	
	# getting extension lengths #
	my $start_x_len = $Sstart_full - $Sstart_fullx;
	my $end_x_len = $Send_fullx - $Send_full;			
	my $proto_len = $Send_full - $Sstart_full + 1; 		# inclusive 
	
	## 5' end ##
	my $FiveP_pam_len = abs($$pam_r[0]) - abs($$pam_r[1]) + 1;
	die " LOGIC ERROR: 5' pam length is < 1!\n" unless $FiveP_pam_len > 0;

	my $FiveP_start = $start_x_len + $$pam_r[0];
	my $FiveP_pam;
	if($FiveP_start < 0){		# not enough extension, returning ""
		$FiveP_pam = "";
		}
	else{
		$FiveP_pam = substr $proto_seq_full_x, $FiveP_start, $FiveP_pam_len;
		}	


	## 3' end #
	my $ThreeP_pam_len = abs($$pam_r[3]) - abs($$pam_r[2]) + 1;
	die " LOGIC ERROR: 3' pam length is < 1!\n" unless $ThreeP_pam_len > 0;

	my $ThreeP_start = $proto_len + $start_x_len + $$pam_r[2] -1;
	my $ThreeP_pam;
	if($ThreeP_start < 0){		# not enough extension, returning ""
		$ThreeP_pam = "";
		}
	else{
		$ThreeP_pam = substr $proto_seq_full_x, $ThreeP_start, $ThreeP_pam_len;
		}
		#print Dumper $proto_len, $ThreeP_start, $ThreeP_pam_len, $ThreeP_pam;
	
		#print Dumper $FiveP_pam, $ThreeP_pam; exit;
	return $FiveP_pam, $ThreeP_pam; 
	}

sub ext_lower_case{
# making extensions to protospacer lower case #
## substr extensions, and full length; lower case extensions; stitching together ##
	my ($proto_seq_full_x, $Sstart_fullx, $Send_fullx, 
		$Sstart_full, $Send_full) = @_;
	
	# getting extension lengths #
	my $start_x_len = $Sstart_full - $Sstart_fullx;
	my $end_x_len = $Send_fullx - $Send_full;			
	my $proto_len = $Send_full - $Sstart_full + 1; 		# inclusive 
	
	# substringing extensions and full #
	my $startx_seq = substr($proto_seq_full_x, 0, $start_x_len);
	my $proto_seq = substr($proto_seq_full_x, $start_x_len, $proto_len); 
	my $endx_seq = substr($proto_seq_full_x, $start_x_len + $proto_len, $end_x_len);

	# lower case extensions #
	$startx_seq =~ tr/[A-Z]/[a-z]/;
	$endx_seq =~ tr/[A-Z]/[a-z]/;
	

		#print Dumper $startx_seq, $proto_seq, $endx_seq; 
	return join("", $startx_seq, $proto_seq, $endx_seq);
	}

sub extend_beyond_full{
# extending beyond full length #
	my ($Sstart, $Send, $strand, $Slen, $extend) = @_;
	
	$Sstart -= $extend;
	$Sstart = 1 if $Sstart < 1;
	$Send += $extend;
	$Send = $Slen if $Send > $Slen;
			
	return $Sstart, $Send;
	}

sub call_blastdbcmd{
	my ($subject_id, $Sstart, $Send, $strand, $db, $invert) = @_;
	
	# flipping strand in order to get rev-comp of spacer hit (if needed) #
	if($invert && $revcomp_b){
		if($strand eq "plus"){ $strand = "minus"; }
		else{ $strand = "plus"; }
		}
	
	# calling #
	my $cmd = "blastdbcmd -db $db -entry '$subject_id' -range '$Sstart-$Send' -strand $strand |";
		#print Dumper $cmd; exit;
	open PIPE, $cmd or die $!;
	my $seq;
	while(<PIPE>){
		chomp;
		die " ERROR: >1 sequence returned by blastdbcmd!\n"
			if /^>/ && $seq;
		next if /^>/;
		$seq .= $_;
		}
	return $seq;
	}

sub extend_full_len{
# extending spacer hit to full length if needed #
	my ($Sstart, $Send, $strand, $Qstart, $Qend, $Qlen, $Slen) = @_;
	
	# missing lengths in spacer blast (if partial hit) #
	my $Q5px_len = $Qstart - 1;			# extension off 5' end of spacer (+ strand); number of bp added 
	my $Q3px_len = $Qlen - $Qend;		# extension off 3' end of spacer (+ strand); number of bp added
	($Q5px_len, $Q3px_len) = ($Q3px_len, $Q5px_len) if $strand eq "minus"; 	# applyin the extension for other strand
	
	# changing subject start-end to full length (as much as allowed) #
	my ($Sstart_full, $Send_full);
	$Sstart_full = $Sstart - $Q5px_len;
	$Sstart_full = 1 if $Sstart_full < 1;			# border issues
	$Send_full = $Send + $Q3px_len; 
	$Send_full = $Slen if $Send_full > $Slen;		# border issues
	
	# getting actual subject extension (acounting for border issues) #
	my $S5px_len = $Sstart - $Sstart_full;		# number of bp added
	my $S3px_len = $Send_full - $Send;			# number of bp added

	# making gaps based on difference between intended extend and actual extend #
	#my $5px_gap = "-" x ($Q5px_len - $S5px_len);
	#my $3px_gap = "-" x ($Q3px_len - $S3px_len);
	#map{$_ = "" unless $_} ($5px_gap, $3px_gap);
	#($5px_gap, $3px_gap) = ($3px_gap, $5px_gap) if $strand eq "minus"; 		# gap applied to the other end
	
	# getting qstart/qend full (as must as subject extension would allow) #	
	($S5px_len, $S3px_len) = ($S3px_len, $S5px_len) if $strand eq "minus"; 	# applying the extension for other strand
	my $Qstart_full = $Qstart - $S5px_len;		
	my $Qend_full = $Qend + $S3px_len;
			
	# sanity check #
	die " LOGIC ERROR $!" unless $Qstart_full > 0 && $Qend_full > 0;
	die " LOGIC ERROR $!" unless $Sstart_full > 0 && $Send_full > 0;
	
	return $Sstart_full, $Send_full, $Qstart_full, $Qend_full; #$5px_gap, $3px_gap;
	}
	
sub check_values{
# checking values to make sure they exist #
	my ($row, $fields_r) = @_;
	
	#print Dumper $row;
	#print Dumper $fields_r; 
	
	my @chk = ("query id", "subject id", "s. start", "s. end",
				"q. start", "q. end", 
				"evalue", "query length", "subject length", "BTOP");
	map{die " ERROR: '$_' not found!\n" unless $$row[$fields_r->{$_}]} @chk;
	}

sub load_fasta{
# loading fasta file as a hash #
	my $fasta_in = shift;
	open IN, $fasta_in or die $!;
	my (%fasta, $tmpkey);
	while(<IN>){
		chomp;
 		s/#.+//;
 		next if  /^\s*$/;	
 		if(/>.+/){
 			$_ =~ s/^>//;
 			$fasta{$_} = "";
 			$tmpkey = $_;	# changing key
 			}
 		else{$fasta{$tmpkey} .= $_; }
		}
	close IN;
		#print Dumper %fasta; exit;
	return \%fasta;
	} #end load_fasta



