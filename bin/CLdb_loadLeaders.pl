#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_loadLeaders.pl -- adding/updating leader entries in to CRISPR_db

=head1 SYNOPSIS

CLdb_loadLeaders.pl [flags] leader_aligned.fasta

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -keep <int>-<int>

Positions in the aligment to keep as conserved leader region (indexed by 1). 
The alignment start and end will be used if the range start
or end is not given (eg. '50-' means keep position 50 to the end). [-]

=item -gap  <bool>

Leave gaps in leader sequence entered into the DB? [TRUE] 

=item -verbose  <bool>

Verbose output. [FALSE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_loadLeaders.pl

=head1 DESCRIPTION

After determining the actual leader region by aligning
approximate leader regions written by CLdb_getLeaderRegions.pl,
automatically trim and load leader sequences & update
loci start-end to include the leader region.

=head2 -keep

The alignment should have the orientation whether the array is downstream
(to the right in the alignment). Therefore, trimming of the non-conserved 
regions of the alignment are likely going to be on the left (upstream) side.
For instance, '-keep 100-' will remove the 1st 99 positions from the alignemnt.

=head1 EXAMPLES

=head2 Triming 49 positions off of the upstream unconserved end of the alignment

CLdb_loadLeaders.pl -d CLdb.sqlite -t 50 test_leader_Ib.fna test_leader_Ib_aln.fna

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
use List::Util qw/min max/;

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
use CLdb::seq qw/revcomp/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $rm_gap);
my $keep = '-';
GetOptions(
	   "database=s" => \$database_file,
	   "keep=s" => \$keep,
	   "gap" => \$rm_gap,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
die " ERROR: provide the aligned leader fasta file!\n"
	unless $ARGV[0];
## keep
my $keep_r = check_keep($keep);


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);

# load alignment #
my $fasta_aln_r = load_fasta($ARGV[0]);

# determining trim start-end #
trim_se($fasta_aln_r, $keep_r, $rm_gap);

# updating / loading_db #
load_leader($dbh, $fasta_aln_r);

# updating loci table #
my $loci_tbl_r = get_loci_table($dbh);

# reset locus start-end #
reset_locus_start_end($loci_tbl_r, $fasta_aln_r);

# adjust locus_start-end to include leaders if needed #
adjust_locus_start_end($loci_tbl_r, $fasta_aln_r);

# update locus table #
update_loci_table($dbh, $loci_tbl_r, $fasta_aln_r);

# disconnect to db #
$dbh->commit;	
$dbh->disconnect();
exit;


#--- Subroutines ---#
sub check_keep{
  my $keep = shift;
  die "ERROR: '-keep' should be a range\n"
    unless $keep =~ /^[^-]*-[^-]*$/;

  my @keep = split /-/, $keep;
  map{ $_ = 1 unless defined $_ and $_ ne ''} @keep[0..1];
  map{ $_ = $_ -1 } @keep;

  if($keep[1] > 0){
    die "ERROR: -keep range must be positive start>end\n"
      unless $keep[0] <= $keep[1];
  }
  else{ $keep[1] = 'end'; }

  return \@keep;
}


sub update_loci_table{
# updating loci table with possibly new locus start-end
  my ($dbh, $loci_tbl_r, $fasta_aln_r) = @_;
  
  # preparing sql #	
  my $cmd = "UPDATE Loci SET 
	locus_start=?, locus_end=? where locus_id = ?";
  $cmd =~ s/[\r\n]//g;
  my $sql = $dbh->prepare($cmd);
  
  # updating #
  my $cnt = 0;
  foreach my $locus_id (keys %$fasta_aln_r){
    
    $sql->bind_param( 1, ${$loci_tbl_r->{$locus_id}}[0] );		# updating locus start
    $sql->bind_param( 2, ${$loci_tbl_r->{$locus_id}}[1] );		# updating locus end
    
    $sql->bind_param(3, $locus_id);
    $sql->execute( );
    
    if($DBI::err){
      print STDERR "ERROR: $DBI::errstr for updating $locus_id \n";
    }
    else{ $cnt++; }
  }
  
  print STDERR "...Number of entries updated in loci table: $cnt\n"
    unless $verbose;
}

sub adjust_locus_start_end{
# expanding locus-start-end if leader extends out beyond locus start-end
  my ($loci_tbl_r, $fasta_aln_r) = @_;
  
  foreach my $locus_id (keys %$fasta_aln_r){	
    die " ERROR: $locus_id does not exist in locus table!\n"	
      unless exists $loci_tbl_r->{$locus_id};
    
    # flipping start-end of locus & leader if needed #
    my @flip_bool = (0,0);
    if( $fasta_aln_r->{$locus_id}{"meta"}{"start"} >   # leader 			
	$fasta_aln_r->{$locus_id}{"meta"}{"end"}){
      ($fasta_aln_r->{$locus_id}{"meta"}{"start"}, 			
       $fasta_aln_r->{$locus_id}{"meta"}{"end"}) = 
	 flip_se($fasta_aln_r->{$locus_id}{"meta"}{"start"},
		 $fasta_aln_r->{$locus_id}{"meta"}{"end"});
      $flip_bool[0] = 1;
			}
    
    if ( ${$loci_tbl_r->{$locus_id}}[0] >					# locus
	 ${$loci_tbl_r->{$locus_id}}[1]){
      
      @{$loci_tbl_r->{$locus_id}}[0..1] =
			 	flip_se(@{$loci_tbl_r->{$locus_id}}[0..1]);
      $flip_bool[1] = 1;
    }
    
    # checking whether leader falls outside of locus #
    ## if yes, updating locus table ##
    if($fasta_aln_r->{$locus_id}{"meta"}{"start"} < 			# leader_start < locus_start
       ${$loci_tbl_r->{$locus_id}}[0] ){
      my $se = join(" -> ", ${$loci_tbl_r->{$locus_id}}[0], $fasta_aln_r->{$locus_id}{"meta"}{"start"});
      
      ${$loci_tbl_r->{$locus_id}}[0] = $fasta_aln_r->{$locus_id}{"meta"}{"start"};
		   	
      
      print STDERR "...Locus_start for $locus_id extended to include leader region ($se)\n";
    }
    elsif( $fasta_aln_r->{$locus_id}{"meta"}{"end"} > 			# leader_end > locus_end
	   ${$loci_tbl_r->{$locus_id}}[1]){
      my $se = join(" -> ", ${$loci_tbl_r->{$locus_id}}[1], $fasta_aln_r->{$locus_id}{"meta"}{"end"});
      
      ${$loci_tbl_r->{$locus_id}}[1] = $fasta_aln_r->{$locus_id}{"meta"}{"end"};
      print STDERR "...Locus_end for $locus_id extended to include leader region ($se)\n";
    }
    
    # flipping back start-end if necessary #
    if($flip_bool[0]){
      ($fasta_aln_r->{$locus_id}{"meta"}{"start"},
       $fasta_aln_r->{$locus_id}{"meta"}{"end"} ) =
	 flip_se($fasta_aln_r->{$locus_id}{"meta"}{"start"},
		 $fasta_aln_r->{$locus_id}{"meta"}{"end"});
    }
    if($flip_bool[1]){
      @{$loci_tbl_r->{$locus_id}}[0..1] =
	flip_se(@{$loci_tbl_r->{$locus_id}}[0..1]);	
    }
  }	
  
  #print Dumper $loci_tbl_r; exit;	
}

sub reset_locus_start_end{
# reseting locus start-end to just array + CAS #
## should handle cases where CAS or array start-end encompasses the other
  my ($loci_tbl_r, $fasta_aln_r) = @_;
	
  foreach my $locus_id (keys %$fasta_aln_r){			# each locus of interest
    # checking for operon or array start-end #
    die " ERROR: no array_start or cas_start found in loci table for cli.$locus_id!"
      unless ${$loci_tbl_r->{$locus_id}}[2] || ${$loci_tbl_r->{$locus_id}}[4];
    
    # if just array or operon (CAS), give locus those values #
    if(! ${$loci_tbl_r->{$locus_id}}[2] ){			# if no operon, using array
      my $min = min(@{$loci_tbl_r->{$locus_id}}[4..5]);
			my $max = max(@{$loci_tbl_r->{$locus_id}}[4..5]);
      if ( ${$loci_tbl_r->{$locus_id}}[0] <			# locus start < end
	   ${$loci_tbl_r->{$locus_id}}[1]){
	${$loci_tbl_r->{$locus_id}}[0] = $min;
	${$loci_tbl_r->{$locus_id}}[1] = $max;
      }
      else{											# locus start > end
	${$loci_tbl_r->{$locus_id}}[0] = $max;
	${$loci_tbl_r->{$locus_id}}[1] = $min;
      }
    }
    elsif(! ${$loci_tbl_r->{$locus_id}}[4] ){			# if no array start'; using operon
      my $min = min(@{$loci_tbl_r->{$locus_id}}[2..3]);
      my $max = max(@{$loci_tbl_r->{$locus_id}}[2..3]);
      if ( ${$loci_tbl_r->{$locus_id}}[0] <			# locus start < end
	   ${$loci_tbl_r->{$locus_id}}[1]){
	${$loci_tbl_r->{$locus_id}}[0] = $min;
	${$loci_tbl_r->{$locus_id}}[1] = $max;
      }
      else{											# locus start > end
	${$loci_tbl_r->{$locus_id}}[0] = $max;
	${$loci_tbl_r->{$locus_id}}[1] = $min;
      }
    }
    else{		 # both array and CAS operon
      my $min = min(@{$loci_tbl_r->{$locus_id}}[2..5]);
      my $max = max(@{$loci_tbl_r->{$locus_id}}[2..5]);
      
      # adjusting using operon & array #
      if ( ${$loci_tbl_r->{$locus_id}}[0] <			# locus start < end
	   ${$loci_tbl_r->{$locus_id}}[1]){
	${$loci_tbl_r->{$locus_id}}[0] = $min;
	${$loci_tbl_r->{$locus_id}}[1] = $max;				
      }
      else{											# locus start > end
	${$loci_tbl_r->{$locus_id}}[0] = $max;
	${$loci_tbl_r->{$locus_id}}[1] = $min;
      }
    }	
  }
  #print Dumper "here", $loci_tbl_r; exit;
}

sub flip_se{
  my ($start, $end) = @_;
  return $end, $start;
}

sub get_loci_table{
# getting loci table for updating #
  my ($dbh) = @_;
  
  my $cmd = "SELECT locus_id, locus_start, locus_end, CAS_start, CAS_end, array_start, array_end from loci";
  my $ret = $dbh->selectall_arrayref($cmd);
  die "ERROR: not entries in loci table!\n" unless defined $ret;
  
  my %loci_tbl;
  foreach my $row (@$ret){
    $loci_tbl{$$row[0]} = [@$row[1..$#$row]];
  }
  
  return \%loci_tbl;
}

sub load_leader{
# loading leader info into db #
  my ($dbh, $fasta_aln_r) = @_;

  # cmd 
  my $cmd = "INSERT INTO Leaders(Locus_ID, Leader_start, Leader_end, Leader_sequence) values (?,?,?,?)";
  my $sql = $dbh->prepare($cmd);

  # update db
  my $cnt = 0;
  foreach my $locus (keys %$fasta_aln_r){  
    $sql->execute( $locus, 									# locus_ID
		   $fasta_aln_r->{$locus}{"meta"}{"start"}, 	# start
		   $fasta_aln_r->{$locus}{"meta"}{"end"}, 		# end
		   $fasta_aln_r->{$locus}{"seq"});					# sequence
    if($DBI::err){
      print STDERR "ERROR: $DBI::errstr for $locus\n";
    }
    else{ $cnt++; }
  }
  $dbh->commit;
  
  print STDERR "...Number entries added/updated in leaders table: $cnt\n";
}

sub trim_se{
# trimming alignment to just required sequences #
## modifying start-end locations accordingly ##
  my ($fasta_aln_r, $keep_r, $rm_gap) = @_;

  foreach my $locus (keys %$fasta_aln_r){
    #print Dumper $fasta_aln_r->{$locus}; exit;
    # trimming alignment sequence (strand doesn't matter for sequence)
    my $seq_len = length $fasta_aln_r->{$locus}{seq};
    $keep_r->[1] = $seq_len if $keep_r->[1] eq 'end';
    ## splitting into 3 parts: upstream, keep, downstream
    my $up_seq = substr($fasta_aln_r->{$locus}{seq}, 
			0, $keep_r->[0]); # no inclusive of -keep start  
    my $keep_seq = substr($fasta_aln_r->{$locus}{seq}, $keep_r->[0], 
			  $keep_r->[1] - $keep_r->[0] + 1); # inclusive of -keep start & end  
    my $down_seq = substr($fasta_aln_r->{$locus}{seq}, 
			  $keep_r->[1], $seq_len - $keep_r->[1]); # not inclusive of -keep end  
    $fasta_aln_r->{$locus}{seq} = $keep_seq; # portion to keep

    ## getting number of gaps
    #print Dumper $up_seq, $keep_seq, $down_seq; exit;
    #my %gap_cnt = (up => 0, keep => 0, down => 0);
    my %gap_cnt;
    $gap_cnt{up} = count_gaps($up_seq);
    $gap_cnt{keep} = count_gaps($keep_seq);
    $gap_cnt{down} = count_gaps($down_seq);

    # editing sequence region start-end
    ## start-end set for + strand
    if($fasta_aln_r->{$locus}{meta}{strand} eq '+'){
      ## start
      $fasta_aln_r->{$locus}{meta}{start} = 
	$fasta_aln_r->{$locus}{meta}{start} + length($up_seq) - $gap_cnt{up};
      ## end
      $fasta_aln_r->{$locus}{meta}{end} = 
	$fasta_aln_r->{$locus}{meta}{end} - length($down_seq) - $gap_cnt{down};
    }
    elsif($fasta_aln_r->{$locus}{meta}{strand} eq '-'){   # - strand, start-end on + strand flipped
      ## start (relates to downstream on + strand)
      $fasta_aln_r->{$locus}{meta}{start} = 
	$fasta_aln_r->{$locus}{meta}{start} + length($down_seq) - $gap_cnt{down};
      ## end (relates to upstream on + strand)
      $fasta_aln_r->{$locus}{meta}{end} = 
	$fasta_aln_r->{$locus}{meta}{end} - length($up_seq) - $gap_cnt{up};
    }
    else{
      die "ERROR: cannot identify strand for locus: $locus\n";
    }    
  }
}


sub trim_se_OLD{
# modifying raw fasta to relect alignment #
## trimming sequence and modifying start-end locations accordingly ##
## trimming from leader end most distant from array ##
  my ($fasta_aln_r, $keep_r, $rm_gap) = @_;

  # debug
  my ($ori_r, $trim_length); 
  #print Dumper $keep_r; 

#print Dumper $fasta_aln_r; exit;
	
  foreach my $locus (keys %$fasta_aln_r){
    die " LOGIC ERROR: $locus not found in ori hash!\n"
      unless exists $ori_r->{$locus};
    
    # getting seq info #
    my $seq_len = length $fasta_aln_r->{$locus}{seq};
    
    # trimming alignment sequences (removing non-conserved alignment positions) #
    if($ori_r->{$locus} eq "norm"){
      #print Dumper 'here';
      if($fasta_aln_r->{$locus}{meta}{strand} eq '+'){		# trimming from 5' end (+ strand)
	# trim seq #
	my $trimmed_seq = substr( $fasta_aln_r->{$locus}{seq}, 0, $trim_length);
	$fasta_aln_r->{$locus}{seq} = substr( $fasta_aln_r->{$locus}{seq}, $trim_length, $seq_len);
	
	# changing start position (moving up toward 3') #
	my $ngap = count_gaps($trimmed_seq);			# counting gaps of what remains
	$fasta_aln_r->{$locus}{meta}{start} += ($trim_length - $ngap);
      }
      elsif( $fasta_aln_r->{$locus}{meta}{strand} eq '-' ){		# trimming from 3' end (+ strand)
	# trim seq #
	my $trimmed_seq = substr($fasta_aln_r->{$locus}{seq}, $seq_len - $trim_length );
	$fasta_aln_r->{$locus}{seq} = substr($fasta_aln_r->{$locus}{seq}, 0, $seq_len - $trim_length);
	
	# changing end position (moving back toward 5') #
	my $ngap = count_gaps($trimmed_seq);			# counting gaps of what remains
	$fasta_aln_r->{$locus}{meta}{end} -= ($trim_length - $ngap);
      }
      else{ die $!; }
    }
    elsif($ori_r->{$locus} eq 'rev' || $ori_r->{$locus} eq 'rev-comp'){
      
      if($fasta_aln_r->{$locus}{"meta"}{"loc"} eq "start"){		# trimming from 3' end (+ strand)
	# trim seq #
	my $trimmed_seq = substr($fasta_aln_r->{$locus}{"seq"}, $seq_len - $trim_length );
	$fasta_aln_r->{$locus}{"seq"} = substr($fasta_aln_r->{$locus}{"seq"}, 0, $seq_len - $trim_length);
	
	# changing end position (moving back toward 5') #
	my $ngap = count_gaps($trimmed_seq);			# counting gaps of what remains				
	$fasta_aln_r->{$locus}{"meta"}{"end"} -= ($trim_length - $ngap);			
      }
      elsif( $fasta_aln_r->{$locus}{"meta"}{"loc"} eq "end" ){		# trimming from 5' end (+ strand)
	# trim seq #
	my $trimmed_seq = substr( $fasta_aln_r->{$locus}{"seq"}, 0, $trim_length);
	$fasta_aln_r->{$locus}{"seq"} = substr( $fasta_aln_r->{$locus}{"seq"}, $trim_length , $seq_len);
	
	# changing start position (moving up toward 3') #
	my $ngap = count_gaps($trimmed_seq);			# counting gaps of what remains				
	$fasta_aln_r->{$locus}{"meta"}{"start"} += ($trim_length - $ngap);
      }		
      else{ die $!; }
    }
    else{ die $!; }
    
    # removing gaps if specified #
    $fasta_aln_r->{$locus}{"seq"} =~ s/-//g if $rm_gap;
  }
  
  print Dumper %$fasta_aln_r; exit;
}

sub count_gaps{
# counting gaps that fall in trimmed alignment range #
  my ($seq) = @_;
  
  my $ngap = 0;
  $ngap++ while $seq =~ /-/g;
  return $ngap;
}

sub get_array_se{
# getting the array start-end from loci table #
  my ($dbh, $ori_r) = @_;
  
  my $cmd = "SELECT array_start, array_end FROM loci where locus_id = ?";
  my $sql = $dbh->prepare($cmd);
  
  my %array_se;
  foreach my $locus (keys %$ori_r){
    (my $tmp = $locus) =~ s/cli\.//;
    $sql->execute($tmp);
    my $ret = $sql->fetchall_arrayref();
    die " ERROR: no array start-end info in Loci table for $locus!\n"
      unless $$ret[0];
    
    $array_se{$locus} = $ret;
  }
  
  #print Dumper %array_se; exit;
  return \%array_se; 
}	

sub load_fasta{
# loading fasta #
  my $fasta_in = shift;
  open IN, $fasta_in or die $!;	
  my (%fasta, $tmpkey);
  while(<IN>){
    chomp;
    s/#.+//;
    next if  /^\s*$/;	
    
    if(/>.+/){
      
      s/>_R_/>/;				# if reversed alignment (mafft add '_R_' if reversed)
      my @line = split /[>|]/;
      die "ERROR: leader sequence: '$_' does not have CLdb_getLeaderRegions.pl formatting!\n"
	unless scalar @line == 6;      

      # die if mafft revcomp (leaders should align with each other without rev-comp
      die "ERROR: leader regions should alignment without the need to reverse-complement\n"
	if $line[1] =~ /^_R_/;
 
      # loading metadata
      my @meta = qw/NA locus_id scaffold start end strand/;
      foreach my $i (1..$#line){
	$fasta{$line[1]}{"meta"}{$meta[$i]} = $line[$i];
      }
      
      $tmpkey = $line[1];	# changing key
    }
    else{$fasta{$tmpkey}{"seq"} .= $_; }
  }
  close IN;
  
  #print Dumper %fasta; exit;
  return \%fasta;
} #end load_fasta

sub update_db{
# updating any entries with CLI identifiers #
  my ($dbh, $loci_r, $header_r) = @_;
  
  # entries ordered #
  my @keys = keys %$header_r;
  @keys = grep(!/^locus_id$/i, @keys);
  my @values = @$header_r{@keys};
  
  # setting up update for all columns #
  my @set;
  foreach my $key (@keys){
    (my $tmp = $key) =~ s/$/ = ?/;
    push @set, $tmp;
  }

  # preparing sql #	
  my $cmd = join(" ", "UPDATE Loci SET", join(",", @set), "where locus_id = ?");
  my $sql = $dbh->prepare($cmd);
  
  # updating #
  foreach my $entry (keys %$loci_r){
    next if $entry eq "new_entry";
    my $row = $loci_r->{$entry};
    
    $sql->execute( (@$row[@values], $$row[$header_r->{"locus_id"}]) );
    
    if($DBI::err){
      print STDERR "ERROR: $DBI::errstr in: ", join("\t", @$row[@values]), "\n";
    }
  }
  
  $dbh->commit;	
  
  print STDERR "...Number of entries updated in loci table: ", (scalar keys %$loci_r) -1, "\n"
    unless $verbose;
}


sub check_for_loci_table{
  my ($table_list_r) = @_;
  die " ERROR: loci table not found in database!\n"
    unless grep(/^loci$/i, @$table_list_r);
}

sub list_tables{
  my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
  return [keys %$all];
}

