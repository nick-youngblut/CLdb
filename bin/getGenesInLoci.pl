#!/usr/bin/env perl

=pod

=head1 NAME

getGenesInLoci.pl -- getting all genes in CRISPR Loci; output is a table of gene information

=head1 SYNOPSIS

getGenesInLoci.pl [flags] 

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -subtype  <char>

Refine query to specific a subtype(s) (>1 argument allowed).

=item -taxon_id  <char>

Refine query to specific a taxon_id(s) (>1 argument allowed).

=item -taxon_name  <char>

Refine query to specific a taxon_name(s) (>1 argument allowed).

=item -exist  <bool>

Do not write any genes already existing in the Genes table (unless '-a' used). [TRUE]

=item -all  <bool>

All genes defined by the Loci table are written. 
Any existing entries written in place of new entry.
Not compatible with '-c' [FALSE]

=item -conflict  <bool>

Just write out genes that are conflicting
(new and old versions). Not compatible with '-a' [FALSE]

=item -quiet  <bool>

Turn off all warnings. [FALSE]

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message.

=back

=head2 For more information:

CLdb_perldoc getGenesInLoci.pl

=head1 DESCRIPTION

Get all CDS features from >= genbanks that fall
into CRISPR loci. 

The output can be piped directly into
loadGenes.pl or the aliases (or other values)
can be edited first.

Genbank files must be in $HOME/genbank/

=head2 Existing genes in Genes table ('-e', '-a', '-c')

By default, only new genes (ie. not found already in
Genes table) are written. This is to
preserve the existing alias values, which may 
have been manually editted. 

All gene entries (new & existing) can be written
using '-a'. To view conflicting gene entries,
use '-c'. For '-c', existing entries will have 'existing_entry'
in the 'Sequence' column.

=head2 In_CAS

Genes are designated as falling in the CAS operon
('In_CAS' field in DB) if they fall inside
the designated operon range (designated in Loci table)
but not inside the CRISPR array range (also 
designated in Loci table), UNLESS the array spans both
sides of the CAS operon (CAS_start-CAS_end).

=head2 WARNING:

CDS features in genbanks must have db_xref tags with 
the fig|peg ID (example: "fig|6666666.40253.peg.2362");
otherwise, the CDS will not be written to the output table.
Extra information can come before 'fig' (eg. 'ITEP:' or 'SEED:'),
but it must end with a colon.

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

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use CLdb::utilities qw/
			file_exists 
			connect2db
			lineBreaks2unix
			get_file_path/;
use CLdb::query qw/
		    table_exists
		    join_query_opts/;
use CLdb::genbank::genbank_get_region qw/genbank_get_region/;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $quiet);
my ($all_genes, $existing, $conflicting);
my (@subtype, @taxon_id, @taxon_name);
my $query = "";
GetOptions(
	   "database=s" => \$database_file,
	   "subtype=s{,}" => \@subtype,
	   "taxon_id=s{,}" => \@taxon_id,
	   "taxon_name=s{,}" => \@taxon_name,
	   "query=s" => \$query,
	   "all" => \$all_genes,   # get all genes (including existing?; overwriting conflicts w/ existing entry). [FALSE]
	   "existing" => \$existing, 	# check for existing, if yes, do not write (unless '-a'). [TRUE]
	   "conflicting" => \$conflicting,     # write out both entries for conflicts? [FALSE]
	   "verbose" => \$verbose,	       # TRUE
	   "quiet" => \$quiet, 		       # turn off warnings
	   "help|?" => \&pod2usage # Help
	   );

#--- I/O error & defaults ---#
file_exists($database_file, "database");
my $db_path = get_file_path($database_file);
my $genbank_path = "$db_path/genbank";
$genbank_path = File::Spec->rel2abs($genbank_path);


#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# check for loci table #
table_exists($dbh, "loci");

# joining query options (for table join) #
# joining query options (for table join) #
my $join_sql = "";
$join_sql .= join_query_opts(\@subtype, "subtype");
$join_sql .= join_query_opts(\@taxon_id, "taxon_id");
$join_sql .= join_query_opts(\@taxon_name, "taxon_name");

# getting loci start - end from loci table (& genbank file names)
my $loci_se_r = get_loci_start_end($dbh, $query, $join_sql);

# getting CDS info for CDS features in loci regions #
my $loci_tbl_r = call_genbank_get_region($loci_se_r, $genbank_path);

# determine if CDS are already in gene table #
## if exists & '-a', existing entry written ## 
check_exists_in_gene_table($dbh, $loci_tbl_r, $all_genes, $conflicting) unless $existing; 

# determine if CDS are in operons #
check_in_CAS($loci_se_r, $loci_tbl_r);

# write locus-gene table #
write_loci_tbl($loci_tbl_r);

# disconnect #
$dbh->disconnect();
exit;


### Subroutines
sub check_exists_in_gene_table{
# determine if CDS are already in gene table #
## if exists & '-a', existing entry written ## 
  my ($dbh, $loci_tbl_r, $all_genes, $conflicting) = @_;
  
  # query gene table #
  my $query = "SELECT Locus_ID,Gene_Start,Gene_End,Gene_ID,Gene_Alias,Gene_Length__AA,In_CAS from Genes";
  my $genes_r = $dbh->selectall_arrayref($query);	
  
  # loading hash w/ entries: locus=>gene_id #
  my %exists;
  map{ $exists{$$_[0]}{$$_[3]} = [@$_[1..$#$_]] } @$genes_r;
  
  # checking for existing #
  my %status;
  foreach my $locus (sort keys %$loci_tbl_r){
    foreach my $feat (keys %{$loci_tbl_r->{$locus}}){
      # gene ID
      my $gene_id;
      if (! exists $loci_tbl_r->{$locus}{$feat}{db_xref}){
	printf STDERR "WARNING: 'db_xref' not find out %s -> %s. Skipping.\n",
	  $locus, $feat;
	next;
      }
      else{
	$gene_id = $loci_tbl_r->{$locus}{$feat}{db_xref};
      }
      
      # existing gene ID?
      if(exists $exists{$locus}{$gene_id} ){		# existing entry for gene
	$status{'existing'}++;
	
	# writing just conflicting entries 
	if($conflicting){				
	  $loci_tbl_r->{$locus}{"N$feat"} = $exists{$locus}{$gene_id};
	}
	# removing existing gene
	else{		
	  delete $loci_tbl_r->{$locus}{$feat};				
	  
	  # replacing existing 
	  if($all_genes){					
	    $loci_tbl_r->{$locus}{$feat} = $exists{$locus}{$gene_id};
	    $status{'replacing'}++;
	  }
	  else{ $status{'keeping'}++; }		# keeping existing entries 
	}
      }
      else{			# entry does not exists 
	if ($conflicting){		# just writing conflicting
	  delete $loci_tbl_r->{$locus}{$feat};
	}
	else{ $status{'adding'}++; }
      }
    }
  }
  
  # status #
  print STDERR "\n### gene entry new/existing/conflicting report ###\n";
  print STDERR "Number of already existing entries in CLdb: ", $status{'existing'}, "\n"
    if exists $status{'existing'};
  print STDERR "Number of entries to be replaced (existing in CLdb but written to output table): ", $status{'replacing'}, "\n"
    if exists $status{'replacing'};
  print STDERR "Number of entries to keep as they were (existing in CLdb & NOT written to output table): ", $status{'keeping'}, "\n"
    if exists $status{'keeping'};
  print STDERR "Number of new entries (no in CLdb & written to output table): ", $status{'adding'}, "\n"
    if exists $status{'adding'};
  print STDERR "###------------------------------------------###\n\n";
  #print Dumper %$loci_tbl_r; exit;
}


sub write_loci_tbl{
  # writing out loci table for manual editting of aliases #
  my ($loci_tbl_r) = @_;
  
  # header
  print join("\t", qw/Locus_ID Gene_Id Gene_start Gene_end Gene_length__AA In_CAS Gene_Alias Sequence/), "\n";

  # tag_ID ordering
  #my @tags = qw/df_xref start end In_CAS product translation/;
  
  # writing each features
  foreach my $loci (keys %$loci_tbl_r){
    foreach my $feature (keys %{$loci_tbl_r->{$loci}}){
      # checking db_xref for 'fig|peg'
      unless ($loci_tbl_r->{$loci}{$feature}{db_xref} =~ /fig\|.+peg/){		
	print STDERR " WARNING: Locus$loci -> $feature does not have a FIG-PEG ID in a db_xref tag!\n"
	  unless $quiet; 
      }
      
      # length already found? #
      my ($len_AA, $seq);
      if($loci_tbl_r->{$loci}{$feature}{translation} =~ /^\d+$/){   # existing entry; no sequence
	$len_AA = $loci_tbl_r->{$loci}{$feature}{translation};
	if($conflicting){ $seq = "existing_entry"; }
	else{ $seq = ""; }
      }
      else{ 
	$len_AA = length $loci_tbl_r->{$loci}{$feature}{translation}; 
	$seq = $loci_tbl_r->{$loci}{$feature}{translation};
      }
      
      # writing line #
      map{ $_ = "" unless exists $loci_tbl_r->{$loci}{$feature}{$_} }
	qw/start end In_CAS product translation db_xref/;
      
      print join("\t", $loci,
		 $loci_tbl_r->{$loci}{$feature}{db_xref},
		 $loci_tbl_r->{$loci}{$feature}{start},
		 $loci_tbl_r->{$loci}{$feature}{end},
		 $len_AA,
		 $loci_tbl_r->{$loci}{$feature}{In_CAS},
		 $loci_tbl_r->{$loci}{$feature}{product},
		 $loci_tbl_r->{$loci}{$feature}{translation}
		), "\n";
    }
  }
}


sub check_in_CAS{
  # checing whether CDS features are in the designated operon locations #
  my ($loci_se_r, $loci_tbl_r) = @_;
  
  foreach my $locus (keys %$loci_tbl_r){
    foreach my $feature (keys %{$loci_tbl_r->{$locus}}){
      #print Dumper "$locus -> $feature";
      die " LOGIC ERROR: $!\n" unless 
	exists $loci_se_r->{$locus};
      
      # defining start - end #
      ## f = feature; o = operon; c = crispr_array ##
      my $f_start = int $loci_tbl_r->{$locus}{$feature}{start};
      my $f_end = int $loci_tbl_r->{$locus}{$feature}{end};
      
      my $o_start = ${$loci_se_r->{$locus}}[3];
      my $o_end = ${$loci_se_r->{$locus}}[4];			
      
      my $c_start = ${$loci_se_r->{$locus}}[5];
      my $c_end = ${$loci_se_r->{$locus}}[6];
      
      # skipping if no CAS_start || CAS_end #
      next if ! defined $o_start || $o_start =~ /^$/
	|| ! defined $o_end || $o_end =~ /^$/;
      
      # making all on the + strand #
      ($f_start, $f_end) = set_to_pos_strand($f_start, $f_end);
      ($o_start, $o_end) = set_to_pos_strand($o_start, $o_end);
      ($c_start, $c_end) = set_to_pos_strand($c_start, $c_end) if $c_start && $c_end;		
      
      # determining location #
      ## gene must fall in operon, without falling into crispr array ##
      if($f_start >= $o_start && $f_end <= $o_end){	# gene in CAS
	# check for overlap w/ CRISPR array if present
	if($c_start && $c_end){					
	  # if not in crispr array	   
	  if( ($f_start < $c_start && $f_end < $c_end) ||
	      ($f_start > $c_start && $f_end > $c_end) ){	
	    $loci_tbl_r->{$locus}{$feature}{In_CAS} = "yes";
	  }
	  # if inside of CRISPR array span
	  else{    
	    # in crispr array, but array spans operon, so calling 'yes'
	    if( $c_start < $o_start && $c_end > $o_end ){    
	      $loci_tbl_r->{$locus}{$feature}{In_CAS} = "yes";
	    }
	    else{
	      # in crispr array & array does not span CAS operon, so not defined as in operon
	      $loci_tbl_r->{$locus}{$feature}{In_CAS} = "no"; 
	    }
	  }	
	}
	# in operon, no CRISPR array
	else{ $loci_tbl_r->{$locus}{$feature}{In_CAS} = "yes"; }	
      }
      # not in operon
      else{ $loci_tbl_r->{$locus}{$feature}{In_CAS} = "no"; }		
    }
  }
  #print Dumper %$loci_tbl_r; exit;
}

sub set_to_pos_strand{
# setting all start-end so start is <= end #
  my ($start, $end) = @_;
  return $start, $end if $start !~ /^\d+$/ || $end !~ /^\d+$/;
  if ($start > $end){ return $end, $start; }
  else{ return $start, $end; }
}

sub call_genbank_get_region{
# calling genbank_get_region to get CDS info from genbank files #
  my ($loci_se_r, $genbank_path) = @_;
  
  my %loci_tbl;
  foreach my $locus (keys %$loci_se_r){
    my $genbank_file = join("/", $genbank_path, ${$loci_se_r->{$locus}}[2]);
    die " ERROR: $genbank_file not found!\n"
      unless -e $genbank_file;
    
    
    my $start = ${$loci_se_r->{$locus}}[0];
    my $end = ${$loci_se_r->{$locus}}[1];
		my $scaffold = ${$loci_se_r->{$locus}}[7];
    my $CAS_status = ${$loci_se_r->{$locus}}[8];
    my $array_status = ${$loci_se_r->{$locus}}[9];
    print STDERR join("\n ",
		      "...Getting features in:",
		      "file:\t\t\t$genbank_file",
		      "scaffold:\t\t$scaffold",
		      "region:\t\t$start-$end",
		      "CAS_status:\t\t$CAS_status",
		      "Array_status:\t\t$array_status"), "\n";

    # return: {featID : feat_tag_ID : tag_value}
    my $ret_r;
    if($start <= $end){
      $ret_r = genbank_get_region($scaffold, $start, $end, $genbank_file);
    }
    else{
      $ret_r = genbank_get_region($scaffold, $end, $start, $genbank_file);
    }

         
    my @col_sel = qw/start end db_xref product translation/;
    
    ## loading loci table -- {locusID : featureID : feat_tag_ID : tag_value}
    foreach my $feat (keys %$ret_r){
      foreach my $colOfInt (@col_sel){
	if(exists $ret_r->{$feat}{$colOfInt}){
	  $loci_tbl{$locus}{$feat}{$colOfInt} = $ret_r->{$feat}{$colOfInt};
	}
	else{
	  print STDERR " WARNING: '$colOfInt' tag not found in ".
	    "feature: $feat for locus ID '$locus'! Skipping feature\n";			    
	}
      }
    }

    printf STDERR " Number of features found: %i\n", scalar keys %{$loci_tbl{$locus}};
  }
		
  # sanity check #
  unless (%loci_tbl){
    $dbh->disconnect();
    die "\nNo CDS found in any of the specified loci regions! Nothing to add to CLdb\n";
  }
		
  #print Dumper %loci_tbl; exit;
  return \%loci_tbl;
}

sub get_loci_start_end{
# getting locus start and locus end from loci table in db #
	my ($dbh, $query, $join_sql) = @_;
	
	my %loci_se;
	my $cmd = "SELECT Locus_id, Locus_start, Locus_end, 
			Genbank_File, CAS_start, CAS_end,
			Array_Start, Array_End, scaffold,
			CAS_Status, Array_Status
			from loci 
			where (CAS_start is not null or CAS_end is not null) 
			$join_sql $query";
	$cmd =~ s/[\t\n]+/ /g;

	my $loci_se_r = $dbh->selectall_arrayref($cmd);	
	die "ERROR: no locus start-end values found!\n" unless $loci_se_r;	

	foreach my $row (@$loci_se_r){
		$loci_se{$$row[0]} = [@$row[1..$#$row]];
		}
		
		#print Dumper %loci_se; exit;
	return \%loci_se;
	}



