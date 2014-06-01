#!/usr/bin/env perl

=pod

=head1 NAME

CLdb_checkArrayCoords.pl -- checking that array start-end in the loci table match with the Spacers & DRs tables

=head1 SYNOPSIS

CLdb_checkArrayCoords.pl [flags]

=head2 Required flags

=over

=item -database  <char>

CLdb database.

=back

=head2 Optional flags

=over

=item -prefix  <char>

Output file prefix

=item -verbose  <bool>

Verbose output. [TRUE]

=item -help  <bool>

This help message

=back

=head2 For more information:

perldoc CLdb_checkArrayCoords.pl

=head1 DESCRIPTION

Checking to make sure that array_start &
array_end in the loci table match
the start-end positions of each array
as designated by the array files.

The array_start & array_end values
in the loci table are assumed to be correct.

If an array file (& the corresponding entries in
the Spacers and DRs table) is wrong,
a new array file is writen.

=head1 EXAMPLES

=head2 Usage:

CLdb_checkArrayCoords.pl -d CLdb.sqlite 

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
use Set::IntervalTree;
use List::Util qw/min max/;

# CLdb #
use FindBin;
use lib "$FindBin::RealBin/../lib";
use lib "$FindBin::RealBin/../lib/perl5/";
use CLdb::query qw/
		    table_exists
		    n_entries/;
use CLdb::utilities qw/
			file_exists 
			connect2db
			get_file_path
		      /;
	

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file);
my $prefix = "coord__";
GetOptions(
	   "prefix=s" => \$prefix,
	   "database=s" => \$database_file,
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );


#--- I/O error & defaults ---#
file_exists($database_file, "database");
my $db_path = get_file_path($database_file);

#--- MAIN ---#
# connect 2 db #
my $dbh = connect2db($database_file);

# database metadata #
my $table_list_r = list_tables();
check_for_loci_table($table_list_r);
my $column_list_r = list_columns("loci", 1);


# getting input loci table #
my ($loci_r) = get_loci_table($dbh, $column_list_r);

# getting input spacer & DR table info
my $spacers_col_r = list_columns("spacers", 1);
my $spacers_r = get_array_table($dbh, $spacers_col_r, 'spacers', 'spacer_id');
my $DRs_col_r = list_columns('DRs', 1);
my $DRs_r = get_array_table($dbh, $DRs_col_r, 'DRs', 'dr_id');
my $DR_min_max_r = get_min_max_pos($DRs_r);

# checking to see if DR_min_max falls within array-start-end for locus_id
my $chk_r = check_array_start_end($loci_r, $column_list_r, $DR_min_max_r);

# correcting array files
correct_array_files( $loci_r, $column_list_r, $chk_r, $prefix, $db_path);

# disconnect to db #
$dbh->disconnect();
exit;

### Subroutines
sub correct_array_files{
  my ($loci_r, $column_list_r, $chk_r, $prefix, $db_path) = @_;

  foreach my $taxon (keys %$chk_r){
    foreach my $scaf (keys %{$chk_r->{$taxon}}){
      foreach my $locus_id (keys %{$chk_r->{$taxon}{$scaf}}){
	
	# loading arary_file 
	my $array_file = ${$loci_r->{$taxon}{$scaf}{$locus_id}}[ $column_list_r->{array_file} ];
	open IN, "$db_path/array/$array_file" or die $!;
	print STDERR "Writing corrected position file to: '$db_path/array/$array_file'\n";

	  # output file
	(my $outfile = $array_file) =~ s/^/$prefix/;
	open OUT, ">$outfile" or die $!;
	
	my $diff = $chk_r->{$taxon}{$scaf}{$locus_id}{coord_diff};
	
	while(<IN>){
	  chomp;
	  next if /^\s*$/;
	  
	  my @l = split /\t/;
	  $l[0] = $l[0] - $diff;
	  $l[3] = $l[3] - $diff if defined $l[3];

	  print OUT join("\t", @l), "\n";
	}
	close IN;
	close OUT;
      }
    }
  }
}

sub check_array_start_end{
  my ($loci_r, $column_list_r, $DR_min_max_r) = @_;

  # print Dumper %$loci_r; exit;
  my %chk;
  foreach my $taxon (keys %$loci_r){
    my %dup_array_chk;
    foreach my $scaf (keys %{$loci_r->{$taxon}}){
      foreach my $locus_id (keys %{$loci_r->{$taxon}{$scaf}}){
	my $array_start = ${$loci_r->{$taxon}{$scaf}{$locus_id}}[ $column_list_r->{array_start} ];
	my $array_end = ${$loci_r->{$taxon}{$scaf}{$locus_id}}[ $column_list_r->{array_end} ];

	next unless $array_start and $array_end;
	
	# flipping to + strand if needed
	($array_end, $array_start) = ($array_start, $array_end)
	  if $array_start > $array_end;

	# checking positions
	print STDERR "WARNING: cannot find DR-min-max for locus_id: $locus_id\n"
	  unless exists $DR_min_max_r->{$locus_id};
	
	unless( $array_start <= $DR_min_max_r->{$locus_id}{min} &&
		$array_end >= $DR_min_max_r->{$locus_id}{max} ){
	  print STDERR "WARNING: $taxon => $locus_id & array positions conflict: ",
	    join(",", $array_start, $array_end, 
		 $DR_min_max_r->{$locus_id}{min}, $DR_min_max_r->{$locus_id}{max}
		), "\n";
	  
	  # keeping track of array errors
	  $chk{$taxon}{$scaf}{$locus_id}{coord_diff} = $DR_min_max_r->{$locus_id}{min} - $array_start;
	}
	
	# checking length
	my $loci_array_len = $array_end - $array_start;
	my $tbl_array_len = $DR_min_max_r->{$locus_id}{max} - $DR_min_max_r->{$locus_id}{min};
	unless( abs($loci_array_len - $tbl_array_len) <= 100 ){
	  print STDERR "WARNING: $taxon => $locus_id & array lengths conflict: $loci_array_len\t$tbl_array_len\n";
	
	  $chk{$taxon}{$scaf}{$locus_id}{len_diff} = $loci_array_len - $tbl_array_len;
	}
      }
    }
  }


  #print Dumper %chk; exit;
  return \%chk;
}

sub get_min_max_pos{
  my ($tbl) = @_;

  # start_pos = 2; end_pos = 3
  my %min_max;
  foreach my $locus_id (keys %$tbl){
    my @starts;
    map{ push @starts, ${$tbl->{$locus_id}{$_}}[2] } keys %{$tbl->{$locus_id}};   
    $min_max{$locus_id}{min} = min @starts;
    
    my @ends;
    map{ push @ends, ${$tbl->{$locus_id}{$_}}[3] } keys %{$tbl->{$locus_id}};   
    $min_max{$locus_id}{max} = max @starts;
  }

 # print Dumper %min_max; exit;
  return \%min_max;
}

sub get_array_table{
# getting values from array table to check for errors
  my ($dbh, $column_list_r, $tbl, $id) = @_;

  # getting loci table #
  my $cmd = join(" ", "SELECT", 
		 join(",", sort{$column_list_r->{$a}<=> $column_list_r->{$b}} keys %$column_list_r), 
		 "from $tbl");
  my $rows_r = $dbh->selectall_arrayref($cmd);
  
  die " ERROR: no entries in loci table!\n" unless @$rows_r;

  # hash: taxon_name=>locus_id=>scaffold=>entry # 
  my %tbl;
  foreach my $row (@$rows_r){
	  # uninitialized to blank #
	  map{$_ = "" unless $_ }  @$row;
	  
	  # loading hash #
	  my $locus_id = $$row[$column_list_r->{"locus_id"}];
	  my $element_id = $$row[$column_list_r->{$id}];
	  $tbl{$locus_id}{$element_id} = $row;
	}

  #print Dumper %tbl; exit;
  return \%tbl;
  
}

sub get_loci_table{
# getting whole loci table #
	my ($dbh, $column_list_r) = @_;
	
	# getting loci table #
	my $cmd = join(" ", "SELECT", 
			join(",", sort{$column_list_r->{$a}<=> $column_list_r->{$b}} keys %$column_list_r), 
			"from loci");
	my $loci_r = $dbh->selectall_arrayref($cmd);
	
	die " ERROR: no entries in loci table!\n" unless @$loci_r;
	
	# hash: taxon_name=>locus_id=>scaffold=>entry # 
	my %loci;
	foreach my $row (@$loci_r){
		# uninitialized to blank #
		map{$_ = "" unless $_ }  @$row;
				
		# loading hash #
		my $taxon_name = $$row[$column_list_r->{"taxon_name"}];
		my $locus_id = $$row[$column_list_r->{"locus_id"}];
		my $scaffold = $$row[$column_list_r->{"scaffold"}];
		$loci{$taxon_name}{$scaffold}{$locus_id} = $row;
		}
		#print Dumper %loci; exit;
	return \%loci;
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

sub list_columns{
	my ($column_list, $silent_ret) = @_;

	my $all = $dbh->selectall_arrayref("pragma table_info($column_list)");
		#print Dumper $$all[0]; exit;
	my %tmp;
	for my $i (0..$#$all){ 
		$$all[$i][1] =~ tr/A-Z/a-z/;		# lower case for matching
		$tmp{$$all[$i][1]} = $i; 
		}
	if($silent_ret){ return \%tmp; }
	else{  print "Columns:\n", join(",\n", keys %tmp), "\n\n";  exit; }
	}


