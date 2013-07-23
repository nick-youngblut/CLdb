#!/usr/bin/env perl

### modules
use strict;
use warnings;
use Pod::Usage;
use Data::Dumper;
use Getopt::Long;
use File::Spec;
use DBI;

### args/flags
pod2usage("$0: No files given.") if ((@ARGV == 0) && (-t STDIN));

my ($verbose, $database_file, $subject_in);
my ($taxon_id, $taxon_name, $delimiter, $just_hit_scafs);
GetOptions(
	   "database=s" => \$database_file,
	   "subject=s" => \$subject_in,			# subject fasta
	   "taxon_id=s" => \$taxon_id, 			# subject taxon_ID	(SUBJECT__[12])
	   "taxon_name=s" => \$taxon_name, 		# taxon_name		(SUBJECT__[12])
	   "delimiter=s" => \$delimiter,		# used if taxon_(name|ID) in subject of blast hit [TRUE]
	   "hits" => \$just_hit_scafs, 			# just loading scaffolds with hits
	   "verbose" => \$verbose,
	   "help|?" => \&pod2usage # Help
	   );

### I/O error & defaults
die " ERROR: provide a database file name!\n"
	unless $database_file;
die " ERROR: cannot find database file!\n"
	unless -e $database_file;
die " ERROR: either provide -taxon_name, -taxon_id, or -delimiter with -taxon_name SUBJECT or -taxon_id SUBJECT\n"
	unless ($delimiter && $taxon_id) || ($delimiter && $taxon_name) || $taxon_id || $taxon_name;
print STDERR " WARNING: no subject fasta provided! Subject sequence information will not be added to CLdb!\n"
	unless $subject_in;

$delimiter = qr/$delimiter/ if $delimiter;

### MAIN
# connect 2 db #
my %attr = (RaiseError => 0, PrintError=>0, AutoCommit=>0);
my $dbh = DBI->connect("dbi:SQLite:dbname=$database_file", '','', \%attr) 
	or die " Can't connect to $database_file!\n";

# checking for existence of blast_hits table #
my $table_list_r = list_tables($dbh);
check_for_blast_table($table_list_r);

# getting taxon names & IDs from database #
my ($taxon_info_r, $taxa_r) = get_db_taxon_info($dbh);

# loading blast hit table into db #
my $scafs_r;
($scafs_r, $taxon_id, $taxon_name) = get_blast_hits($dbh, $taxon_info_r, $taxa_r, $taxon_id, $taxon_name, $delimiter);

# loading subject fasta #
load_subject_fasta($dbh, $subject_in, $scafs_r, $taxon_id, $taxon_name) if $subject_in;

# disconnect #
$dbh->disconnect();
exit;


### Subroutines 
sub load_subject_fasta{
# loading fasta of subject with hits #
	my ($dbh, $subject_in, $scafs_r, $taxon_id, $taxon_name) = @_;
	
	# preparing sql #
	my $cmd = "INSERT INTO blast_subject(Taxon_id, Taxon_name, Scaffold_name, Scaffold_sequence) values (?,?,?,?)";
	my $sql = $dbh->prepare($cmd);
	
	# reading in fasta #
	open IN, $subject_in or die $!;
	my ($name, $seq);
	my $entry_cnt = 0;
	while(<IN>){
		chomp;
 		 s/#.+//;
 		next if  /^\s*$/;	
 		
 		if(/>.+/ || eof(IN) ){
 		 
 		 	# loading last sequence #
 			if($name){
	 			if(! $just_hit_scafs){ # keeping just hit scaffolds & exists scaff
	 				$entry_cnt = add_entry($sql, $taxon_id, $taxon_name, $name, $seq, $entry_cnt)
	 					if exists $scafs_r->{$name};	 			
	 				}		
		 		else{				# adding all scaffolds 
	 				$entry_cnt = add_entry($sql, $taxon_id, $taxon_name, $name, $seq, $entry_cnt)
		 			}
	 			}
	 		
	 		# getting new sequence
 			$_ =~ s/^>| .+//g;
 			$name = $_;
 			$seq = "";
 			}
 		else{$seq .= $_; }
		}
	close IN;
	
	$dbh->commit;
	
	print STDERR "...Number of blast subject entries added/updated in database: $entry_cnt\n";
	}
	
sub add_entry{
# adding gentry for blast subject #
	my ($sql, $taxon_id, $taxon_name, $name, $seq, $entry_cnt) = @_;
	
	$sql->execute( ($taxon_id, $taxon_name, $name, $seq) );
	if($DBI::err){
		print STDERR "ERROR: $DBI::errstr in: ", join("\t", $taxon_id, $taxon_name, $name, $seq), "\n";
		}
	else{ $entry_cnt++; }

	return $entry_cnt;
	}

sub get_blast_hits{
# getting blast hits to a genome #
	my ($dbh, $taxon_info_r, $taxa_r, $taxon_id, $taxon_name, $delimiter) = @_;

	# loading entry #
	my $cmd = "INSERT INTO blast_hits(Spacer_group,Subject,pident,length,mismatch,gapopen,qstart,qend,sstart,send,evalue,bitscore,CRISPR_array,Taxon_ID,Taxon_name) values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)";
	my $sql = $dbh->prepare($cmd);

	# data stream of blast hits #
	my $entry_cnt = 0;
	my %scafs;
	while(<>){
		chomp;
		next if /^\s*$/;
		my @line = split /\t/;
		die " ERROR: blast hit table should have 12-14 columns!\n"
			unless $#line >= 11 && $#line <= 13;
		
		# removing group from groupID #
		$line[0] =~ s/^Group//;
		
		# parsing subject if needed #
		if($delimiter){
			my @t_s = split /$delimiter/, $line[1], 2;
			die " ERROR: the delimiter did not split the string: \"$line[1]\"\n"
				unless scalar @t_s == 2;
			
			# SUBJECT #
			if($taxon_id && $taxon_id =~ /^SUBJECT/){
				my @tmp = split /__/, $taxon_id;
				if($tmp[1] && $tmp[1] == 2){
					$taxon_id = $t_s[1];
					$line[1] = $t_s[0];		# subject
					}
				else{
					$taxon_id = $t_s[0];	
					$line[1] = $t_s[1];		# subject
					}
				}
			elsif($taxon_name && $taxon_name eq "SUBJECT"){
				my @tmp = split /__/, $taxon_name;
				if($tmp[1] && $tmp[1] == 2){				
					$taxon_name = $t_s[1];	
					$line[1] = $t_s[0];		# subject
					}
				else{
					$taxon_name = $t_s[0];
					$line[1] = $t_s[1];		# subject
					}			
				}
			else{ die " ERROR: delimiter used, but neither taxon_id nor taxon_name == \"SUBJECT\"\n"; }
			}
		
		# taxon ID #
		$taxon_id = "" unless $taxon_id;
		$taxon_name = "" unless $taxon_name;
	
		# check for existence of subjects in locus table #
		if( $line[13] && ! exists $taxa_r->{"taxon_id"}{$line[13]} ){ 
			print STDERR " WARNING: taxon_id -> \"$line[13]\" not found in Loci table! The Genome is not in CLdb!\n";		
			}
		if( $line[14] && ! exists $taxa_r->{"taxon_name"}{$line[14]} ){ 
			print STDERR " WARNING: taxon_name -> \"$line[14]\" not found in Loci table! The Genome is not in CLdb!\n";
			}
		
		# getting date of loading #
		
		
		# loading db #
		$sql->execute( @line, $taxon_id, $taxon_name );
		if($DBI::err){
			print STDERR "ERROR: $DBI::errstr in: ", join("\t", @line), "\n";
			}
		else{ $entry_cnt++; }
		
		# loading scafs #
		$scafs{$line[1]} = 1 unless $just_hit_scafs;
		}
		
	$dbh->commit;	
	
	print STDERR "...Number of BLAST hit entries added/updated in database: $entry_cnt\n";
	
	return \%scafs, $taxon_id, $taxon_name;
	}

sub get_db_taxon_info{
	my ($dbh) = @_;
	
	my $query = "SELECT taxon_id, taxon_name from loci";
	my $ret = $dbh->selectall_arrayref($query);
	
	my %taxa;
	foreach my $row (@$ret){
		$taxa{"taxon_id"}{$$row[0]} = 1;
		$taxa{"taxon_name"}{$$row[1]} = 1;
		}
		#print Dumper @$ret; exit;
	return $ret, \%taxa;		# taxon_id, taxon_name
	}

sub check_for_blast_table{
	my ($table_list_r) = @_;
	die " ERROR: loci table not found in database!\n"
		unless grep(/^blast_hits$/i, @$table_list_r);
	}

sub list_tables{
	my $dbh = shift;
	my $all = $dbh->selectall_hashref("SELECT tbl_name FROM sqlite_master", 'tbl_name');
	return [keys %$all];
	}

__END__

=pod

=head1 NAME

CLdb_loadBlastHits.pl -- loading spacer BLASTn hits (& subject sequences)

=head1 SYNOPSIS

CLdb_loadBlastHits.pl [flags] < blast_table.txt

=head2 Required flags

=over

=item -d 	CLdb database.

=back

=head2 Optional flags

=over

=item -subject

Subject used for blasting.

=item -taxon_id

Taxon_ID of subject (FIG ID)

=item -taxon_name

Taxon_name of subject (taxonomic classification)

=item -delimiter

Character(s) delimiting taxon_id/taxon_name & scaffold of subject.
Example: column2 in BLAST file = "E.coli__scaffold1"

=item -hits

Load just scaffolds with spacer blast hits? [TRUE]

=item -v	Verbose output. [TRUE]

=item -h	This help message

=back

=head2 For more information:

perldoc CLdb_loadBlastHits.pl

=head1 DESCRIPTION

Loading spacer BLAST results into CLdb for further analysis.

Either a Taxon_ID or or Taxon_name must be provided.
If the taxon_id or taxon_name is in the 'subject' column
of the blast table, provide "SUBJECT" for -taxon_id or -taxon_name.
If the ID/name is in the 2nd part of the subject name,
provide "SUBJECT_2".

=head1 EXAMPLES

=head2 Usage:

CLdb_loadBlastHits.pl -da CLdb.sqlite -subject 1.F.A.1A.3_peg_merged.fna < spacers_1.F.A.1A.3_blast.txt -taxon_name "1.F.A.1A.3"

=head1 AUTHOR

Nick Youngblut <nyoungb2@illinois.edu>

=head1 AVAILABILITY

sharchaea.life.uiuc.edu:/home/git/CLdb/

=head1 COPYRIGHT

Copyright 2010, 2011
This software is licensed under the terms of the GPLv3

=cut

