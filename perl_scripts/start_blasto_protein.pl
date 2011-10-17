#!/usr/bin/perl
########################################################################
# Script name  :    start_blasto_protein.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
#
#  BLASTS LIBRARY OF PROTEIN SEQUENCES AGAINST OTHOLOG DATABASE
########################################################################
## INCLUDES
####################
use Getopt::Std;
use File::Basename;
use Bio::Perl;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Data::Dumper;
use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
#use strict;
use warnings;

require 'lib/Blaststuff.pm';
require 'lib/Iostuff2.pm';

####################
# VARIABLES
####################
my %options;
getopt( q(tdwics), \%options );

my $exitcode = $? >> 8; # exit code
my $exitsignal = $? & 127; # signal that caused exit




############## PARAMETERS START ##############
#BEGIN{
	my %options_hash = %{ Iostuff::read_options("options.txt") };
	my $root_directory = $options_hash{root_directory} || die "Please specify root directory\n";
$ENV{BLASTDATADIR} = $options_hash{blastdb} || $root_directory."/db/";
$ENV{BLASTDB} = $options_hash{blastdb} || $root_directory."/db/"; 
	my $tmp_dir = "/tmp/";
	my %abb_taxa_hash =(
AT=> q(1),
CE => q(1),
DM=> q(1),
HS=> q(1),
SC=> q(1),
SP=> q(1),
EC=> q(1));



##### FOR SPLITTED ANALYSES ONLY
my $start_entry = $options{c} || 0;
my $stop_entry = $options{s} || 1000000000;

###### USER PARAMETERS
my $incrementor = $options{i} || "";
### E-Value
my $e_value      = $options_hash{e_value} || q(1e-10);
### NAME OF TAXON
my $current_taxon = $options{t} || die qq(please specify Taxon..\n);

#my $taxa_counter = $options{z} || "0";

my $taxa_shortcut = $options{w} || die qq(please specify shortcut for Taxon\n);

### Orthology Database Options
my $database_type = $options_hash{orthology_database_type} || die qq(Please specify type of orthology database\n);
my $database_file = $options_hash{orthology_database_file} || "OG.db";

my $blast_database = $options_hash{blastdb} || $root_directory."/db/";

### Read orthologous group information
my %cluster_list =();
my $kog_file ="";
if($database_type eq "o"){
	%cluster_list = %{ Iostuff::read_orthomcl_list($blast_database."/OrthoMCL_list.txt")};
	$protein_template_source = $blast_database."/seqs_orthomcl-2.fasta";
}
if($database_type eq "k"){
	%cluster_list = %{ Iostuff::read_kog_list($blast_database."/KOG_list.txt")};
	$kog_file = $blast_database."/KOG_list.txt";
	$protein_template_source = $blast_database."/kyva";
}
if($database_type eq "m"){
	$kog_file = $blast_database."/oma.txt";
}

### DIRECTORY CONTAINING THE FASTA FILES
my $fasta_directory = $options_hash{fasta_directory} || $root_directory."/est_libraries/";
### OUTPUT DIRECTORIES
my $output_directory    = $root_directory.$options_hash{project_name};
my $annotation_file = $output_directory."/annotations/".$current_taxon."".$incrementor.".txt";


### COUNTS NUMBER OF FASTA ENTRIES
my $no_current_fasta_entry      = 1;

my %cluster_groups;

my $translation = $options_hash{translation};

### WEIGHTED BLASTO-ANALYSIS
my $weight = 1;
my $taxa_in_db = 7;


### TEMPORARY FILES
#my ($template_nucleotide_fh, $template_nucleotide_file) = tempfile(UNLINK => 1);
my ($template_nucleotide_fh, $template_nucleotide_file) = tempfile();
die "Could not create temporary filehandle: $! \n" if !(-e $template_nucleotide_file);
#my ($temp_alignment_fh_in, $temp_alignment_file_in) = tempfile(UNLINK => 1);
my ($temp_alignment_fh_in, $temp_alignment_file_in) = tempfile();
die "Could not create temporary filehandle: $! \n" if !(-e $temp_alignment_file_in);
#my ($temp_blast_query_fh, $temp_blast_query_file) = tempfile(UNLINK => 1);
my ($temp_blast_query_fh, $temp_blast_query_file) = tempfile();
die "Could not create temporary filehandle: $! \n" if !(-e $temp_blast_query_file);
my ($temp_blast_output_fh, $temp_blast_output_file) = tempfile();
die "Could not create temporary filehandle: $! \n" if !(-e $temp_blast_output_file);

#$temp_blast_query_file = "blastquery";

### SETTING BLAST PARAMETERS
my @blast_params = (
        program  => q(blastp),
		outfile => $temp_blast_output_file
    );
     print qq(Current Taxon is $current_taxon\n);

  my $blast_factory = Bio::Tools::Run::StandAloneBlast->new(@blast_params);
### SET DATABASE
    $blast_factory->database($database_file);
### SET E-Value
    $blast_factory->e($e_value);
# Turn filtering off (to avoid unhandable error messages)
   $blast_factory->F("F");

### FOR EACH FASTA FILE
    my $inner_counter = 0;
    my $fasta_file    = $fasta_directory."/" . $current_taxon . ".fa";
#}  
    
############## PARAMETERS END ##############
  
####################
# MAIN
####################
    
MAIN: {
	
	Iostuff::print_start();	
    
    ### IF WHOLE FILE HAS TO BE PROCESSED    
      if(!defined $options{c} && !defined $options{s}){
        $temp_blast_query_file = $fasta_file;
    }
    else{
#### OPEN FASTA FILE
    my $input_fasta_sequence = Bio::SeqIO->new(
        -file   => $fasta_file,
        -format => "fasta"
    );
### FLAG FOR CORRECT PARSING THE BLAST OUTPUT
    while ( my $current_query_seq = $input_fasta_sequence->next_seq ) {
	 my $cluster_group;   
### BLAST EST SEQUENCE FROM A USER-DEFINED RANGE OF THE EST LIBRARY (DEFAULT: all)
	if($no_current_fasta_entry < $start_entry ){$no_current_fasta_entry++;next;}
	if($no_current_fasta_entry >= $stop_entry){last;}
        #print "\tSelect Query Sequence No. ".$no_current_fasta_entry++;
		write_sequence(">>$temp_blast_query_file", 'fasta', $current_query_seq) || warn "couldnt write BLAST query file\n";
	$no_current_fasta_entry++;
	}
    }
    ### FLAG FOR CORRECT PARSING THE BLAST OUTPUT
my $parse_check;

### BLAST SEQUENCE
#        my $blast_report = $blast_factory->blastall($current_query_seq);
    print "\tPerforming BLAST Search\n";
    print "\tFile: $fasta_file\n";
    chomp(my $num_seqs = `grep ">" -c $temp_blast_query_file`);
    print "\t#Seqs: $num_seqs\n";
    print "\tEstimated time: ".($num_seqs/1.94)." seconds\n";
    
#    print "cat $temp_blast_query_file | blastall -p blastp -d /c1/scratch/fabian/projects/OrthoSelect_geobio/evaluation/db/OG.db  -o $temp_blast_output_file\n";
#print "cat $temp_blast_query_file | blastall -p blastp -e $e_value -d $root_directory/db/$database_file  -o $temp_blast_output_file\n";
`cat $temp_blast_query_file | blastall -p blastp -e $e_value -d $root_directory/db/$database_file  -o $temp_blast_output_file`;
#    print "cat blast_query | blastall -p blastx -d /c1/scratch/fabian/projects/OrthoSelect_geobio/evaluation/db/OG.db  -o out\n";

### HANDLE BLAST ERRORS
	if (my $signal = $? & 127) {warn "External command BLAST exited on signal $signal\n";}
	elsif (our $! = $? >> 8) {$! = $exitcode;warn "External command BLAST existed with error: ".$!."\n";}
	        my %sequence_hash = ();
        my $id = "";
print "\tReading sequence file into memory...\n";

#strace  time perl perl_scripts/kog_binning.pl -t Chaetoderma_nitidulum  -w Chaetoderma_nitidulum -m 10 
### Slurping File
my $fh = new IO::Handle;
open $fh, '<', $fasta_file;
my $text = "";
sysread $fh, $text, -s $fh;
$fh->flush;
close $fh;
print "\tAll Sequences read into memory\n";        
my @whole_lines = split(/\n/,$text);
my $curr_sequence ="";
my $debug_counter = 0;
my $debug_counter2 = 0;
my $curr_header = "";
my $curr_domain = "";
foreach(@whole_lines){
    	my $line = $_;
	if($line =~ />(\w*)/){
            $id = $1;
            if(!defined $id){
                warn "There was a problem with header $line\n";
                next;
            }
          #  $sequence_hash{$id} = ">$id";
        }
        else{
      #  print "\tremembering $id\n";
            $sequence_hash{$id} .= $line;
        }
}
         my $blast_report = new Bio::SearchIO(
                -format => 'blast',
#                -file   => 'out',
                -file   => $temp_blast_output_file
                );
                ### PARSE BLAST RESULTS
    while ( my $result = $blast_report->next_result) {
        my $accession_number = $result->query_name();
#	$accession_number =~ /(\w*)\|(\w*)/;
	$accession_number =~ /(\w*)/;
	$accession_number = $1;
        ( my $id = $taxa_shortcut . "_" . $accession_number ) =~ s/\s//;
      #  print "For result $id\n";
              #exit;
        if ( !$blast_report ) { print qq(\t--> No BLAST HITS\n); next;}
### PARSE BLAST-RESULT
my %cluster_groups = ();
    if($options_hash{orthology_database_type} eq "o"){
	    $parse_check = Blaststuff::parse_blast_orthomcl( $blast_report_file, $current_taxon, $id, $output_directory, $incrementor,\%cluster_list, $weight, $taxa_in_db);
    	}
	else{
		$parse_check = Blaststuff::parse_blast_kog({ blast_report=>\$result, best_ortho_groups => \%cluster_groups,current_taxon=>$current_taxon, current_id=>$id, output_directory=>$output_directory, incrementor=>$incrementor,cluster_list_ref=>\%cluster_list, weight=>$weight, taxa_in_db=>$taxa_in_db, kog_file=>$kog_file});

#		$parse_check = Blaststuff::parse_blast_kog( $blast_report_file, $current_taxon, $id, $output_directory, $incrementor,\%cluster_list, $weight, $taxa_in_db,$kog_file);
    	}
   # 	print "\t\tback from parsing...with ".keys(%cluster_groups)."\n";
### NOT HITS FOR BLAST SEARCH	
	if ( !  $parse_check) {
		print qq(\t-->No BLAST HIT (reference not defined)\n); 
		next;
	}
     else {
  #   	 %cluster_groups = %{$parse_check};
   #   	if ( keys %cluster_groups == 0) {
    #      print qq(\t--> No BLAST HIT (cluster_groups == 0)\n); 
     #     next;
	#	}
  ### PROCESS EACH BEST OG    
      	foreach(keys %cluster_groups){
        	$cluster_group = $_;
#			if(!exists $useful_kogs{$cluster_group}){
					#print "\t\tKog $cluster_group not usefull\n";
	#				next;
		#		}
#		print "\tHiiiiiiit\n";
### ANNOTATION FOR OG
	my $annotation_arg = $cluster_groups{$cluster_group}{"annotation"};
            my $seqobj = Bio::Seq->new(
                    -seq        => $sequence_hash{$accession_number},,
                       -display_id => $id,
                       -alphabet   => q(protein)
                   );

### CREATE DIRECTORY FOR ORTHOLOGOUS GROUP
		    if ( !-e "$output_directory/basis_hits/$cluster_group" ) {
			    mkdir("$output_directory/basis_hits/$cluster_group", 0777 )
			    or die q(could not create hit-directory);
		    }
## SAVE PROTEIN SEQUENCE
	my $prot_file = "$output_directory/basis_hits/$cluster_group/".$cluster_group."_".$current_taxon."_prot_hits".$incrementor.".fasta";
	print "Saving to $prot_file\n";
#	next;
#if (sysopen(LOGFILE, $prot_file, O_EXCL|O_CREAT|O_WRONLY)) {
	 open (FILE, ">>$prot_file");    # open the file
	 flock FILE, 2;                 # try to lock the file
	write_sequence( ">>$prot_file",'fasta', $seqobj) || warn "couldnt write ".$seqobj->id."\n";
	flock FILE, 8;                 # unlock the file
        close(FILE);                   # close the file
	print "DIED  - did not really save\n" if ! -e $prot_file;

#}
## LOG ANNOTATIONS
	open $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
	flock $annotation_fh, 2;                 # try to lock the file
	print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";

	flock $annotation_fh, 8;
	close $annotation_fh or warn "Could not close: $!";


#	open my $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
#	print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";
#	close $annotation_fh or warn "Could not close: $!";
	}
        }
    }
}

	
####################
# END
####################





