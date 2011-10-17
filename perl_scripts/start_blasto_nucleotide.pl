#!/usr/bin/perl
########################################################################
# Script name  :    start_blasto_nucleotide.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
#
# BLASTS EST LIBRARY AGAINST ORTHOLOG DATABASE AND TRANSLATES ESTS
########################################################################
## INCLUDES
####################
use Getopt::Std;
use File::Basename;
use Fcntl qw(:DEFAULT);
use Bio::Perl;
use Pod::Usage;
use Bio::AlignIO;
use Bio::SearchIO;
use Bio::Tools::BPlite;
use Bio::Tools::Genewise;
use Bio::LocatableSeq;
use Bio::Tools::Run::StandAloneBlast;
use Data::Dumper;
use Pod::Usage;
use File::Temp qw/ tempfile tempdir /;
use strict;
use warnings;
require 'lib/Blaststuff.pm';
require 'lib/Iostuff2.pm';

####################
# VARIABLES
####################
my %options;
getopt( q(tdwicsm), \%options );

my $exitcode = $? >> 8; # exit code
my $exitsignal = $? & 127; # signal that caused exit




############## PARAMETERS START ##############

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
### SHORTCUT FORM OF TAXON
my $taxa_shortcut = $options{w} || die qq(please specify shortcut for Taxon\n);
### Minimum length of Hit
my $minimum_length_of_hit =  $options{m} || "10";

### Orthology Database Options
my $database_type = $options_hash{orthology_database_type} || die qq(Please specify type of orthology database\n);
my $database_file = $options_hash{orthology_database_file} || "OG.db";

my $blast_directory = $root_directory."/db/";
my $blast_database = $root_directory."/db/OG.db";

my $protein_template_source ="";
### Read orthologous group information
my %cluster_list =();
my $kog_file ="";
if($database_type eq "o"){
	%cluster_list = %{ Iostuff::read_orthomcl_list($blast_directory."/OrthoMCL_list.txt")};
	$protein_template_source = $blast_directory."/seqs_orthomcl-2.fasta";
}
if($database_type eq "k"){
	%cluster_list = %{ Iostuff::read_kog_list($blast_directory."/KOG_list.txt")};
	$kog_file = $blast_directory."/KOG_list.txt";
	$protein_template_source = $blast_directory."/kyva";
}
if($database_type eq "m"){
	$kog_file = $blast_directory."/oma.txt";
}

### DIRECTORY CONTAINING THE FASTA FILES
my $fasta_directory = $options_hash{fasta_directory} || $root_directory."/est_libraries/";
### OUTPUT DIRECTORIES
my $output_directory    = $root_directory.$options_hash{project_name};
my $annotation_file = $output_directory."/annotations/".$current_taxon."".$incrementor.".txt";


### COUNTS NUMBER OF FASTA ENTRIES
my $no_current_fasta_entry = 1;

my %cluster_groups =();
### METHOD FOR TRANSLATION
my $translation = $options_hash{translation};

### WEIGHTED BLASTO-ANALYSIS
my $weight = 1;
my $taxa_in_db = 7;


### PARAMETERS FOR FILE LOCKING
my $lock_file = 2;
my $unlock_file = 8;

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

### SETTING BLAST PARAMETERS
my @blast_params = (
        program  => q(blastx),
		outfile => $temp_blast_output_file
    );
 #    print qq(Current Taxon is $current_taxon\n);

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

############## PARAMETERS END ##############


####################
# MAIN
####################
    
MAIN: {
	
#	Iostuff::print_start();	
       if(! -e $fasta_file){
	warn "\t\tThere was a problem: The file $fasta_file does not exists.Aborting....\n";
	exit;
    }
 
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
    my $parse_check;
   
### BLAST SEQUENCE
#        my $blast_report = $blast_factory->blastall($current_query_seq);
    print "\tTaxon: $current_taxon\n";
    print "\tFile: $fasta_file\n";
    chomp(my $num_seqs = `grep ">" -c $temp_blast_query_file`);
    print "\t#Seqs: $num_seqs\n";
    print "\tEstimated time for Blast Search: ";
    printf "%d seconds (or less)\n",($num_seqs/1);
    print "\tPerforming BLAST Search...";
    
   # print "cat $temp_blast_query_file | blastall -p blastx -e $e_value -d $blast_database\n";
   `cat $temp_blast_query_file | blastall -p blastx -e $e_value -d $blast_database  -o $temp_blast_output_file`;
    print "done\n";
#    `cat $temp_blast_output_file`;
### HANDLE BLAST ERRORS
	if (my $signal = $? & 127) {warn "External command BLAST exited on signal $signal\n";}
	elsif (our $! = $? >> 8) {$! = $exitcode;warn "External command BLAST existed with error: ".$!."\n";}
	#next;
        my %sequence_hash = ();
        my $id = "";
### Slurping File
	print "\tReading sequence file into memory...";
	my $fh = new IO::Handle;
	open $fh, '<', $fasta_file;
	my $text = "";
	sysread $fh, $text, -s $fh;
	$fh->flush;
	close $fh;
	print "done\n";        
	my @whole_lines = split(/\n/,$text);
	my $curr_sequence ="";
	my $debug_counter = 0;
	my $debug_counter2 = 0;
	my $curr_header = "";
	my $curr_domain = "";
	
	print "\tParsing BLAST results and saving sequences in appropriate folders. This can take some time....";
	
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
	            $sequence_hash{$id} .= $line;
		}
	}
         my $blast_report = new Bio::SearchIO(
                -format => 'blast',
#                -file   => 'out',
                -file   => $temp_blast_output_file
                );
      #   print Dumper %sequence_hash;
      #   exit;
         
### PARSE BLAST RESULTS
    while ( my $result = $blast_report->next_result) {
        my $accession_number = $result->query_name();
	$accession_number =~ /(\w*)\|?(\w*)/;
	$accession_number = $1;
        ( my $id = $taxa_shortcut . "_" . $accession_number ) =~ s/\s//;
	#print "For result $id\n";
              #exit;
        if ( !$blast_report ) { print qq(\t--> No BLAST HITS\n); next;}
### PARSE BLAST-RESULT
my %cluster_groups = ();
	if($options_hash{orthology_database_type} eq "o"){
		$parse_check = Blaststuff::parse_blast_orthomcl( $temp_alignment_file_in, $current_taxon, $id, $output_directory, $incrementor,\%cluster_list, $weight, $taxa_in_db);
    	}
	else{
		$parse_check = Blaststuff::parse_blast_kog({ blast_report=>\$result, best_ortho_groups => \%cluster_groups, current_taxon=>$current_taxon, current_id=>$id, output_directory=>$output_directory, incrementor=>$incrementor,cluster_list_ref=>\%cluster_list, weight=>$weight, taxa_in_db=>$taxa_in_db, kog_file=>$kog_file});
    	}
### NO HITS FOR BLAST SEARCH
	if(!defined $parse_check){
	#	print qq(\t--> No BLAST HIT\n); 
		next;
	}
#	%cluster_groups = %{$parse_check};
	if ( keys %cluster_groups == 0) {
         # print qq(\t--> No BLAST HIT\n); 
          next;
	}
	else{
	      
### PROCESS EACH BEST OG
	foreach(keys %cluster_groups){
		my $cluster_group = $_;
### DEFAULT SCORE FOR TRANSLATION		
		my ($estscan_score, $genewise_score) = 1000; 
### ANNOTATION FOR OG	
	my $annotation_arg = $cluster_groups{$cluster_group}{"annotation"};
### ACCESSION NUMBER OF BEST HIT OF OG
my $cluster_template = $cluster_groups{$cluster_group}{"template_acc"};
	$cluster_template =~ /.*\|(.*)/;
	$cluster_template = $1 if $database_type eq "o";
#	$cluster_template =~ /\|//;
### SAVE EST SEQUENCE TEMPORARY (FOR TRANSLATION)
#print "Should use $accession_number";
        if(! exists $sequence_hash{$accession_number}){
            warn "There was a problem finding $accession_number in the fasta file\n";
        }
	my $seqobj = Bio::Seq->new(
                    -seq        => $sequence_hash{$accession_number},
                       -display_id => $id,
                       -alphabet   => q(dna)
                   );
	write_sequence( ">$template_nucleotide_file",'fasta', $seqobj) || warn "couldnt write ".$seqobj->id."\n";

### GET PROTEIN SEQUENCE FOR BEST HIT FROM OG
	my $record = Iostuff::getfastarecord($protein_template_source,$cluster_template);
	# print " get $protein_template_source from $cluster_template is $record \n";
	#exit;	 
	if (!$record eq ""){
            $record =~ s/\s//g;
### SAVE BEST HIT FROM OG AS TEMPLATE/REFERENCE SEQUENCE
	my $template_object = Bio::Seq->new(
                    -seq        => $record,
                    -display_id => "protein_template",
                    -alphabet   => q(protein)
                   );
	write_sequence( ">$temp_alignment_file_in",'fasta', $template_object) || warn "couldnt write ".$template_object->id."\n";

#### TRANSLATION using GENEWISE
	my $genewise_rec = `genewise $temp_alignment_file_in $template_nucleotide_file -pep -silent -quiet`;
	$genewise_rec =~ s/.*Mak.*\n// if $genewise_rec =~ /Making/;
	$genewise_rec =~ s/.*>/>/s if not $genewise_rec =~ /^>/;
	$genewise_rec =~ s/\/\///g;
	    
#### TRANSLATION using BIOPERL (6-FRAME)
	my %translated_hash = %{Iostuff::translate_bioperl($seqobj)};
#### TRANSLATION using ESTSCAN
	my $estscan_rec = `ESTScan  $template_nucleotide_file -t - `;

#	goto SKIP_COMPARING;
#### COMPARE SCORES FROM TRANSLATION AGAINST REFERENCE SEQUENCE
	my $FHE;
	#GENEWISE vs. protein template
	unless($genewise_rec eq "" || ! defined $genewise_rec){
		open($FHE,">$temp_alignment_file_in");
		print {$FHE} $genewise_rec."\n>protein_template\n".$record."\n";
		close($FHE);
		my $genewise_result = Iostuff::test_translation_quality($temp_alignment_file_in, $template_nucleotide_file, $minimum_length_of_hit);
		$genewise_score = defined $genewise_result? $genewise_result : 1000;
	}
	#BIOPERL vs. protein template
	my %translate_score_hash =();
	foreach(keys %translated_hash){
		my $frame = $_;
		my $counter =1;
		open($FHE,">$temp_alignment_file_in");
		print {$FHE} ">nucleotide_template\n".$translated_hash{$frame}{seq}."\n>protein_template\n".$record."\n";
		close($FHE);	
		my $resuelt = Iostuff::test_translation_quality($temp_alignment_file_in, $template_nucleotide_file, $minimum_length_of_hit);
		$translated_hash{$frame}{score} = defined $resuelt? $resuelt : 1000;
	#	print "score $frame: ".$translated_hash{$frame}{score}."\n"; 
	}
	my @sorted_translation_scores = sort { $translated_hash{$b}{score} <=> $translated_hash{$a}{score} } keys %translated_hash;
	my $translate_score = $sorted_translation_scores[-1];
	
	#ESTSCan vs. protein template

	unless($estscan_rec eq ""){
		open($FHE,">$temp_alignment_file_in");
		print {$FHE} $estscan_rec."\n>protein_template\n".$record."\n";
		close($FHE);	
		my $estscan_result = Iostuff::test_translation_quality($temp_alignment_file_in, $template_nucleotide_file, $minimum_length_of_hit);
		$estscan_score = defined $estscan_result? $estscan_result : 1000;
	}
### IF ALL TRANSLATIONS ARE TOO BAD
	$genewise_score = 1000 if ! defined $genewise_score;
	if($genewise_score > 1 && $translated_hash{$translate_score}{score} > 1 && $estscan_score > 1){
		$annotation_arg .= "\tNo Translation\t".$genewise_score."\t".$translated_hash{$translate_score}{score}."\t".$estscan_score."\t";
	}
	else{
		my $final_record ="";
## ESTSCAN WAS BEST
	if($genewise_score >= $estscan_score && $translated_hash{$translate_score}{score} >= $estscan_score){
		$estscan_rec =~ s/>.*//g;
		$estscan_rec =~ s/\s//g;
		$final_record = $estscan_rec;
		$annotation_arg .= "\tESTScan\t".$genewise_score."\t".$translated_hash{$translate_score}{score}."\t[".$estscan_score."]\t";

	}
## GENEWISE WAS BEST
	if($estscan_score > $genewise_score && $translated_hash{$translate_score}{score} >= $genewise_score ){
		$genewise_rec =~ s/>.*//g;
		$genewise_rec =~ s/\s//g;
		$final_record = $genewise_rec;	
		$annotation_arg .= "\tGeneWise\t[".$genewise_score."]\t".$translated_hash{$translate_score}{score}."\t".$estscan_score."\t";
	}
## BIOPERL WAS BEST
	if($genewise_score > $translated_hash{$translate_score}{score} && $estscan_score > $translated_hash{$translate_score}{score}){
		$final_record = $translated_hash{$translate_score}{seq};
		$annotation_arg .= "\tBioPerl\t".$genewise_score."\t[".$translated_hash{$translate_score}{score}."]\t".$estscan_score."\t";
	}
	SKIP_COMPARING:
		$estscan_rec =~ s/>.*//g;
		$estscan_rec =~ s/\s//g;
	my $final_protein = Bio::Seq->new(
                    -seq        => $estscan_rec,
                       -display_id => $id,
                       -alphabet   => q(protein)
                   );
		$genewise_rec =~ s/>.*//g;
		$genewise_rec =~ s/\s//g;

	my $final_protein2 = Bio::Seq->new(
                    -seq        => $genewise_rec,
                       -display_id => $id."E",
                       -alphabet   => q(protein)
                   );
	### CREATE DIRECTORY FOR ORTHOLOGOUS GROUP
		    if ( !-e "$output_directory/basis_hits/$cluster_group" ) {
                     #   print "making dir: $output_directory/basis_hits/$cluster_group\n";
			    mkdir("$output_directory/basis_hits/$cluster_group", 0777 )
			    or die q(could not create hit-directory);
		    }
## SAVE EST SEQUENCE in appropriate OG-DIRECTORY
	my $nucl_file = "$output_directory/basis_hits/$cluster_group/".$cluster_group."_".$current_taxon."_nucl_hits".$incrementor.".fasta";
## SAVE TRANSLATED EST SEQUENCE in appropriate OG-DIRECTORY
	my $prot_file = "$output_directory/basis_hits/$cluster_group/".$cluster_group."_".$current_taxon."_prot_hits".$incrementor.".fasta";
	#print "Saving to $prot_file\n";

	open (FILE, ">>$nucl_file");    # open the file
	flock FILE, 2;                 # try to lock the file
	write_sequence( ">>$nucl_file",'fasta', $seqobj) || warn "couldnt write ".$seqobj->id."\n";
	flock FILE, 8;                 # unlock the file
        close(FILE);                   # close the file
	
	open (FILE, ">>$prot_file");    # open the file
	flock FILE, 2;                 # try to lock the file
 	write_sequence( ">>$prot_file",'fasta', $final_protein) || warn "couldnt write ".$final_protein->id."\n" ;
 #	write_sequence( ">>$prot_file",'fasta', $final_protein2) || warn "couldnt write ".$final_protein2->id."\n" ;
	flock FILE, 8;                 # unlock the file
        close(FILE);                   # close the file
	
	print "\tDIED  - did not really save\n" if ! -e $prot_file;
}
## LOG ANNOTATIONS
		$annotation_arg =~ s/\t1000/\tNN/g;
		my $annotation_fh;
		open $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
		flock $annotation_fh, 2;                 # try to lock the file
		print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";
		flock $annotation_fh, 8;
		close $annotation_fh or warn "Could not close: $!";
	}
	else{
		$annotation_arg .= "\tNo Translation possible (E.g.Could not get protein template)";
		my $annotation_fh ="";
		open  $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
		flock $annotation_fh, 2;                 # try to lock the file
		print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";
		flock $annotation_fh, 8;
		close $annotation_fh or warn "Could not close: $!";
	}
	 }
        }
    }

}




####################
# END
####################

