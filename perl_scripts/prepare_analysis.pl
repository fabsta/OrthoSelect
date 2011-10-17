#! /usr/bin/perl
########################################################################
# Script name  :    prepare_analysis.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
#
# PREPARES THE SYSTEM FOR SUBSEQUENT ANALYSIS
########################################################################
## INCLUDES
####################
use Getopt::Std;
use File::Basename;
use strict;
use warnings;
my $LINE_LENGTH = 64;

require 'lib/Iostuff2.pm';


############## PARAMETERS START ##############
my%options=();
getopt("to",\%options);

### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };


my $orthology_database_type = $options_hash{orthology_database_type} || die "please specify type of orthology database\n";
my $orthology_database_name = "OG.db";

my $fasta_directory = $options_hash{fasta_directory} || "est_libraries";
my $fasta_ok = Iostuff::check_directory_exists($fasta_directory);
	die "Please make sure $fasta_directory exists / is writeable" if !$fasta_ok;
my $blast_directory = $options_hash{blastdb} || "db";
my $blast_ok = Iostuff::create_single_directory($blast_directory);
	
	die "Please make sure $blast_directory exists / is writeable" if !$blast_ok;
	
my @fasta_files = <$fasta_directory/*.fa>;

my $taxa_file = $options_hash{taxa_list} || die "please specify name for file containing taxa name\n";

my @required_programs = (
		"ESTScan",
		"blastall",
		"muscle",
		"t_coffee",
		"Gblocks",
#		"wget",
		"clustalw",
		#"noisy",
		"genewise"
		);

############## PARAMETERS END ##############

####################
# MAIN
####################
MAIN:{


### CHECK IF REQUIRED PROGRAMS ARE PRESENT
	print "=" x $LINE_LENGTH, "\n";
	print "\tChecking Installed Programs \n";
	print "=" x $LINE_LENGTH, "\n";

my @req_prog_errors = @{ Iostuff::check_required_programs(\@required_programs,1)};
	if (scalar @req_prog_errors){
		print "\nInitial tests are unable to find the following required programs:\n";
		foreach (@req_prog_errors)	{
			print $_,"\t";
		}
	die "\n\nPlease make sure these programs can be found in your PATH or retry to automatically download them\n";
	}
	
### CHECK IF FASTA FILES CONTAIN SEQUENCES
	print "=" x $LINE_LENGTH, "\n";
	print "\tChecking Correct Fasta Format \n..............\n";
	print "=" x $LINE_LENGTH, "\n";

foreach(@fasta_files){
	my $current_fasta_file = $_;
	my $number_of_sequences = `grep -c ">" $current_fasta_file`;
	warn "No Sequences in File $current_fasta_file or not in FASTA format\n" if $number_of_sequences == 0;

	Iostuff::check_fasta_format($current_fasta_file);

}


unless(-e "$blast_directory/OG.db.phr"){
	### DOWNLOAD ORTHOLOG DATABASE
	if($orthology_database_type eq "k"){
		my $db_download_ok = Iostuff::download_kog_database($blast_directory,$orthology_database_name);
		print "....Sucessfully set up kog database\n" if $db_download_ok;
		my $db_blast_ok = Iostuff::check_blast_db(qq($blast_directory/$orthology_database_name));
		print "...Sucessfully tested kog database\n" if $db_blast_ok;
		
	}
	if($orthology_database_type eq "o"){
		my $db_download_ok = Iostuff::download_orthomcl_database($blast_directory, $orthology_database_name);
		print "...Sucessfully set up orthomcl database\n" if $db_download_ok;
		my $db_blast_ok = Iostuff::check_blast_db(qq($blast_directory/$orthology_database_name));
		print "...Sucessfully tested orthomcl database\n" if $db_blast_ok;
	}
}
	

### CREATE TAXA FILE
### TRANSLATE ALL EST SEQUENCES
	print "=" x $LINE_LENGTH, "\n";
	print "\tGenerating file containing all taxa in the study \n..............\n";
	print "=" x $LINE_LENGTH, "\n";
	system("perl perl_scripts/create_taxa_file.pl -o $taxa_file ");
#	Iostuff::print_start();
#	print "Analysis finished\nYou can start the orthology assignment now by typing e.g.: perl start_orthology_assignment_single.pl\n";

### TRANSLATE ALL EST SEQUENCES
	print "=" x $LINE_LENGTH, "\n";
	print "\tTranslating all EST sequences \n..............\n";
	print "=" x $LINE_LENGTH, "\n";
	my @grep_protein_data = `grep \" p\" $taxa_file`;
	if(!@grep_protein_data){
		system("perl perl_scripts/translate_ests.pl -i est_libraries -o $options_hash{translated_sequences_folder}");
	}
	Iostuff::print_start();
	print "Analysis finished\nYou can start the orthology assignment now by typing e.g.: perl start_orthology_assignment_single.pl\n";
}


####################
# END
####################
