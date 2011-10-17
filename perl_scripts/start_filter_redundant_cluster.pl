#! /usr/bin/perl -w
########################################################################
# Script name  :    gene_selection.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <Fab.Schreiber@gmail.com>
# 
# Selection of "best" orthologs using a computer cluster 
#
# 
########################################################################
## INCLUDES
####################
use File::Basename;
use Getopt::Std;
use File::Temp qw/ tempfile/;
use warnings;
use strict;
require 'lib/Iostuff2.pm';

############## PARAMETERS START ##############

my %options;
getopt( q(tw), \%options );


### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

my $root_directory = $options_hash{root_directory} || die "Please specify root directory\n";


my $input_directory = $options_hash{distance_calculation_input} || die "please specify input directory for analysis\n";
	$input_directory = $root_directory."/".$options_hash{project_name}."/".$input_directory;
my @directories_to_analyse = <$input_directory/*>;
my $output_directory = $options_hash{distance_calculation_output} || die "please specify input directory for analysis\n";
	$output_directory = $root_directory."/".$options_hash{project_name}."/".$output_directory;

	
	
if(scalar(@directories_to_analyse) == 0){
	die "No Folder in directory $input_directory\nMake sure the right folder is selected\n";
}
my $alignment_method = $options_hash{alignment_method} || die "please specify alignment method for analysis\n";
my $different_alignment_method = $options_hash{different_alignment_method};
if($alignment_method eq "o" and $different_alignment_method eq ""){
	die "Please specify usage of the other alignment method\n";
}
my $distance_matrix_type = $options_hash{distance_matrix_type} || die "please specify distance matrix type for analysis\n";
my $translated_sequences_folder = $options_hash{translated_sequences_folder} || die "please specify taxon database for analysis\n";



### CREATE TEMPORARY FILE
	    my ($file_handle,$temp_bash_script) =  tempfile(UNLINK => 1);
	    unless($temp_bash_script){
		    print "Could not create temporary filehandle: $! \n";
	    }	
############## PARAMETERS END ##############

$temp_bash_script = "bash_script.sh";
####################
# MAIN
####################
MAIN: {
### CREATE OUTPUT DIRECTORY
	Iostuff::create_single_directory($output_directory);
	
	foreach my $current_directory (@directories_to_analyse){
		my $command = "perl perl_scripts/select_best_ortholog.pl -d $current_directory -o $output_directory -m $alignment_method -t $distance_matrix_type  -x $translated_sequences_folder -a $different_alignment_method";
		Iostuff::write_bash_script($temp_bash_script,$command, basename($current_directory));
		#print $command."\n";
		system("qsub -o /dev/null -e /dev/null $temp_bash_script");
		}
	Iostuff::print_start();
print "Analysis finished\n";

}


####################
# END
####################

