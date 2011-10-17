#! /usr/bin/perl -w
########################################################################
# Script name  :    start_filter_redundant_single.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <Fab.Schreiber@gmail.com>
# 
# Selection of "best" orthologs using on a single machine
#
# 
########################################################################
## INCLUDES
####################
use Getopt::Std;
use File::Basename;
use warnings;
use strict;
require 'lib/Iostuff2.pm';

############## PARAMETERS START ##############

my %options;
getopt( q(t), \%options );


### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };


my $project_name = $options_hash{project_name} || die "please specify project_name for analysis\n";

my $input_directory = $options_hash{distance_calculation_input} || die "please specify input directory for analysis\n";

my $root_directory = $options_hash{root_directory} || die "please specify root directory for analysis\n";
my @directories_array = <$root_directory/$project_name/$input_directory/*>;
die "Specified directory does not exist or is empty\n" if (scalar(@directories_array) == 0);

my $output_directory = $options_hash{distance_calculation_output} || die "please specify output directory for analysis\n";
   $output_directory = "$root_directory/$project_name/$output_directory";
my $alignment_method = $options_hash{alignment_method} || die "please specify alignment method for analysis\n";
my $different_alignment_method = $options_hash{different_alignment_method};
my $translated_sequences_folder = $options_hash{translated_sequences_folder} || die "please specify taxon database for analysis\n";

my $distance_matrix_type = $options_hash{distance_matrix_type} || die "please specify distance matrix type for analysis\n";


############## PARAMETERS END ##############


####################
# MAIN
####################

MAIN: {
	### CREATE OUTPUT DIRECTORY
	Iostuff::create_single_directory($output_directory);

	foreach my $current_directory (@directories_array){
	#	next if !$current_directory =~ /KOG0009/;
	#	print ("nice perl perl_scripts/select_best_ortholog.pl -d ".$current_directory." -o $output_directory -m $alignment_method -t $distance_matrix_type -a $different_alignment_method");
	#	exit;
		system("nice perl perl_scripts/select_best_ortholog.pl -d ".$current_directory." -o $output_directory -m $alignment_method -t $distance_matrix_type -x $translated_sequences_folder -a $different_alignment_method");
		#exit;
	}
Iostuff::print_start();
print "Analysis finished\n";
}


#nice perl perl_scripts/select_best_ortholog.pl -d /c1/scratch/fabian/projects/OrthoSelect_geobio/OrthoSelect/merged_with_tricho//new/selected_genes/KOG0841 -o /c1/scratch/fabian/projects/OrthoSelect_geobio/OrthoSelect/merged_with_tricho//new/reduced_genes -m m -t g -a
####################
# END
####################
