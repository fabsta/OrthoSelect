#! /usr/bin/perl
########################################################################
# Script name  :    start_orthology_assignment_single.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
# Start assignment of ESTs to ortholog groups using a single computer
# 
########################################################################
## INCLUDES
####################
use warnings;
use strict;
use File::Basename;
use Getopt::Std;

require 'lib/Iostuff2.pm';

############## PARAMETERS START ##############
my %options;
getopt( q(tw), \%options );


### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

my $minimum_length_of_hit = $options_hash{minimum_length_of_hit} || die "please specify minimum length of hit\n";


### CREATE DIRECTORIES
my $root_directory = $options_hash{root_directory} || die "please specify the root directory\n";
my $project_directory = $root_directory."/".$options_hash{project_name} || die "please specify the root directory\n";


my $fasta_directory = $options_hash{fasta_directory} || $root_directory."/est_libraries";

############## PARAMETERS END ##############

### CREATE DIRECTORIES
Iostuff::create_project_directories($project_directory);


####################
# MAIN
####################

MAIN: {
	Iostuff::print_start();	
## READ FILE CONTAINING TAXA TO ANALYSE
	Iostuff::check_taxa_list_syntax($options_hash{taxa_list});
	my $taxa = `cat $options_hash{taxa_list}`;
	my @taxa_entries = split(/\n/,$taxa);
	foreach (@taxa_entries){
		/"(\w*)"\s?"(\w*)"\s?(p|n)?/;
		my $current_taxon = $1;  # Name of Taxon
		my $current_taxa_short = $2; # Shortcut for Taxon
		my $type_of_sequence = $3; # Nucleotide or protein data
		print "\nAnalysing $current_taxon ESTs\n";
		my $script_version  = (defined $type_of_sequence && $type_of_sequence eq "p") ? "start_blasto_protein.pl" : "start_blasto_nucleotide.pl";
		my $command = " perl_scripts/$script_version -t $current_taxon  -w $current_taxa_short -m $minimum_length_of_hit";
		system("perl $command"); # start analysis
	}
	Iostuff::print_start();
	print "Analysis finished\nThe next step is the selection of genes\n";
}

####################
# END
####################
