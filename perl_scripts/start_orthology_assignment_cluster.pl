#! /usr/bin/perl
########################################################################
# Script name  :    start_orthology_assignment_cluster.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
# Start assignment of ESTs to ortholog groups using a computer cluster
# 
########################################################################
## INCLUDES
####################
use warnings;
use strict;
use File::Basename;
use Getopt::Std;
use File::Temp qw/ tempfile tempdir /;
require 'lib/Iostuff2.pm';

############## PARAMETERS START ##############
my %options;
getopt( q(tw), \%options );

### READ OPTIONS FOR ANALYS FROM options_file
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

my $minimum_length_of_hit = $options_hash{minimum_length_of_hit} || die "please specify minimum length of hit\n";

### CREATE DIRECTORIES
my $root_directory = $options_hash{root_directory} || die "please specify the root directory\n";
my $project_directory = $root_directory."/".$options_hash{project_name} || die "please specify the root directory\n";
my $fasta_directory = $options_hash{fasta_directory} || $root_directory."/est_libraries";
### CREATE TEMPORARY FILE
#my ($file_handle,$temp_bash_script) =  tempfile(UNLINK => 1);
#unless($temp_bash_script){
#    print "Could not create temporary filehandle: $! \n";
#}	

my $temp_bash_script= "tmp_script"; 
unless($temp_bash_script){
    print "Could not create temporary filehandle: $! \n";
}	


############## PARAMETERS END ##############

Iostuff::create_project_directories($project_directory);		    


####################
# MAIN
####################	
	my $a =0;    
MAIN: {		    
	Iostuff::print_start();	
	Iostuff::check_taxa_list_syntax($options_hash{taxa_list});
	my $taxa = `cat $options_hash{taxa_list}`;
	my @taxa_entries = split(/\n/,$taxa);
	foreach (@taxa_entries){
		/"(\w*)"\s?"(\w*)"\s?(p|n)?/;
		my $current_taxon = $1;  # Name of Taxon
		my $current_taxa_short = $2; # Shortcut for Taxon
		my $type_of_sequence = $3; # Nucleotide or protein data
		
		my $number_of_entries = `grep -c ">" $fasta_directory/$current_taxon.fa `;
		my $number_of_threads = $options_hash{no_threads};

		for(my $i=0; $i* $number_of_threads < $number_of_entries; $i++){
			my $start_sequence = $i * $number_of_threads;
			my $stop_sequence = ($i+1) * $number_of_threads;
			my $script_version  = (defined $type_of_sequence && $type_of_sequence eq "p") ? "start_blasto_protein.pl" : "start_blasto_nucleotide.pl";
			my $command = "perl perl_scripts/$script_version -t $current_taxon  -w $current_taxa_short -c $start_sequence -s $stop_sequence -m $minimum_length_of_hit ";
			Iostuff::write_bash_script($temp_bash_script,$command, $current_taxon);
			#exit;
			system("qsub  -o /dev/null -e /dev/null $temp_bash_script");
		}
	}
	Iostuff::print_start();
	print "Analysis finished\nThe next step is the selection of genes\n";
}

####################
# END
####################

