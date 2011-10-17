#!/usr/bin/perl
########################################################################
# Script name  :    gene_selection.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <Fab.Schreiber@gmail.com>
# 
# What this script does 
#
# 
########################################################################
## INCLUDES
####################
use Getopt::Std;
use Bio::Perl;
use File::Basename;
use Data::Dumper;

require 'lib/Iostuff2.pm';
############## PARAMETERS START ##############
%options=();
getopt("dt",\%options);
my $LINE_LENGTH = 64;

### READ OPTIONS FOR ANALYSIS
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };



### READ FILE CONTAINING TAXA OF INTEREST
my $taxa_file = $options_hash{taxa_list} || die "please specify location of file containing the taxa (e.g. taxa_list.txt)\n";
## Check Syntax of File
Iostuff::check_taxa_list_syntax($options_hash{taxa_list});

my %taxon_hash = %{  Iostuff::read_taxa($options_hash{taxa_list}) };
my $orthology_database_type = $options_hash{orthology_database_type} || die "please specify type of orthology database\n";
my $annotation_file = "db/KOG_list.txt";

my $statistic_directory = $options_hash{statistic_directory} || die "please specify location of hits\n";

### Directory containing the hits
my $hits_directory = $options{d} || $options_hash{root_directory}."/".$options_hash{project_name}."/".$statistic_directory;

my @directory_array =<$hits_directory/*>; 

my %result_hash = ();


############## PARAMETERS END ##############

####################
# MAIN
####################

MAIN: {
	print "OG\tAnnotation\t";
	foreach my $taxon(sort keys %taxon_hash){
		print $taxon."\t";
	}
	print "\n";
	if($statistic_directory =~ /basis_hits/){
		&stat_for_basis();
	}
	elsif($statistic_directory =~ /selected_genes/){
		&stat_for_basis();
	}
	else{
		&stat_for_others();
	}
	
}



####################
sub stat_for_basis(){
####################
	foreach my $current_directory (@directory_array){
		my $OG = basename($current_directory);
			my $taxa_in_folder = `grep ">" $current_directory/*prot_hits.fasta`;
			my %taxa_count = ();
			foreach(keys %taxon_hash){
				$taxa_count{$_}=1 if($taxa_in_folder =~ /$taxon_hash{$_}/);				
				#$count_taxa++ if($taxa_in_folder =~ /$_/);
				#print $_." with $count_taxa \n";
			}
				
	my $taxa_for_alignment_counter = 0;
		
	if($orthology_database_type eq "o"){
		print $OG."\tNN\t";
	}
	else{
		chomp(my $annot = `grep $OG $annotation_file`);
		$annot =~ m/\[\w*\]\sKOG\d*\s(.*)/;
		print $OG."\t".$1."\t";
	}
		foreach(sort keys %taxon_hash){
			if(exists $taxa_count{$_}){print "1\t"; $taxa_for_alignment_counter++;$result_hash{$_}++;}
			else{print "0\t";}
		}
			print $taxa_for_alignment_counter."|".(keys %taxon_hash)."\n";	
	}
	print "Overall\tNN\t";
	foreach(sort keys %taxon_hash){
		print $result_hash{$_}."|".(scalar @directory_array)."\t";		
	}
	print "NN\t\n";

}


sub stat_for_others(){
	foreach my $current_directory (@directory_array){
		my $OG = basename($current_directory);
	#	print $OG."\ngrep \">\" $current_directory/*final.fasta\n\n";
		my $taxa_in_folder = `grep ">" $current_directory/*final.fasta`;
			my %taxa_count = ();
			foreach(keys %taxon_hash){
				$taxa_count{$_}=1 if($taxa_in_folder =~ /$taxon_hash{$_}/);				
				
				$count_taxa++ if($taxa_in_folder =~ /$taxon_hash{$_}/);
		#		print $_." with $count_taxa \n" if $_ =~ /Trich/;
			}
			
			
	#		foreach(sort keys %taxa_count){
	#			print $_."--> ".$taxa_count{$_}."\n";
	#		}
#			print Dumper %taxa_count;	
#			exit if $OG eq "";
	my $taxa_for_alignment_counter = 0;
	
	if($orthology_database_type eq "o"){
		print $OG."\tNN\t";
	}
	else{
		chomp(my $annot = `grep $OG $annotation_file`);
		$annot =~ m/\[\w*\]\sKOG\d*\s(.*)/;
		print $OG."\t".$1."\t";
	}
	foreach(sort keys %taxon_hash){
		if(exists $taxa_count{$_}){print "1\t"; $taxa_for_alignment_counter++;$result_hash{$_}++;}
		else{print "0\t";}
	}
	print $taxa_for_alignment_counter."|".(keys %taxon_hash)."\n";	
		#exit if $OG eq "KOG0009";

	}
	print "Overall\tNN\t";
	foreach(sort keys %taxon_hash){
		print $result_hash{$_}."|".(scalar @directory_array)."\t";		
	}

	print "NN\t\n";

}

####################
# END
####################






