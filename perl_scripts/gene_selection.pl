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
use Getopt::Std;
use Bio::Perl;
use File::Basename;
use Data::Dumper;

require 'lib/Iostuff2.pm';

%options=();
getopt("dt",\%options);
my $LINE_LENGTH = 64;

### READ OPTIONS FOR ANALYSIS
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

#Iostuff::print_options(\%options_hash, "gene_selection_option,required_taxa_file");


### GENE SELECTION OPTION
my $gene_selection_option = $options_hash{gene_selection_option} || die "please specify gene selection strategy\n";

### READ FILE CONTAINING TAXA OF INTEREST
my $required_taxa_file = $options_hash{required_taxa_file} || die "Please specify location of file containing the required taxa (e.g. required_taxa_list.txt)\n";
## Check Syntax of File
Iostuff::check_taxa_list_syntax($options_hash{required_taxa_file});

my %required_taxa_hash = %{ Iostuff::read_required_taxa($required_taxa_file)};
#print Dumper(%required_taxa_hash);
### READ FILE CONTAINING FUNCTIONAL CLASSES
my $required_fun_classes_file = $options_hash{blastdb}."/".$options_hash{fun_file} || die "please specify location of functional annotation file (e.g. fun.txt)\n";
my %required_fun_classes_hash = %{ Iostuff::read_fun_class_file($required_fun_classes_file)};

### Directory containing the hits
my $hits_directory = $options{d} || $options_hash{root_directory}.$options_hash{project_name}."/basis_hits";

my @directory_array =<$hits_directory/*>; 

my %taxa_count =();

my $name_for_taxa_set = $options_hash{gene_selection_directory} || die "please specify output directory for selected genes (e.g. selected_genes) \n";

my $output_directory = $options_hash{root_directory}.$options_hash{project_name}."/".$name_for_taxa_set;

#my $output_directory = $options_hash{root_directory}.$options_hash{project_name}."/test";



### CREATE OUTPUT DIRECTORY
Iostuff::create_single_directory($output_directory);

my $log_file = $options_hash{root_directory}.$options_hash{project_name}."/stats/".$name_for_taxa_set.".txt";


MAIN: {
	if($gene_selection_option eq "s"){
		Iostuff::check_taxa_list_syntax($options_hash{required_taxa_file});
		my %required_taxa_hash = %{ Iostuff::read_required_taxa($required_taxa_file)};
		
		&select_genes_taxalist(\@directory_array,\%required_taxa_hash,$output_directory);
	}
	if($gene_selection_option eq "m"){
		my %required_taxa_hash = %{ Iostuff::read_required_taxa_monophylum($required_taxa_file)};
	#	exit;
	#print Dumper %required_taxa_hash;
		&select_genes_taxalist_monophylum(\@directory_array,\%required_taxa_hash,$output_directory);
	}
	
	
	#if($gene_selection_option eq "f"){
	#	&select_genes_fun_class(\@directory_array,\%required_fun_classes_hash,$output_directory);
	#}
	Iostuff::print_start();
	print "Analysis finished\nFilter redundant sequences from each OG by typing: perl start_filter_redundant_single.pl\n";

}



####################################################################
sub select_genes_taxalist($$$){
####################################################################
my @directory_array = @{ ($tmp = shift) };
my %required_taxa_hash = %{ ($tmp2 = shift) };
my $output_directory = shift;
	print "=" x $LINE_LENGTH, "\n";
	print "\tSelecting Genes according to given taxa \n";
	print "=" x $LINE_LENGTH, "\n";

foreach my $current_directory (@directory_array){
	my $taxa_in_folder = `ls $current_directory/*prot*`;
	print "\t".basename($current_directory).":\t";

	#print $taxa_in_folder;
	my $count_taxa =0;
	foreach(keys %required_taxa_hash){
		$count_taxa++ if($taxa_in_folder =~ /$_/);
	#	print $_." with $count_taxa \n";
	}
	
## COPY DIRECTORY IF CONTAINED ALL TAXA 
	#print "$count_taxa ".(keys %required_taxa_hash)."\n";
	if($count_taxa == keys %required_taxa_hash){
		Iostuff::copy_folder($current_directory, $output_directory);
		print "will be selected\n";
	}
	else{
		print "will be ignored\n";
		
	}
	}

}

####################################################################
sub select_genes_taxalist_monophylum($$$){
####################################################################
my @directory_array = @{ ($tmp = shift) };
my %required_taxa_hash = %{ ($tmp2 = shift) };
my $output_directory = shift;
	print "=" x $LINE_LENGTH, "\n";
	print "\tSelecting Genes according to given monophyla \n";
	print "=" x $LINE_LENGTH, "\n";

#print Dumper %required_taxa_hash ;

foreach my $current_directory (@directory_array){
#	my $taxa_in_folder = `ls $current_directory/*prot*`;
	print "\n\n".basename($current_directory).":\n\n";
	my $taxa_in_folder = `grep ">" $current_directory/*prot*`;
	my $monophylum_counter = 0;
#	print $taxa_in_folder;
	my $count_taxa =0;
monophylum_loop:
	foreach my $monophylum (keys %required_taxa_hash){
#monophylum_loop: 
	print $monophylum.": ";
		foreach my $monophylum_species (keys %{$required_taxa_hash{$monophylum}}){
			if($taxa_in_folder =~ /$monophylum_species/){
				print $monophylum_species."\t";
				$monophylum_counter++; 
			#	print $monophylum_counter."\n";
		next monophylum_loop;
			}
		}
		#print "\n";
	}
	
## COPY DIRECTORY IF CONTAINED ALL TAXA 
	if($monophylum_counter == keys %required_taxa_hash){
		Iostuff::copy_folder($current_directory, $output_directory);
		print "\t".basename($current_directory)." will be selected\n";
	
		#`echo "basename($current_directory) >> einsle"`;
	}

}
}


####################################################################
sub select_genes_fun_class($$$){
####################################################################
my @directory_array = @{ ($tmp = shift) };
my %mol_functions_hash = %{ ($tmp2 = shift) };
my $output_directory = shift;

foreach my $current_directory (@directory_array){
	my $current_kog = basename($current_directory);
	
	## COPY DIRECTORY IF CONTAINED ALL TAXA
	my $fun_class = Iostuff::get_fun_class4kog($current_kog);
	if(exists  $mol_functions_hash{$fun_class}){
	#	print $fun_class."\n";
		Iostuff::copy_folder($current_directory, $output_directory);
	}
}

}











=head1 NAME

gene_selection.pl - Implements gene selection strategies used in phylogenomic analysis

=head1 SYNOPSIS


  gene_selection.pl
	[-t Configuration file]
	

=head1 DESCRIPTION

This script implements gene selection strategies used in phylogenomic analysis.

On the one hand genes can be selected to be further included in the study

based on the 

=item existence of sequences of specified taxa 

=item affilitation of the gene to a specified functional class of molecules


=head2 An Example 

C<perl gene_selection.pl -t config.txt>

Read 'config.txt' and performs the gene_selection according to the parameters

set in the 'config.txt' file

=head1 SEE ALSO

L<perlpodspec>

=head1 COPYRIGHT

Copyright 2008 Fabian Schreiber <fschrei@gwdg.de>.

Permission is granted to copy, distribute and/or modify this 
document under the terms of the GNU Free Documentation 
License, Version 1.2 or any later version published by the 
Free Software Foundation; with no Invariant Sections, with 
no Front-Cover Texts, and with no Back-Cover Texts.

=head1 AUTHOR

Fabian Schreiber <fschrei@gwdg.de>


=cut

