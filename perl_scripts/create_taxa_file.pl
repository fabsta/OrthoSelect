#! /usr/bin/perl
########################################################################
# Script name  :    nph-blast.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# . 
#
# Script to create the file containing all taxa present in the study
########################################################################
## INCLUDES
####################
use Getopt::Std;
use File::Basename;
use strict;
use warnings;
require 'lib/Iostuff2.pm';

####################
# VARIABLES
####################

my%options=();
getopt("to",\%options);

### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

# Directory with EST libraries
my $db_directory = $options_hash{fasta_directory}  || "est_libraries";

my @fasta_files = <$db_directory/*.fa>;
# Name of file containing taxa
my $taxa_file = $options{o} ||"taxa_list.txt";


####################
# MAIN
####################

Main:{
	print "\tCreating file $taxa_file that contains all taxa used for the study....";
open my $TAXA_FILE_FH,">",$taxa_file || die "Couldn't open '$taxa_file': $!";
foreach(@fasta_files){
	my $fasta_file = $_;
	chop(my $number_of_fasta_entries = `grep -c ">" $fasta_file`);
	(my $fasta_file_short = $fasta_file) =~ s/\.fa//g;

	my @phylip_split = split(/_/, basename($fasta_file_short));
	my $first_part = substr($phylip_split[0],0,7);
	my $second_part = substr($phylip_split[1],0,2);
	print {$TAXA_FILE_FH} "\"".(basename($fasta_file_short))."\" \"".$first_part."_".$second_part."\"\n" or warn "Could not write: $!";

}
close($TAXA_FILE_FH)|| warn "Couldn't close '$taxa_file': $!";
	print "done\n";
}

####################
# END
####################
