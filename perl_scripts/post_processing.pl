#!/usr/bin/perl -w
########################################################################
# Script name  :    gblock_alignments.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# This script is a simple client interface for displaying the blast 
# search form, running the blast search, and generating the search 
# result. 
#
# 
########################################################################
use Getopt::Std;
use Bio::Perl;
use Bio::SeqIO;
use File::Copy;
use Cwd;
use File::Basename;
use strict;
use warnings;

require 'lib/Iostuff2.pm';


my %options=();
getopt("dt",\%options);


### READ OPTIONS FOR ANALYS FROM 'options.txt'
my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };

my $post_alignment_filter = $options_hash{post_alignment_filter} || "";
my $post_alignment_threshold = $options_hash{post_alignment_threshold} || "50";
my $root_directory = $options_hash{root_directory} || die "Please specify root directory in options.txt\n";
my $perform_aliscore = $options_hash{aliscore} || "n";
my $input_directory = $options_hash{post_process_directory} || die "Please choose input directory for post processing analysis\n";
$input_directory = $options_hash{root_directory}."/".$options_hash{project_name}."/".$input_directory;
my @directories_to_analyse = <$input_directory/*>;
if(scalar(@directories_to_analyse) == 0){
	die "No Folder in directory\nMake sure the right folder is selected\n";
}

my $post_processing_method = $options_hash{post_method} || die "Please choose a post processing method\n";

MAIN: {
	if($post_processing_method eq "g"){
	#	&perform_gblocks();
	}
	
	if($post_processing_method eq "n"){
	#	&perform_noisy();
	}
	
	if($post_alignment_filter eq "t"){
	
	}
	if($perform_aliscore eq "y"){
		&perform_aliscore();
	}

	Iostuff::print_start();
	print "Analysis finished\n";
}

####################
####################
sub perform_noisy{
foreach(@directories_to_analyse){
	my $current_directory = $_;
	my $kog = basename($current_directory);
	
	print "\nAccessing directory ".$kog."\n";
### NAME OF ALIGNMENT	
	my $alignment_file = $current_directory."/".$kog."_final.fasta";
	Iostuff::check_fasta_format($alignment_file);
### NAME OF Noisy ALIGNMENT
	my $alignment_gblocked_file = $current_directory."/".$kog."_final_noisy.fasta";
#	my $seq_type = "--seqtype P";
	my $arg = qq( $alignment_file);
	
#### GBLOCK Alignment
#	print("noisy $arg");

	system("noisy $arg");

#### RENAME GBLOCKED Alignment
#	print("\nmv  $kog* $current_directory/");
	system("mv $kog* $current_directory/");
#next;
	if($post_alignment_filter eq "t"){
		
		    my $in  = Bio::SeqIO->new(-file => "inputfilename" ,
                           -format => 'Fasta');
		    my $out = Bio::SeqIO->new(-file => ">outputfilename" ,
                           -format => 'Fasta');

		    while ( my $seq = $in->next_seq() ) {
			    $out->write_seq($seq);
		    }
		
	}

#	rename($alignment_file."-gb", $alignment_gblocked_file);
	}
	return 1;

}

####################
sub perform_gblocks{
####################	
foreach(@directories_to_analyse){
	my $current_directory = $_;
	my $kog = basename($current_directory);

	print "\nAccessing directory ".$kog."\n";
### NAME OF ALIGNMENT	
	my $alignment_file = $current_directory."/".$kog."_final.fasta";
	Iostuff::check_fasta_format($alignment_file);
### NAME OF GBLOCKED ALIGNMENT
	my $alignment_gblocked_file = $current_directory."/".$kog."_final_gblocked.fasta";
### COUNT NUMBER OF SEQUENCES IN FILE  	 
	chop( my $number_of_sequences = `grep -c ">" $alignment_file` ); 
#	print $number_of_sequences." present\n";
#### GBLOCKS PARAMETER   
	my $b1 = ($options_hash{gblocks_b1}) ? $options_hash{gblocks_b1}  : (int($number_of_sequences/2)) +1;
	my $b2 = ($options_hash{gblocks_b2}) ? $options_hash{gblocks_b2} : (int($number_of_sequences * 0.65));
	$b2 =2 if $number_of_sequences ==2;
	my $b3 = ($options_hash{gblocks_b3}) ? $options_hash{gblocks_b3} :  q(10);
	my $b4 = ($options_hash{gblocks_b4}) ? $options_hash{gblocks_b4} : q(5);
	my $b5 = ($options_hash{gblocks_b5}) ? $options_hash{gblocks_b5} : q(a);
	my $arg = qq( $alignment_file -t=p -b1=$b1 -b2=$b2 -b3=$b3 -b4=$b4 -b5=$b5 -p=n);
#### GBLOCK Alignment
#	print("Gblocks $arg");
#next;
	system("Gblocks $arg");

#### RENAME GBLOCKED Alignment
	rename($alignment_file."-gb", $alignment_gblocked_file);
	unless(!-e $alignment_gblocked_file){
		if($post_alignment_filter eq "t"){
			my $alignment_gblocked_filtered_file = $current_directory."/".$kog."_final_gblocked_filtered.fasta";
			my $annotation_file = $current_directory."/LOG_rejected_sequences.txt";
			my $in  = Bio::SeqIO->new(-file => $alignment_gblocked_file ,
			           -format => 'Fasta');
			my $out = Bio::SeqIO->new(-file => ">$alignment_gblocked_filtered_file" ,
                           -format => 'Fasta');

			while ( my $seq = $in->next_seq() ) {
			    ## Avoid too short sequences
				next if $seq->length < 10;
			    ## count gaps in sequence
				my @gaps = split(/-/, $seq->seq);
			    #print "length is ".$seq->length." with ".scalar(@gaps)." gaps and content ";
			    
				next if $seq->length == 0;
				my $character_content = 100 - int(100 * scalar(@gaps)) / $seq->length;
			  #  print $character_content."\n";
				if($character_content >= $post_alignment_threshold){
					$out->write_seq($seq);
				}
				else{
					my $annotation_arg = "Id: ".$seq->id."\tcharacter_content: $character_content \%";
				    	open  my $annotation_fh, '>>', $annotation_file or warn "Could not open: $!";
					print {$annotation_fh} $annotation_arg."\n" or warn "Could not write: $!";
					close $annotation_fh or warn "Could not close: $!";
				}
			}
		    ### CHECK IF SEQUENCES HAVE BEEN REMOVED
			if(!-e $annotation_file){
				unlink($alignment_gblocked_filtered_file);
			}
		}
	}
	
}

	### PERFORM ALISCORE
sub perform_aliscore(){
	    my $cwd = getcwd();
	    chdir($cwd) || die "Could not change to root directory: \n";
	foreach(@directories_to_analyse){
		my $current_directory = $_;
		my $current_working_dir = $cwd;
		chdir($cwd) || die "Could not change to root directory: \n";
		my $kog = basename($current_directory);
		my $aliscore_dir = $current_directory."/aliscore";
		print "\nAccessing directory ".$kog."\n";
	### NAME OF ALIGNMENT	
		my $alignment_file = $current_directory."/".$kog."_final.fasta";
		my $alignment_file_aliscore = $aliscore_dir."/".$kog."_final.fasta";
		my $aliscore_alignment_final = $current_directory."/".$kog."_final_aliscore.fasta";
		my $alicut_script_source = "perl_scripts/alicut.pl";
		my $alicut_script_destination = $aliscore_dir."/alicut.pl";
	
		Iostuff::create_single_directory($aliscore_dir);
		### Copy alignment to aliscore directory
		copy($alignment_file,$alignment_file_aliscore) or die "Copy failed: $!";
		print "\tMasking random sequence similarities (RSS) using Aliscore....";
		my $command = "perl perl_scripts/aliscore.pl -i $alignment_file_aliscore";
		`$command`;
		print "done\n";
		my $check_aliscore = `ls $aliscore_dir/*List*`;
		chomp($check_aliscore);
		if(!$check_aliscore){
			warn "\tThere was a problem with Aliscore.Skipping..\n";
			next;
		}
		#print "$alicut_script_source\t$alicut_script_destination\n";
		copy($alicut_script_source,$alicut_script_destination) or die "Copy failed: $!";
		my $current_cwd = $aliscore_dir;
		if(chdir($current_cwd)){
			my $t = `pwd`;
			#print "\tin $t\n";
			#print "cat \"\\n\" | perl alicut.pl";
			print "\tRemoving RSS from alignment....done\n";
			`cat \"\\n\" | perl alicut.pl`;
			#print "done\n";
		}
		else{
			warn "\tCould not change directory: $!\n\t\tCould not remove masked RRS from alignment\n";
		}
#		print "cd $current_working_dir\n";
		chdir($current_working_dir) || warn "Could not go back to working directory\n";
		`cp $aliscore_dir/ALICUT_*.fasta $aliscore_alignment_final`;
		`rm -rf $aliscore_dir`;
		#exit;
	}
}

	return 1;
}





=head1 NAME

gblock_alignments.pl - Computes conserved blocks by using Gblocks

=head1 SYNOPSIS


  start_blast.pl
	[-t Configuration file]
	

=head1 DESCRIPTION

Computes conserved blocks useful for phylogenetic analysis using the Software Gblocks

=head2 An Example 

C<perl gblock_alignment.pl -t config.txt>

Read 'config.txt' and performs the gblocking of alignments according to the parameters

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
