	#!/usr/bin/perl 
########################################################################
# Script name  :    select_best_ortholog.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <Fab.Schreiber@gmail.com>
# 
# 
# SELECTS SEQUENCES MOST PROBABLE BEING ORTHOLOG FROM SET OF PROTEIN SEQUENCES
# 
########################################################################
## INCLUDES
####################
use Getopt::Std;
use Bio::Perl;
use Bio::SearchIO;
use Data::Dumper;
use File::Basename;
use File::Temp qw/ tempfile /;
use strict;

require 'lib/Blaststuff.pm';
require 'lib/Iostuff2.pm';

############## PARAMETERS START ##############
my %options=();
getopt("domxat",\%options);

my $root_folder = "/c1/scratch/fabian/projects/OrthoSelect_geobio/evaluation";
my $input_directory = $options{d} || die "please specify input directory for analysis\n";
my $output_directory = $options{o} || die "please specify output directory for analysis\n";
my $alignment_method = $options{m} || die "please specify alignment method for analysis\n";
my $translated_sequences_folder = $options{x} || die "please specify taxon database for analysis\n";
my $different_alignment_method = $options{a};
if($alignment_method eq "o" and $different_alignment_method eq ""){
	die "Please specify usage of the other alignment method\n";
}
my $distance_matrix_method = $options{t} || die "please specify method for distance matrix calculation\n";

my $orthologous_group = basename($input_directory);

### VARIABLES FOR TEMPORARY FILES
#  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($concated_nucleotide_sequences_fh,$concated_nucleotide_sequences) =  tempfile();
	die "Could not create temporary filehandle: $! \n" if !(-e $concated_nucleotide_sequences);
	$concated_nucleotide_sequences = $output_directory."/".$orthologous_group."/".$orthologous_group."_nucl_hits.fasta";
#  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($concated_protein_sequences_fh,$concated_protein_sequences) =  tempfile();
	die "Could not create temporary filehandle: $! \n" if !(-e $concated_protein_sequences);
	$concated_protein_sequences = $output_directory."/".$orthologous_group."/".$orthologous_group."_prot_hits.fasta";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($concated_protein_sequences_cleaned_fh,$concated_protein_sequences_cleaned) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $concated_protein_sequences_cleaned);
#	my $concated_protein_sequences_cleaned = $output_directory."/".$orthologous_group."/".$orthologous_group."_prot_hits.fasta_cleaned";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($concated_nucleotide_sequences_reduced_fh,$concated_nucleotide_sequences_reduced) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $concated_nucleotide_sequences_reduced);
#	my $concated_nucleotide_sequences = $output_directory."/".$orthologous_group."/".$orthologous_group."_nucl_hits.fasta";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($concated_estscan_protein_sequences_fh,$concated_estscan_protein_sequences) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $concated_estscan_protein_sequences);
#	my $concated_estscan_protein_sequences = "$output_directory/".$orthologous_group."/".$orthologous_group."_prot_hits_estscan.fasta";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($tmp_alignment_best_hits_fh,$tmp_alignment_best_hits) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $tmp_alignment_best_hits);
	#my $tmp_alignment_best_hits =  "$output_directory/".$orthologous_group."/".$orthologous_group."_best_hmm_hits.fasta";
##  DISTANCE MATRIX 
	my ($fasta_file_dst_fh,$fasta_file_dst) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $fasta_file_dst);
##  ALIGNMENT OF ALL HITS 
	my $first_alignment_file = $output_directory."/".$orthologous_group."/".$orthologous_group."_prot_hits_aligned.fasta";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($fasta_file_reduced_fh,$fasta_file_reduced) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $fasta_file_reduced);
##  TEMPORARY FILE
	my ($tmp_file_fh,$tmp_file) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $tmp_file);
##  FILE CONTAINING ALIGNED BEST SEQUENCE / TAXON
	my $final_alignment_file = $output_directory."/".$orthologous_group."/".$orthologous_group."_final.fasta";
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($final_alignment_file_hmm_fh,$final_alignment_file_hmm) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $final_alignment_file_hmm);
##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
	my ($final_alignment_file_estscan_fh,$final_alignment_file_estscan) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $final_alignment_file_estscan);

	my $final_alignment_file_phylip = $output_directory."/".$orthologous_group."/".$orthologous_group."_final.phy";
   
## CONTAINER FOR BEST SEQUENCE / TAXON 
   my @best_seqs=();

############## PARAMETERS END ##############

### CREATE OUTPUT DIRECTORY
Iostuff::create_single_directory($output_directory);

####################
# MAIN
####################

MAIN: {
	print "=" x 64, "\n";
	print "\tAccessing orthologous group ".$orthologous_group."\n";
	print "=" x 64, "\n";
	my $protein_data = 0;
 ### EXIT IF ONLY ONE SEQUENCE EXISTS
	(my $number_of_seqs = `grep ">"  $input_directory/*_prot_hits*.fasta | wc -l`) =~ s/\n//g;
	print "\tThere exist $number_of_seqs sequences for this orthologous group\n";
	if($number_of_seqs < 2){
		die "\t\tGroup contains less than 2 sequences. Skipped\n" if $number_of_seqs < 2;
	}
   ### CREATE SUBDIRECTORY 
	Iostuff::create_single_directory($output_directory."/".$orthologous_group);
 ### CONCATENATE ALL SINGLE NUCLEOTIDE AND PROTEIN SEQUENCES
	system("cat $input_directory/*_prot_hits*.fasta > $concated_protein_sequences");
	print "$input_directory/*_nucl_hits*.fasta\n";
	my @nucleotide_files = `cat $input_directory/*_nucl_hits*.fasta`;
	if(!@nucleotide_files){
		$protein_data = 1;
		print "\nonly protein data...\n";
	}
	else{
		system("cat $input_directory/*_nucl_hits*.fasta > $concated_nucleotide_sequences");
	}   
	### TEMPORARY ACTION
	open(my $I,"<","$concated_protein_sequences") || die "Could not open \n";
	open(my $O,">",$tmp_file) || die "Could not open \n";
	while(<$I>){
		if(/^>/){
			my $tp = $_;
			$tp =~ s/\|/_/g;
			print $O $tp;
		}
		else{
			print $O $_;
		}
	}
	close($I);
	close($O);
	`mv $tmp_file $concated_protein_sequences`;
   
	if(!$protein_data){
		open($I,"<","$concated_nucleotide_sequences") || die "Could not open \n";
		open($O,">",$tmp_file) || die "Could not open \n";
   
		while(<$I>){
			if(/^>/){
				my $tp = $_;
				$tp =~ s/\|/_/g;
				print $O $tp;
			}
			else{
				print $O $_;
		       }
		};
		close($I);
		close($O);
		`mv $tmp_file $concated_nucleotide_sequences`;
	}
	if(-e $final_alignment_file){
		die "\t\tAlignment file already exists.Nothing to be done.Exiting..\n";
	}
#if(-e $final_alignment_file){   goto FRAMEDP;}
### REMOVE DUPLICATE AND TOO SIMILAR SEQUENCES
	#unless (-e $concated_protein_sequences_cleaned){
#		print "\t\tRemoving of duplicate entries\n";
#		system("perl perl_scripts/remove_duplicates.pl -i $concated_protein_sequences -o $concated_protein_sequences_cleaned");
#		system("cp $concated_protein_sequences_cleaned $concated_protein_sequences");
	#}	
   
      ##REPLACE * by X in ALIGNMENTS
     # `perl -i -pe 's/\*/X/g' $concated_protein_sequences`;
	print "\tCalculating Distance matrix\n";
	if($distance_matrix_method eq "a"){  
		$fasta_file_dst = Iostuff::calculate_clustalw_distance($concated_protein_sequences);
	}
	else{
		&calculate_pairwise_distances($concated_protein_sequences, $fasta_file_dst);
	}
#### COMPUTE DISTANCE MATRIX
	if(!-e $fasta_file_dst){
		die "Could not find distance matrix\n";
	}
#### Read all taxa and distances in hash

	my %taxa_distances_container_hash = %{&read_distance_matrix($fasta_file_dst)};
	
	print "\t\tSelect sequences most probable being ortholog\n";
	foreach my $current_taxon (sort keys %taxa_distances_container_hash){
		push(@best_seqs ,&select_lowest_distance(\%taxa_distances_container_hash, $current_taxon));
	}
### WRITE BEST SEQUENCES/TAXON TO FILE
	   Iostuff::write_array_of_acc_to_file(\@best_seqs, $concated_protein_sequences,$fasta_file_reduced );
	    print "\t\tCompute multiple sequence alignment\n";
	if($alignment_method eq "t"){   
	#### DO T_COFFEE ALIGNMENT

		system("nice t_coffee -in $fasta_file_reduced -outfile $final_alignment_file -output fasta -quiet");
	}
	if($alignment_method eq "m"){   	
		#### DO MUSCLE ALIGNMENT
		system("nice muscle -in $fasta_file_reduced -out $final_alignment_file -stable -quiet");
	}
	if($alignment_method eq "o"){   	
		#### DO DIFFERENT ALIGNMENT
		system("nice $different_alignment_method");
	}
	if($alignment_method eq "p"){   	
		#### DO PRANK ALIGNMENT
		$fasta_file_reduced = "test.fa";
		print "nice prank -noxml -notree -quiet -d=$fasta_file_reduced -o=$final_alignment_file\n";
		system("nice prank -noxml -notree -quiet -d=$fasta_file_reduced -o=$final_alignment_file");
		system("nice prank -convert -f 12 -d=$final_alignment_file -o=$final_alignment_file_phylip");
	}
	## PROTEIN OR NUCLEOTIDE DATA?
	if($protein_data){
		print "\texit because it is protein data\n";
		exit();
	}
	print "\t\tCreate new nucleotide file from sequences in the alignment\n";
	FRAMEDP:
		my $make_reduced_nucleotide_file_check = &make_reduced_nucleotide_file($final_alignment_file,$concated_nucleotide_sequences, $concated_nucleotide_sequences_reduced);
		if(!$make_reduced_nucleotide_file_check){
			die "Could not create reduced nucleotide file. Aborting...\n";
		}
		print "\t\tRe-translate nucleotide sequences using ESTScan with different matrices\n";
		my $make_estscan_translations_check = &make_estscan_translations($final_alignment_file,$concated_estscan_protein_sequences, $concated_nucleotide_sequences_reduced);
		if(!$make_estscan_translations_check){
			die "Could not translate sequences. Aborting...\n";
		}
		print "\t\tRemoving of duplicate entries\n";
		system("perl perl_scripts/remove_duplicates.pl -i $concated_estscan_protein_sequences -o $concated_protein_sequences_cleaned");
		system("cp $concated_protein_sequences_cleaned $concated_protein_sequences");

		
		print "\t\tCalculate distance matrix\n";
		if($distance_matrix_method eq "a"){  
			$fasta_file_dst = Iostuff::calculate_clustalw_distance($concated_protein_sequences);
		}
		else{
			&calculate_pairwise_distances($concated_protein_sequences, $fasta_file_dst);
		}
#### COMPUTE DISTANCE MATRIX
		if(!-e $fasta_file_dst){
			die "Could not find distance matrix\n";
		}
#### Read all taxa and distances in hash
		%taxa_distances_container_hash = %{&read_distance_matrix($fasta_file_dst)};
		print "\t\tSelect sequences most probable being ortholog\n";
	
		foreach my $current_taxon (sort keys %taxa_distances_container_hash){
			push(@best_seqs ,&select_lowest_distance(\%taxa_distances_container_hash, $current_taxon));
		}
	   ### WRITE BEST SEQUENCES/TAXON TO FILE
		Iostuff::write_array_of_acc_to_file(\@best_seqs, $concated_protein_sequences,$fasta_file_reduced );
		print "\t\tCompute multiple sequence alignment\n";
		if($alignment_method eq "t"){   
	#### DO T_COFFEE ALIGNMENT
			system("nice t_coffee -in $fasta_file_reduced -outfile $final_alignment_file_estscan -output fasta -quiet");
		}
		if($alignment_method eq "m"){   	
		#### DO MUSCLE ALIGNMENT
			system("nice muscle -in $fasta_file_reduced -out $final_alignment_file_estscan -stable -quiet");
		}
		if($alignment_method eq "o"){   	
		#### DO DIFFERENT ALIGNMENT
			system("nice $different_alignment_method");
		}
		if($alignment_method eq "p"){   	
		#### DO PRANK ALIGNMENT
			$fasta_file_reduced = "test.fa";
			print "nice prank -noxml -notree -quiet -d=$fasta_file_reduced -o=$final_alignment_file_estscan\n";
			system("nice prank -noxml -notree -quiet -d=$fasta_file_reduced -o=$final_alignment_file_estscan");
			system("nice prank -convert -f 12 -d=$final_alignment_file -o=$final_alignment_file_phylip");
	}
	HMMADD:
		my $existing_translations = `grep ">" $translated_sequences_folder/* | wc -l`;
		chomp($existing_translations);
		if(!$existing_translations){
			die "\tWarning: There are no translations in directory $translated_sequences_folder.\n\tTranslate all ESTs by typing:\"perl perl_scripts/prepare_analysis.pl\"\n";
		}

		my $add_taxa_per_hmm_check = &add_taxa_per_hmm("$output_directory/".$orthologous_group."/", $orthologous_group, $translated_sequences_folder,$tmp_alignment_best_hits);
		print "\t\tComputing final multiple sequence alignment\n";
	
		if($alignment_method eq "t"){   
	#### DO T_COFFEE ALIGNMENT
			system("nice t_coffee -in $tmp_alignment_best_hits -outfile $final_alignment_file_hmm -output fasta -quiet");
		}
		if($alignment_method eq "m"){   	
		#### DO MUSCLE ALIGNMENT
			system("nice muscle -in $tmp_alignment_best_hits -out $final_alignment_file_hmm -stable -quiet");
		}
		if($alignment_method eq "o"){   	
		#### DO DIFFERENT ALIGNMENT
			system("nice $different_alignment_method");
		}
		if($alignment_method eq "p"){   	
		#### DO PRANK ALIGNMENT
			$fasta_file_reduced = "test.fa";
			print "nice prank -noxml -notree -quiet -d=$tmp_alignment_best_hits -o=$final_alignment_file_hmm\n";
			system("nice prank -noxml -notree -quiet -d=$tmp_alignment_best_hits -o=$final_alignment_file_hmm");
			system("nice prank -convert -f 12 -d=$final_alignment_file -o=$final_alignment_file_phylip");
		}
	
}

####################
# END MAIN
####################


sub read_distance_matrix(){
	my $file = shift;
	my %taxa_distances_container_hash = ();
	my $no_seq_for_taxon =0;
	my $distances_converted ="";
	my ($taxon, $accession_number, $distances);
	my @splitted = ();
	
	open my $DISTANCE_MATRIX,'<',$file || die "Couldn't open '$file': $!";
	while(my $current_line = <$DISTANCE_MATRIX>){
	## SKIP FIRST LINE OF FILE
	#matches:     100	
		next if $current_line =~ /^\s{4}\d*\n/;
	#matches: Acropor_mi|GS01PD01.b1     0.442  0.554  0.518  0.434  0.446  0.524  0.584  0.460	
#		if($current_line =~ /^(\w*)\|(\w*)(\.\w\d)?\s+(.*)/){
		if($current_line =~ /^(\w*_\w*)_(\w*)(\.\w\d)?\s+(.*)/){
		
			$taxa_distances_container_hash{$taxon}{$no_seq_for_taxon} =$distances_converted unless $distances_converted eq "";
			$distances_converted ="";
			$taxon = $1;
			$accession_number = $2;
			$distances = $4;
			$distances_converted = $accession_number."|";
			$no_seq_for_taxon++;
			@splitted = split(/\s{2,11}/,$distances);
	## Seperate distances by '|'		
			foreach (@splitted){
				next if $_ eq "";
			# Skip identifier	
				next if $_ =~ /\.w/;
				$_=~ s/\s//g;
				$distances_converted .=$_."|";
			}
		}
		else{
	#matches:      0.442  0.554  0.518  0.434  0.446  0.524  0.584  0.460	
			@splitted = split(/\s{2,11}/,$current_line);
	## Seperate distances by '|'		
			foreach(@splitted){
				next if $_ eq ""; 
				$_=~ s/\s//g; 
				$distances_converted .=$_."|";
			}
		}
	}
	$taxa_distances_container_hash{$taxon}{$no_seq_for_taxon} =$distances_converted;
			
	close($DISTANCE_MATRIX)|| warn "Couldn't close '$file': $!";
    return \%taxa_distances_container_hash;
}


sub select_lowest_distance(){
	my %taxon_sequence_container_hash = %{(my $tm = shift)};
	my $current_taxon = shift;
	my @positions_for_current_taxon = keys %{$taxon_sequence_container_hash{$current_taxon}};
	my @list_of_selected_seqs =();
	my %count_best_sequence_hash =();
	# Checking distance for taxon 'taxon_i' ...
	foreach my $taxon_i (keys %taxon_sequence_container_hash){
		next if ($taxon_i eq $current_taxon);
		my $current_lowest_distance = 100; 
		my $favorite_seq_of_taxon = "";
		my $favorite_position = 0;
	# for every sequence from Taxon 'taxon_i' ...
	foreach my $seq (keys %{$taxon_sequence_container_hash{$taxon_i}}){
	# to every other sequence from taxon '$pos_i' 
		my @arr = split(/\|/, $taxon_sequence_container_hash{$taxon_i}{$seq});
			foreach my $pos_i (@positions_for_current_taxon){
			## Remember sequence with lowest distance to taxon '$pos_i'
				if ($arr[$pos_i] < $current_lowest_distance){
					$current_lowest_distance = $arr[$pos_i]; 
					$favorite_seq_of_taxon = $arr[0];
					$favorite_position = $pos_i;
				}
			}
		}
	$count_best_sequence_hash{$favorite_position}++;
	}
	my @sorted = sort { $count_best_sequence_hash{$a} < $count_best_sequence_hash{$b} } keys %count_best_sequence_hash;
	my @tmp = split(/\|/,$taxon_sequence_container_hash{$current_taxon}{$sorted[0]});
	return $current_taxon."_".$tmp[0];
}





sub get_median(){
	my $file =shift;
	my @seq_lengths_array =();
	my $seq_file = Bio::SeqIO->new(-file => $file ,
                           -format => 'Fasta');

    while ( my $seq = $seq_file->next_seq() ) {
	    my $temp_seq = $seq->seq;
	    $temp_seq =~ s/\-//g;
	
	    push(@seq_lengths_array, length($temp_seq	));
	    
    }
    @seq_lengths_array = sort { $a <=> $b } @seq_lengths_array;	
    my $mid = int((scalar(@seq_lengths_array))/2);
    
    my $median = (scalar(@seq_lengths_array) % 2) ? $seq_lengths_array[$mid] : ($seq_lengths_array[$mid] + $seq_lengths_array[$mid-1])/2;

    return $median; 
}



sub calculate_pairwise_distances(){
	my $first_alignment_file =shift;
	my $distance_file = shift;	
	my %distance_hashoh = ();	
## Temporary file	
	my ($file_handle,$temp_alignment_file_in) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $temp_alignment_file_in);
	my ($file_handle_out,$temp_alignment_file_out) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $temp_alignment_file_out);

	
	my $alignment_object = Bio::SeqIO->new(-file => $first_alignment_file ,
                           -format => 'Fasta');
	open my $DISTANCE_OUT,'>',$distance_file || die "Couldn't open '$distance_file': $!";

	(my $no_seqs = `grep ">" $first_alignment_file | wc -l`) =~ s/\s//g;;
	print {$DISTANCE_OUT} "    ".$no_seqs."\n";
	    
	    my $no_of_current_source_seq=-0;
	    my $no_of_current_dest_seq = 0;
	    my $general_counter = 1;
	    
    while ( my $seq_source = $alignment_object->next_seq() ) {
	    my $no_of_current_seq=1;
	    print "\t\tCalculate distance for Sequence ".$general_counter++."/$no_seqs\n\t\t";
	    $no_of_current_source_seq++;
	    $no_of_current_dest_seq = 0;
	    print {$DISTANCE_OUT} $seq_source->id."        ";
		my $alignment_object2 = Bio::SeqIO->new(-file => $first_alignment_file ,
                           -format => 'Fasta');
		    while ( my $seq_destination = $alignment_object2->next_seq() ) {
			    my $number_of_matches =0;
			    my $number_of_non_gaps =0;
			    $no_of_current_dest_seq++;
			 ## DISTANCE BETWEEN SAME SEQUENCES   
			    if($no_of_current_source_seq == $no_of_current_dest_seq){

				    printf {$DISTANCE_OUT} ("0.000");
				    print(".");
			   
				    if($no_of_current_seq++ % 8){
					    print {$DISTANCE_OUT} "  ";
				    }
				    else{
					    print {$DISTANCE_OUT} "\n           ";
				    }
				    next;
			    }
		my $calculated_distance = 0;
	my $copied_distance =0;
	if($no_of_current_source_seq > $no_of_current_dest_seq){
		$calculated_distance = $distance_hashoh{$no_of_current_dest_seq}{$no_of_current_source_seq};
		$copied_distance =1;
	}
	else{

			    # Save sequences in file
			    
			    my $seq_out = Bio::SeqIO->new('-file' => ">$temp_alignment_file_in",
                                       '-format' => 'Fasta');
                                       $seq_out->write_seq($seq_source, $seq_destination);
			       system("muscle -in $temp_alignment_file_in -out $temp_alignment_file_out -quiet");
			die "Could not create  alignment output file: $! \n" if !(-e $temp_alignment_file_out);

			my $pairwise_alignment_object = Bio::SeqIO->new(-file => $temp_alignment_file_out ,
                          -format => 'Fasta');
			my $seq_1 = $pairwise_alignment_object->next_seq();
			my $seq_2 = $pairwise_alignment_object->next_seq();
			next if !$seq_1;
			next if !$seq_2;
			
			    #### ACCESS NUCLEOTIDWISE
			    for(my $i=1; $i <= $seq_1->length; $i++){ 
			## Skip columns containg gap ('-')
				    next if($seq_1->subseq($i,$i) eq "-" && $seq_2->subseq($i,$i));
				    $number_of_matches++ if (($seq_1->subseq($i,$i) eq $seq_2->subseq($i,$i))); 
			    }
			$calculated_distance = 1-($number_of_matches / ($seq_source->length));
			$distance_hashoh{$no_of_current_source_seq}{$no_of_current_dest_seq} = $calculated_distance;
	}

			    printf {$DISTANCE_OUT} ("%1.3f",$calculated_distance);
			    printf  (".") if !$copied_distance;
			    #printf ("%1.3f",$calculated_distance);
			    
			    if($no_of_current_seq++ % 8){
				    print {$DISTANCE_OUT} "  ";
			    }
			    else{
				    print {$DISTANCE_OUT} "\n           ";
			    }
		    
		    }
		    print {$DISTANCE_OUT} "\n";
		    print "\n";
   }
    	close($DISTANCE_OUT)|| warn "Couldn't close '$distance_file': $!";
}


sub make_reduced_nucleotide_file(){
my $final_alignment_file = shift;
my $concated_nucleotide_sequences = shift;
my $concated_nucleotide_sequences_reduced = shift;
	### GREP HEADERS FROM FINAL
	my @final_alignment_headers = `grep ">" $final_alignment_file`;
	if(scalar(@final_alignment_headers) == 0){
		warn "Could not create alternative translations using ESTScan with different matrices\n";
	}
	else{
		my %final_alignment_header_hash = ();
		foreach(@final_alignment_headers){
			chomp(my $header = $_);
			$header =~ s/>//g;
			$header =~ s/\|/_/g;
			$header =~ s/E$//g;
			$final_alignment_header_hash{$header} = 1;
		}
	#	print "NUCL in: $concated_nucleotide_sequences\nOut:$concated_nucleotide_sequences_reduced\n";
		my $concated_nucleotide_sequences_input  = Bio::SeqIO->new(-file => $concated_nucleotide_sequences ,
                           -format => 'Fasta');
		my $concated_nucleotide_sequences_reduced_output  = Bio::SeqIO->new(-file => ">".$concated_nucleotide_sequences_reduced ,
                           -format => 'Fasta');
		while(my $current_seq = $concated_nucleotide_sequences_input->next_seq()){
			#print $current_seq->id."\n";
			if(exists $final_alignment_header_hash{$current_seq->id}){
				$concated_nucleotide_sequences_reduced_output->write_seq($current_seq);
			}
		}
		### CHECK CREATED NUCL FILE
		chomp(my $no_headers_new_nucl_file = `grep ">" $concated_nucleotide_sequences_reduced`);
		if(!$no_headers_new_nucl_file){
			warn "\tNew Nucleotide file does not contain sequences\n";
			return 0;
		}
	return 1;
	
}


sub make_estscan_translations(){
my $final_alignment_file = shift;
my $concated_estscan_protein_sequences = shift;
my $concated_nucleotide_sequences_reduced = shift;
unlink($concated_estscan_protein_sequences);
			 my @est1 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/At.smat `;
			 my @est2 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/Dm.smat `;
			 my @est3 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/Dr.smat `;
			 my @est4 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/Mm.smat `;
			 my @est5 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/Rn.smat `;
			 my @est6 =  `ESTScan  $concated_nucleotide_sequences_reduced -t - -M \$ESTSCANDIR/Hs.smat `;
			open my $OUTEST, ">", $concated_estscan_protein_sequences or die "Could not open $concated_estscan_protein_sequences $!\n";
			foreach(@est1){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."AT\n";}else{print $OUTEST $_;}}
			foreach(@est2){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."DM\n";}else{print $OUTEST $_;}}
			foreach(@est3){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."DR\n";}else{print $OUTEST $_;}}
			foreach(@est4){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."MM\n";}else{print $OUTEST $_;}}
			foreach(@est5){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."RN\n";}else{print $OUTEST $_;}}
			foreach(@est6){if(/>(\w*_\w*_\w*)/){print $OUTEST ">".$1."HS\n";}else{print $OUTEST $_;}}
		close $OUTEST or die "Could not close $concated_estscan_protein_sequences\n";
		`cat $final_alignment_file >> $concated_estscan_protein_sequences`;
	
		if(! -e $concated_estscan_protein_sequences){
			warn "There was a problem translating the sequences using ESTScan\n";
			return 0;
		}
		return 1;
}



sub add_taxa_per_hmm(){
	my $curr_folder = shift;
	my $kog = shift;
	my $translated_sequences_folder = shift;
	my $tmp_alignment_best_hits = shift;
	my ($hmmsearch_out_fh,$hmmsearch_out) =  tempfile();
	die "Could not create temporary filehandle: $! \n" if !(-e $hmmsearch_out);
	my ($hmmsearch_out_fh2,$hmmsearch_out2) =  tempfile(UNLINK => 1);
	die "Could not create temporary filehandle: $! \n" if !(-e $hmmsearch_out2);
	my %taxa_hash = %{Iostuff::read_taxa("taxa_list.txt")};
	my $final_alignment = "$curr_folder/".$kog."_final.fasta";
	my $hmm_alignment = "$curr_folder/".$kog."_final.fasta";
	my $collected_hmm_hits = "$curr_folder/".$kog."_collected_hmm_hits.fasta";
	my $out_alignment = "$curr_folder/".$kog."_final.fasta";
	my $kog_hmm = "$curr_folder/".$kog."_final.hmm";
#goto HMMPARSE;

#	unlink($collected_hmm_hits);
	unlink($tmp_alignment_best_hits);
	foreach(sort keys %taxa_hash){
		my $missing_taxon = $_;
   			### BUILD HMM
   		if(! -e $kog_hmm){
			print "\t\t\tBuilding HMM...\n";
   			`hmmbuild $kog_hmm $hmm_alignment` ;
   		}
		
   			### CHECK IF BUILD WAS SUCCESS
   		if(! -e $kog_hmm){
   			print "\t\tThere was a problem building the HMM.Skipping...\n";
   			next;
		}			
		my $taxon_database = "$translated_sequences_folder/$missing_taxon.fa";
    		### HMMSEARCH3
    		##  FILE CONTAINING BEST SEQUENCE / TAXON (not aligned) 
			print "\t\t\tSearching translated EST sequences from $missing_taxon for additional hits...\n";
#		print "\t\thmmsearch $kog_hmm $taxon_database > $hmmsearch_out\n";
		`hmmsearch $kog_hmm $taxon_database > $hmmsearch_out`;
		if(! -e $hmmsearch_out && ! -s $hmmsearch_out){
			print "\tCould not perfom hmmersearch. Call was hmmsearch $kog_hmm $taxon_database > $hmmsearch_out\n";
			return();
		}
	#	`cp $hmmsearch_out $missing_taxon.lol`;
		my @hmmsearch_output = `grep \"E =\" $hmmsearch_out`;
		if(!@hmmsearch_output){
		;
		#	print "\tno hit\n";
		}
	#	print "\tchecking output with ".scalar(@hmmsearch_output)." lines\n";
		foreach(@hmmsearch_output){
			/(.*):\sdomain.*E\s=\s(.*)/;
				my ($acc_number_to_search, $score) = ($1,$2);
				if(!  $acc_number_to_search || ! defined $score){
					print "\t\tCould not find accession number/score  in HMM output ($_)\n";
					next;
				}
		#	print "\tget record $acc_number_to_search from $taxon_database\n";
			my $record = Iostuff::getfastarecord($taxon_database,$acc_number_to_search);
			$record =~ s/\s//g;

			if($record eq ""){
				print "\t\tCould not get record ".$acc_number_to_search." from $taxon_database\n";
				next;
			}
	    ## NEW SEQ OBJECT
			my $seq = Bio::Seq->new( -seq => $record,
                                 -id  => $acc_number_to_search);
             #                    print "\twriting $record to $collected_hmm_hits\n";
			write_sequence(">>$collected_hmm_hits",'fasta',$seq);
			last;
		}
	#	open my $HMMSEARCH_OUTPUT, "<", $hmmsearch_out or die "Could not open $hmmsearch_out for HMMSEARCH\n";
	#	while(<$HMMSEARCH_OUTPUT>){
	#		print $_."\n";
	#	    next if /^#/;
	#	    next if /^$/;
	#	    if(/Query:\s+(\w+_(\w+))/){
	#		my   $query_sequence = $1;
	#		my   $species = $2;
	#		print "query seq: $query_sequence\tspecies: $species\n";
	#		my $tmp_line = <$HMMSEARCH_OUTPUT>;
	#		my $tmp_line1 = <$HMMSEARCH_OUTPUT>;
	#		my $tmp_line2 = <$HMMSEARCH_OUTPUT>;
	#		my $tmp_line3 = <$HMMSEARCH_OUTPUT>;
	#		my $tmp_line4 = <$HMMSEARCH_OUTPUT>;
	#		$tmp_line4 =~ /\s+(\d+(\.\d+e-?\d+)?)/;
	#		print "\ttmp_line4: $tmp_line4\n";
	#		my @splitted = split /\s/, $tmp_line4;
	#		my $acc_number_to_search = $splitted[-1];
	#		if(! defined $acc_number_to_search){
	#			print "\t\tCould not find accession number in HMM output\n";
	#			next;
	#		}
	#		my $score = $1;
	#		print "\tget record $acc_number_to_search from $taxon_database\n";
	#		my $record = Iostuff::getfastarecord($taxon_database,$acc_number_to_search);
	#		$record =~ s/\s//g;

	#		if($record eq ""){
	#			print "\t\tCould not get record ".$acc_number_to_search." from $taxon_database\n";
	#			next;
	#		}
	    ## NEW SEQ OBJECT
	#		my $seq = Bio::Seq->new( -seq => $record,
     #                            -id  => $acc_number_to_search);
      #                           print "\twriting $record\n";
	#		write_sequence(">>$collected_hmm_hits",'fasta',$seq);
	#	}
	 #   }
    }
    if(! -e $collected_hmm_hits || ! -s $collected_hmm_hits){
    	warn "\tProblem with collecting the hits from HMMER search detected.Aborting.....\n";
    	return();
    }
    #### KEEP ONLY THE GOOD SEQUENCES - E-Value threshold 1e-10
  #  print "\thmmparse...\n";
HMMPARSE:
	my %acc_to_keep = ();
	my $search_term = "\"\\w*_\\w*_\\w*\"";
#	print "hmmsearch  -E 1e-10 $kog_hmm $collected_hmm_hits\n";
	my @hmm_output = `hmmsearch  -E 1e-10 $kog_hmm $collected_hmm_hits`;

	foreach(@hmm_output){
			if(/(.*):\sdomain.*E\s=\s(.*)/){
				my ($acc_number_to_search, $score) = ($1,$2);
				if(!  $acc_number_to_search || ! defined $score){
					print "\t\tCould not find accession number/score  in HMM output ($_)\n";
					next;
				}
		#	print "\tget record $acc_number_to_search from $taxon_database\n";
				$acc_to_keep{$acc_number_to_search} = 1;
			}
		}
	
	
#	my @hmm_output = `hmmsearch  -E 1e-10 $kog_hmm $collected_hmm_hits`;
# 	foreach(@hmm_output){
#		print $_."\n";
#		if(/^\s+\d+.*(\w*_\w*_\w*)/){
#			print "\t\tmatches first\n";
#			my @tmp_split = split(/\s\s/, $_);
			
#			foreach(@tmp_split){
#				if(/\w*_\w*/){
#					chomp;
#					$acc_to_keep{$_} = 1;
#					print "\tmatches second\n";
#				}
#			}
#		}
#	}
#	print Dumper %acc_to_keep;
#	my $in = Bio::SearchIO->new(-format => 'hmmer',
 #                              -file   => $hmmsearch_out);
#	while( my $result = $in->next_result ) {
#		while( my $hit = $result->next_hit ) {
#			$acc_to_keep{$hit->name()} = 1;
 #       }
  #  }
#if(! -e $collected_hmm_hits && ! -s $collected_hmm_hits){
#	print "\tCould not fetch sequences from sequence pool\n";
#	return();
#}
 	my $all_hmm_hits= Bio::SeqIO->new('-file' => "$collected_hmm_hits",
                                     '-format' => 'Fasta');
    	my $best_hmm_hits= Bio::SeqIO->new('-file' => ">$tmp_alignment_best_hits",
                                     '-format' => 'Fasta');
                          #             print Dumper %acc_to_keep;
	while ( my $seq = $all_hmm_hits->next_seq() ) {
		$best_hmm_hits->write_seq($seq) if exists $acc_to_keep{$seq->id};
		my $ac = $seq->id;
		chomp(my $grep_acc = `grep $ac $collected_hmm_hits`);
	}
## REMOVE TEMPORARY FILES
	unlink($collected_hmm_hits);
	unlink($kog_hmm);
    }
}


END:{
	### REMOVE UNIMPORTANT FILES
	`mv $final_alignment_file_hmm $final_alignment_file`;
}
####################
# END
####################
