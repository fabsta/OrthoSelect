package Blaststuff;
########################################################################
# Script name  :    nph-blast.pl
#
# Date created :    August 2008
#
# Author       :    Fabian Schreiber <faschrei@gwdg.de>
# 
# This package includes methods for parsing blast reports from
# BLAST searches against different ortholog databases
#
# 
########################################################################
## INCLUDES
####################
use warnings;
use Bio::Perl;
use Data::Dumper;
#use Carp;
use File::Temp qw/ tempfile tempdir /;
use strict;


############## PARAMETERS START ##############
my $LINE_LENGTH = 64;
############## PARAMETERS END ##############

########### BLAST PARSER  FOR ORTHOMCL ###########################################
sub parse_blast_orthomcl {    
####################################################################
	### VARIABLES#{{{
   my ($blast_report, $current_taxon, $current_id,$output_directory,
	   $incrementor,$tmp,$weight,$taxa_in_db) = @_;
    my %cluster_list = %{$tmp};
    my $ret_frame     = 0;
    my $max_no_of_hits =10;
    my %existing_tax;
    my $last_frame = "";    #}}}
    my %seqobj_hash;
    my $curr_e = 0;
    my $curr_bit_score = 0;
    my $first_evalue;
    my %existing_taxa;
    my $no_hits = 0;
    my $estwise_template_acc=0;
    my $frame_for_translation=0;
    my %blasto_hash =();
    
    my $fav_ortho_group;
    my $largest_value=0;
    my $largest_bitscore=0;
    my $largest_com_score=0;
    my %best_ortho_groups = ();
    $blast_report = new Bio::SearchIO(
    -format => 'blast',
    -file   => $blast_report
    );
    
### PARSE BLAST RESULTS
    while ( my $result = $blast_report->next_result) {

        while ( my $hit = $result->next_hit) {
		 my $curr_id   = $hit->name;
		 $curr_id =~ /(\w*)/;
		 my $curr_taxa = $1;
		 $hit->description =~ /\|\s*(\w*)\s\|/;
		 my $current_ortho_group = $1;
 
	    while ( my $hsp = $hit->next_hsp ) {
		    my $frame_for_translation_query = ( $hsp->query->frame + 1 ) * $hsp->query->strand;
                    my $frame_for_translation_hit = ( $hsp->hit->frame + 1 ) * $hsp->hit->strand;
                    $curr_e = $hsp->evalue();
		     $curr_taxa = substr($curr_id,0,2);
		  unless($current_ortho_group eq "no_group"){
			    if($curr_e ==0){
				    $curr_e = 1e-200;
				  }

			       			  next if exists $blasto_hash{$current_ortho_group}{$curr_taxa}{"no"};
			       			    $blasto_hash{$current_ortho_group}{"no"}++;
			       			    $blasto_hash{$current_ortho_group}{$curr_taxa}{"no"}++;
			    	 $blasto_hash{$current_ortho_group}{"evalue"} += -log($curr_e) if $blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits;
			    
			    
			    if(not defined $blasto_hash{$current_ortho_group}{"template_acc"}){
					$curr_id =~ s/\s//g;
					    $blasto_hash{$current_ortho_group}{"template_acc"} = $curr_id if not $curr_taxa eq "Dm";
			    }
#### Remember best hits/taxon /best kog			    
			    
			    
			    $blasto_hash{$current_ortho_group}{"f_evalue"}  = $curr_e if not defined $blasto_hash{$current_ortho_group}{"f_evalue"} ;
			    $blasto_hash{$current_ortho_group}{"bit_score"} += $curr_bit_score if $blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits;
			    $blasto_hash{$current_ortho_group}{"no"} ++;
			    $blasto_hash{$current_ortho_group}{"frame"} =$frame_for_translation_query if not defined $blasto_hash{$current_ortho_group}{"frame"};
			
			    unless(exists $blasto_hash{$current_ortho_group}{$curr_taxa} || exists $blasto_hash{$current_ortho_group}{"Y"}){
					   $blasto_hash{$current_ortho_group}{"no_taxa"} ++;
					   if($curr_id =~ /Y/){
						   $blasto_hash{$current_ortho_group}{"Y"} =1;
					   }
					   else{
						   $blasto_hash{$current_ortho_group}{$curr_taxa} =1;
			   		   }
			    }
		    }
		}
	}
    }
    foreach (keys %blasto_hash){
	    my $current_ortho_group = $_;
	### SCORE
	if($blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits){
	    $blasto_hash{$current_ortho_group}{"result"} = $blasto_hash{$current_ortho_group}{"evalue"} / $blasto_hash{$current_ortho_group}{"no"};
	}
	else{
	    $blasto_hash{$current_ortho_group}{"result"} = $blasto_hash{$current_ortho_group}{"evalue"} / $max_no_of_hits;
	}
	
	
	### BITSCORE
	if($blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits){
	    $blasto_hash{$current_ortho_group}{"bitscore_result"} = $blasto_hash{$current_ortho_group}{"bit_score"} / $blasto_hash{$current_ortho_group}{"no"};
	}
	else{
	    $blasto_hash{$current_ortho_group}{"bitscore_result"} = $blasto_hash{$current_ortho_group}{"bit_score"} / $max_no_of_hits;
	}
	
	    ## WEIGHTED SCORE   
	 #   $blasto_hash{$current_ortho_group}{"weighted_result"} = $blasto_hash{$current_ortho_group}{"result"} * ($blasto_hash{$current_ortho_group}{"no_taxa"} / $taxa_in_db)**$weight;
	    
	 #  print $blasto_hash{$current_ortho_group}{"result"}."\t and weighted ".$blasto_hash{$current_ortho_group}{"weighted_result"}."\twith ".$blasto_hash{$current_ortho_group}{"no_taxa"}." hits and bitscore of ".$blasto_hash{$current_ortho_group}{"bitscore_result"}."\n";
	  
	  


	    $fav_ortho_group = $current_ortho_group if not defined $fav_ortho_group;

	### 	   
	$blasto_hash{$current_ortho_group}{"com_score"} =  $blasto_hash{$current_ortho_group}{"result"} * $blasto_hash{$current_ortho_group}{"bitscore_result"};
### UPDATE BEST GROUP INFORMATION
	    
	    if ($largest_com_score < $blasto_hash{$current_ortho_group}{"com_score"} || $blasto_hash{$current_ortho_group}{"result"} == 0){
		    $fav_ortho_group = $current_ortho_group;
		    $largest_value = $blasto_hash{$current_ortho_group}{"result"};
		    $largest_bitscore = $blasto_hash{$current_ortho_group}{"bitscore_result"};
		    $largest_com_score = $blasto_hash{$current_ortho_group}{"com_score"};
		    $frame_for_translation = $blasto_hash{$current_ortho_group}{"frame"};
		    $estwise_template_acc = $blasto_hash{$current_ortho_group}{"template_acc"};
		    
	    }
## ADD HITS ABOVE THRESHOLD above 100
	 if($blasto_hash{$current_ortho_group}{"result"} > 100){
	   $best_ortho_groups{$current_ortho_group}{"template_acc"} = $blasto_hash{$current_ortho_group}{"template_acc"};
}
 	
	    
    }
## Add Best group to array of groups
if(defined $fav_ortho_group){
		  $best_ortho_groups{$fav_ortho_group}{"template_acc"} = $blasto_hash{$fav_ortho_group}{"template_acc"};
  }
######## QUIT IF NO ORTHOLOGOUS GROUP WAS SELECTED    
    if(not defined $fav_ortho_group){return $fav_ortho_group;}
    print "\t--> ";
        foreach(keys %best_ortho_groups){
		$fav_ortho_group = $_;
		print $fav_ortho_group."\t";
######## ANNOTATION
		my $arg =  "$current_id\t$current_taxon\t$fav_ortho_group\tN\tNN\t".$blasto_hash{$fav_ortho_group}{"f_evalue"}."\t".$blasto_hash{$fav_ortho_group}{"template_acc"}."";
		$best_ortho_groups{$fav_ortho_group}{"annotation"} = $arg;
	}
	print "\n";
	return \%best_ortho_groups;
}



########### BLAST PARSER FOR KOG ###########################################
sub parse_blast_kog {    
####################################################################
	### VARIABLES#{{{
#   my ($blast_report, $current_taxon, $current_id,$output_directory,
#	   $incrementor,$cluster_list_ref,$weight,$taxa_in_db, $kog_file) = @_;
   
   my ($arg_ref) = @_;
my $blast_report = ${$arg_ref->{blast_report}};
my $current_taxon = $arg_ref->{current_taxon};
my $current_id = $arg_ref->{current_id};
my $output_directory = $arg_ref->{output_directory};
my $cluster_list_ref = $arg_ref->{cluster_list_ref};
my $incrementor = $arg_ref->{incrementor};
my $weight = $arg_ref->{weight};
my $best_ortho_groups_hashref = $arg_ref->{best_ortho_groups};
my $taxa_in_db = $arg_ref->{taxa_in_db};
my $kog_file = $arg_ref->{kog_file};
   
  	my %cluster_list = %{$cluster_list_ref};
  

    my $ret_frame     = 0;
    my $max_no_of_hits =10;
    my %existing_tax;
    my $last_frame = "";    
    my %seqobj_hash;
    my $curr_e = 0;
    my $curr_bit_score = 0;
    my $first_evalue;
    my %existing_taxa;
    my $no_hits = 0;
    my $estwise_template_acc=0;
    my $frame_for_translation=0;
    my %blasto_hash =();
    my $naja = 0;
    	my %abb_taxa_hash =(
AT=> q(1),
CE => q(1),
DM=> q(1),
HS=> q(1),
SC=> q(1),
SP=> q(1),
EC=> q(1));

    my $fav_ortho_group;
    my $largest_value=0;
    my $largest_bitscore=0;
    my $largest_com_score=0;
  #  $blast_report = new Bio::SearchIO(
  #  -format => 'blast',
  #  -file   => $blast_report
  #  );
### PARSE BLAST RESULTS
  #  while ( my $result = $blast_report->next_result) {
        while ( my $hit = $blast_report->next_hit) {
            ( my $curr_id   = $hit->name )        =~ s/\s//g;
            ( my $curr_taxa = $hit->description ) =~ s/\||\s//g;
                while ( my $hsp = $hit->next_hsp ) {
		    my $frame_for_translation_query = ( $hsp->query->frame + 1 ) * $hsp->query->strand;
                    my $frame_for_translation_hit = ( $hsp->hit->frame + 1 ) * $hsp->hit->strand;
                    $curr_e = $hsp->evalue();
                    $curr_bit_score = $hsp->bits();
				     $curr_taxa = substr($curr_id,0,2);
				     my $taemp = lc $curr_taxa;
				     unless($taemp =~ /[A-Za-z]{2}/){
						  #   print "grep $curr_id $kog_file\n";
						     my $temp_grep = `grep $curr_id $kog_file`;
						     unless(defined $temp_grep){
						     $temp_grep =~ /\s+(.*):\s+(.*)/;
						     unless(defined $2){$curr_id = $2;}
		     		}
		     }
#		     print $curr_e." evalue, $curr_id id, $curr_taxa taxa\n";
		  	if(exists $cluster_list{$curr_id}){
		## Current Identifier does not match known ones	    
			    my $current_ortho_group = $cluster_list{$curr_id}; 
			    if($curr_e ==0){
				    $curr_e = 1e-200;
				  }

			       	next if exists $blasto_hash{$current_ortho_group}{$curr_taxa}{"no"};
			       			    $blasto_hash{$current_ortho_group}{"no"}++;
			       			    $blasto_hash{$current_ortho_group}{$curr_taxa}{"no"}++;
			    	 $blasto_hash{$current_ortho_group}{"evalue"} += -log($curr_e) if $blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits;
			    if(not defined $blasto_hash{$current_ortho_group}{"template_acc"}){
					$curr_id =~ s/\s//g;
					    $blasto_hash{$current_ortho_group}{"template_acc"} = $curr_id if not $curr_taxa eq "Dm";
			    }
			    
			    $blasto_hash{$current_ortho_group}{"f_evalue"}  = $curr_e if not defined $blasto_hash{$current_ortho_group}{"f_evalue"} ;
### Sum up first 5 bit-scores			    
			    $blasto_hash{$current_ortho_group}{"bit_score"} += $curr_bit_score if $blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits;
			    
			    $blasto_hash{$current_ortho_group}{"frame"} =$frame_for_translation_query if not defined $blasto_hash{$current_ortho_group}{"frame"};
			
			    unless(exists $blasto_hash{$current_ortho_group}{$curr_taxa} || exists $blasto_hash{$current_ortho_group}{"Y"}){
					   $blasto_hash{$current_ortho_group}{"no_taxa"} ++;
					   if($curr_id =~ /Y/){
						   $blasto_hash{$current_ortho_group}{"Y"} =1;
					   }
					   else{
						   $blasto_hash{$current_ortho_group}{$curr_taxa} =1;
			   		   }
			    }
		    }
		}
	}
  #  }
### SELECT BEST HIT    
 #   print Dumper %blasto_hash;
#    exit;
    foreach (keys %blasto_hash){
	    my $current_ortho_group = $_;
	
	### SCORE
	if($blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits){
	    $blasto_hash{$current_ortho_group}{"result"} = $blasto_hash{$current_ortho_group}{"evalue"} / $blasto_hash{$current_ortho_group}{"no"};
	}
	else{
	    $blasto_hash{$current_ortho_group}{"result"} = $blasto_hash{$current_ortho_group}{"evalue"} / $max_no_of_hits;
	}
	
	
	### BITSCORE
	if($blasto_hash{$current_ortho_group}{"no"} < $max_no_of_hits){
	    $blasto_hash{$current_ortho_group}{"bitscore_result"} = $blasto_hash{$current_ortho_group}{"bit_score"} / $blasto_hash{$current_ortho_group}{"no"};
	}
	else{
	    $blasto_hash{$current_ortho_group}{"bitscore_result"} = $blasto_hash{$current_ortho_group}{"bit_score"} / $max_no_of_hits;
	}
	
	    ## WEIGHTED SCORE   
	 #   $blasto_hash{$current_ortho_group}{"weighted_result"} = $blasto_hash{$current_ortho_group}{"result"} * ($blasto_hash{$current_ortho_group}{"no_taxa"} / $taxa_in_db)**$weight;
	    
#	    print $blasto_hash{$current_ortho_group}{"result"}."\t and weighted ".$blasto_hash{$current_ortho_group}{"weighted_result"}."\twith ".$blasto_hash{$current_ortho_group}{"no_taxa"}." hits and bitscore of ".$blasto_hash{$current_ortho_group}{"bitscore_result"}."\n";
	  
	  


	    $fav_ortho_group = $current_ortho_group if not defined $fav_ortho_group;

	### 	   
	$blasto_hash{$current_ortho_group}{"com_score"} =  $blasto_hash{$current_ortho_group}{"result"} * $blasto_hash{$current_ortho_group}{"bitscore_result"};
### UPDATE BEST GROUP INFORMATION
	    
	    if ($largest_com_score < $blasto_hash{$current_ortho_group}{"com_score"} || $blasto_hash{$current_ortho_group}{"result"} == 0){
		    $fav_ortho_group = $current_ortho_group;
		    $largest_value = $blasto_hash{$current_ortho_group}{"result"};
		    $largest_bitscore = $blasto_hash{$current_ortho_group}{"bitscore_result"};
		    $largest_com_score = $blasto_hash{$current_ortho_group}{"com_score"};
		    $frame_for_translation = $blasto_hash{$current_ortho_group}{"frame"};
		    $estwise_template_acc = $blasto_hash{$current_ortho_group}{"template_acc"};
		    
	    }
## ADD HITS ABOVE THRESHOLD above 100
	 if($blasto_hash{$current_ortho_group}{"result"} > 100){
						   $$best_ortho_groups_hashref{$current_ortho_group}{"template_acc"} = $blasto_hash{$current_ortho_group}{"template_acc"};
}
    }
## Add Best group to array of groups
if(defined $fav_ortho_group){
						  $$best_ortho_groups_hashref{$fav_ortho_group}{"template_acc"} = $blasto_hash{$fav_ortho_group}{"template_acc"};
						  }

######## QUIT IF NO ORTHOLOGOUS GROUP WAS SELECTED    
    if(not defined $fav_ortho_group){
  #  		print "\tno fav best group\n";
#		    return \%best_ortho_groups;
		    return();
		    }
#    print "preferred orthologous group is ".$fav_ortho_group." with $largest_value from group $fav_ortho_group \n";
#    print "\t--> ";
    foreach(keys %{$best_ortho_groups_hashref}){
 	   $fav_ortho_group = $_;
#    	print $fav_ortho_group."\t";
######## ANNOTATION
		my $annotation_file = $output_directory."/annotations/".$current_taxon."".$incrementor.".txt";
		my $annot = `grep $fav_ortho_group $kog_file`;
	
	    if (defined $annot){
			$annot =~ /\[(.*)\]\s+KOG\d*\s+(.*)/;
			my $cog_class = $1;
			$annot = $2;
			my $arg =  "$current_id\t$current_taxon\t$fav_ortho_group\t$cog_class\t$annot\t".$blasto_hash{$fav_ortho_group}{"f_evalue"}."\t".$blasto_hash{$fav_ortho_group}{"template_acc"}."";
#print "\targ: $arg and \n";
			$$best_ortho_groups_hashref{$fav_ortho_group}{"annotation"} = $arg;
		}
		else{
			my $arg =  "$current_id\t$current_taxon\t$fav_ortho_group\tN\tNN\t".$blasto_hash{$fav_ortho_group}{"f_evalue"}."\t".$blasto_hash{$fav_ortho_group}{"template_acc"}."";
#print "\targ: $arg and \n";
			$$best_ortho_groups_hashref{$fav_ortho_group}{"annotation"} = $arg;
		}
	}
	if(keys(%{$best_ortho_groups_hashref})){
#print "returning ".keys(%{$best_ortho_groups_hashref})."\n";
		return keys(%{$best_ortho_groups_hashref});
	}
	else{
#print "returning ".keys(%{$best_ortho_groups_hashref})."\n";
		return();
	}
return();
}


1;
