package Iostuff;
########################################################################
# Script name  :    nph-blast.pl
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
#use warnings;
use Bio::Perl;
use Data::Dumper;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
#use Carp;
use File::Temp qw/ tempfile tempdir /;
our $LINE_LENGTH = 64;
#use strict;


####################################################################
sub print_start(){
####################################################################
print "=" x $LINE_LENGTH, "\n";
print "-o0o-" x 12, "\n";
print "- " x 30, "\n";
print ">>==<<==" x 8, "\n";

}


####################################################################
sub get_cog4acc(){
####################################################################
	my $acc =shift;
	my $file = shift || "kog_list";
	my $cog_cluster="";
	my $right_cog=0;
	open my $COG_LIST,'<',$file || die "Couldn't open '$file': $!";
	while(<$COG_LIST>){
		if(/\[\w+\]\s+(KOG\d*)\s.*/){
			$cog_cluster =$1;
			#print "found cluster $cog_cluster\n";
			
		}
		if(/\s+(.*):\s+$acc/){
			close($COG_LIST);
				return $cog_cluster || warn "Couldn't close '$file': $!";
		}
	}
	return undef;
	
}


####################################################################
sub get_fun_class4kog($){
####################################################################
	my $kog =shift;
	my $file = shift || "kog_list";
	my $grep  =`grep $kog $file`;
	
	return undef if not defined $grep;
	##matches
	# [N] KOG3905 Dynein light intermediate chain
	$grep =~ /\[(\w)\]\s?KOG\d*\s?.*/;
	my $fun_class = $1;
	return undef if not defined $fun_class;
	return $fun_class;
}



####################################################################
#sub read_taxa($) {
####################################################################
#    my $file = shift;
#    my %taxa_hash;
#    open my $TAXA_FILE,'<' ,$file || die "Couldn't open '$file': $!";
#    while (<$TAXA_FILE>) {
#        chop;
#        ( my $tmpvar1, my $tmpvar2 ) = split( /\s*/, $_ );
#	$tmpvar1=~ s/\s//g;
#	$tmpvar2=~ s/\s//g;
#	$taxa_hash{$tmpvar1} = $tmpvar2;
# }
#     close($TAXA_FILE) || warn "Couldn't close '$file': $!";
#    return \%taxa_hash;
#}

####################################################################
sub read_taxa($){
####################################################################
my $file = shift || "stats/db_list.txt";
my %taxa_hash=();
open my $TAXA_FILE,'<', $file || die "Couldn't open '$file': $!";
while(<$TAXA_FILE>){
	## Matches:
	# "Epydatia_sp" "Epydatia_s" 5
	/"(\w+)"\s+"(\w+)"\s?/;
	$taxa_hash{$1} = $2;
}
close($TAXA_FILE) || warn "Couldn't close '$file': $!";
return \%taxa_hash;
}


####################################################################
sub read_kyva2gb($) {
####################################################################
    my $file = shift;
    my %taxa_hash;
    open my $TAXA_FILE,'<', $file || die "Couldn't open '$file': $!";
    while (<$TAXA_FILE>) {
        #chop;
	#print $_;
	/(\w*)\s*(\w*)/;
        my $tmpvar1=$1;
	my $tmpvar2=$2;
	$tmpvar1=~ s/\s//g;
	$tmpvar2=~ s/\s//g;
	$taxa_hash{$tmpvar1} = $tmpvar2;
	#print $tmpvar1."\t".$tmpvar2."\n";
    }
     close($TAXA_FILE) || warn "Couldn't close '$file': $!";
    return \%taxa_hash;
}


####################################################################
sub read_dbs($) {
####################################################################	
    my $file = shift;
    my @dbs_array;
    open my $DB_FILE, '<',$file || die "Couldn't open '$file': $!";
    while (<$DB_FILE>) {
        chop;
        unless (/#/) {
            push( @dbs_array, $_ );
        }
    }
     close($DB_FILE) || warn "Couldn't close '$file': $!";
    return \@dbs_array;
}

####################################################################
sub read_kog_list($){
####################################################################
	my $file = shift;
	my %kog_list=();
	my $cog_cluster ="";
	open my $COG_LIST,'<',$file || die "Couldn't open '$file': $!";;
	while(<$COG_LIST>){
		#[S] KOG4851 Uncharacterized conserved protein
		if(/\[\w+\]\s+(KOG\d*)\s.*/){
			$cog_cluster =$1;
		}
		#  ath:  At4g03180
		if(/\s+(.*):\s+(.*)/){
		#if(/\s+(.*):\s+([A-Za-z]*.*)/){			
			my $acc = $2;
			$kog_list{$acc} = $cog_cluster;
		#	print "$acc -> $cog_cluster \n";
		}
	}
     close($COG_LIST) || warn "Couldn't close '$file': $!";;
return \%kog_list;
}

####################################################################
sub get_acc4cog(){
####################################################################
	my $file = shift;
	my $cog =shift;
	my $taxon = shift;
	my @cog_cluster =();
	my $right_cog_flag=0;
	my $right_cog = 0;
	open my $COG_LIST,'<',$file || die "Couldn't open '$file': $!";;
	while(<$COG_LIST>){
		if(/\[\w+\]\s+(KOG\d*)\s.*/){
			my $cog_cluster =$1;
			#print "found cluster $cog_cluster\n";
			if($right_cog){
			     close($COG_LIST) || warn "Couldn't close '$file': $!";;
			     return \@cog_cluster;
			}
			print "$cog_cluster $cog\n";
			if($cog_cluster eq $cog){
				$right_cog=1;
				next;
			}
		}
		if($right_cog){
			if(/\s+$taxon:\s+(.*)/){
				my $curr_member=$1;
				print "add $curr_member\n";
				push(@cog_cluster,$curr_member);
			}
		}
	
	}
     close($COG_LIST) || warn "Couldn't close '$file': $!";;

	return 0;
	
}

####################################################################
sub do_alignment() {
####################################################################	
    my $output_directory = shift;
    print "\nDo Protein Alignment\n";
    system("perl clustal.pl $output_directory/prot_hits.fasta\n");

    # 			system("perl clustal_test.pl $output_directoryprot_hits_blast.fasta\n");
    system(
"t_coffee $output_directory/prot_hits.fasta -outfile=$output_directory/prot_alignment_t-coffee.fasta -output=fasta -quiet"
    );
}



	
####################################################################
sub check_sponge_classes() {
####################################################################
    my %hex = (
        "Heteroc_sp" => "1",
        "Oopsaca_sp" => "1",
        "Cratero_me" => "1"
    );
    my %dem = (
        "Ephydat_sp" => "1",
        "Carteri_fo" => "1",
        "Pachydi_gl" => "1",
        "Amphime_qu" => "1",
        "Lubomir_ba" => "1"
    );
    my %cal = ( "Sycon_ra" => "1",
    "Leucett_ch" => "1");
    my %hom = (
        "Oscarel_ca" => "1",
        "Oscarel_lo" => "1"
    );
    my $dem = 0;
    my $hex = 0;
    my $hom = 0;
    my $cal = 0;

    my $taxa_strin = shift;
    my @taxa_array = split( "#", $taxa_strin );
    foreach (@taxa_array) {
        /(.*)\|/;
        ( my $id_to_check = $1 ) =~ s/\s//g;

        # print "check ".$id_to_check."\n";
        $dem++ if exists $dem{$id_to_check};
        $hex++ if exists $hex{$id_to_check};
        $cal++ if exists $cal{$id_to_check};
        $hom++ if exists $hom{$id_to_check};
##print "Dem: $dem\tHex: $hex\tCal: $cal\tHom: $hom\n";
    }
    print "Dem: $dem\tHex: $hex\tCal: $cal\tHom: $hom\n";

    ( $dem && $hex && $cal && $hom ) ? return 1 : return 0;

}

####################################################################
sub find_orthologous_group(){
####################################################################	
		my $o_file = shift;
		print $o_file;
		my %ortho_group_hash =();
		my $current_ortho_group;
		open(O_FILE,$o_file);
		while(<O_FILE>){
			#if(/\[/){print $_;}
			if(/^\[(.*)\]\s*(.*) (.*)/){
				$current_ortho_group = $2;	
				print $current_ortho_group."\n";
			}
			else{
				/.*:\s*(.*)/;
				print "\t".$1."\n";

			}
			
		}
	}

####################################################################
sub translate_est_estscan($){
####################################################################
### CREATE TEMPORARY FILE
	    my ($file_handle,$temp_bash_script) =  tempfile(UNLINK => 1);
	    unless($temp_bash_script){
		    print "Could not create temporary filehandle: $! \n";
	    }	
	    
	my $ESTScan_nucl_file = shift;
	my $ESTScan_prot_file = $ESTScan_nucl_file."_prot";
	my $ESTScan_out_file = $ESTScan_nucl_file."_out";
	system("nice ESTScan -m -100 -d -50 -i -50 -M /gobics/home/fabian/.bin/bioinf/estscan-3.0.2/Hs.smat -p 4 -N 0 -w 60  -S -s 1 $ESTScan_nucl_file -t $ESTScan_prot_file -o $ESTScan_out_file");
	unlink($ESTScan_out_file) || warn "couldnt unlink file $ESTScan_out_file \n";
	

   return $ESTScan_prot_file;
}

####################################################################
sub translate_est_genewise($){
####################################################################
	my $ESTScan_nucl_file = shift;
	my $ESTScan_prot_file = $ESTScan_nucl_file."_prot";
	my $ESTScan_out_file = $ESTScan_nucl_file."_out";
	system("nice ESTScan -m -100 -d -50 -i -50 -M /gobics/home/fabian/.bin/bioinf/estscan-3.0.2/Hs.smat -p 4 -N 0 -w 60  -S -s 1 $ESTScan_nucl_file -t $ESTScan_prot_file -o $ESTScan_out_file");
	unlink($ESTScan_out_file) || warn "couldnt unlink file $ESTScan_out_file \n";
	

   return $ESTScan_prot_file;
}


####################################################################
sub read_options($){
####################################################################
my $file = shift || "options.txt";
my %options_hash=();

open my $OPTIONS,'<',$file || die "Couldn't open '$file': $!";
while(<$OPTIONS>){
	## SKIP COMMENTS
	next if /#/;
	next if /^$/;
	/\s*(\w*)\s?=\s*"(.*)"\s?/;
	next if(!defined $1 || !defined $2);
	$options_hash{$1} = $2;
}
close($OPTIONS)|| warn "Couldn't close '$file': $!";
#print Dumper(%options_hash);
return \%options_hash;
}

####################################################################
sub read_required_taxa($){
####################################################################
my $file = shift || "required_taxa_list.txt";
my %taxa_hash=();
	open my $TAXA,'<',$file || die "Couldn't open '$file': $!";
	while(<$TAXA>){
		## Matches:
		# Epydatia_sp Epydatia_s 5
		/"(\w+)"\s+"(\w+)"\s?/;
		my $key = $1;
		my $value = $2;
		$key=~ s/\s//g;
		$value=~ s/\s//g;
		$taxa_hash{$key} = $value;
	}
	close($TAXA)|| warn "Couldn't close '$file': $!";
	return \%taxa_hash;
}

####################################################################
sub read_required_taxa_monophylum($){
####################################################################
my $file = shift || "required_taxa_list.txt";
my %taxa_hash=();
	open my $TAXA,'<',$file || die "Couldn't open '$file': $!";
	while(<$TAXA>){
		next if /^\n/;
		## Matches:
		#"Monophylum1" = "Acropor_mi","Allomyc_ma","Amphime_qu"
		/"(\w+)"\s?=\s?(.*)/;
		my $monophylum = $1;
		my @monophylum_taxa = split(/,/,$2);
		foreach(@monophylum_taxa){
			/"(\w*)"/;
		#	print " $monophylum  --> $1 \n";
			$taxa_hash{$monophylum}{$1} = 1;
		}
	}
	close($TAXA)|| warn "Couldn't close '$file': $!";
	return \%taxa_hash;
}


####################################################################
sub read_fun_class_file($){
####################################################################
	my $file = shift || "fun.txt";
		if(! -e $file){
		die "File containing functional classification does not exist\n";
	}

	my %func_class_hash=();

	open my $FUNC_CLASS,'<',$file || die "Couldn't open '$file': $!";
	while(<$FUNC_CLASS>){
		## Matches:
		# [J] Translation, ribosomal structure and biogenesis
		$func_class_hash{$1} = 1 if /\s?\[(\w*)\]\s?(.*)/;
	}
	close($FUNC_CLASS)|| warn "Couldn't close '$file': $!";
	return \%func_class_hash;
}




####################################################################
sub create_single_directory($){
####################################################################
my $name = shift;
if ( !-e $name ) {
	mkdir($name, 0777 )
	or die "could not create '$name': $!";
    }
    return 1;
}



####################################################################
sub create_project_directories($){
####################################################################
my $project_name = shift;

if ( !-e $project_name ) {
	mkdir($project_name, 0777 )
    or warn q(could not create hit-directory);
    }
if ( !-e "$project_name/annotations" ) {
			    mkdir("$project_name/annotations", 0777 )
			    or die q(could not create annotation-directory);
		    }
if ( !-e "$project_name/basis_hits" ) {
			    mkdir("$project_name/basis_hits", 0777 )
			    or die q(could not create basis_hits-directory);
		    }
#if ( !-e "$project_name/stats" ) {
#			    mkdir("$project_name/stats", 0777 )
#			    or die q(could not create stats-directory);
#		    }

#if ( !-e "$project_name/logs" ) {
#			    mkdir("$project_name/logs", 0777 )
#			    or die q(could not create logs-directory);
#		    }
	return 1;
}

####################################################################
sub copy_folder($$){
####################################################################
my $source_folder = shift;
my $destination_folder =shift;
system("cp -r $source_folder $destination_folder");

return 1;
}

####################################################################
sub write_array_of_acc_to_file($$$){
####################################################################
	my @best_seqs = @{(my $tm = shift)};
	my $file = shift;
	my $output_file = shift;
	
	    foreach my $key (@best_seqs){
		    $key =~ /\w*_(\w*)/;
		    $key =~ s/\|/_/;
		    
		    my $tst =  $1;
		    $tst =~ s/\s//g;
		#    print "$tst from $file \n";
		    my $record = getfastarecord($file, $tst);
	#	    print "$record was $record\n";
		    next if $record eq 0;
		    $record =~ s/\-|\s//g;
		#    print "saving with $key\n";
		    my $seqobj = Bio::Seq->new(
			   -seq        => $record,
                           -display_id => $key,
                           -alphabet   => 'dna'
                    );
		   write_sequence( ">>$output_file",
                    'fasta', $seqobj);
	    }
    return 1;
}


####################################################################
sub getfastarecord($ $){
####################################################################
my $file = shift;
my $query = shift;
my $flag=0;
my $entry ="";

open my $INPUT,"<",$file || die "Couldn't open '$file': $!";
while(<$INPUT>){
	if(/>/){
		if($flag ==1){
			close($INPUT)|| warn "Couldn't close '$file': $!";
			return $entry;
		}
			if(/$query/){
		#		print "$query matches $_ \n";
				$flag =1;
# 				$entry = ">$query\n";
			}
			else{$flag =0;}
	}
	else{
		if($flag == 1){
			$entry .=$_;
		}
	}
	}
	close($INPUT)|| warn "Couldn't close '$file': $!";
	return $entry if $flag;
	return undef;
}


####################################################################
sub write_bash_script($$$){
####################################################################
my $file = shift;
my $arg = shift;
my $identifier = shift || "";

open my $BASH_SCRIPT,">",$file || die "Couldn't open '$file': $!";


print {$BASH_SCRIPT} "#!/bin/bash\n";
print {$BASH_SCRIPT} "#\$ -S /bin/bash\n";
print {$BASH_SCRIPT} "#\$ -cwd\n";
#print {$BASH_SCRIPT} "#\$ -M your_email\@internet.com\n";
#print {$BASH_SCRIPT} "#\$ -m e\n";
print {$BASH_SCRIPT} "#\$ -N ".$identifier."\n";
print {$BASH_SCRIPT} "#\$ -V\n";
print {$BASH_SCRIPT} $arg."\n";

close($BASH_SCRIPT)|| warn "Couldn't close '$file': $!";
}



####################################################################
sub read_orthomcl_list($){
####################################################################
	my $file = shift;
	my %orthomcl_clusters =();
	open my $ORTHOMCL_LIST,'<',$file || die "Couldn't open '$file': $!";
	while(my $line = <$ORTHOMCL_LIST>){
		# matches OG2_94726: tni|GSTENP00026403001 fru|NEWSINFRUP00000154752
		$line =~ m/(OG\d_\d*):\s?(.*)/;
		my $orthomcl_group = $1;
		my @group_members = split(/\s/, $2);
		foreach my $group_member (@group_members){
			$orthomcl_clusters{$group_member} = $1;
		
		}
	#	print Dumper(%orthomcl_clusters);
	#	exit;
	}
     close($ORTHOMCL_LIST) || warn "Couldn't close '$file': $!";;
return \%orthomcl_clusters;
}



####################################################################
sub check_taxa_list_syntax($){
####################################################################
	my $file = shift;
	if(! -e $file){
		die "File containing relevant taxa (required_taxa_list.txt) does not exist\n";
	}
	my $number_of_lines = `cat $file | wc -l`;
	(my $taxa_string = `cat $file`) =~ s/\n//g;
	# "Apis_mellifera" "Apis_melli" p
	my $count_apostrophe = ($taxa_string =~ tr/"//);
	if((4 * $number_of_lines) != ($count_apostrophe)){
		die "Wrong format!!!\nCheck Syntax of file $file \n";
	}
	#	print "Syntax Check of File $file successful\n";
		return 1;
}


####################################################################
sub check_fasta_format($){
####################################################################
	my $fasta_file = shift;
	
	my $number_of_fasta_blocks_in_format =`perl -pe '/>.*\n(\w|\n)+/' $fasta_file | grep ">" -c`; 
	my $number_of_fasta_headers = `grep ">" $fasta_file -c`;
	
	#print "\tChecking Correct Fasta Format of $fasta_file....";

	
	if($number_of_fasta_blocks_in_format != $number_of_fasta_headers){
		die "Wrong format!!!\nCheck Syntax of file $fasta_file \n";
	}
	#	print " successful\n";
		return 1;
}


####################################################################
sub check_directory_exists{
####################################################################
	my $directory = shift;
	unless(-e $directory && -w $directory){
		warn "Directory $directory does not exists or is not accessible\n";
		return 0;
	}
	return 1;
}


####################################################################
sub check_required_programs{	
####################################################################
	my @progs = @{my $tmp = shift};	
	my $add_vars = shift;
	#my @progs=@_;
	my @errors;
	my $found;
	if($add_vars){
	#	print "$add_vars is add_vars \n";
	unless (exists $ENV{ESTSCANDIR})	{
		die "\nTo ensure the smooth running of the pipeline you need to set the environmental variable ESTSCANDIR\n".
		      "This should be placed in the .cshrc / .bashrc file\n".
			  "For further queries please consult the user guide\n";
	}
	unless (exists $ENV{WISECONFIGDIR}){
		die "\nTo ensure the smooth running of the pipeline you need to set the environmental variable WISECONFIGDIR\n".
		      "This should be placed in the .cshrc / .bashrc file\n".
			  "For further queries please consult the user guide\n";
	}
	}
	#else	{	  
	#	$ENV{PERL5LIB} .= ":".$ENV{ESTSCANDIR}."/MkTables";
	#	$ENV{PATH} .= ":".$ENV{ESTSCANDIR};
	#	$ENV{PATH} .= ":".$ENV{ESTSCANDIR}."/MkTables";
	#}
	
	
	REQ_PROG: foreach my $prog (@progs)	{
		print "Checking $prog ...";
		foreach my $path (split(":","$ENV{'PATH'}"))	{ 
			if (-x "$path/$prog")	{
				#$ENV{PERL5LIB} .= ":$path ";
				#$ENV{PATH} .= ":$path";
				print "successful\n";
				next REQ_PROG;
			
			}
			elsif (-x "$path/MkTables/$prog")	{
				#$ENV{PERL5LIB} .= ":$path/MkTables/";
				#$ENV{PATH} .= ":$path/MkTables";
				print "successful\n";
				next REQ_PROG;
			}
			else	{
				$found=0;
			}
		}
		push @errors, $prog if ($found==0);	
	} 
	return \@errors;
}





####################################################################
sub download_kog_database{
####################################################################
	my $install_directory = shift;
	my $database_name_tmp = shift;
	my $kog_sequence_file = qq($install_directory/kyva);
	my $kog_annotation_file = qq($install_directory/fun.txt);
	my $kog_classification_file = qq($install_directory/KOG_list.txt);
	
	my $database_name = qq($install_directory/$database_name_tmp);
	
	
	print "=" x $LINE_LENGTH, "\n";
	print "\tDownloading and setting up Database $database_name....\n";
	print "=" x $LINE_LENGTH, "\n";
	## Download database in fasta format
	 system("wget ftp://ftp.ncbi.nih.gov/pub/COG/KOG/kyva -O $kog_sequence_file \n");
	 	die "Couldnt download sequence file for KOG\nCheck internet connection" if (!-e $kog_sequence_file);
	 ## Download annotation file
	 system("wget ftp://ftp.ncbi.nih.gov/pub/COG/KOG/fun.txt -O $kog_annotation_file \n");
	 	die "Couldnt download annotation file for KOG\nCheck internet connection" if (!-e $kog_annotation_file);
	 ## Download kog classification file
	 system("wget ftp://ftp.ncbi.nih.gov/pub/COG/KOG/kog -O $kog_classification_file \n");
	 	die "Couldnt download classification file for KOG\nCheck internet connection" if (!-e $kog_classification_file);
	
	
	 ## Format fasta file to blastable database
	print "=" x $LINE_LENGTH, "\n";
	print "Formatting Database $database_name....\n";
	print "=" x $LINE_LENGTH, "\n";
	 system("formatdb -i $kog_sequence_file -n $database_name");
		 die "Couldnt format database" if (!-e "$database_name.phr");
	
	
	return 1;
}

####################################################################
sub download_orthomcl_database{
####################################################################
	my $install_directory = shift;
	my $database_name_tmp = shift;
	my $orthomcl_sequence_file_gzipped = qq($install_directory/seqs_orthomcl-2.fasta.gz);
	my $orthomcl_sequence_file = qq($install_directory/seqs_orthomcl-2.fasta);
	my $orthomcl_annotation_file_gzipped = qq($install_directory/groups_orthomcl-2.txt.gz);
	my $orthomcl_annotation_file_gzipped2 = qq($install_directory/groups_orthomcl-2.txt);
	my $orthomcl_annotation_file = qq($install_directory/OrthoMCL_list.txt);
	my $database_name = qq($install_directory/$database_name_tmp);
	
	print "=" x $LINE_LENGTH, "\n";
	print "\tDownloading and setting up Database $database_name....\n";
	print "=" x $LINE_LENGTH, "\n";
	## Download database in fasta format
	 system("wget http://www.orthomcl.org/common/downloads/2/seqs_orthomcl-2.fasta.gz -O $orthomcl_sequence_file_gzipped \n");
	 	die "Couldnt download sequence file for KOG\nCheck internet connection" if (!-e $orthomcl_sequence_file_gzipped);
	##Extract database
	 system("gunzip $orthomcl_sequence_file_gzipped ");
	## Download annotation file
	 system("wget http://www.orthomcl.org/common/downloads/2/groups_orthomcl-2.txt.gz -O $orthomcl_annotation_file_gzipped \n");
	 	die "Couldnt download annotation file for OrthoMCL\nCheck internet connection" if (!-e $orthomcl_annotation_file_gzipped);
	
	system("gunzip $orthomcl_annotation_file_gzipped ");
	system("mv $orthomcl_annotation_file_gzipped2  $orthomcl_annotation_file");
	
	 ## Format fasta file to blastable database
	print "=" x $LINE_LENGTH, "\n";
	print "Formatting Database $database_name....\n";
	print "=" x $LINE_LENGTH, "\n";
	 system("formatdb -i $orthomcl_sequence_file -n $database_name");
		 die "Couldnt format database" if (!-e "$database_name.phr");
	
	 ## Check if Database is blastable
#	 system("formatdb -i $orthomcl_sequence_file -n $database_name");
	 
	
	
	##Perform short database search check
	
	return 1;
}

####################################################################
sub check_blast_db{
####################################################################
	my $database_name = shift;

	print "=" x $LINE_LENGTH, "\n";
	print "Testing Database $database_name....\n";
	print "=" x $LINE_LENGTH, "\n";
	 ## Check if Database is blastable
#	 print("blastall -p blastx -i test_seq.fa -d $database_name -o test.out");
	 system("blastall -p blastx -i test_seq.fa -d $database_name -o test.out");

	 if(! -e "test.out"){
	 	die "Problem testing database...\nExiting\n";
	 }
	 else{
		 return 1;
	 }
}

####################################################################
sub print_options{
####################################################################
	my %options = %{my $tmp = shift};
	my $selected_options = shift;
	my @selected = split(/,/,$selected_options);
	print "=" x $LINE_LENGTH, "\n";
	print "\tYou selected the following options\n";
	print "=" x $LINE_LENGTH, "\n";

	foreach my $option (@selected){
		print "#\t $option = ".$options{$option}."\n";
	}
	print "\n\nAre the options correct? Type \"y\" for yes:\n	";
	my $user_input = <STDIN>;
	chop($user_input);
	unless($user_input eq "\n" || $user_input eq "y" || $user_input eq "Y"){
		die "\tProgram aborted....\n";
	}


}

####################################################################
sub calculate_clustalw_distance{
####################################################################
my $alignment_file = shift;
(my $clustalw_distance_file = $alignment_file) =~ s/\.fasta/\.dst/;
(my $clustalw_tree_file = $alignment_file) =~ s/\.fasta/\.ph/;
#(my $clustalw_alignment_file = $alignment_file) =~ s/\.fasta/\.aln/;

my $clustalw_options = " -infile=$alignment_file -tree -outputtree=dist -outorder=input ";

unless(-e $alignment_file){
	die "Alignment file $alignment_file does not exist\n";
}

## CALCULATE DISTANCE MATRIX
system("clustalw  $clustalw_options");

## REMOVE UNWANTED FILES
system("rm $clustalw_tree_file");

return $clustalw_distance_file;
}


####################################################################
sub check_shell{
####################################################################
my $shell_command = `ps -p $$`;
if($shell_command =~ /bash/){return "bash";}
if($shell_command =~ //){return "bash";}

return 0;
}

####################################################################
sub test_translation_quality{
####################################################################
my $infile = shift;
my $temp_alignment_file_out = shift;
my $minimum_length_of_hit = shift;
my $number_of_matches =0;
my $number_of_comparisons =0;
die "Could not create temporary filehandle: $! \n" if !(-e $temp_alignment_file_out);

# Get two sequences
   my $str = Bio::SeqIO->new(-file=>$infile, '-format' => 'Fasta')
   || die "Couldnt open file $infile \n";
   my $translated_sequence = $str->next_seq();
   my $template_sequence = $str->next_seq();
   
## Avoid bad / too short query sequences
   return 1000 if $translated_sequence->length < $minimum_length_of_hit;
   # Run bl2seq on them
   my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastp',
                                               'outfile' => $temp_alignment_file_out,
					       #'filter' => 'F'
					       );
   #E-value
 $factory->e(1e-10);
  # Turn filtering off
  $factory->F(F);

   my $bl2seq_report = $factory->bl2seq($translated_sequence, $template_sequence);
## RETURN DEFAULT IF NO HITS
   return 1000 if !$bl2seq_report;
### PARSE BLAST RESULTS
    while ( my $result = $bl2seq_report->next_result) {
        while ( my $hit = $result->next_hit) {
	return $hit->significance;
	}
    }
return 1;
}





####################################################################
sub translate_bioperl{
####################################################################
my $seq = shift;
my %translated_hash = ();

     my $len  = $seq->length();
     for(my $frame=1; $frame < 4; $frame++){
     my $seqobj = Bio::Seq->new(
                        -seq        => $seq->subseq($frame,$seq->length),
                        -display_id => "najaa",
                        -alphabet   => 'dna'
                    );
     $translated_hash{"+$frame"}{seq} = $seqobj->translate('X', 'X', '0', '1', '0', '0', '0', '0')->seq; 
     }
    
    my $rev = $seq->revcom;
     for(my $frame=1; $frame < 4; $frame++){
     my $seqobj = Bio::Seq->new(
                        -seq        => $rev->subseq($frame,$rev->length),
                        -display_id => "najaa",
                        -alphabet   => 'dna'
                    );
     $translated_hash{"-$frame"}{seq} = $seqobj->translate('X', 'X', '0', '1', '0', '0', '0', '0')->seq; 
     }
return \%translated_hash;
}

####################################################################
sub calculate_pairwise_distance{
####################################################################
my $infile = shift;
my $number_of_matches =0;
my $number_of_comparisons =0;
my ($temp_alignment_fh_out, $temp_alignment_file_out) = tempfile(UNLINK => 1);
die "Could not create temporary filehandle: $! \n" if !(-e $temp_alignment_file_out);

# Get two sequences
   my $str = Bio::SeqIO->new(-file=>$infile, '-format' => 'Fasta')
   || die "Couldnt open file $infile \n";
   my $seq3 = $str->next_seq();
   my $seq4 = $str->next_seq();
   
   # Run bl2seq on them
   my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastp',
                                               'outfile' => $temp_alignment_file_out);
  $factory->e(1e-10);

   my $bl2seq_report = $factory->bl2seq($seq3, $seq4);
  return 1 if !$bl2seq_report;
### PARSE BLAST RESULTS
    while ( my $result = $bl2seq_report->next_result) {
        while ( my $hit = $result->next_hit) {
	return $hit->significance;
	}
    }
return 1;
}

1;
