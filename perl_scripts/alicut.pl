#!/usr/bin/perl
use strict                                                                                                   ;
use File::Copy                                                                                               ;
use Tie::File                                                                                                ;
use Fcntl                                                                                                    ;

# updated on 13th febrauary 2009 by patrick kück

printf "\n%68s\n", "------------------------------------------------------------"                            ;
printf "%52s\n"  , "Welcome to ALICUT V2.0 !"                                                                ;
printf "%60s\n"  , "a Perlscript which cuts ALISCORE identified RSS"                                         ;
printf "%56s\n"  , "written by Patrick Kueck (ZFMK, Bonn)"                                                   ;
printf "%68s\n\n", "------------------------------------------------------------"                            ;

REPEAT1:
print
<<start;

	To start ALICUT V2.0, press return.
	To leave ALICUT V2.0, press <q> and return.
	To remain stem positions of identified rss type <s> and return.
	For general information and usage type <help> and return.
	For a more detailed error-report, type <ErrorX>
	(X = Error associated number from 1 to 2).
	
start

printf "\n%68s\n\t", "------------------------------------------------------------"                            ;

chomp ( my $answer = <STDIN> );

$answer =~ /q/i and do { exit; };

$answer =~ /s/i and do { 

print
<<structure_info;

	RSS identified stem positions are remained !
	
structure_info
};

# help informations
$answer =~ /help/i and do {

 print
 <<info;
 
	General Information and Usage:
	-------------------------------
	ALICUT V2.0 removes all ALISCORE identified RSS positions 
	in given FASTA file(s) which are listed in the FASTA file cor-
	responding "List" file(s) of ALISCORE. If structure sequences
	are implemented, ALICUT V2.0 will automatically replace brackets 
	of non rss positions by dots when they are paired with rss 
	identified positions. To remain all stem positions of identified rss 
	within the FASTA file type <s> and return. For using ALICUT V2.0, 
	you need the original ALISCORE FASTA inputfile(s) and "List"
	file(s) in the same folder as ALICUT V2.0.
	
	The "List" outfile(s) must contain the identified RSS positions
	in one single line, separated by whitespace.
	
	ALICUT V2.0 can handle unlimited FASTA files in one single run.
	The sole condition is that no filename, neither the ALISCORE FASTA
	inputnames(s) nor the ALISCORE "List" outfile(s) are changed
	since the use of ALISCORE! ALICUT V2.0 first takes the ALISCORE 
	"List" outfile(s), and then searches for the belonging FASTA file(s). 
	If both are detected, ALICUT V2.0 will cut the RSS identified 
	positions listed in the "List" file(s) out of the associated
	FASTA file and will save the changes in a new FASTA outfile,
	named "ALICUT_FASTAinputname.fas", within a new generated folder 
	named "ALICUT_output". The original FASTA file(s) are moved to 
	the new folder named "/ALICUT_fasta-input_orig". If two "List" files
	were generated from the same FASTA file, maybe the first was
	depending on a NJ tree and the second one on all single comparisons,
	the first outfile will be overwritten by the second.
	
	Additionally ALICUT V2.0 generates an info file "ALICUT_info" within 
	the new folder "/ALICUT_output". In this file information about
	number and percent of removed positions is given for each executed
	FASTA file. If one file was overwritten by the second you can see here
	which one was left and which was rejected. Further ALICUT V2.0 generates
	a structure info file for each single FASTA file which lists all 
	corresponding stem positions, all loop positions and the percentage of
	both structure elements.
	
info
goto REPEAT1;
};


# error reports
$answer =~ /error1/i and do {

 print
 <<error1;

	Error1-Reports:
	-------------------------------
	Error2 => Error2 reports are depending on potential wrong "List" 
	files. ALICUT V2.0 checks each "List" file before it starts and 
	if something seems not ok, ALICUT V2.0 will print a small report 
	on the desktop:
	
	>"Can not find original named ALISCORE listfile in .txt format!"
	ALICUT V2.0 has identified the "List" file, but maybe it is not 
	in .txt format. Check the file extension. Are the "List" file 
	with ALICUT V2.0 in the same folder?
	
	>"File is not a original named ALISCORE listfile!"
	This means that the textfile is not a original ALISCORE "List" 
	outfile. Maybe the filename has changed or the file is just 
	another textfile which can be denied. In last case ignore the 
	error report. Else, change the filename back to the original 
	standard. To find the "List" outfile" of ALISCORE, ALICUT V2.0
	orientates itselfs by the word "List" which must be implemented 
	in the "List" outfile. Upper or lower cases can be ignored.

	>"File has no ALISCORE list format!"
	ALICUT V2.0 has identified the "List" file but maybe the format
	within is wrong. Check the input of the "List" file. Only posi-
	tion numbers, separated by one whitespace, all in the first line,
	are allowed. Further, the line must start and end with a position
	number. If the format alters, ALICUT will hopefully identify it
	before it starts.

error1
goto REPEAT1;
};

$answer =~ /error2/i and do {

 print
 <<error2;

	Error2-Reports:
	-------------------------------
	Error3 => Error3 reports depending most likely on the given FASTA 
	inputfile. 
	
	>"Can not find FASTA file!"
	Maybe the name of the given "List" file does not agree with the 
	associated FASTA file. Normally, if nothing was changed from the 
	ALISCORE output, the FASTA filename is part of the "List" file-
	name: <FASTAfile.fas_List...txt>. Alternatively the FASTA file
	is not within the same folder.
	
	>"File FASTA file is empty!"
	Well, what do I have to say about this error?!
	
	>"taxon name missing in FASTA file!"
	ALICUT V2.0 can not find min. one taxonname within FASTA file.
	
	>"Wrong sequence signs in sequence of FASTA file!"
	ALICUT V2.0 scans every FASTA sequence for forbidden signs 
	within it, e.g. "=", "2", the execution of the depending 
	FASTA file will be interrupted.
	
	>"No equal sequence lengths in FASTA file!"
	ALICUT V2.0 indicates if sequences have no equal length indicates 
	that this FASTA file does not contain aligned sequences. So the 
	execution of this file will be interrupted.
	
	>"To short Sequence lengths in FASTA file!"
	If the sequences were to short to cut all demand positions, 
	denotes that the given "List" file does not depent to the given 
	FASTA file. Maybe you have to check if both filenames are really 
	correct to each other.
	
	>"Multiple taxa of name in FASTA file"
	If there are 2 equal sequence names within the same FASTA file,
	ALICUT V2.0.pl will skip this file.
	
	>"Sequence missing in FASTA file!"
	ALICUT V2.0 can not find min. one sequence within FASTA file.
	
error2
exit;
};

# Define global variables
my @cut_positions   = () ;
my @fasta_files     = () ;
my $structurestring = () ;
my $length_seq      = () ;
my $file_struc      = () ;
my $j               = 0  ;


## Open new folder for original inputfiles
#mkdir   "ALICUT_fasta-input_orig", 0755 
#or warn "Error1: Can not generate folder \"ALICUT_fasta-input_orig\": $!\n"  ;

# Open new folder for ALICUT outputs
#mkdir   "ALICUT_output", 0755 
#or warn "Error1: Can not generate folder \"ALICUT_output\": $!\n" ;

# Starting point for next commands
READING:

# Read IN of all .txt files within the same folder as ALICUT and handle it
foreach my $file ( <*.txt> ) {
	
	# Set counter +1
	$j++;
	
	# Define local variables for each read in of .txt file
	my ( %seen, @inputfile, @zahlenreihe, @fasta_cut ) = ();
	
	
	# Evaluation of each .txt file for ALISCORE-list outfile
	if  ( $file =~ /List/i ) { open IN, "<$file" or die "Can not find original named ALISCORE listfile in .txt format" }
	else                     { warn "\tError1: $file is not a original named ALISCORE listfile!\n"        ; next READING }
	

	# Read in of the ALISCORE-list outfile
	chomp  ( my $line = <IN> ) ;
	
	unless ( $line =~ /^(\d+ )+\d+$|^\d+$/ ) { warn "\tError1: $file has no ALISCORE list format\n" ; next READING }

	
	# Total number of randomized identified positions
	my $number_cut_positions     = ( @cut_positions = split " ", $line ) ;
	
	
	# For opening the right fastafile for each listfile the last signs of the listfile must be removed from "filename.fas_List_random.txt" to "filename.fas"
	my $wo                       = index  ( $file, "List" )          ;
	my $file_fasta               = substr ( $file, 0, ( $wo - 1 ) )  ;
	   $file_struc               = substr ( $file, 0, ( $wo - 4 ) )  ;
	
	
	#checking line feeds and making shure that just \n remains
	(open IN , $file_fasta ) or warn "\tError2: $file_fasta can not be opened!\n" and next READING ;
	
	(tie ( my @data, 'Tie::File', $file_fasta )) ;
	
	warn "\tError2: $file_fasta is empty!\n" and next READING if 0 == @data ;
	
	map { s/\r\n/\n/g } @data ;
	map { s/\r/\n/g   } @data ;

	untie @data ;
	

	# Read in of the original ALISCORE fasta infile which belongs to the listfile
	open INfas, "<$file_fasta" or warn "\tError2: Can not find $file_fasta!\n" and next READING ;
	
	
	# Read IN the list associated FASTA file to array @inputfile
	chomp ( @inputfile = <INfas> ) ;
	warn "\tError2: File $file_fasta is empty!\n" if 0 == @inputfile and next READING ;
	
	
	# Handle the FASTA file in the way that sequencename and sequence alternate in each line
	@inputfile                   = fas_bearbeiten ( @inputfile ) ;
	my $number_characters_before = $inputfile[1] ;
	
	
	# Generate a hash: key=>taxon, value => sequenz
	my %sequence                 = @inputfile ;
	my @values                   = values %sequence ;
	
	
	# Determine basepositions before und after cut. Output of cuttings as total number and in percent
	   $length_seq               = length $values[0] ;
	my $number_characters_after  = $length_seq-$number_cut_positions ;
	
	
	# Evaluation of sequences
	for my $key ( keys %sequence ){
		
		
		# if whitespace are between ">" and the next sign within sequence name, delete these whitespaces
		$key =~ s/^\>\s*/\>/g ;
		
		# if whitespaces between last sign and newline in sequence name, delete these whitespaces
		$key =~ s/\s*$//g ;
		
		my %FASTA = () ;
		
		warn "\tError2: Taxon name missing in $file_fasta!\n"                                 and next READING if $key =~ /^\>$/                                                            ;
		warn "\tError2: Sequence missing in $file_fasta!\n"                                   and next READING if $sequence{$key} =~ /^\n$|^$/                                              ;
		warn "\tError2: Multiple taxa of $key in $file_fasta\n"                               and next READING if defined $FASTA{$key}                                                      ;
		warn "\tError2: $key in $file_fasta is not in FASTA format\n"                         and next READING if	$key !~ /^\>/                                                             ;
		warn "\tError2: Sequence of $file_fasta involves forbidden signs in $key\n"           and next READING if $sequence{$key} =~ /\d|\:|\;|\,|\*|\>|\<|\%|\$|\§|\"|\!|\^|\_|\+|\#|[|]/  ;
		warn "\tError2: Sequences of $file_fasta have no equal length\n"                      and next READING if length $sequence{$key} != $length_seq                                     ;             
		warn "\tError2: Sequence lengths in $file_fasta are to short to cut all positions\n"  and next READING if $length_seq < $cut_positions[ $#cut_positions ]                           ;
		
		# Structure tagging, handling and dot allocation
		if ( $sequence{$key} =~ /.*\(.*\).*/ ){
			
			   $structurestring =  $sequence{$key} ; 
			   $structurestring =~ s/-/./g ;
			   $sequence{$key}  =  &structure_handling ( $structurestring ) ;
		}
	}
	
	my $percent_left = ( $number_characters_after / $length_seq ) * 100 ;
	
	
	# Generate a array of numbers from 1 till last position to cut
	for ( 0..($length_seq-1), my $i=0 ) { push @zahlenreihe, $i; $i++ }
	
	
	# The ALISCORE-list file numbered the first seqposition with 1, the first position in perl is numbered 0 !
	# Therefore every position of the list-file must be substrated with -1 to get the identical positions !
	for my $value ( @cut_positions ) { $value = $value - 1 }
	
	
	# Read IN of cut positions to a hash
	for ( @cut_positions ){ $seen{$_} = 1 }
	
	
	# Compare array of numbers with @cut_positions
	for ( @zahlenreihe ){ unless ( $seen{$_} ){ push @fasta_cut, $_ } }
	
	
	open OUT, ">ALICUT_$file_fasta" ;
	
	
	# Assume uncut positions to $final and print them out to ALICUT_$file_fasta
	for ( keys %sequence ){
		
		my @bases = split "", $sequence{$_}          ;
		my @final = map { $bases[$_] } @fasta_cut    ;
		my $final = $_."\n".( join "", @final )."\n" ;
		
		print OUT "$final" ;
	}
	
	
	# Readout of extra infos to ALICUT_info
	open  OUTinfo, ">>ALICUT_info.txt"                                                                ;
	print OUTinfo  "\nUsed ALISCORE \"List\": $file\n"                                                ;
	print OUTinfo  "Used FASTA file: $file_fasta\n"                                                   ;
	print OUTinfo  "Basenumber before cut:\t$length_seq\n"                                            ;
	print OUTinfo  "Basenumber after cut:\t$number_characters_after\n"                                ;
	print OUTinfo  "Percent left:\t$percent_left%\n\n"                                                ;
	print          "\tDone  : $file cut to ALICUT_$file_fasta\n"                                      ;

	push @fasta_files, $file_fasta                                                                    ;

}


# Count the number of right handled FASTA files
my $N_bearbeitet = @fasta_files ;

close IN       ;
close INfas    ;
close OUTinfo  ;


# Copy all files to their belonging folders and delete the originals
#foreach (@fasta_files) { move    ( "$_", "ALICUT_fasta-input_orig/$_" ) or warn   "\nError4: Can not move \"$_\" !\n" }
#move    ( "ALICUT_info", "ALICUT_output/ALICUT_info" ) or warn   "Error4: Can not move \"ALICUT_info\"\n"  ;
#move    ( "ALICUT_Struc_info_${file_struc}", "ALICUT_output/ALICUT_Struc_info_${file_struc}" ) or warn   "Error4: Can not move \"ALICUT_Struc_info_${file_struc}\"\n"  ;


# Print OUT number of right handled FASTA files in relation to total number of files
printf "\n%68s\n",   "------------------------------------------------------------" ;
printf "%42s\n",     "$N_bearbeitet FASTA file(s) correctly handled!"               ;
printf "%57s\n",     "Further infos are printed out in Alicut_info.txt!"            ;
printf "\n%63s\n",   "ALICUT V2.0 Finished! Thank you for using it! Good bye!"      ;
printf "%68s\n",     "------------------------------------------------------------" ;


# set timer
my ( $user, $system, $cuser, $csystem ) = times ;

print <<TIME;

			***  time used: $user sec  ***

TIME
#

exit ;


sub fas_bearbeiten{
	
	my @infile = @_                   ;
	
	grep  s/(\>.*)/$1\t/,     @infile ;
	grep  s/ //g,             @infile ;
	grep  s/\n//g,            @infile ;
	grep  s/\t/\n/g,          @infile ;
	grep  s/\>/\n\>/g,        @infile ;
	my $string = join "",     @infile ;
	@infile    = split "\n",  $string ;
	shift                     @infile ;
	return                    @infile ;
}


sub structure_handling{

	my @pair_infos            =  ()                           ;
	my @forward               =  ()                           ;
	my @structurestring       =  ()                           ;
	my @loops                 =  ()                           ;
	my @pairs                 =  ()                           ;
	my %structure_of_position =  ()                           ;
	my %seen_struc            =  ()                           ;
	my @structures            =  split "", $structurestring   ;
	
	
	# Stem assignment
	my  $i = 0                                                                                                         	                  ;
	CHECKING:
	for ( @structures ){ $i++                                                                                                             ;
		
		SWITCH:
		$structure_of_position{$i} = $_                                                                                                   ;
		
		if ( $_  =~ /\(/ ){ push @forward, $i                                                                          and next CHECKING  }
		if ( $_  =~ /\)/ ){ my $pair_1 = pop @forward; push @pairs, ( $pair_1, $i ); push @pair_infos, ( $pair_1.":".$i ); next CHECKING  }
		if ( $_  =~ /\./ ){ push @loops,   $i                                                                          and next CHECKING  }
	}
	
	@pair_infos  =  reverse @pair_infos                                                                                                   ;
	
	
	
	
	# Generate listfiles for structure_info file
	my $pairlist =  join "\n\t\t\t\t\t", @pair_infos   ;
	my $looplist =  join "\n\t\t\t\t\t", @loops        ;
	
	
	# Number and proportion of stem and loop positions for structure info file
	my $N_total  =  @structures                        ;
	my $N_stems  =  @pair_infos                        ;
	my $N_loops  =  $N_total - ( $N_stems * 2 )        ;
	my $P_loops  =  ( $N_loops / $N_total ) * 100      ;
	my $P_stems  =  100 - $P_loops                     ;

	
	# Open structure info outfile
	open OUTstruc, ">ALICUT_Struc_info_${file_struc}txt"                             ;
	
	# Print out
	print OUTstruc "\nOriginal structure information identified in $file_struc:\n\n"  ;
	print OUTstruc "- Number of characters:\t\t\t$N_total\n"                          ;
	print OUTstruc "- Number of single loop characters:\t$N_loops [$P_stems %]\n"     ;
	print OUTstruc "- Number of paired stem characters:\t$N_stems [$P_loops %]\n"     ;
	print OUTstruc "\n- Paired stem positions:\t\t$pairlist\n\n"                      ;
	print OUTstruc "\n- Loop positions:\t\t\t$looplist\n"                             ;

	close OUTstruc;
	
	if  ( $answer =~ /s/i ){
		
		my @cut_positions2 = ();
		
		# Remain rss identified stem positions within the MSA
		for ( @pairs ){ $seen_struc{$_} = 1                                           }
		for ( @cut_positions ){ unless ( $seen_struc{$_} ){ push @cut_positions2, $_  } }
		@cut_positions = @cut_positions2                                              ;
	}
	
	else{
		
		my %pair = @pairs;
		
		# Replace paired structure positions of rss identified positions by dots
		for my $bp_for ( keys %pair ){
			
			for my $rss ( @cut_positions ){
				
				if ( $bp_for        == $rss ){ $structure_of_position{$pair{$bp_for}}  = "." ; last }
				if ( $pair{$bp_for} == $rss ){ $structure_of_position{$bp_for}         = "." ; last }
			}
		}
	}
	
	for    ( my $k=1; $k<=$length_seq; $k++ ){ push @structurestring, $structure_of_position{$k}   }
	my     $structure_string_neu = join "", @structurestring                                       ;
	return $structure_string_neu                                                                   ;
	
}	
	
	
	
	
	
	
	
	

