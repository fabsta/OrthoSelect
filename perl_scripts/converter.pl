#!/usr/bin/perl
### LIBRARIES{{{
use Getopt::Std;
use Bio::Perl;
use Data::Dumper;
use Pod::Usage;
use File::Copy;
use Bio::SeqIO;
use File::Basename;
use strict;
use Cwd;

#use warnings;
### VARIABLES#{{{
my %options;
getopt( "iots", \%options );
MAIN: {

	my $dir = getcwd;
	my $output_directory = "$dir/est_libraries";
	my @input_files = glob("$output_directory/*");
	if(!@input_files){
		print "\tWarning: There are no input files in directory \"est_libraries\".Exiting...\n";
	}
	my $taxon_counter = 1;
	foreach(@input_files){
		print "\tConverting headers of file $_\n";
		convert_fasta_file({sequence_file => $_, taxon_counter => $taxon_counter++});
	}
	exit;
	my $file_containing_conversions = "taxa_conversion.txt";
	if(-e $file_containing_conversions){
		my %taxa_conversion_hash = %{&read_taxa($file_containing_conversions)};
		foreach my $file (keys %taxa_conversion_hash){
			&convert_fasta_file($file,$taxa_conversion_hash{$file}{taxa_name},$taxa_conversion_hash{$file}{source},$output_directory);
		}
	}
	else{
		my $input = $options{i} || die "Please specify input file\n";
		my $output = $options{t} || die "Please specify input file\n";
		my $source = $options{s} || die "Please specify source of sequences (JGI (t), JGI-ESTs(e) or TBestDB (d) or dbEST (n))";
		&convert_fasta_file($input,$output,$source,$output_directory);
	}
}
	
####################################################################
sub convert_fasta_file(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $fasta_file = $arg_ref->{sequence_file};
my $taxon_counter = $arg_ref->{taxon_counter};
#### OTHER VARIABLES
#my $assignment_outfile = $assignment_file."_converted";
my $outfile = $fasta_file."_converted";


if( -e $outfile && -s $outfile){
#if(-e $assignment_outfile &&  -e $outfile ){
	print "\t$outfile\n";
	print "\t[CONVERT SEQUENCES] Converted files already exist. Skipping\n";
	return 1;
}

	my $conversion_file = "conversion.txt";
	### Slurping File
	my $fh = new IO::Handle;
	open $fh, '<', $fasta_file;
	my $text = q{};
	sysread $fh, $text, -s $fh;
	$fh->flush;
	close $fh;
#	print "\t[CONVERT SEQUENCES] All Sequences read into memory\n";
	
	my %converter_hash = ();
	my @whole_lines = split(/\n/,$text);
	my $curr_sequence = q{};
	my $curr_header = q{};
	my $curr_domain = q{};
	#print "\t[CONVERT SEQUENCES] Processing file line by line\n";
	
	open my $OUTFILE_FH, '>', $outfile or warn "Can't open ".$outfile.": $!";
	open my $CONVERSION_FH, '>>', $conversion_file or warn "Can't open ".$conversion_file.": $!";
	my $sequence_counter = 1;
	foreach(@whole_lines){
		if(/>(.*)/){
			$curr_header = $1;
			###
			# Conversion to a format handable by OrthoSelect
			# taxon1_sp
			###
			
#			print {$OUTFILE_FH} ">Tax".$taxon_counter."S".$sequence_counter."_sp\n";
#			print {$CONVERSION_FH} "$curr_header-->Tax".$taxon_counter."S".$sequence_counter."_sp\n";
#			$converter_hash{$curr_header} = "Tax".$taxon_counter."S".$sequence_counter."_sp";
			print {$OUTFILE_FH} ">seq".$sequence_counter."\n";
			print {$CONVERSION_FH} "$curr_header-->seq".$sequence_counter."\n";
			$converter_hash{$curr_header} = "seq".$sequence_counter;
			$sequence_counter++;
		}
		else{
			print {$OUTFILE_FH} $_."\n";
		}
	}
#	print "\t[CONVERT SEQUENCES] Done converting headers of fasta file\n";
	#exit;
	close $OUTFILE_FH or warn "Can't close ".$outfile.": $!";
	close $CONVERSION_FH or warn "Can't close ".$conversion_file.": $!";

	
	###
	# Copy converted file to old file
	###
	if(-e $outfile && -s $outfile){
		move($outfile, $fasta_file)
	}
	else{
		print "\tCould not convert sequences. If you think that this is a bug, sent an email with your sequences attached to the author \n";
	}#	exit;
	return 1;
}

	
####################################################################
sub reconvert_fasta_file(){
####################################################################
#### PARAMETER VARIABLES
my ($arg_ref) = @_;
my $fasta_file = $arg_ref->{sequence_file};
my $assignment_file = $arg_ref->{all_predictions};
#### OTHER VARIABLES
my $assignment_outfile = $assignment_file."_converted";
my $outfile = $fasta_file."_converted";

}


sub read_taxa(){
my $file = shift;
die "Could not read taxa_conversion.txt\nCheck if it exists!\n" if !-e $file;
my %taxa_conversion_hash = ();

	 open (FILE, "$file") or die "Could not open file\n";
	

while(<FILE>){
	/"(.*)"="(.*)"="(.*)"/;
	my $file = $1;
	my $taxon =$2;
	my $source =$3;
	die "Please check format of taxa_conversion.txt" if not defined $file;
	$taxa_conversion_hash{$file}{taxa_name} =$taxon;
	$taxa_conversion_hash{$file}{source} =$source;
}
        close(FILE)or die "Could not close file\n";
return \%taxa_conversion_hash;
}
