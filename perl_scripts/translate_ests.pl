#!/usr/bin/perl
### LIBRARIES{{{
use Getopt::Std;
use Bio::Perl;
use Data::Dumper;
use Pod::Usage;
use Bio::SeqIO;
use File::Basename;
use strict;
require 'lib/Iostuff2.pm';
use warnings;
### VARIABLES#{{{
my %options;
getopt( "io", \%options );

my $options_file = $options{t} || "options.txt";
my %options_hash = %{ Iostuff::read_options($options_file) };
Iostuff::check_taxa_list_syntax($options_hash{taxa_list});
my %taxa_shortcuts = %{Iostuff::read_taxa($options_hash{taxa_list})};
MAIN: {

	my $input_directory = $options{i} || die "please specify output directory for analysis\n";
	my $output_directory = $options{o} || die "please specify output directory for analysis\n";
	
	if(! -e $output_directory){
		Iostuff::create_single_directory($output_directory);	
	}
	else{
		my $existing_translations = `grep ">" $output_directory/* | wc -l`;
		chomp($existing_translations);
		if($existing_translations){
			print "\tWarning: Translations already exists in $output_directory.\n\tAborting translation...\n";
			exit;
		}
	}
	my @folders = <$input_directory/*>;
	
	print "\tStart Translation of EST sequences...(This may take several minutes)....\n";
	
	foreach my $file(@folders){
		my $taxon = basename($file);
		$taxon =~ s/\.fa//g;
		print "\tCurrent Taxon: ".$taxon."\n";
		my $output_file = "$output_directory/$taxon.fa";
		my @est1 =  `ESTScan  $file -t - -M \$ESTSCANDIR/At.smat `;
		my @est2 =  `ESTScan  $file -t - -M \$ESTSCANDIR/Dm.smat `;
		my @est3 =  `ESTScan  $file -t - -M \$ESTSCANDIR/Dr.smat `;
		my @est4 =  `ESTScan  $file -t - -M \$ESTSCANDIR/Mm.smat `;
		my @est5 =  `ESTScan  $file -t - -M \$ESTSCANDIR/Rn.smat `;
		my @est6 =  `ESTScan  $file -t - -M \$ESTSCANDIR/Hs.smat `;
			 
	### ADD UNIQUE IDENTIFIER TO FASTA HEADER.
	### IMPORTANT FOR "REMOVE DUPLICATE"-METHOD THAT WOULD REMOVE
	### SEQUENCES WITH IDENTICAL FASTA HEADERS
		open my $OUTEST, ">", $output_file or die "Could not open $output_file $!\n";
		foreach(@est1){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."AT\n";}else{print $OUTEST $_;}}
		foreach(@est2){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."DM\n";}else{print $OUTEST $_;}}
		foreach(@est3){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."DR\n";}else{print $OUTEST $_;}}
		foreach(@est4){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."MM\n";}else{print $OUTEST $_;}}
		foreach(@est5){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."RN\n";}else{print $OUTEST $_;}}
		foreach(@est6){if(/>(\w*_?\w*_?\w*)/){print $OUTEST ">".$taxa_shortcuts{$taxon}."_".$1."HS\n";}else{print $OUTEST $_;}}
		close $OUTEST or die "Could not close $output_file\n";
		}
#	exit;
	print "\n\tTranslation finished\n";
	}
