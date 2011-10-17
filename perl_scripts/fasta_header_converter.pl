#!/usr/bin/perl
### LIBRARIES{{{
use Getopt::Std;
use Bio::Perl;
use Data::Dumper;
use Pod::Usage;
use Bio::SeqIO;
use File::Basename;
use strict;
#use warnings;
### VARIABLES#{{{
my %options;
getopt( "iots", \%options );
MAIN: {
	
	my $output_directory = "est_libraries";
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
	
	sub convert_fasta_file{
		my ($input, $output,$source, $output_directory) = @_;
		#$input = "/".$input;
		my $in  = Bio::SeqIO->new(-file => $input ,
                           -format => 'Fasta');
    my $out = Bio::SeqIO->new(-file => ">$output_directory/$output.fa" ,
                           -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
#	    print $seq->id."\n";
## JGI-Transcript
	if($source eq "t"){
	    $seq->id =~ /jgi\|.*\|(\d*)\|/;
	 #   print $1."\n";
	    my $accession = $1;
	    if(!defined $accession){ 
		    die "File not in JGI-Format.\nCant change fasta header\n";
	    }
	    $seq->id($accession."|".$output);
	    $seq->description("");
	}
## JGI-EST
	if($source eq "e"){
	    $seq->id =~ /(\w*[:\d]*)/;
	    #print $1."\n";
	    my $accession = $1;
	    
	    if(!defined $accession){ 
		    die "File not in JGI-Format.\nCant change fasta header\n";
	    }
	    $accession =~ s/:/_/g;
	    $seq->id($accession."|".$output);
	    $seq->description("");
	   # print $seq->id."\n";
	}
## TBestDB	
	if($source eq "d"){
	    $seq->description =~ /Id\s?:\s?(\w*)/;
	#    print $seq->id."\t desc: ".$seq->description."\n";
	    my $accession = $1;
	    if(!defined $accession){ 
		    die "File not in TBestDB-Format.\nCant change fasta header\n";
	    }
	    $seq->id($accession."|".$output);
	    $seq->description("");
	}
## dbEST
	if($source eq "n"){
		$seq->id =~ /gi\|(\w*)\|gb/;
		my $accession = $1;
		if(!defined $accession){ 
		    die "File not in dbEST-Format.\nCant change fasta header\n";
		}
		$seq->id($accession."|".$output);
		$seq->description("");
	}
	$out->write_seq($seq);
    }
    
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
