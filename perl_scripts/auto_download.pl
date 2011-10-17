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
# 
########################################################################
use Getopt::Std;
use File::Basename;
use LWP::Simple;
use strict;
use warnings;
my $LINE_LENGTH = 64;

#require 'lib/Iostuff2.pm';

my %options=();
getopt("top",\%options);

my $env_path = $ENV{PATH};
my $estscandir_path = $ENV{ESTSCANDIR};
my $wiseconfigdir_path = $ENV{WISECONFIGDIR};
my  $curr_dir = `pwd`;
$curr_dir =~ s/\n//;
my $install_bioperl = $options{p} || "";
my $program_dir = "programs/";
		
# Create directory to store the programs
system("mkdir $program_dir");
my %download_hash =(
	"bioperl" => "ftp://bioperl.org/pub/bioperl/DIST/bioperl-1.5.1.zip",
	"Btlib" => "http://downloads.sourceforge.net/estscan/BTLib-0.17.tar.gz?modtime=1175013253&big_mirror=0",
	"ESTScan" => "http://downloads.sourceforge.net/estscan/ESTScan2-2.1.tar.gz?modtime=1166457155&big_mirror=0",
	"blastall" => "ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/blast-2.2.18-universal-macosx.tar.gz",
	"muscle" => "http://www.drive5.com/muscle/downloads3.7/muscle3.7_src.tar.gz",
	"clustalw" => "ftp://ftp.ebi.ac.uk/pub/software/clustalw2/2.0.9/clustalw-2.0.9-src.tar.gz",
	"t_coffee" => "http://www.tcoffee.org/Packages/T-COFFEE_distribution.tar.gz",
	"Gblocks" => "http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Universal_0.91b.tar.Z",
	"noisy" => "http://www.bioinf.uni-leipzig.de/Software/noisy/Noisy-1.5.7.tar.gz",
	"genewise" => "ftp://ftp.ebi.ac.uk/pub/software/unix/wise2/wise2.2.0.tar.gz"
	);
my %install_hash =(
	"bioperl" => "unzip bioperl-1.5.1.zip -d $program_dir| cd $program_dir/bioperl-1.5.1| perl Makefile.PL | make | sudo make install",
	"Btlib" => "tar -zxvf BTLib-0.17.tar.gz -C $program_dir| cd $program_dir/BTLib-0.17| perl Makefile.PL | make | sudo make install",
	"ESTScan" => "tar -zxvf ESTScan2-2.1.tar.gz -C $program_dir|cd $program_dir/ESTScan2-2.1|perl Makefile.PL|make|sudo make install| export ESTSCANDIR $curr_dir/$program_dir/ESTScan2-2.1/",
	"blastall" => "tar -zxvf blast-2.2.18-universal-macosx.tar.gz -C $program_dir| export PATH $curr_dir/$program_dir/blast-2.2.18/bin",
	"muscle" => "mkdir $program_dir/muscle3.7_src|tar -zxvf muscle3.7_src.tar.gz -C $program_dir/muscle3.7_src| cd $program_dir/muscle3.7_src|perl -i -pe 's/-lm -static/-lm/g' Makefile|make | export PATH $curr_dir/$program_dir/muscle3.7_src",
	"clustalw" => "tar -zxvf clustalw-2.0.9-src.tar.gz  -C $program_dir| cd $program_dir/clustalw-2.0.9| ./configure| make| sudo make install| sudo mv /usr/local/bin/clustalw2 /usr/local/bin/clustalw",
	"t_coffee" => "tar -zxvf T-COFFEE_distribution.tar.gz -C $program_dir| cd $program_dir/T-COFFEE_distribution_Version_7.04/|./install|mkdir \$HOME/.t_coffee|mkdir \$HOME/.t_coffee/cache|mkdir \$HOME/.t_coffee/methods|mkdir \$HOME/.t_coffee/mcoffee|mkdir \$HOME/.t_coffee/tmp,
	#export PATH $curr_dir/$program_dir/T-COFFEE_distribution_Version_6.30/bin",
	"Gblocks" => "gzip -d Gblocks_Universal_0.91b.tar.Z | tar xvf Gblocks_Universal_0.91b.tar -C $program_dir| export PATH $curr_dir/$program_dir/Gblocks_0.91b",
	"noisy" => "tar -zxvf Noisy-1.5.7.tar.gz -C $program_dir| cd $program_dir/Noisy-1.5.7| ./configure | make | sudo make install |",
	"genewise" => "tar -zxvf wise2.2.0.tar.gz -C $program_dir| cd $program_dir/wise2.2.0/src/| make all| export PATH $curr_dir/$program_dir/wise2.2.0/src/bin/| export WISECONFIGDIR $curr_dir/$program_dir/wise2.2.0/wisecfg/"
	);
my %progname_hash =(
	"bioperl" => "bioperl-1.5.1.zip",
	"Btlib" => "BTLib-0.17.tar.gz",
	"ESTScan" => "ESTScan2-2.1.tar.gz",
	"blastall" => "blast-2.2.18-universal-macosx.tar.gz",
	"muscle" => "muscle3.7_src.tar.gz",

	"clustalw" => "clustalw-2.0.9-src.tar.gz",
	"t_coffee" => "T-COFFEE_distribution.tar.gz",
	"Gblocks" => "Gblocks_Universal_0.91b.tar.Z",
	"noisy" => "Noisy-1.5.7.tar.gz",
	"genewise" => "wise2.2.0.tar.gz",
	);

;
my $bash_ok = &check_shell();
	die "Cant determine type of shell used\nAborting...\n" if !$bash_ok;

print "You are using a $bash_ok shell\n";	


## LINUX AS OPERATING SYSTEM
if(defined $options{o} && $options{o} eq "linux"){
	$download_hash{"blastall"} = "ftp://ftp.ncbi.nih.gov/blast/executables/LATEST/blast-2.2.18-ia32-linux.tar.gz";
	$install_hash{"blastall"} = "tar -zxvf blast-2.2.18-ia32-linux.tar.gz -C $program_dir| export PATH $curr_dir/$program_dir/blast-2.2.18/bin",
	$progname_hash{"blastall"} = "blast-2.2.18-ia32-linux.tar.gz";
	
	$download_hash{"Gblocks"} = "http://molevol.cmima.csic.es/castresana/Gblocks/Gblocks_Linux_0.91b.tar.Z";
	$install_hash{"Gblocks"} = "gzip -d Gblocks_Linux_0.91b.tar.Z | tar xvf Gblocks_Linux_0.91b.tar -C $program_dir| export PATH $curr_dir/$program_dir/Gblocks_0.91b",
	$progname_hash{"Gblocks"} = "Gblocks_Linux_0.91b.tar.Z";
}
my @required_programs = (
		"ESTScan",
		"blastall",
		"muscle",
		"t_coffee",
		"Gblocks",
		"clustalw",
		#"noisy",
		"genewise",
		);

MAIN:{


### CHECK IF REQUIRED PROGRAMS ARE PRESENT
	print "=" x $LINE_LENGTH, "\n";
	print "\tChecking Installed Programs \n";
	print "=" x $LINE_LENGTH, "\n";
	system("source \$HOME/.profile");
	my @req_prog_errors = @{ &check_required_programs(\@required_programs,0)};
	if (scalar @req_prog_errors){
		print "\nInitial tests are unable to find the following required programs:\n";
		foreach (@req_prog_errors)	{
			print $_,"\n";
			if(/ESTScan/){push(@req_prog_errors,"Btlib");}
		}
		
	print "\n\nThe missing programs will be installed\n";
	}
		#my $bioperl_installed = `perl -e 'use Bio::Perl'`;
			push(@req_prog_errors,"bioperl") if $install_bioperl eq "t";
		#	print $bioperl_installed;
		#	print "huhuuuu" if $bioperl_installed eq "";
		#	print "install it" if $bioperl_installed ne "";
			#print "naja, ist doch gut, oda?\n\n\n";exit;
		#exit;
		
	#	foreach (@req_prog_errors)	{
		foreach my $prog(sort(@req_prog_errors))	{
	#	foreach my $prog(sort keys %download_hash)	{
			print "Downloading...".$prog," (This can take some time..)\n";
			next if $prog =~ /wget/;
			# get the data
			my ($content_type, $document_length, $modified_time, $expires, $server) = head($download_hash{$prog});
			if(not defined $document_length){
				warn "Could not download ".$progname_hash{$prog}." from ".$download_hash{$prog}." : Server can be down\nTry installing manually or to a later point\n";
				next
			}
			print "The file is ".($document_length/ 1000000)." Mb big\n";
			#print "This may take some time....\n" if ($document_length/ 1000000) > 10;
			#next;
			unless(is_success(getstore($download_hash{$prog}, $progname_hash{$prog}))){
				&write_path_variables($env_path, $estscandir_path, $wiseconfigdir_path);
				print "Could not download ".$progname_hash{$prog}." from ".$download_hash{$prog}." : Server can be down\nTry installing manually or to a later point\n";
			#	next;
			}
#			my $content = get($download_hash{$prog});
#			die "Couldn't get it!" unless defined $content;
			#system("wget ".$download_hash{$prog});
			print "Installing...".$prog,"\n";
			
		#	sleep(5);
			#print("wget ".$download_hash{$_});
	#INSTALL PROGRAMM
		my @install_commands = split(/\|/,$install_hash{$prog});
		foreach(@install_commands){
			print $_."\n";
#next;
			if(/cd\s+/){
				$_ =~ s/cd\s+//;
				$_ =~ s/\s//;
				
#				print "change to $_ \n";
		##check if directory exists
				die "Directory $_ does not exist \n" if (!-e $_);		
				chdir($_) or die "Couldnt change directory\n";
			}
			elsif(/export/){
				/\s*export\s*(\w*)\s*(.*)/;
				my $variable = $1;
				my $env_command = $2;
#				print "PATH before: ".$ENV{$variable}."\n";
#				print "will be changed to \n";
#				print "\$ENV{$variable}=$env_command \n";
				if($variable =~ /PATH/){
					$env_path .= ":$env_command";
				}
				else{
					$estscandir_path = $env_command if $variable =~ /ESTSCANDIR/;
					$wiseconfigdir_path = $env_command if $variable =~ /WISECONFIGDIR/;
				}
				
			#	system("echo $_ > \$HOME/.profile");
			}
			else{
			#	print $_."\n";
				system($_);
			}
#			or die "couldnt do $_\n";
		
		}
		chdir($curr_dir) or die "Couldnt change directory\n";
		#exit;
			print "Installed $prog\n\n";
		## REMOVE DOWNLOADED PACKAGE
			system("rm ".$progname_hash{$prog});
		}
#	die "\n\nPlease make sure these programs can be found in your PATH\n";
	#}
#	system("cat \$HOME/.profile");
#	system("source ".$ENV{'HOME'}."/.profile");
	system(". ~/.profile");
	
&print_start();
print "Download and Installation finished\nRe-read the your profile file by typing: source ~/.profile\n\n";
#`source \$HOME/.profile`;
}


END: {
 &write_path_variables($env_path, $estscandir_path, $wiseconfigdir_path);

}


####################################################################
sub write_path_variables{
####################################################################
my ($env_path, $estscandir_path, $wiseconfigdir_path) = @_;
## PATH
	print "save Path ";
#	print("echo export PATH=$env_path >> \$HOME/.profile\n");
	system("echo export PATH=$env_path >> \$HOME/.profile");
	
## ESTSCANDIR
unless(!defined $estscandir_path ){
#	print("echo export ESTSCANDIR=$estscandir_path >> \$HOME/.profile\n");
	system("echo export ESTSCANDIR=$estscandir_path >> \$HOME/.profile");
}
	## WISECONFIGDIR
unless(! defined $wiseconfigdir_path){
#	print("echo export WISECONFIGDIR=$wiseconfigdir_path >> \$HOME/.profile\n");
	system("echo export WISECONFIGDIR=$wiseconfigdir_path >> \$HOME/.profile");
}
#	print 	"source ".$ENV{'HOME'}."/.profile\n";
	#system("source ".$ENV{'HOME'}."/.profile");
#		`source \$HOME/.profile` || print "couldnt backtick\n";
#		exec("source \")
	system(". ~/.profile");

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
sub print_start(){
####################################################################
print "=" x $LINE_LENGTH, "\n";
print "-o0o-" x 12, "\n";
print "- " x 30, "\n";
print ">>==<<==" x 8, "\n";

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


=head1 NAME

prepare_analysis.pl - The script checks whether everything is set for the analysis 

=head1 SYNOPSIS


  prepare_analysis.pl
	[-t Configuration file]
	[-o Output file containing all taxa present]
	

=head1 DESCRIPTION

The script performs some simple checks whether everything is set
for the analysis to start.
	

=head2 An Example 

C<perl prepare_analysis.pl -t config.txt -o taxa.txt>

Reads 'config.txt', checks fasta files in 'fasta_directory' if valid fasta format,
checks if the required external programs are available and writes all taxa and
a shortened version of the taxa name in phylip-format (10 letters) to file 'taxa.txt'

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
