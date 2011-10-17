
#!/usr/bin/perl
### LIBRARIES{{{
use Getopt::Std;
use File::Basename;
use Bio::Perl;
use Data::Dumper;
use Pod::Usage;
use Bio::SearchIO;
use Bio::Tools::BPlite;
use Bio::Tools::Run::StandAloneBlast;
use Bio::AlignIO;
#use Tie::Hash::Regex;
use strict;
use warnings;
#use fastastuff;
my %options;
getopt( q(iost), \%options );



my $infile = $options{i} || die "Please enter input file\n";
my $outfile = $options{o} ||  "out_95.fa";
my $seq_type = $options{s} || "p";
my $threshold = $options{t} || "95";

my $initial_number_of_files = 0;
my %save_seq_hash = ();

my $in  = Bio::SeqIO->new(-file => $infile ,
                           -format => 'Fasta');
FILE_ITERATION:
    while ( my $seq = $in->next_seq() ) {
  		if(!defined $seq->id || !$seq->seq()){
  	#		warn "\t\t".$seq->id." has no sequence\n";
  			next FILE_ITERATION;
  		}
    	next FILE_ITERATION if length($seq->seq) < 100;
  #  	print "Saving ".$seq->id."\n";
  	$seq->id =~ /(\w*_\w*)_(\w*)/;
  		my ($taxon,$id) = ($1,$2);
  		if(!defined $taxon || !defined $id){
  			warn "\t\tProblem with ".$seq->id."\n";
  		}
  #  	if(exists $save_seq_hash{$taxon}{$id}){
  # 		warn "\t\tEntry ".$seq->id." already exists\n";
  #  		next FILE_ITERATION;
  #  	}
    #	print Dumper %save_seq_hash;
    	my $bl2seq_query = Bio::Seq->new(-seq        => $seq->seq, -display_id => "query", -alphabet   => 'protein');

    	## SEARCH FOR EQUAL MATCHES
    	foreach my $curr_id(keys %{$save_seq_hash{$taxon}}){
    #		print "\t\t$save_seq_hash{$taxon}{$curr_id}\n";
    		if(exists $save_seq_hash{$taxon}{$curr_id} && $save_seq_hash{$taxon}{$curr_id} eq $seq->seq){
	    #		warn "\t\tEntry ".$seq->id." is the same as ".$taxon."_$curr_id\n";
	    	#	print "\t\t$save_seq_hash{$taxon}{$_}\n\t\t".$seq->seq."\n";
    			next FILE_ITERATION;
    		}
    	my $bl2seq_target = Bio::Seq->new(-seq        =>$save_seq_hash{$taxon}{$curr_id} , -display_id => "query", -alphabet   => 'protein');
 		my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastp',
                                               'outfile' => "naja");
	if($seq_type eq "n"){
		my $factory = Bio::Tools::Run::StandAloneBlast->new('program' => 'blastn',
                                               'outfile' => "naja");
	}
		   #E-value
			 $factory->e(1e-10);
		  # Turn filtering off
		  $factory->F("F");
	#		print "\t\tBl2seq: ".$seq->id." (Query) --> $taxon _".$curr_id." (Target)\n";
		   my $bl2seq_report = $factory->bl2seq($bl2seq_query, $bl2seq_target);
			## RETURN DEFAULT IF NO HITS
		   next if !$bl2seq_report;
			### PARSE BLAST RESULTS
		    while ( my $result = $bl2seq_report->next_result) {
		        while ( my $hit = $result->next_hit) {
#					print $hit->significance."\t".$hit->identity."\n";
		#			print "E-Value: ".$hit->significance."\t";
					#print "Percent: ".$hit->bits."\n";
				#	print Dumper $hit;
					if ( my $hsp = $hit->next_hsp ) {
		#				print "\t\t".$hsp->num_identical."/".$hsp->length." --> (".($hsp->num_identical / $hsp->length * 100)." %)\n";;
#						print "\t\t".$hsp->frac_identical{total};
#						my $tmp = $hsp->frac_identical;
	#					print $tmp->total
				if((($hsp->num_identical / $hsp->length) * 100) > $threshold){
	#				print "\t".$bl2seq_query->seq."\n\t".$bl2seq_target->seq."\n";
		#			print "\tSelect longer sequence: ";
					if(length($bl2seq_query->seq) > length($bl2seq_target->seq)){
			#			print length($bl2seq_query->seq)." > ".length($bl2seq_target->seq)."\n";
						delete $save_seq_hash{$taxon}{$curr_id};
					}
					else{
				#		print length($bl2seq_query->seq)." < ".length($bl2seq_target->seq)."\n";
						next FILE_ITERATION;
					}
				}
			}

    
		#			exit;
				}
		    }
}
 	  	$save_seq_hash{$taxon}{$id} = $seq->seq;
	  	$initial_number_of_files++;
#	  	print Dumper %save_seq_hash;
#	    $out->write_seq($seq);
    }
    my %ids_to_keep = ();
    my $end_counter = 0;
    foreach my $tax(sort keys %save_seq_hash){
  #  print "\t\t$tax\n";
#	print Dumper $save_seq_hash{$tax};
	    foreach(keys %{$save_seq_hash{$tax}}){
	#		print $tax."_".$_."\n";
			$ids_to_keep{$tax."_".$_} = 1;
	    	$end_counter++;
	}
}

 $in  = Bio::SeqIO->new(-file => $infile ,
                           -format => 'Fasta');
 my $out = Bio::SeqIO->new(-file => ">$outfile" ,
                           -format => 'Fasta');
    while ( my $seq = $in->next_seq() ) {
 	   $out->write_seq($seq) if exists $ids_to_keep{$seq->id};
    }
unlink("naja");
	#  	print Dumper %save_seq_hash;
#print "Reduction of  $end_counter/$initial_number_of_files\n";


