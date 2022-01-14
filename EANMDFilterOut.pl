#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2022
# Filter EANMD outCombined file v1.20 2022/01/14
# hukaining@gmail.com

use strict;
use warnings;
use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
use re 'eval';

our $opfn="filteredOut";
my $verbose;
our $startcodon="ATG";

#GetOptions("i=s"=>\$seqfilename,"o=s" => \$opfn,"verbose"=>\$verbose)
#or die("Error in command line arguments\n perl Searchhelitron -o <outputfile> (inputfile)\n");
GetOptions("o=s" => \$opfn, "s=s"=>\$startcodon, "verbose"  => \$verbose)
or die("[-]Error in command line arguments
  Usage: perl EANMDFilterOut.pl [options]  <input EANMD outCombined file>
    options:
    [-o string|outprefix Default: filteredOut]
    [-s string|Start_codon to remain. default: \"ATG\"]
    Note: Filter the records in EANMD outCombined file by start codon, MXE and CDS length. v1.20 2022/01/14.\n");


if (not @ARGV) {
	die ("[-] Error: Not find a input EANMD outCombined file.\n");
}
our $loadingstarttime=time();
print "Start loading EANMD outCombined file: @ARGV \n";
print "Start codon: $startcodon\n";
open OUT, "> $opfn.filter.$startcodon.outCombined.txt" or die ("[-] Error: Can't open or create $opfn.filter.$startcodon.outCombined.txt\n");
print OUT "QueryCol1\tSEUSDSCoordinates\tQueryCol3\tTranscript_id\tStrand\tExons\tStart_exon\tStop_exon\tSE_exon_Number\tSE(US)_Pos\tSE_length\tOri_CDS_length\tOri_Star_codon_to_exon_end_seq_len\trm/add_SE_start_to_end_seq_len";
print OUT "\tSEseq\tOri_CDSexons_seq\trm/add_SE_CDSexons_seq\tOri_last_junction_pos\tOri_last_dj\tOri_NMD\tStart_codon\tOri_AA\trm/add_SE_AA";
print OUT "\tAA_len+1\tOri_AA_1st_stop_pos\tOri_AA_stop_pos\tSEed_AA_1st_stop_pos\tSEed_AA_stop_pos\tFrame_shift_flag\tNew_1st_stop_pos_dj\tNMD_flag\tNMD_in/ex_flag\tsource\tSEupstreamCDS\tSEupstreamAApos\tUSexonNumber\tDSexonNumber\tExonNumbersBetweenUSDS\n";
our $count1=0;
our $count2=0;

while(defined(our $line = <>)){
    $count1++;
	my @tmp=split("\t",$line);
    my $Pos=$tmp[9];
    my $Start_codon=$tmp[20];
    my $source=$tmp[32];
    my $CDSlength=$tmp[11];
    my $oriAAlen=$tmp[24];
    my $inUSDSexonnumbers=$tmp[37];
        $inUSDSexonnumbers =~ s/\R//g;
    if ($Pos ne "5UTR" && $Start_codon !~ m/$startcodon/ig){ ### Rule ###2021.12.09 add i
        next;
    }elsif($Pos ne "5UTR" && $source eq "USSEDS" && $inUSDSexonnumbers >1){
        next;
    }elsif($Pos ne "5UTR" && $source eq "USDS" && $inUSDSexonnumbers >0){
        next;
    }elsif($Pos ne "5UTR" && $Pos ne "Start_codon" && (($oriAAlen-1)*3) ne $CDSlength){ ### Rule. remove annotation error 2022.01.14.
        next;
    }else{
        $count2++;
        print OUT "$line";
    }

}




our $loadingendtime=time();
print "Done!\nOutPut file: $opfn.filter.$startcodon.outCombined.txt\n";
print "Input $count1 record(s).\nRemain $count2 record(s).\n";
printf "Using %g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;