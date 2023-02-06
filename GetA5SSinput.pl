#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2023
# EANMD convert rMATS results to an EANMD AS_events input list. GetA5SSinput v1.100 2023/02/06
# hukaining@gmail.com
#
#use 5.0100;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
# use List::Util qw/ min max sum /;

our $opfn = "";
# our $opnmdinput="";
our $FDR=0.05;
our $DeltaPSIcutoff=0.15;
our $compares1=2;
our $compares2=2;
my $verbose;

GetOptions("o=s" => \$opfn,"i=f"=>\$DeltaPSIcutoff,"f=f"=>\$FDR,"c1=i"=>\$compares1,"c2=i"=>\$compares2,"verbose"=>\$verbose)
	or die("[-]Error in command line arguments\nConvert A5SS rMATS to EANMD input, GetA5SSinput v1.10 2023/02/06\nUsage: perl GetA5SSinput.pl [options] <input rmats A5SS result>\n
  options:\n
	 [-o output prefix | default: rMATS_filtered.out]\n
     [-i float|delta PSI cutoff [0-1.0]. default: 0.15]\n
	 [-f float|FDR cutoff [0-1.0]. default: 0.05]\n
	 [-c1 int|The first sample numbers.  default: 2]\n
	 [-c2 int|The second sample numbers. default: 2]\n
  Note: convert the rMATS results to an EANMD AS_events input list. \n");

if ($opfn eq ""){
    $opfn="rMATS_filtered.out";
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}else{
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.txt\n");
# print "Min depth of averge read counts: $mindepth\ndelta PSI cutoff: $DeltaPSIcutoff\nFDR cutoff: $FDR\nThe first sample numbers: $compares1\nThe second sample numbers: $compares2\n"; #Min count of inclusion events: $mincount\nThe US and DS fold change cutoff: $MXEfold\n";
print "Delta PSI cutoff: $DeltaPSIcutoff\nFDR cutoff: $FDR\nThe first sample numbers: $compares1\nThe second sample numbers: $compares2\n"; #Min count of inclusion events: $mincount\nThe US and DS fold change cutoff: $MXEfold\n";
 
# print OUT "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tID";
print OUT "ID\tGeneID\tgeneSymbol\tchr\tstrand\tlongExonStart_0base\tlongExonEnd\tshortES\tshortEE\tflankingES\tflankingEE\tID";
# longExonStart_0base	longExonEnd	shortES	shortEE	flankingES	flankingEE
open OUTEAINPUT, "> $opfn.EANMD.input.txt" or die ("[-] Error: Can't open or creat $opfn.EANMD.input.txt\n");

for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tIJC_SAMPLE_1_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tSJC_SAMPLE_1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tIJC_SAMPLE_2_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tSJC_SAMPLE_2_$i";
}
print OUT "\tIncFormLen\tSkipFormLen\tPValue\tFDR";
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tIncLevel1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tIncLevel2_$i";
}
print OUT "\tIncLevelDifference";
# across_short_boundary_count	long_to_flanking_count	exclusive_to_long_count	short_to_flanking_count
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tacross_short_boundary_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tacross_short_boundary_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tlong_to_flanking_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tlong_to_flanking_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\texclusive_to_long_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\texclusive_to_long_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tshort_to_flanking_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tshort_to_flanking_count2_$i";
}
# print OUT "\tMin_average_counts\tMinJCofInclusion\tUp_Down_JC_fold\n";
print OUT "\n";

our $count1 = 0;
our $count2 = 0;
LINE: while(our $row = <>){	
#while(our $seq = <SEQFILENAME>){
	#chomp $row;
	# if ($row =~ m/^\#/) {next;}
	
	my @col =split(/\t/,$row);
	my $tmpID = $col[0];
    if ($col[0] eq "ID") {next;}
        $count1++;
	my $tmpGeneID = $col[1];
	my $tmpgeneSymbol = $col[2];
	my $tmpchr = $col[3];
	my $tmpstrand = $col[4];
    my $tmpSEstart = $col[5];
    my $tmpSEend = $col[6];
    my $tmpUSstart = $col[7];
    my $tmpUSend = $col[8];
    my $tmpDSstart = $col[9];
    my $tmpDSend = $col[10];
    
    my @IJC_S1 = split(/,/,$col[12]); # split IJC_sampl_1
    my @SJC_S1 = split(/,/,$col[13]); # split SJC_sample_1

    my @IJC_S2 = split(/,/,$col[14]);
    my @SJC_S2 = split(/,/,$col[15]);

    my $IncFormLen = $col[16];
    my $SkipFormLen = $col[17];
    my $tmppvalue = $col[18];
    my $tmpFDR = $col[19]; #FDR
 
    my @IncLevel1 = split(/,/,$col[20]);
    my @IncLevel2 = split(/,/,$col[21]);
    my $tmpPSI = $col[22]; 

    my @US2SE = split(/,/,$col[23]); #individual Up Stream Junction Count.
    my @SE2DS = split(/,/,$col[24]); #individual Down Stream Junction Count.
    my @SEcount = split(/,/,$col[25]); #individual SE Count.
    my @US2DS = split(/,/,$col[25]); #individual Up Stream to Down Stream Junction Count (Skipped Junction Count).

### Start to filter
    ## 1. FDR cutoff $FDR.

    if ($tmpFDR > $FDR){
        next;
    }

    # 2. Min average read count has to greater than $mindepth .
    # # 2. Max average read count has to greater than $mindepth .


    # my $meanS1reads = ((sum @IJC_S1) + (sum @SJC_S1))/$compares1;
    # my $meanS2reads = ((sum @IJC_S2) + (sum @SJC_S2))/$compares2;
    # my $minavgreads = min($meanS1reads, $meanS2reads);
    # my $maxavgreads = max($meanS1reads, $meanS2reads);
    # if ($minavgreads < $mindepth){
    #     next;
    # }
    # if ($maxavgreads < $mindepth){
    #     next;
    # }
    
    ### 3. PSI cutoff. $DeltaPSIcutoff

    if (abs($tmpPSI)<$DeltaPSIcutoff){
        next;
    }

    # ## 4. Inclusion events' min US|DS JC greater equal to $mincount.
    # # 4. Inclusion events' max US|DS JC greater equal to $mincount.
    # my $tmpmincount=0;
    # if ($tmpPSI>0){
    #     # $tmpmincount = min(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
    #     $tmpmincount = max(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
    # }else{
    #     # $tmpmincount = min(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]));
    #     $tmpmincount = max(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]));
    # }
    
    # if ($tmpmincount<$mincount){
    #     next;
    # }

    # ### 5. Inclusion events' min (US|DS JC)/max(US|DS JC) greater equal to $MXEfold.
    # my $tmpMXEfold=0;
    # if ($tmpPSI>0){
    #     $tmpMXEfold = min(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]))/max(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
    # }else{
    #     $tmpMXEfold = min(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]))/max(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]))
    # }

    # if ($tmpMXEfold< $MXEfold){
    #     next;
    # }

    ########## Output
    $count2++;
    # my $res1=join("\t",@col[0..11],@IJC_S1,@SJC_S1,@IJC_S2,@SJC_S2,@col[16..19],@IncLevel1,@IncLevel2,$tmpPSI,@US2SE,@SE2DS,@SEcount,@US2DS,$minavgreads,$tmpmincount,$tmpMXEfold);
    my $res1=join("\t",@col[0..11],@IJC_S1,@SJC_S1,@IJC_S2,@SJC_S2,@col[16..19],@IncLevel1,@IncLevel2,$tmpPSI,@US2SE,@SE2DS,@SEcount,@US2DS);

    # my $resIJCS1=join("\t",@IJC_S1);
    # my $resSJCS1=join("\t",@SJC_S1);
    # my $resIJCS2=join("\t",@IJC_S2);
    # my $resSJCS2=join("\t",@SJC_S2);
    # my $res2=join();

    # print OUT "$res1\t$resIJCS1\t$resSJCS1\t$resIJCS2\t$resSJCS2\t\n";
    print OUT "$res1\n";
    $tmpgeneSymbol =~ s/"//g; ### 2021-09-09 rm double quotes.
    $tmpGeneID =~ s/"//g;
    
  
    ## A5SS 
    # ES: Exon Start; EE: Exon End.
    my $A5SSES;
    my $A5SSEE; 

    my $flankingES;
    my $flankingEE;

    my $shortES;
    my $shortEE;

    if ($tmpstrand eq "+"){ #A5SS diff to A3SS; SE: Long; US: Short; DS: Flanking;

        # $A5SSES = $tmpSEstart;
        # $A5SSEE = $tmpUSstart;

        # $flankingES = $tmpDSstart;
        # $flankingEE = $tmpDSend;

        # $shortES = $tmpUSstart;
        # $shortEE = $tmpUSend;

        $A5SSES = $tmpUSend;
        $A5SSEE = $tmpSEend;

        $flankingES = $tmpDSstart;
        $flankingEE = $tmpDSend;

        $shortES = $tmpUSstart;
        $shortEE = $tmpUSend;
    
    }else{
        
        # $A5SSES = $tmpUSend;
        # $A5SSEE = $tmpSEend;

        # $flankingES = $tmpUSstart;
        # $flankingEE = $tmpUSend;

        # $shortES = $tmpDSstart;
        # $shortEE = $tmpDSend;

        $A5SSES = $tmpSEstart;
        $A5SSEE = $tmpUSstart;

        $flankingES = $tmpUSstart;
        $flankingEE = $tmpUSend;

        $shortES = $tmpDSstart;
        $shortEE = $tmpDSend;
    }
    # print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$tmpSEstart:$tmpSEend:$tmpstrand"."@"."$tmpchr:$tmpUSstart:$tmpUSend:$tmpstrand"."@"."$tmpchr:$tmpDSstart:$tmpDSend:$tmpstrand\t$tmpGeneID\n";
    # print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$A5SSES:$A5SSEE:$tmpstrand"."@"."$tmpchr:$flankingES:$flankingEE:$tmpstrand"."@"."$tmpchr:$shortES:$shortEE:$tmpstrand\t$tmpGeneID\n";
    print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$A5SSES:$A5SSEE:$tmpstrand"."@"."$tmpchr:$shortES:$shortEE:$tmpstrand"."@"."$tmpchr:$flankingES:$flankingEE:$tmpstrand\t$tmpGeneID\n";
	
  


	
}

close OUT;
print "All $count1 AS events. Rest $count2 AS events.\n";
# print " Rest $count2 AS events.\n";