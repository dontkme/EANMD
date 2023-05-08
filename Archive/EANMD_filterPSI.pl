#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2022
# EANMD Filter the rMATS results and output an EANMD AS_events input list. EANMD_fileterPSI v1.250 2023/01/10
# hukaining@gmail.com
#
#use 5.0100;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max sum /;

sub mean {
    return sum(@_)/@_;
}

use re 'eval';

our $opfn="";
# our $opnmdinput="";
my $verbose;
our $FDR=0.05;
our $DeltaPSIcutoff=0.15;
our $MXEfold=0.05;
# our $Ainfo;
our $mindepth=20;
our $mincount=2;
# our $filterpercent=0.3;
# our $samplenum=4;
our $compares1=2;
our $compares2=2;



GetOptions("o=s" => \$opfn,"d=i"=>\$mindepth,"m=i"=>\$mincount,"i=f"=>\$DeltaPSIcutoff,"f=f"=>\$FDR,"mf=f"=>\$MXEfold,"c1=i"=>\$compares1,"c2=i"=>\$compares2,"verbose"=>\$verbose)
	or die("[-]Error in command line arguments\n    Filter PSI v1.25 2023/01/10\nUsage: perl EANMD_filterPSI.pl [options] <input rmats result>\n
  options:\n
	 [-o output prefix. default: rMATS_filtered.out]\n
	 [-d int|min depth of average read counts. default: 20]\n
     [-m int|min count of inclusion events' UP|Downstream junction. default: 2]\n
	 [-i float|delta PSI cutoff [0-1.0]. default: 0.15]\n
	 [-f float|FDR cutoff [0-1.0]. default: 0.05]\n
	 [-c1 int|The first sample numbers.  default: 2]\n
	 [-c2 int|The second sample numbers. default: 2]\n
	 [-mf float|The US and DS fold change cutoff. default: 0.05]\n
  Note: Filter the rMATS results and output an EANMD AS_events input list. \n");
	#open IN,"chrA09.snp.vcf";
if ($opfn eq ""){
    $opfn="rMATS_filtered.out";
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}else{
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.txt\n");
open OUTEAINPUT, "> $opfn.EANMD.input.txt" or die ("[-] Error: Can't open or creat $opfn.EANMD.input.txt\n");
# print "Input $samplenum samlpes.";
# if ($compares1<=$samplenum and $compares2<=$samplenum){
#   print "Compare sample$compares1 to sample$compares2.\n";
#   }else{
#     die ("[-]Error compare Sample ID.\nPlease confirm to compare sample$compares1 to sample$compares2?\n");
#   } 
print "Min depth of averge read counts: $mindepth\ndelta PSI cutoff: $DeltaPSIcutoff\nFDR cutoff: $FDR\nThe first sample numbers: $compares1\nThe second sample numbers: $compares2\nMin count of inclusion events: $mincount\nThe US and DS fold change cutoff: $MXEfold\n";
  
print OUT "ID\tGeneID\tgeneSymbol\tchr\tstrand\texonStart_0base\texonEnd\tupstreamES\tupstreamEE\tdownstreamES\tdownstreamEE\tID";
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
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tupstream_to_target_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tupstream_to_target_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\ttarget_to_downstream_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\ttarget_to_downstream_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\ttarget_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\ttarget_count2_$i";
}
for (my $i=1;$i<=$compares1;$i++){
    print OUT "\tupstream_to_downstream_count1_$i";
}
for (my $i=1;$i<=$compares2;$i++){
    print OUT "\tupstream_to_downstream_count2_$i";
}
print OUT "\tMin_average_counts\tMinJCofInclusion\tUp_Down_JC_fold\n";


our $count1=0;
our $count2=0;

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

    #### Start to filter
    ### 1. FDR cutoff $FDR.

    if ($tmpFDR > $FDR){
        next;
    }

    # 2. Min average read count has to greater than $mindepth .
    ## 2. Max average read count has to greater than $mindepth .


    my $meanS1reads = ((sum @IJC_S1) + (sum @SJC_S1))/$compares1;
    my $meanS2reads = ((sum @IJC_S2) + (sum @SJC_S2))/$compares2;
    my $minavgreads = min($meanS1reads, $meanS2reads);
    my $maxavgreads = max($meanS1reads, $meanS2reads);
    if ($minavgreads < $mindepth){
        next;
    }
    # if ($maxavgreads < $mindepth){
    #     next;
    # }
    
    ### 3. PSI cutoff. $DeltaPSIcutoff

    if (abs($tmpPSI)<$DeltaPSIcutoff){
        next;
    }

    ## 4. Inclusion events' min US|DS JC greater equal to $mincount.
    # 4. Inclusion events' max US|DS JC greater equal to $mincount.
    my $tmpmincount=0;
    if ($tmpPSI>0){
        # $tmpmincount = min(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
        $tmpmincount = max(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
    }else{
        # $tmpmincount = min(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]));
        $tmpmincount = max(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]));
    }
    
    if ($tmpmincount<$mincount){
        next;
    }

    ### 5. Inclusion events' min (US|DS JC)/max(US|DS JC) greater equal to $MXEfold.
    my $tmpMXEfold=0;
    if ($tmpPSI>0){
        $tmpMXEfold = min(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]))/max(mean(@US2SE[0..$compares1-1]),mean(@SE2DS[0..$compares1-1]));
    }else{
        $tmpMXEfold = min(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]))/max(mean(@US2SE[$compares1..($compares1+$compares2-1)]),mean(@SE2DS[$compares1..($compares1+$compares2-1)]))
    }

    if ($tmpMXEfold< $MXEfold){
        next;
    }

    ########## Output
    $count2++;
    my $res1=join("\t",@col[0..11],@IJC_S1,@SJC_S1,@IJC_S2,@SJC_S2,@col[16..19],@IncLevel1,@IncLevel2,$tmpPSI,@US2SE,@SE2DS,@SEcount,@US2DS,$minavgreads,$tmpmincount,$tmpMXEfold);
    # my $resIJCS1=join("\t",@IJC_S1);
    # my $resSJCS1=join("\t",@SJC_S1);
    # my $resIJCS2=join("\t",@IJC_S2);
    # my $resSJCS2=join("\t",@SJC_S2);
    # my $res2=join();

    # print OUT "$res1\t$resIJCS1\t$resSJCS1\t$resIJCS2\t$resSJCS2\t\n";
    print OUT "$res1\n";
    $tmpgeneSymbol =~ s/"//g; ### 2021-09-09 rm double quotes.
    $tmpGeneID =~ s/"//g;
    print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$tmpSEstart:$tmpSEend:$tmpstrand"."@"."$tmpchr:$tmpUSstart:$tmpUSend:$tmpstrand"."@"."$tmpchr:$tmpDSstart:$tmpDSend:$tmpstrand\t$tmpGeneID\n";
	
  


	
}

close OUT;
print "All $count1 AS events. Rest $count2 AS events.\n";