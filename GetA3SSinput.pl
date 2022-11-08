#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2022
# EANMD convert rMATS results to an EANMD AS_events input list. GetA3SSinput v1.000 2022/11/08
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
my $verbose;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose)
	or die("[-]Error in command line arguments\nConvert A3SS rMATS to EANMD input, GetA3SSinput v1.00 2022/11/08\nUsage: perl EANMD_filterPSI.pl [options] <input rmats result>\n
  options:\n
	 [-o output prefix. default: rMATS_filtered.out]\n
  Note: convert the rMATS results to an EANMD AS_events input list. \n");

if ($opfn eq ""){
    $opfn="rMATS_filtered.out";
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}else{
    print "Output file: $opfn.txt $opfn.EANMD.input.txt\n";
}

open OUT, "> $opfn.txt" or die ("[-] Error: Can't open or creat $opfn.txt\n");
open OUTEAINPUT, "> $opfn.EANMD.input.txt" or die ("[-] Error: Can't open or creat $opfn.EANMD.input.txt\n");
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

    #### Start to filter
    ### 1. FDR cutoff $FDR.

    # if ($tmpFDR > $FDR){
    #     next;
    # }

    # 2. Min average read count has to greater than $mindepth .
    ## 2. Max average read count has to greater than $mindepth .


    # my $meanS1reads = ((sum @IJC_S1) + (sum @SJC_S1))/$compares1;
    # my $meanS2reads = ((sum @IJC_S2) + (sum @SJC_S2))/$compares2;
    # my $minavgreads = min($meanS1reads, $meanS2reads);
    # my $maxavgreads = max($meanS1reads, $meanS2reads);
    # if ($minavgreads < $mindepth){
    #     next;
    # }
    # # if ($maxavgreads < $mindepth){
    # #     next;
    # # }
    
    # ### 3. PSI cutoff. $DeltaPSIcutoff

    # if (abs($tmpPSI)<$DeltaPSIcutoff){
    #     next;
    # }

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
    # # my $resIJCS1=join("\t",@IJC_S1);
    # # my $resSJCS1=join("\t",@SJC_S1);
    # # my $resIJCS2=join("\t",@IJC_S2);
    # # my $resSJCS2=join("\t",@SJC_S2);
    # # my $res2=join();

    # # print OUT "$res1\t$resIJCS1\t$resSJCS1\t$resIJCS2\t$resSJCS2\t\n";
    # print OUT "$res1\n";
    $tmpgeneSymbol =~ s/"//g; ### 2021-09-09 rm double quotes.
    $tmpGeneID =~ s/"//g;
    ## A3SS 
    
    my $A3SSES;
    my $A3SSEE; 

    my $flankingES;
    my $flankingEE;

    my $shortES;
    my $shortEE;

    if ($tmpstrand eq "+"){

        $A3SSES = $tmpSEstart;
        $A3SSEE = $tmpUSstart;

        $flankingES = $tmpDSstart;
        $flankingEE = $tmpDSend;

        $shortES = $tmpUSstart;
        $shortEE = $tmpUSend;
    
    }else{
        
        $A3SSES = $tmpUSend;
        $A3SSEE = $tmpSEend;

        $flankingES = $tmpUSstart;
        $flankingEE = $tmpUSend;

        $shortES = $tmpDSstart;
        $shortEE = $tmpDSend;
    }
    # print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$tmpSEstart:$tmpSEend:$tmpstrand"."@"."$tmpchr:$tmpUSstart:$tmpUSend:$tmpstrand"."@"."$tmpchr:$tmpDSstart:$tmpDSend:$tmpstrand\t$tmpGeneID\n";
    print OUTEAINPUT "$tmpgeneSymbol\t$tmpchr:$A3SSES:$A3SSEE:$tmpstrand"."@"."$tmpchr:$flankingES:$flankingEE:$tmpstrand"."@"."$tmpchr:$shortES:$shortEE:$tmpstrand\t$tmpGeneID\n";
	
  


	
}

close OUT;
print "All $count1 AS events. Rest $count2 AS events.\n";
# print " Rest $count2 AS events.\n";