#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2022
# Transformation NMD unique flag results into GTF v1.100 2022/03/01
# hukaining@gmail.com

use strict;
use warnings;
#use 5.010;
use Getopt::Long;
use Pod::Usage;
use Time::HiRes 'time';
#use Parallel::ForkManager;
#our $MAX_processes=2;
#my $pm=Parallel::ForkManager->new($MAX_processes);
use re 'eval';
our $opfn="TransOut";
my $verbose;
our $as_no=2;
our $uf_no=23;
our $feature="CDS";
our $upstreml=250;
our $downstreml=250;
our $exonSL=50;
our $exonEL=50;
our $annot="";
our $genome="";
our $sepchr="_chr";
our %Chrid2seq;

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=i"=>\$as_no,"I=i"=>\$uf_no,"e=i"=>\$exonSL,"E=i"=>\$exonEL,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "a=s"=>\$annot,"f=s"=>\$feature,"g=s"=>\$genome, "s=s"=>\$sepchr)
or die("[-]Error in command line arguments
  Usage: perl TransFlag2GTF.pl [options] <input AS flag file>
    options:
    [-o string|outprefix Default: TransOut]
    [-a string|Annoation filter]
    [-g string|Genome FASTA file]
    [-i int|AS events column Number. Default: 2]
    [-I int|UniqFlag column Number. Default: 23]
    [-e int|Exon start length to exact. Default: 50]
    [-E int|Exon end length to exact. Default: 50]
    [-u int|Exon start upstream length Default: 250]
    [-d int|Exon end downstream length Default: 250]
    [-s string|Delimiter of gene name and chromosome. Default: \"_chr\"]
    [-f string|Specify feature type in GTF annotation. Default: CDS]
	 
    Note: Transformation NMD unique flag results into GTF v1.100 2022/03/01\n");



###################sub TRseq##########
sub TRseq($)
{
	my ($pinseq) = @_;
	my $pinseqtr = reverse $pinseq;
	$pinseqtr =~ tr/ACGTacgt/TGCAtgca/;
	return $pinseqtr;
}
##################TRseq End#############


# if ($genome && $annot){
#     print "Input genome: $genome\nInput GTF annotation: $annot\n";
if ($genome){
    print "Input genome: $genome\n";
    print "Upstream length: $upstreml\nDownstream length: $downstreml\n";
    print "Exon start length: $exonSL\nExon end length: $exonEL\n";

     #### Print ouputs
     ### US
    open OUTUSUP, "> $opfn.US.u$upstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.US.u$upstreml.fa \n");
    open OUTUSDOWN, "> $opfn.US.d$downstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.US.d$downstreml.fa \n");
    open OUTUSEXUP, "> $opfn.US.u$upstreml.es$exonSL.fa" or die ("[-] Error: Can't open or create file: $opfn.US.u$upstreml.es$exonSL.fa \n");
    open OUTUSEXDOWN, "> $opfn.US.d$downstreml.ed$exonEL.fa" or die ("[-] Error: Can't open or create file: $opfn.US.d$downstreml.ed$exonEL.fa \n");
     ### SE 
    open OUTSEUP, "> $opfn.SE.u$upstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.SE.u$upstreml.fa \n");
    open OUTSEDOWN, "> $opfn.SE.d$downstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.SE.d$downstreml.fa \n");
    open OUTSEEXUP, "> $opfn.SE.u$upstreml.es$exonSL.fa" or die ("[-] Error: Can't open or create file: $opfn.SE.u$upstreml.es$exonSL.fa \n");
    open OUTSEEXDOWN, "> $opfn.SE.d$downstreml.ed$exonEL.fa" or die ("[-] Error: Can't open or create file: $opfn.SE.d$downstreml.ed$exonEL.fa \n");

    ### DS 
    open OUTDSUP, "> $opfn.DS.u$upstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.DS.u$upstreml.fa \n");
    open OUTDSDOWN, "> $opfn.DS.d$downstreml.fa" or die ("[-] Error: Can't open or create file: $opfn.DS.d$downstreml.fa \n");
    open OUTDSEXUP, "> $opfn.DS.u$upstreml.es$exonSL.fa" or die ("[-] Error: Can't open or create file: $opfn.DS.u$upstreml.es$exonSL.fa \n");
    open OUTDSEXDOWN, "> $opfn.DS.d$downstreml.ed$exonEL.fa" or die ("[-] Error: Can't open or create file: $opfn.DS.d$downstreml.ed$exonEL.fa \n");
    

    ################
    # Loading Genome.
    ################
    # if (not $ARGV[0]) {
    #     die ("[-] Error: Not find a input Genome FASTA file.\n");
    # }
    our $loadingstarttime=time();
    # print @ARGV;
    # if (@ARGV eq "") {
    #     die ("[-] Error: Not find a input genome FASTA file.\n");
    #     }
    open GENOME, "< $genome" or die ("[-] Error: Can't open genome file: $genome.\n");
    print "Start loading genomeic sequences.\n";

    our $Chri=0;
    our @Chrname=();
    our @Chrseq=();
    # our %Chrid2seq;
    #@ARGV = qw#''  Not_Find_a_File#;
    #say @ARGV;
    #say $0;
    while(defined(our $seq = <GENOME>)){

        if ($seq =~ m/^.*>/){
        # $seq=~ m/^.*>([a-zA-Z0-9_.-]+) ?.*/;
        $seq=~ m/^\s*>\s*(\S+)\s*(.*)/;
        print "$1\n";
        $Chrname[$Chri]= $1;
        $Chri++;
        }else{
            $seq =~ s/\s//;
            $seq =~ tr/MRWSYKVHDBmrwsykvhdb/CGTGCGGCGGCGTGCGGCGG/;#1.31add snp replace
            #$Chrseq[$Chri-1] .=$seq;
            $Chrid2seq{$Chrname[$Chri-1]} .=$seq;
        }
    }

        # for (our $i=0;$i<$Chri-1;$i++){
        #     $Chrid2seq{$Chrname[$i]}=$Chrseq[$i];
        # }

    close GENOME;
    our $loadingendtime=time();
    print "$Chri Sequences\n";
    print "Finished loading!\n";
    printf "%g Sec %g Min\n",$loadingendtime-$loadingstarttime,($loadingendtime-$loadingstarttime)/60;

   



    # our $starttime=time();
    # print "Running. Please wait for a minite.\n";
}

#####################################
#Start main 
#####################################
print "GTF out Feature: $feature\n";
print "AS events column Number: $as_no\nUniqFlag column Number: $uf_no\n";
print "Delimiter of gene name and chromosome: $sepchr\n";

if (not $ARGV[0]) {
    die ("[-] Error: Not find a input AS_events-UniqueFlag file.\n");
}

our $starttime=time();
print "Running. Please wait for a minite.\n";
open OUTGTF, "> $opfn.gtf" or die ("[-] Error: Can't open or create file: $opfn.gtf \n");
print OUTGTF "#\n";

our $annotcount=0;
our @tmp;

while(defined(our $inrow = <>)){

    if ($inrow =~ m/^\#/) {next;}
    if ($inrow =~ m/^""/) {next;}

    $annotcount++;
    if ($annotcount % 500 == 0){
        print "Dealed with $annotcount results.\n";
    }

    @tmp = split(/\t/, $inrow);
    my $ASevents=$tmp[$as_no-1];
    my $flag=$tmp[$uf_no-1];
    $ASevents =~ s/"//g;
    $flag =~ s/"//g;
    $flag =~ s/\s//g;  ### Output value.
    my @AS_detail=split(/@/, $ASevents);
    # print "@AS_detail\t$flag\n";

    ########
    # split AS_events
    ########

    my @AS_name_a = split(/$sepchr/,$AS_detail[0]);
    # print "$AS_name_a[0]\n";
    my $Gene_name = $AS_name_a[0];
    my @SE=split(/:/, $AS_detail[0]);
    my $SES=$SE[1];
    my $SEE=$SE[2];
    # print "@SE\n";
    my @US=split(/:/, $AS_detail[1]);
    my @DS=split(/:/, $AS_detail[2]);
    # print "@DS\n";
    my $chr=$US[0];
    my $MP=$US[3];

    my $USS=$US[1];
    my $USE=$US[2];
    my $DSS=$DS[1];
    my $DSE=$DS[2];

    #fix bed 0-base:
    $SES++;
    $USS++;
    $DSS++;

    # print "$chr\n";
    # print "$MP\n";
    my $outlineUS = join("\t", $chr, "EANMD",$feature, $USS,$USE,".",$MP,".","gene_id \"$Gene_name\"; transcript_id \"$ASevents\"; gene_name \"$Gene_name\"; gene_source \"US\"; gene_biotype \"$flag\";\n");
    print OUTGTF "$outlineUS";
    my $outlineSE = join("\t", $chr, "EANMD",$feature, $SES,$SEE,".",$MP,".","gene_id \"$Gene_name\"; transcript_id \"$ASevents\"; gene_name \"$Gene_name\"; gene_source \"SE\"; gene_biotype \"$flag\";\n");
    print OUTGTF "$outlineSE";
    my $outlineDS = join("\t", $chr, "EANMD",$feature, $DSS,$DSE,".",$MP,".","gene_id \"$Gene_name\"; transcript_id \"$ASevents\"; gene_name \"$Gene_name\"; gene_source \"DS\"; gene_biotype \"$flag\";\n");
    print OUTGTF "$outlineDS";



    ###### Print sequences.
    if ($genome){

        my $seqname=$ASevents;

        if ($MP eq "+"){
            # + SE
            my $SEUPseq=substr($Chrid2seq{$chr}, $SES-$upstreml-1,$upstreml);
            my $SEDownseq=substr($Chrid2seq{$chr}, $SEE,$downstreml);
            my $SEUPEXseq=substr($Chrid2seq{$chr}, $SES-$upstreml-1,$upstreml+$exonSL);
            my $SEDownEXseq=substr($Chrid2seq{$chr}, $SEE-$exonEL,$downstreml+$exonEL);

            print OUTSEUP "\>$seqname|SEu$upstreml $chr:",$SES-$upstreml,"-",$SES-1,"\n$SEUPseq\n";
            print OUTSEEXUP "\>$seqname|SEu$upstreml.es$exonSL $chr:",$SES-$upstreml,"-",$SES-1+$exonSL,"\n$SEUPEXseq\n";
            print OUTSEDOWN "\>$seqname|SEd$downstreml $chr:",$SEE+1,"-",$SEE+$downstreml,"\n$SEDownseq\n";
            print OUTSEEXDOWN "\>$seqname|SEd$downstreml.ed$exonEL $chr:",$SEE+1-$exonEL,"-",$SEE+$downstreml,"\n$SEDownEXseq\n";

            # + US
            my $USUPseq=substr($Chrid2seq{$chr}, $USS-$upstreml-1,$upstreml);
            my $USDownseq=substr($Chrid2seq{$chr}, $USE,$downstreml);
            my $USUPEXseq=substr($Chrid2seq{$chr}, $USS-$upstreml-1,$upstreml+$exonSL);
            my $USDownEXseq=substr($Chrid2seq{$chr}, $USE-$exonEL,$downstreml+$exonEL);

            print OUTUSUP "\>$seqname|USu$upstreml $chr:",$USS-$upstreml,"-",$USS-1,"\n$USUPseq\n";
            print OUTUSEXUP "\>$seqname|USu$upstreml.es$exonSL $chr:",$USS-$upstreml,"-",$USS-1+$exonSL,"\n$USUPEXseq\n";
            print OUTUSDOWN "\>$seqname|USd$downstreml $chr:",$USE+1,"-",$USE+$downstreml,"\n$USDownseq\n";
            print OUTUSEXDOWN "\>$seqname|USd$downstreml.ed$exonEL $chr:",$USE+1-$exonEL,"-",$USE+$downstreml,"\n$USDownEXseq\n";

             # + DS
            my $DSUPseq=substr($Chrid2seq{$chr}, $DSS-$upstreml-1,$upstreml);
            my $DSDownseq=substr($Chrid2seq{$chr}, $DSE,$downstreml);
            my $DSUPEXseq=substr($Chrid2seq{$chr}, $DSS-$upstreml-1,$upstreml+$exonSL);
            my $DSDownEXseq=substr($Chrid2seq{$chr}, $DSE-$exonEL,$downstreml+$exonEL);

            print OUTDSUP "\>$seqname|DSu$upstreml $chr:",$DSS-$upstreml,"-",$DSS-1,"\n$DSUPseq\n";
            print OUTDSEXUP "\>$seqname|DSu$upstreml.es$exonSL $chr:",$DSS-$upstreml,"-",$DSS-1+$exonSL,"\n$DSUPEXseq\n";
            print OUTDSDOWN "\>$seqname|DSd$downstreml $chr:",$DSE+1,"-",$DSE+$downstreml,"\n$DSDownseq\n";
            print OUTDSEXDOWN "\>$seqname|DSd$downstreml.ed$exonEL $chr:",$DSE+1-$exonEL,"-",$DSE+$downstreml,"\n$DSDownEXseq\n";



        }elsif($MP eq "-"){

            # + SE switch up down
            # my $SEUPseq=TRseq(substr($Chrid2seq{$chr}, $SES-$upstreml-1,$upstreml));
            # my $SEDownseq=TRseq(substr($Chrid2seq{$chr}, $SEE,$downstreml));
            # my $SEUPEXseq=TRseq(substr($Chrid2seq{$chr}, $SES-$upstreml-1,$upstreml+$exonSL));
            # my $SEDownEXseq=TRseq(substr($Chrid2seq{$chr}, $SEE-$exonEL,$downstreml+$exonEL));

            # print OUTSEUP "\>$seqname|SEu$upstreml $chr:",$SES-$upstreml,"-",$SES-1,"\n$SEUPseq\n";
            # print OUTSEEXUP "\>$seqname|SEu$upstreml.es$exonSL $chr:",$SES-$upstreml,"-",$SES-1+$exonSL,"\n$SEUPEXseq\n";
            # print OUTSEDOWN "\>$seqname|SEd$downstreml $chr:",$SEE+1,"-",$SEE+$downstreml,"\n$SEDownseq\n";
            # print OUTSEEXDOWN "\>$seqname|SEd$downstreml.ed$exonEL $chr:",$SEE+1-$exonEL,"-",$SEE+$downstreml,"\n$SEDownEXseq\n";

            # + SE switch up down
            my $SEDownseq=TRseq(substr($Chrid2seq{$chr}, $SES-$downstreml-1,$downstreml));
            my $SEupseq=TRseq(substr($Chrid2seq{$chr}, $SEE,$upstreml));
            my $SEDownEXseq=TRseq(substr($Chrid2seq{$chr}, $SES-$downstreml-1,$downstreml+$exonEL));
            my $SEUPEXseq=TRseq(substr($Chrid2seq{$chr}, $SEE-$exonSL,$downstreml+$exonSL));

            print OUTSEDOWN "\>$seqname|SEd$downstreml $chr:",$SES-$downstreml,"-",$SES-1,"\n$SEDownseq\n";
            print OUTSEEXDOWN "\>$seqname|SEd$downstreml.ed$exonEL $chr:",$SES-$downstreml,"-",$SES-1+$exonEL,"\n$SEDownEXseq\n";
            print OUTSEUP "\>$seqname|SEu$upstreml $chr:",$SEE+1,"-",$SEE+$upstreml,"\n$SEupseq\n";
            print OUTSEEXUP "\>$seqname|SEu$upstreml.es$exonSL $chr:",$SEE+1-$exonSL,"-",$SEE+$upstreml,"\n$SEUPEXseq\n";

            # + UP DS switch up down ,switch USS <-> DSS, USE <-> DSE

            my $tmpS=$DSS;
            my $tmpE=$DSE;
            $DSS=$USS;
            $DSE=$USE;
            $USS=$tmpS;
            $USE=$tmpE;

            my $USDownseq=TRseq(substr($Chrid2seq{$chr}, $USS-$downstreml-1,$downstreml));
            my $USupseq=TRseq(substr($Chrid2seq{$chr}, $USE,$upstreml));
            my $USDownEXseq=TRseq(substr($Chrid2seq{$chr}, $USS-$downstreml-1,$downstreml+$exonEL));
            my $USUPEXseq=TRseq(substr($Chrid2seq{$chr}, $USE-$exonSL,$downstreml+$exonSL));

            print OUTUSDOWN "\>$seqname|USd$downstreml $chr:",$USS-$downstreml,"-",$USS-1,"\n$USDownseq\n";
            print OUTUSEXDOWN "\>$seqname|USd$downstreml.ed$exonEL $chr:",$USS-$downstreml,"-",$USS-1+$exonEL,"\n$USDownEXseq\n";
            print OUTUSUP "\>$seqname|USu$upstreml $chr:",$USE+1,"-",$USE+$upstreml,"\n$USupseq\n";
            print OUTUSEXUP "\>$seqname|USu$upstreml.es$exonSL $chr:",$USE+1-$exonSL,"-",$USE+$upstreml,"\n$USUPEXseq\n";


            #### DS
            my $DSDownseq=TRseq(substr($Chrid2seq{$chr}, $DSS-$downstreml-1,$downstreml));
            my $DSupseq=TRseq(substr($Chrid2seq{$chr}, $DSE,$upstreml));
            my $DSDownEXseq=TRseq(substr($Chrid2seq{$chr}, $DSS-$downstreml-1,$downstreml+$exonEL));
            my $DSUPEXseq=TRseq(substr($Chrid2seq{$chr}, $DSE-$exonSL,$downstreml+$exonSL));

            print OUTDSDOWN "\>$seqname|DSd$downstreml $chr:",$DSS-$downstreml,"-",$DSS-1,"\n$DSDownseq\n";
            print OUTDSEXDOWN "\>$seqname|DSd$downstreml.ed$exonEL $chr:",$DSS-$downstreml,"-",$DSS-1+$exonEL,"\n$DSDownEXseq\n";
            print OUTDSUP "\>$seqname|DSu$upstreml $chr:",$DSE+1,"-",$DSE+$upstreml,"\n$DSupseq\n";
            print OUTDSEXUP "\>$seqname|DSu$upstreml.es$exonSL $chr:",$DSE+1-$exonSL,"-",$DSE+$upstreml,"\n$DSUPEXseq\n";






        }

    }


}



######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
print "$annotcount records.\n";
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;