#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2022
# Transformation NMD unique flag results into GTF v1.000 2022/03/01
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

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=i"=>\$as_no,"I=i"=>\$uf_no,"e=i"=>\$exonSL,"E=i"=>\$exonEL,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "a=s"=>\$annot,"f=s"=>\$feature,"g=s"=>\$genome, "s=s"=>\$sepchr)
or die("[-]Error in command line arguments
  Usage: perl TransFlag2GTF.pl [options] <input AS flag file>
    options:
    [-o string|outprefix Default: TransOut]
    #[-a string|GTF annoation file]
    [-g string|Genome FASTA file]
    [-i int|AS events column Number. Default: 2]
    [-I int|UniqFlag column Number. Default: 23]
    [-e int|Exon start length to exact. Default: 50]
    [-E int|Exon end length to exact. Default: 50]
    [-u int|Exon start upstream length Default: 250]
    [-d int|Exon end downstream length Default: 250]
    [-s string|Delimiter of gene name and chromosome. Default: \"_chr\"]
    [-f string|Specify feature type in GTF annotation. Default: CDS]
	 
    Note: Transformation NMD unique flag results into GTF v1.000 2022/03/01\n");






# if ($genome && $annot){
#     print "Input genome: $genome\nInput GTF annotation: $annot\n";
if ($genome){
    print "Input genome: $genome\n";
    print "Upstream length: $upstreml\nDownstream length: $downstreml\n";
    print "Exon start length: $exonSL\nExon end length: $exonEL\n";

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
    our %Chrid2seq;
    #@ARGV = qw#''  Not_Find_a_File#;
    #say @ARGV;
    #say $0;
    while(defined(our $seq = <GENOME>)){

        if ($seq =~ m/^.*>/) {
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
    if ($annotcount % 100 == 0){
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

    $annotcount++;
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


}



######################################
#End  main
######################################
our $endtime=time();
#say $starttime;
#say $endtime;
print "$annotcount records.\n";
printf "Done! %g Sec %g Min\n",$endtime-$starttime,($endtime-$starttime)/60;