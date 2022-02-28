#!/usr/bin/env perl

#AUTHORS
# Kaining Hu (c) 2022
# Transformation NMD unique flag results into GTF v0.1000 2022/02/28
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
our $annot="";
our $genome="";
our $sepchr="_chr";

GetOptions("o=s" => \$opfn,"verbose"=>\$verbose,"i=i"=>\$as_no,"I=i"=>\$uf_no,"u=i"=>\$upstreml,"d=i"=>\$downstreml, "a=s"=>\$annot,"f=s"=>\$feature,"g=s"=>\$genome, "s=s"=>\$sepchr)
or die("[-]Error in command line arguments
  Usage: perl TransFlag2GTF.pl [options] <input AS flag file>
    options:
    [-o string|outprefix Default: TransOut]
    [-a string|GTF annoation file]
    [-g string|Genome FASTA file]
    [-i int|AS events column Number. Default: 2]
    [-I int|UniqFlag column Number. Default: 23]
    [-s string|Delimiter of gene name and chromosome. Default: \"_chr\"]
    [-f string|Specify feature type in GFF annotation.default: CDS]
    [-u int|upstream length Default: 250]
    [-d int|downstream length Default: 250]
	 
    Note: Transformation NMD unique flag results into GTF v0.1000 2022/02/28\n");



if ($genome && $annot){
    print "Input genome: $genome\nInput GTF annotation: $annot\n";
}
