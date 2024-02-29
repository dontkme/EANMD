#!/usr/bin/perl

#AUTHORS
# Kaining Hu (c) 2022
# Convert SUPPA2 ioe file to EANMD AS_events input list. TransIOE2IRinput.pl v0.100 2024/02/27
# hukaining@gmail.com
#
use 5.0100;
use strict;
use warnings;
use Getopt::Long;
use Pod::Usage;
use re 'eval';

our $opfn="";

# our $exlen = 100;

# GetOptions("o=s" => \$opfn,"l=i"=>\$exlen)
GetOptions("o=s" => \$opfn)
	or die("[-]Error in command line arguments\n    onvert SUPPA2 ioe file to EANMD AS_events input list v0.100 2024/02/27\n
    Usage: perl TransIOE2IRinput.pl [options] <input IR ioe>\n
  Options:\n
	 [-o output prefix. default: trans2AS_events.out]\n
     
  Note: Convert SUPPA2 ioe file to EANMD AS_events input list. \n");


# if ($exlen<0) { 
#     die ("[-] Error: -l option can't use negtive value\n");
# }

# print "Extend exon length: $exlen\n";
if ($opfn eq ""){
    $opfn="trans2AS_events.out";
    print "Output file: $opfn.EANMD.input.txt\n";
}else{
    print "Output file: $opfn.EANMD.input.txt\n";
}


open OUTEAINPUT, "> $opfn.EANMD.input.txt" or die ("[-] Error: Can't open or creat $opfn.EANMD.input.txt\n");


our $count1=0;

## Main 
LINE: while(our $row = <>){

    my @col =split(/\t/,$row);
    my $tmpID = $col[0];

    if ($col[0] eq "seqname") {next;}
    $count1++;

    my $gene_id=$col[1];
    my $event_id=$col[2];

    my @coods=split(/[:]/, $event_id);

    my $chr=$coods[1];
    my $Strand=$coods[5]; # ENSMUSG00000033793.12;RI:chr1:5093363:5093452-5095615:5095728:+
    my @cood1 =split(/-/,$coods[3]);
    my $USstart = $coods[2]-1;
    my $USend= $cood1[0];
    my $SEstart=$cood1[0];
    # my @cood2=split(/-/,$coods[3]);
    my $SEend=$cood1[1]-1;
    my $DSstart = $cood1[1]-1;
    my $DSend =$coods[4];
    # my ( $SEend, $DSend) =split(/-/,$coods[3]);
    # my $SEend=split(/\t/,$coods[3])[1];
    # my $DSstart=split(/\t/,$coods[3])[2];
    # say $chr;
    # $DSstart=$DSstart-1;
    # my $USstart=$USend-$exlen-1; # BED Base 0;
    # my $DSend=$DSstart+$exlen;

    my $AS="$chr:$SEstart:$SEend:$Strand".'@'."$chr:$USstart:$USend:$Strand".'@'."$chr:$DSstart:$DSend:$Strand";
    # say $AS;
    print OUTEAINPUT "$gene_id\t$AS\t$gene_id\n";

    # say @coods;

}
close OUTEAINPUT;
print "All $count1 AS events.\n";