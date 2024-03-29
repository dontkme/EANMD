#!/usr/bin/env Rscript
# combinefile <- commandArgs(trailingOnly = TRUE)
# # print(c("combinefile: ", combinefile))
# print(combinefile)
###### EANMDflagcount_reverse.R v1.04
##### Written by Kaining Hu 2022-08-23
library(getopt)

spec <- matrix(

  c("Output",  "o", 1, "character", "Output prefix",
    #"Rank", "r", 1, "character",  "flag rank",
    "Input",  "i", 2, "character",  "Input combined file",
    "help",   "h", 0, "logical",  "Detail Help"),
  byrow = TRUE, ncol = 5
  )


opt <- getopt(spec = spec)
if (is.null(opt$Output)) {
  opt$Output <- opt$Input
}
# print(opt)


if(!is.null(opt$help) || is.null(opt$Input) || is.null(opt$Output) ) {
  # ... 
  cat(paste(getopt(spec = spec, usage = T), "\n"))
  
  quit()
}



# TFNMD <-read.csv(opt$Input,sep="\t",header = T,stringsAsFactors = F,comment.char = "")
TFNMD <-read.table(opt$Input,sep="\t",header = T,stringsAsFactors = F)
#TFNMD <-read.table(opt$Input,sep="\t",header = T,stringsAsFactors = F,comment.char = "")
attach(TFNMD)
 # TFNMD$key <- paste(TFNMD$X.QueryCol1,TFNMD$SEUSDSCoordinates,sep="_")
 TFNMD$key <- paste(TFNMD$QueryCol1,TFNMD$SEUSDSCoordinates,sep="_")
 # sortedTFNMD <- TFNMD[order(TFNMD$X.QueryCol1,TFNMD$SEUSDSCoordinates,TFNMD$NMD_in.ex_flag),]
 sortedTFNMD <- TFNMD[order(TFNMD$QueryCol1,TFNMD$SEUSDSCoordinates,TFNMD$NMD_in.ex_flag),]
 keytable2 <- table(sortedTFNMD$key, sortedTFNMD$NMD_in.ex_flag)
 # head(keytable2)
 keytable2outname <- paste(opt$Output, "Key2NMDex_in_table.txt", sep = ".")
 # print(keytable2outname)
write.table(keytable2, file = keytable2outname, sep = "\t", col.names = NA)
SEUS_Pos_table2 <- table(sortedTFNMD$key, sortedTFNMD$SE.US._Pos)
 # head(SEUS_Pos_table2)
SEUS_Pos_table2outname <- paste(opt$Output, "Key2SEUS_POS_table.txt", sep = ".")
 # print(SEUS_Pos_table2outname)
write.table(SEUS_Pos_table2, file = SEUS_Pos_table2outname, sep = "\t", col.names = NA)
Framekeytable2 <- table(sortedTFNMD$key, sortedTFNMD$Frame_shift_flag)
 # head(Framekeytable2)
Framekeytable2outname <- paste(opt$Output, "Key2Frame_flag_table.txt", sep = ".")
 # print(Framekeytable2outname)
write.table(Framekeytable2,file = Framekeytable2outname,sep = "\t", col.names = NA)
# mergedfile<- merge(merge(keytable2,SEUS_Pos_table2,by="key"),Framekeytable2 ,by="key")
# head(mergedfile)
  keytable2 <- as.data.frame.matrix(keytable2)
  SEUS_Pos_table2 <- as.data.frame.matrix(SEUS_Pos_table2)
  Framekeytable2 <- as.data.frame.matrix(Framekeytable2)
  ma <- merge(keytable2, Framekeytable2, by = 0)
  row.names(ma) <- ma$Row.names
  ma$Row.names <- NULL
  mall <- merge(ma, SEUS_Pos_table2, by = 0)
  mallname <- paste(opt$Output, "FinalUniqNMDflag.txt", sep = ".")
  names(mall) <- make.names(names(mall)) # make.names of safe column names
  mallnrow <- nrow(mall) # count number of rows
  print(mallnrow)
  for (i in 1:mallnrow) {  #### Reverse ordered 2021-07-28
    if (!is.null(mall$No_info[i]) && mall$No_info[i]>0 && mall$No_info[i]>(mall$NMD_ex[i] + mall$NMD_in[i])) {   #### Add number compare 2021-08-23
      mall$Finalflag[i] <- "No_info"
    }else if (!is.null(mall[i,3]) && mall[i,3] > 0 && mall[i,3] > (mall$NMD_ex[i] + mall$NMD_in[i])) {
      mall$Finalflag[i]="No_NMD"
    }else if (!is.null(mall$V1.x[i]) && mall$V1.x[i]>0 && mall$V1.x[i]>(mall$NMD_ex[i] + mall$NMD_in[i]) ){
      mall$Finalflag[i]="Not_detected"
    }else if(!is.null(mall$Same.stop_codon.Need.check[i]) && mall$Same.stop_codon.Need.check[i]>0){
      mall$Finalflag[i]="Same_stop_codon_Need_check" #Add need_check 2021.08.13
    }else if(!is.null(mall$No.stop_codon[i]) && mall$No.stop_codon[i]>0){
       mall$Finalflag[i]="No_stop_codon" #Add No_stop_codon 2021.08.13
    }else if(!is.null(mall$Start_codon[i]) &&  mall$Start_codon[i]>0){
       mall$Finalflag[i]="Start_codon" #Add Start_codon 2021.08.16
    }else if (!is.null(mall$NMD_ex[i]) && !is.null(mall$NMD_in[i]) && mall$NMD_ex[i]>0 && mall$NMD_in[i]>0){
      mall$Finalflag[i]="NMD_ex_in"
    }else if(!is.null(mall$NMD_ex[i]) && !is.null(mall$NMD_in[i]) && mall$NMD_ex[i] >0 && mall$NMD_in[i]==0){
      mall$Finalflag[i]="NMD_ex"
    }else if(!is.null(mall$NMD_ex[i]) && !is.null(mall$NMD_in[i]) && mall$NMD_ex[i]==0 && mall$NMD_in[i]>0){
      mall$Finalflag[i]="NMD_in"
    # }else if(mall$NMD_ex[i]==0 && mall$NMD_in[i]==0 && mall$Start_codon[i]>0){
    }
    
    # if (mall$NMD_ex[i]>0 && mall$NMD_in[i]>0){
    #   mall$Finalflag[i]="NMD_ex_in"
    # }else if(mall$NMD_ex[i] >0 && mall$NMD_in[i]==0){
    #   mall$Finalflag[i]="NMD_ex"
    # }else if(mall$NMD_ex[i]==0 && mall$NMD_in[i]>0){
    #   mall$Finalflag[i]="NMD_in"
    # # }else if(mall$NMD_ex[i]==0 && mall$NMD_in[i]==0 && mall$Start_codon[i]>0){
    # }else if( mall$Start_codon[i]>0){
    #   mall$Finalflag[i]="Start_codon"
    # }else if(mall$X5UTR[i]>0 ){
    #   mall$Finalflag[i]="5UTR"
    # }else if((mall$Upstream.stop_codon[i] + mall$Downstream.stop_codon[i])>0){
    #   mall$Finalflag[i]="ORF_changing"
    # }else if(mall$ORF.preserving[i]>0){
    #   mall$Finalflag[i]="ORF_preserving"
    # }else if(mall$No.stop_codon[i]>0){
    #   mall$Finalflag[i]="No_stop_codon"
    # }else if(mall$Same.stop_codon.Need.check[i]>0){
    #   mall$Finalflag[i]="Same_stop_codon_Need_check"
    # }else if(mall$X3UTR[i]>0){
    #   mall$Finalflag[i]="3UTR"
    # }else if(mall$Same.stop_codon[i]>0){
    #   mall$Finalflag[i]="Same.stop_codon"
    # }else{
    #   mall$Finalflag[i]="Other"
    # }
  }
  
  write.table(mall,file=mallname,sep = "\t",col.names=NA)
  print(table(mall$Finalflag))
  piedata <- as.data.frame(table(mall$Finalflag))
  # piedata <- as.data.frame.matrix(table(mall$Finalflag))
  piedataoutname <- paste(opt$Output,"FinalFlagPieData.txt",sep=".")
  write.table(piedata,file=piedataoutname,sep = "\t",col.names=NA)
  # print(piedata)
  # library(RColorBrewer)
  # myPalette <- brewer.pal(9, "Purples") 
  
  
