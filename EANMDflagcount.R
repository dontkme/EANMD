#!/usr/bin/env Rscript
# combinefile <- commandArgs(trailingOnly = TRUE)
# # print(c("combinefile: ", combinefile))
# print(combinefile)
###### EANMDflagcount.R v1.40
##### Written by Kaining Hu 2024-09-03
library(getopt)
library(dplyr)

spec <- matrix(

  c("Output",  "o", 1, "character", "Output prefix",
    #"Rank", "r", 1, "character",  "flag rank",
    "Input",  "i", 2, "character",  "Input combined file",
    "help",   "h", 0, "logical",  "Detail Help"),
  byrow = TRUE, ncol = 5
  )


opt <- getopt(spec=spec)
if(is.null(opt$Output)){
  opt$Output <- opt$Input
}
# print(opt)


if (!is.null(opt$help) || is.null(opt$Input) || is.null(opt$Output)) {
  # ... 
  cat(paste(getopt(spec = spec, usage = T), "\n"))

  quit()
}



# TFNMD <-read.csv(opt$Input,sep="\t",header = T,stringsAsFactors = F,comment.char = "")
TFNMD <- read.table(opt$Input, sep = "\t", header = T, stringsAsFactors = F)
#TFNMD <-read.table(opt$Input,sep="\t",header = T,stringsAsFactors = F,comment.char = "")
attach(TFNMD)
 # TFNMD$key <- paste(TFNMD$X.QueryCol1,TFNMD$SEUSDSCoordinates,sep="_")
 TFNMD$AS_events <- paste(TFNMD$QueryCol1, TFNMD$SEUSDSCoordinates, sep = "_")
 TFNMD$AS_events_T <- paste(TFNMD$AS_events, TFNMD$Transcript_id, TFNMD$SE_exon_Number, TFNMD$source ,sep = "@")
 # sortedTFNMD <- TFNMD[order(TFNMD$X.QueryCol1,TFNMD$SEUSDSCoordinates,TFNMD$NMD_in.ex_flag),]
 sortedTFNMD <- TFNMD[order(TFNMD$QueryCol1, TFNMD$SEUSDSCoordinates, TFNMD$NMD_in.ex_flag),]



 keytable2 <- table(sortedTFNMD$AS_events, sortedTFNMD$NMD_in.ex_flag)
#  keytable2$Sum <- rowSums(keytable2[, 2:ncol(keytable2)])
keytable2 <- as.data.frame.matrix(keytable2)
keytable2$AS_events <- row.names(keytable2)
# print(keytable2)
keytable2 <- keytable2 %>%
  rowwise() %>%
  mutate(
    SumTrans = sum(c_across(where(is.numeric)))
  )

if (!is.null(keytable2$NMD_in)) {

  keytable2 <- keytable2 %>%
  rowwise() %>%
  mutate(
    NMD_in_P = round(NMD_in / SumTrans, 6)
  )
}

if (!is.null(keytable2$NMD_ex)) {

  keytable2 <- keytable2 %>%
  rowwise() %>%
  mutate(
    NMD_ex_P = round(NMD_ex / SumTrans, 6)
  )
}

if (!is.null(keytable2$NMD_in) && !is.null(keytable2$NMD_ex)) {

  keytable2 <- keytable2 %>%
  rowwise() %>%
  mutate(
    NMD_P = round((NMD_in + NMD_ex) / SumTrans, 6)
  )
}
# keytable2

 # head(keytable2)
 keytable2outname <- paste(opt$Output, "Key2NMDex_in_table.txt", sep = ".")
 # print(keytable2outname)
write.table(keytable2, file = keytable2outname, sep = "\t", col.names = NA)

sortedTFNMD <- sortedTFNMD %>% left_join(keytable2, by = "AS_events")


# sortedTFNMD <- sortedTFNMD %>% mutate(SE_Pos_P = round(SE_exon_Number/Exons, 8)) %>% mutate(NMD_Score = round(SE_Pos_P * NMD_P, 8))  ### Add MMD_score transformed_value = NMD_P * (1 / (1 + exp(-k * (SE_position - 0.5))))
# sortedTFNMD <- sortedTFNMD %>% mutate(SE_Pos_P = round(SE_exon_Number/Exons, 8)) %>% mutate(NMD_Score = round(NMD_P * (1 / (1 + exp(-10 * (SE_Pos_P - 0.5)))), 8)) #sigmoid 
#sortedTFNMD <- sortedTFNMD %>% mutate(SE_Pos_P = round(SE_exon_Number/Stop_exon, 8)) %>% mutate(NMD_Score = ifelse(SE_Pos_P<=1, round(NMD_P * (1 / (1 + exp(-10 * (SE_Pos_P - 0.3)))), 8), round(NMD_P * (1 / (1 + exp(5 * (SE_Pos_P - 1.3)))), 8))) #sigmoid Use stop_exon
sortedTFNMD <- sortedTFNMD %>% rowwise() %>% mutate(SE_Pos_P = round((SE_exon_Number - Start_exon+1)/(Stop_exon - Start_exon+1), 8)) %>% mutate(NMD_Score = ifelse(SE_Pos_P<=1, round(NMD_P * (1 / (1 + exp(-10 * (SE_Pos_P - 0.25)))), 8), round(NMD_P * (1 / (1 + exp(5 * (SE_Pos_P - 1.3)))), 8))) %>% #sigmoid Use stop_exon CDS + Exons Exon position
  # mutate(NMD_Score = ifelse(source == "USDS", NMD_Score * 1.5, NMD_Score)) %>% # Buff USDS NMD_score
  # mutate(NMD_Score = ifelse(!(SE_length  %in% c("Null", "", "NA")), ifelse(SE_length %% 3 != 0, NMD_Score * 2, NMD_Score), NMD_Score)) %>%  # Buff frame shift
  mutate(NMD_Score = ifelse(!(SEed_AA_1st_stop_pos  %in% c("Null", "", "-")), ifelse( as.numeric(SEed_AA_1st_stop_pos) * 3 > 50, NMD_Score * 2, NMD_Score), NMD_Score)) # Buff new stop condon longer than 50. # 2024.09.19

Summary.sortedTFNMD <- sortedTFNMD %>%
group_by(AS_events) %>%
  summarise(
    MaxNMD_Score = max(NMD_Score, na.rm = TRUE),
    MinNMD_Score = min(NMD_Score, na.rm = TRUE),
    MeanNMD_Score = mean(NMD_Score, na.rm = TRUE),
    SDNMD_Score = sd(NMD_Score, na.rm = TRUE))

sortedTFNMD <- sortedTFNMD %>% left_join(Summary.sortedTFNMD, , by = "AS_events")

write.table(sortedTFNMD, file = paste(opt$Output, "AS_events_NMD_P.txt", sep = "."), sep = "\t", col.names = NA)




NewExons <- sortedTFNMD %>% filter(Ori_CDSexons_seq != "") %>% select(AS_events_T, rm.add_SE_CDSexons_seq, QueryCol3)
OriExons <- sortedTFNMD %>% filter(Ori_CDSexons_seq != "") %>% select(Transcript_id, Ori_CDSexons_seq, QueryCol3) %>% distinct(Transcript_id, .keep_all = T)
OriExons2 <- OriExons
colnames(OriExons2) <- colnames(NewExons)
TX2GeneList <- rbind(
  NewExons,
  OriExons2
  )


TX2GeneList<-TX2GeneList[,c(1,3,2)]
colnames(TX2GeneList) <- c("TXNAME", "GENEID", "SEQUENCE")
write.table(TX2GeneList, file = paste(opt$Output, "rm.add_SE_output_sequences_OriTX2Gene.txt", sep = "."), sep = "\t", row.names = F)
# Specify the output FASTA file path
output_fasta <- paste(opt$Output, "rm.add_SE_output_sequences.fa", sep = ".")

# Write to FASTA format
write_fasta <- function(df, output_file) {
  # Open the connection to the file
  file_conn <- file(output_file, open = "w")
  
  # Loop over each row of the data frame and write in FASTA format
  for (i in 1:nrow(df)) {
    # Write the FASTA header (with a '>' character)
    writeLines(paste0(">", df$AS_events_T[i]), con = file_conn)
    
    # Write the corresponding sequence
    writeLines(df$rm.add_SE_CDSexons_seq[i], con = file_conn)
  }
  
  # Close the file connection
  close(file_conn)
}

# Call the function to write the data to the FASTA file
write_fasta(NewExons, output_fasta)

# Check the file
cat("FASTA file saved as:", output_fasta)

# Write to FASTA format
write_Ori_fasta <- function(df, output_file) {
  # Open the connection to the file
  file_conn <- file(output_file, open = "w")
  
  # Loop over each row of the data frame and write in FASTA format
  for (i in 1:nrow(df)) {
    # Write the FASTA header (with a '>' character)
    writeLines(paste0(">", df$Transcript_id[i]), con = file_conn)
    
    # Write the corresponding sequence
    writeLines(df$Ori_CDSexons_seq[i], con = file_conn)
  }
  
  # Close the file connection
  close(file_conn)
}
output_Ori_fasta <- paste(opt$Output, "Ori_sequences.fa", sep = ".")
# Call the function to write the data to the FASTA file
write_Ori_fasta(OriExons, output_Ori_fasta)




 SEUS_Pos_table2<- table(sortedTFNMD$AS_events, sortedTFNMD$SE.US._Pos)
 # head(SEUS_Pos_table2)
 SEUS_Pos_table2outname <- paste(opt$Output, "Key2SEUS_POS_table.txt", sep = ".")
 # print(SEUS_Pos_table2outname)
  write.table(SEUS_Pos_table2, file = SEUS_Pos_table2outname,sep = "\t", col.names=NA)
 Framekeytable2<- table(sortedTFNMD$AS_events, sortedTFNMD$Frame_shift_flag)
 # head(Framekeytable2)
 Framekeytable2outname <- paste(opt$Output, "Key2Frame_flag_table.txt", sep = ".")
 # print(Framekeytable2outname)
  write.table(Framekeytable2, file = Framekeytable2outname,sep = "\t", col.names=NA)
# mergedfile<- merge(merge(keytable2,SEUS_Pos_table2,by="key"),Framekeytable2 ,by="key")
# head(mergedfile)
  keytable2 <- as.data.frame(keytable2)
  SEUS_Pos_table2 <- as.data.frame.matrix(SEUS_Pos_table2)
  SEUS_Pos_table2$AS_events <- row.names(SEUS_Pos_table2)
  Framekeytable2 <- as.data.frame.matrix(Framekeytable2)
  Framekeytable2$AS_events <- row.names(Framekeytable2)
  # ma<- merge(keytable2, Framekeytable2,by = 0)
  ma <- keytable2 %>% left_join(Framekeytable2, by = "AS_events")
  row.names(ma) <- ma$Row.names
  ma$Row.names <- NULL
  # mall <- merge(ma, SEUS_Pos_table2, by = 0)
  mall <- ma %>% left_join(SEUS_Pos_table2, by = "AS_events")
  mallname <- paste(opt$Output, "FinalUniqNMDflag.txt", sep = ".")
  names(mall) <- make.names(names(mall)) # make.names of safe column names
  mallnrow <- nrow(mall) # count number of rows
  print(mallnrow)
  for (i in 1:mallnrow) {
    if (!is.null(mall$Start_codon[i]) && mall$Start_codon[i] > 0) {    ### 2021-10-20 Start codon as 1st order.
      mall$Finalflag[i] = "Start_codon"
    }else if(!is.null(mall$NMD_ex[i]) && !is.null(mall$NMD_in[i]) && mall$NMD_ex[i] > 0 && mall$NMD_in[i] > 0) {
      mall$Finalflag[i] = "NMD_ex_in"
    }else if ( !is.null(mall$NMD_in[i]) && mall$NMD_in[i] > 0 && (is.null(mall$NMD_ex[i])|| mall$NMD_ex[i] == 0) ) {
      mall$Finalflag[i] = "NMD_in"
    }else if (!is.null(mall$NMD_ex[i]) && mall$NMD_ex[i]>0 && (mall$NMD_in[i] == 0 || is.null(mall$NMD_in[i]))) {
      mall$Finalflag[i] = "NMD_ex"
    }else if (!is.null(mall$X5UTR[i]) && mall$X5UTR[i] > 0) {
      mall$Finalflag[i] = "5UTR"
    }else if ((mall$Upstream.stop_codon[i] + mall$Downstream.stop_codon[i])>0) {
      mall$Finalflag[i] = "ORF_changing"
    }else if (!is.null(mall$ORF.preserving[i]) && mall$ORF.preserving[i] > 0) {
      mall$Finalflag[i] = "ORF_preserving"
    }else if (!is.null(mall$No.stop_codon[i]) && mall$No.stop_codon[i] > 0) {
      mall$Finalflag[i] = "No_stop_codon"
    }else if (!is.null(mall$Same.stop_codon.Need.check[i]) && mall$Same.stop_codon.Need.check[i] > 0) {
      mall$Finalflag[i] = "Same_stop_codon_diff_from_annot"
    # }else if(mall$Stop_codon[i]>0){
    #   mall$Finalflag[i]="Stop_codon"
    }else if(!is.null(mall$X3UTR[i]) && mall$X3UTR[i]>0){
      mall$Finalflag[i] = "3UTR"
    # }else if(mall$Same.stop_codon[i]>0){
    #   mall$Finalflag[i]="Same.stop_codon" 
    }else{
      mall$Finalflag[i] = "Other"
    }
  }
  
  mall <- mall %>% left_join(Summary.sortedTFNMD, by = "AS_events") ## 2024.09.19 Add NMD_score summary
  write.table(mall, file = mallname, sep = "\t",col.names = NA)
  print(table(mall$Finalflag))
  piedata <- as.data.frame(table(mall$Finalflag))
  # piedata <- as.data.frame.matrix(table(mall$Finalflag))
  piedataoutname <- paste(opt$Output, "FinalFlagPieData.txt", sep = ".")
  write.table(piedata,file = piedataoutname, sep = "\t", col.names = NA)
  # print(piedata)
  # library(RColorBrewer)
  # myPalette <- brewer.pal(9, "Purples") 
  

