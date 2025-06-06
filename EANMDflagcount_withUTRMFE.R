#!/usr/bin/env Rscript

# combinefile <- commandArgs(trailingOnly = TRUE)
# # print(c("combinefile: ", combinefile))
# print(combinefile)
###### EANMDflagcount.R v1.51
##### Written by Kaining Hu 2025-01-09
options(warn = -1)
library(getopt)
library(dplyr)
library(stringr)

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



 # 转换函数
convert_to_suppa <- function(coord) {
  # Step 1: 按照 @ 分割
  parts <- str_split(coord, "@")[[1]]
  
  # Step 2: 分别提取三部分的起始和终止坐标
  part1 <- str_split(parts[1], ":")[[1]]
  part2 <- str_split(parts[2], ":")[[1]]
  part3 <- str_split(parts[3], ":")[[1]]
  
  # 提取染色体
  chr <- part1[1]
  
  # 提取起始坐标 (part2的结束和part1的开始)
  start_exon_end <- part2[3]
  alt_exon_start <- part1[2]
  alt_exon_end <- part1[3]
  end_exon_start <- part3[2]
  
  # 提取链方向
  strand <- part1[4]
  
  # Step 3: 组合成 SUPPA 格式
  suppa_format <- paste0(chr, ":", start_exon_end, "-", alt_exon_start, ":", alt_exon_end, "-", end_exon_start, ":", strand)
  
  return(suppa_format)
}

TFNMD <- TFNMD %>%
  mutate(AS.SUPPA = sapply(SEUSDSCoordinates, convert_to_suppa))
 # sortedTFNMD <- TFNMD[order(TFNMD$X.QueryCol1,TFNMD$SEUSDSCoordinates,TFNMD$NMD_in.ex_flag),]
 sortedTFNMD <- TFNMD[order(TFNMD$QueryCol1, TFNMD$SEUSDSCoordinates, TFNMD$NMD_in.ex_flag),]

sortedTFNMD.SUPPA <- sortedTFNMD %>% select(AS_events, AS.SUPPA) %>%
  distinct(AS_events, .keep_all = TRUE)


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
sortedTFNMD <- sortedTFNMD %>% rowwise() %>% mutate(SE_Pos_P = round((SE_exon_Number - Start_exon+1)/(Stop_exon - Start_exon+1), 8)) %>% 
                            mutate(SE_Pos_FStart = round((SEupstreamCDS)/(Ori_CDS_length), 8)) %>%
                            mutate(SE_Pos_FEnd = round((SEupstreamCDS + SE_length)/(Ori_CDS_length), 8)) %>%
                            mutate(SE_Pos_FMid = round((SEupstreamCDS + SE_length/2)/(Ori_CDS_length), 8)) %>% ## 20241014 add SE fraction
                            mutate(SE_Pos_2LE = Exons - SE_exon_Number) %>%
                            mutate(Min_stop_pos = pmin(as.numeric(Ori_AA_1st_stop_pos), as.numeric(SEed_AA_1st_stop_pos), na.rm=TRUE)) %>%
                            # mutate(MinStop_Score = ifelse((Min_stop_pos * 3 ) <= 200, (Min_stop_pos * 3)/200 *0.5 , 1)) %>% ## 20241016 1st stop pos <200 nt efficiency # 2024.11.20 remove
                            mutate(MaxDJ = pmax(as.numeric(Ori_last_dj), as.numeric(New_1st_stop_pos_dj), na.rm=TRUE)) %>%
  # mutate(NMD_Score = ifelse(SE_Pos_P<=1, 
  #                           round(NMD_P * (1 / (1 + exp(-10 * (SE_Pos_P - 0.251)))), 8), 
  #                           round(NMD_P * (1 / (1 + exp(5 * (SE_Pos_P - 1.249)))), 8))) %>% #sigmoid Use stop_exon CDS + Exons Exon position
    # mutate(NMD_Score = ifelse(SE_Pos_P<=1, 
    #                         round( (1 / (1 + exp(-10 * (SE_Pos_P - 0.251)))), 8), 
    #                         round( (1 / (1 + exp(5 * (SE_Pos_P - 1.249)))), 8))) %>%  # 2024.11.20 remove
  # mutate(NMD_Score = ifelse(source == "USDS", NMD_Score * 1.5, NMD_Score)) %>% # Buff USDS NMD_score
  # mutate(NMD_Score = ifelse(!(SE_length  %in% c("Null", "", "NA")), ifelse(SE_length %% 3 != 0, NMD_Score * 2, NMD_Score), NMD_Score)) %>%  # Buff frame shift
  # mutate(NMD_Score = ifelse(!(SEed_AA_1st_stop_pos  %in% c("Null", "", "-")), ifelse( as.numeric(SEed_AA_1st_stop_pos) * 3 < 51, NMD_Score * 0.25, NMD_Score), NMD_Score)) %>% # Buff new stop condon longer than 50. # 2024.09.19
  # mutate(NMD_Score = ifelse((SE_exon_Number - Start_exon + 1) <= 2, NMD_Score * 0.25, NMD_Score)) %>%
  # mutate(NMD_Score = ifelse((!(SEupstreamCDS  %in% c("Null", "", "-",NA)) & SEupstreamCDS <= 51), NMD_Score * 0.25, NMD_Score))  %>% # 2024.11.05 add UTR-3 sequence and RNA # 2024.11.20 remove
  mutate(SE_2Start = SE_exon_Number - Start_exon, # 2024.11.25 rm +1 mean intron number(s).
         Min_stop_pos_F = Min_stop_pos * 3/Ori_Star_codon_to_exon_end_seq_len, ## 2024.11.21 fix bug
         Min_stop_pos_F.MaxDJ = Min_stop_pos_F * MaxDJ,
         Min_stop_Pos_to_End_nt = Ori_Star_codon_to_exon_end_seq_len - Min_stop_pos * 3, # 2024.11.21 add Min_stop_Pos_to_End_nt
         # Min_stop_pos_F.MaxDJ.absMaxDJ = Min_stop_pos_F*MaxDJ*abs(MaxDJ)
         SEmod3 = SE_length %% 3,
         LE_length = Ori_Star_codon_to_exon_end_seq_len - Ori_last_junction_pos,
         UTR3_length = Ori_Star_codon_to_exon_end_seq_len - Ori_CDS_length,
         UTR3_introns = Exons - Stop_exon,
         UTR3_by_LE_length = UTR3_length/LE_length, # 2024.11.01 add new 3 features.
         UTR3_F_FL = UTR3_length/Ori_Star_codon_to_exon_end_seq_len,
         UTR3_DJ = UTR3_length - LE_length,
         New_UTR3_length = rm.add_SE_start_to_end_seq_len - as.numeric(SEed_AA_1st_stop_pos) * 3 + 3, # 2024.11.04 # 2024.11.06 add stop-codon to 3'UTR length
          # pmin(as.numeric(Ori_AA_1st_stop_pos), as.numeric(SEed_AA_1st_stop_pos), na.rm=TRUE))
         ) %>%
    mutate(Max_UTR3_length = pmax(as.numeric(UTR3_length), as.numeric(New_UTR3_length), na.rm=TRUE), 
            Min_UTR3_length = pmin(as.numeric(UTR3_length), as.numeric(New_UTR3_length), na.rm=TRUE), 
            Ori_UTR3_seq = str_sub(Ori_CDSexons_seq, as.numeric(Ori_CDS_length)+1, -1),
            New_UTR3_seq = ifelse(!(SEed_AA_1st_stop_pos  %in% c("Null", "", "-",NA)), str_sub(rm.add_SE_CDSexons_seq, as.numeric(SEed_AA_1st_stop_pos) * 3 -2, -1), NA)
            ) %>% # 2024.11.06 add stop-codon to 3'UTR length
    mutate(Ori_UTR3_seq = ifelse(Ori_UTR3_seq == "", "0", Ori_UTR3_seq),
           New_UTR3_seq = ifelse(New_UTR3_seq == "", "0", New_UTR3_seq)
          ) %>% 
    mutate(
          Max_UTR3_seq = ifelse(!(Max_UTR3_length  %in% c( "", NA)) & (Max_UTR3_length == UTR3_length), Ori_UTR3_seq,
          ifelse((Max_UTR3_length == New_UTR3_length), New_UTR3_seq, NA)),
          Min_UTR3_seq = ifelse(!(Min_UTR3_length  %in% c( "", NA)) & (Min_UTR3_length == UTR3_length), Ori_UTR3_seq,
          ifelse((Min_UTR3_length == New_UTR3_length), New_UTR3_seq, NA)) 
    ) # 2024.11.06 add Max and Min UTR3 seq

    All_Ori_UTR3_seq <- sortedTFNMD$Ori_UTR3_seq
    All_New_UTR3_seq <- sortedTFNMD$New_UTR3_seq # 2024.11.07 single files.
    All_Ori_UTR3_seq.file.raw <- paste(opt$Output, "temp.All_Ori_UTR3_seq.raw.txt", sep = ".") 
    All_New_UTR3_seq.file.raw <- paste(opt$Output, "temp.All_New_UTR3_seq.raw.txt", sep = ".") 
    write.table(All_Ori_UTR3_seq, file = All_Ori_UTR3_seq.file.raw, col.names = F, row.names = F, quote = F)
    write.table(All_New_UTR3_seq, file = All_New_UTR3_seq.file.raw, col.names = F, row.names = F, quote = F)

    sortedTFNMD_LT30K <- sortedTFNMD %>%
    mutate(
        Ori_UTR3_seq = ifelse(nchar(Ori_UTR3_seq) > 30000, 
                                  str_sub(Ori_UTR3_seq, -30000, -1), 
                                  Ori_UTR3_seq),
        New_UTR3_seq = ifelse(nchar(New_UTR3_seq) > 30000, 
                                  str_sub(New_UTR3_seq, -30000, -1), 
                                  New_UTR3_seq)
    )

    All_Ori_UTR3_seq <- sortedTFNMD_LT30K$Ori_UTR3_seq
    All_New_UTR3_seq <- sortedTFNMD_LT30K$New_UTR3_seq ### 2025-01-09 add max 30K as RNAfold input.

    # Define a function to get MFE from RNAfold
# get_mfe <- function(sequence) {
#   # Save the sequence to a temporary file
#   temp_file <- tempfile()
#   writeLines(sequence, temp_file)
  
#   # Call RNAfold on the sequence file, capture output
#   output <- system(paste("RNAfold -j4 --noPS -i ", temp_file), intern = TRUE)
  
#   # Clean up temp file
#   unlink(temp_file)
  
#   # Parse the output for MFE value
#   # MFE is in the last part of the second line in the format: "( ... )   -X.XX kcal/mol"
#   mfe_line <- output[2]  
#   mfe_value <- as.numeric(sub(".*\\((.*)\\).*", "\\1", mfe_line)) # Extract MFE value
  
#   return(mfe_value)
# }

# # Example usage
# sequence <- "GCGCUUCGCCGAGCGC"
# mfe_value <- get_mfe(sequence)
# print(mfe_value)

# sortedTFNMD <- sortedTFNMD %>% rowwise() %>% 
#     mutate(Ori_UTR3_seq_MFE = get_mfe(Ori_UTR3_seq),
#           New_UTR3_seq_MFE = get_mfe(New_UTR3_seq)
#     ) ## Add MFE use RNAfold

# Function to run RNAfold in batch mode and extract MFE values
# get_mfe_batch <- function(sequences) {
#   # Remove any NA or empty sequences
#   sequences <- sequences[sequences != "" & !is.na(sequences)]
  
#   # Save sequences to a temporary file in FASTA format
#   temp_input_file <- tempfile(fileext = ".fa")
#   temp_output_file <- tempfile(fileext = ".out")
  
#   writeLines(paste0(">seq", seq_along(sequences), "\n", sequences), temp_input_file)
  
#   # Run RNAfold and capture the output in a single call
#   system(paste("RNAfold -j4 --noPS < ", temp_input_file, " > ", temp_output_file))
  
#   # Read output and parse MFE values
#   output <- readLines(temp_output_file)
  
#   # Check for expected output structure
#   mfe_lines <- output[seq(2, length(output), by = 2)]
#   if (length(mfe_lines) != length(sequences)) {
#     mfe_values <- rep(NA, length(sequences)) # Use NA if there's a mismatch
#   } else {
#     mfe_values <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", mfe_lines))
#   }
  
#   # Clean up temporary files
#   unlink(c(temp_input_file, temp_output_file))
  
#   # Return MFE values along with NA for empty or missing sequences
#   return(mfe_values)
# }

# Now applying get_mfe_batch to the dataframe using mutate
# sortedTFNMD <- sortedTFNMD %>% ungroup() %>%
#   mutate(
#     Ori_UTR3_seq_MFE = if_else(Ori_UTR3_seq != "" & !is.na(Ori_UTR3_seq), get_mfe_batch(Ori_UTR3_seq), NA_real_),
#     New_UTR3_seq_MFE = if_else(New_UTR3_seq != "" & !is.na(New_UTR3_seq), get_mfe_batch(New_UTR3_seq), NA_real_)
#   )# 2024.11.07 

# Summary.sortedTFNMD <- sortedTFNMD %>% ungroup() %>%
# group_by(AS_events) %>%
#   summarise(
#     MaxNMD_Score = max(NMD_Score, na.rm = TRUE),
#     MinNMD_Score = min(NMD_Score, na.rm = TRUE),
#     MeanNMD_Score = mean(NMD_Score, na.rm = TRUE),
#     SDNMD_Score = sd(NMD_Score, na.rm = TRUE))

# sortedTFNMD <- sortedTFNMD %>% left_join(Summary.sortedTFNMD, , by = "AS_events") # 2024.11.20 remove NMD_score summary

write.table(sortedTFNMD, file = paste(opt$Output, "AS_events_NMD_P.txt", sep = "."), sep = "\t", row.names = F, quote = F)




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
    # }else if(!is.null((mall$NMD_ex[i]) && !is.null(mall$NMD_in[i])) && (mall$NMD_ex[i] > 0 && mall$NMD_in[i] > 0)) {
    }else if(!is.null(mall$NMD_ex[i]) && !is.null(mall$NMD_in[i]) && (mall$NMD_ex[i] > 0 && mall$NMD_in[i] > 0)) { # 2025-03-06. claude suggested
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
  
  # mall <- mall %>% left_join(Summary.sortedTFNMD, by = "AS_events") ## 2024.09.19 Add NMD_score summary. ## 2024.11.20 remove Summary.sortedTFNMD
  mall <- mall %>% left_join(sortedTFNMD.SUPPA, by = "AS_events") # 2024.09.20 Add AS.SUPPA

  write.table(mall, file = mallname, sep = "\t",col.names = NA)
  print(table(mall$Finalflag))
  piedata <- as.data.frame(table(mall$Finalflag))
  # piedata <- as.data.frame.matrix(table(mall$Finalflag))
  piedataoutname <- paste(opt$Output, "FinalFlagPieData.txt", sep = ".")
  write.table(piedata,file = piedataoutname, sep = "\t", col.names = NA)
  # print(piedata)
  # library(RColorBrewer)
  # myPalette <- brewer.pal(9, "Purples") 
  
cat("Start calculate MFE by RNAfold, it may take minutes. ")
start_time <- Sys.time()

  #  All_Ori_UTR3_seq <- sortedTFNMD$Ori_UTR3_seq
  # All_New_UTR3_seq <- sortedTFNMD$New_UTR3_seq 
  All_Ori_UTR3_seq.file <- paste(opt$Output, "temp.All_Ori_UTR3_seq.txt", sep = ".") 
  All_New_UTR3_seq.file <- paste(opt$Output, "temp.All_New_UTR3_seq.txt", sep = ".") 
  write.table(All_Ori_UTR3_seq, file = All_Ori_UTR3_seq.file, col.names = F, row.names = F, quote = F)
  write.table(All_New_UTR3_seq, file = All_New_UTR3_seq.file, col.names = F, row.names = F, quote = F)
  # get_mfe_batch <- function(sequences) {
  # Remove any NA or empty sequences
  # sequences <- sequences[sequences != "" & !is.na(sequences)]
  All_Ori_UTR3_seq.MFE <- paste(opt$Output, "temp.All_Ori_UTR3_seq.MFE.txt", sep = ".")
  All_New_UTR3_seq.MFE <- paste(opt$Output, "temp.All_New_UTR3_seq.MFE.txt", sep = ".") 
  # Save sequences to a temporary file in FASTA format
  # temp_input_file <- tempfile(fileext = ".fa")
  # temp_output_file <- tempfile(fileext = ".out")
  
  # writeLines(paste0(">seq", seq_along(sequences), "\n", sequences), temp_input_file)
  
  # Run RNAfold and capture the output in a single call
  system(paste("RNAfold -j7 --noPS < ", All_Ori_UTR3_seq.file, " > ", All_Ori_UTR3_seq.MFE)) # CPU threads: 7
  system(paste("RNAfold -j7 --noPS < ", All_New_UTR3_seq.file, " > ", All_New_UTR3_seq.MFE))
  
  # Read output and parse MFE values
  output1 <- readLines(All_Ori_UTR3_seq.MFE)
  output2 <- readLines(All_New_UTR3_seq.MFE)
  
  # Check for expected output structure
  mfe_lines1 <- output1[seq(2, length(output1), by = 2)]
  mfe_lines2 <- output2[seq(2, length(output2), by = 2)]
  # if (length(mfe_lines) != length(sequences)) {
  #   mfe_values <- rep(NA, length(sequences)) # Use NA if there's a mismatch
  # } else {
  mfe_values1 <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", mfe_lines1))
  mfe_values2 <- as.numeric(gsub(".*\\((.*)\\).*", "\\1", mfe_lines2))

  sortedTFNMD$Ori_UTR3_seq_MFE <- mfe_values1
  sortedTFNMD$New_UTR3_seq_MFE <- mfe_values2
  # sortedTFNMD <- sortedTFNMD %>% rowwise() %>% 
  #       mutate(Ori_UTR3_seq_MFE_per_nt = Ori_UTR3_seq_MFE/UTR3_length,
  #             New_UTR3_seq_MFE_per_nt = New_UTR3_seq_MFE/New_UTR3_length
  #       )

  sortedTFNMD <- sortedTFNMD %>% 
    rowwise() %>% 
    mutate(
        Ori_UTR3_seq_MFE_per_nt = if_else(UTR3_length > 0, Ori_UTR3_seq_MFE / pmin(UTR3_length, 30000), NA_real_), # 2025-01-09, add pmin 30K
        New_UTR3_seq_MFE_per_nt = if_else(New_UTR3_length > 0, New_UTR3_seq_MFE / pmin(New_UTR3_length, 30000), NA_real_) # 2025-01-09, add pmin 30K
    )
  # }
  MEFoutfile <- paste(opt$Output, "AS_events_NMD_P.MFE.txt", sep = ".")
  # write.table(sortedTFNMD, file = paste(opt$Output, "AS_events_NMD_P.MFE.txt", sep = "."), sep = "\t", row.names = F, quote = F)
  write.table(sortedTFNMD, file = MEFoutfile, sep = "\t", row.names = F, quote = F)
  # Clean up temporary files
  # unlink(c(temp_input_file, temp_output_file))
  unlink(c(All_Ori_UTR3_seq.MFE, All_New_UTR3_seq.MFE, All_Ori_UTR3_seq.MFE, All_New_UTR3_seq.MFE))
  # Return MFE values along with NA for empty or missing sequences
  # return(mfe_values)
# }
end_time <- Sys.time()

run_time <- end_time - start_time
cat("Runing time: ", run_time) # 2024.11.07 separate RNAfold part. 

cat("Predict NMD efficiency ") # 2024.11.22 add xbt prediction
library(xgboost)

xbt.input <- read.table(file= MEFoutfile, header = T, sep = "\t")
# or sortedTFNMD
# xbt.input <- sortedTFNMD # Input
bst_NMD_efficiency_0 <- xgb.load("NMD_efficiency.xgboost.model") # load model

xbt.input.select<- xbt.input %>% dplyr::select(AS_events_T, Min_stop_pos_F, Min_stop_pos_F.MaxDJ, SE_length, SE_Pos_2LE, UTR3_length, MaxDJ, Start_exon, Min_stop_Pos_to_End_nt, Ori_UTR3_seq_MFE_per_nt)
dtest <- xgb.DMatrix(data = as.matrix(xbt.input.select[, -1]))

# predictions <- predict(bst_NMD_efficiency, as.matrix(xbt.input.select[, -1]))
predictions <- predict(bst_NMD_efficiency_0, dtest)

xbt.input.pred <- xbt.input
xbt.input.pred$NMD_Score <- predictions
xbt.input.pred <- xbt.input.pred %>% mutate(NMD_Score = ifelse(SE_Pos_FStart > 0, NMD_Score, NA))

used_best_threshold <- 0.268355399370193 # Threshold

xbt.input.pred <- xbt.input.pred %>% mutate(NMD_or = ifelse(NMD_P > 0, "NMD", "NonNMD"))
xbt.input.pred <- xbt.input.pred %>% mutate(NMD_Pre = ifelse(NMD_Score > used_best_threshold, "NMD", "NonNMD"))

xbtoutfile <- paste(opt$Output, "AS_events_NMD_P.MFE.xbtPred.txt", sep = ".")

Summary.xbt.input.pred <- xbt.input.pred %>% ungroup() %>%
group_by(AS_events) %>%
  summarise(
    MaxNMD_Score = max(NMD_Score, na.rm = TRUE),
    MinNMD_Score = min(NMD_Score, na.rm = TRUE),
    MeanNMD_Score = mean(NMD_Score, na.rm = TRUE),
    SDNMD_Score = sd(NMD_Score, na.rm = TRUE),
    MedianNMD_Score = median(NMD_Score, na.rm = TRUE))

xbt.input.pred <- xbt.input.pred %>% left_join(Summary.xbt.input.pred, , by = "AS_events") # 2024.11.20 remove NMD_score summary
write.table(xbt.input.pred, file = xbtoutfile, sep = "\t", row.names = F, quote = F)

mall <- mall %>% left_join(Summary.xbt.input.pred, by = "AS_events") ## 2024.09.19 Add NMD_score summary. ## 2024.11.20 remove Summary.sortedTFNMD . 2024.11.22 Summary.xbtoutfile
mall <- mall %>% mutate(NMD_Pre_Flag = ifelse(MedianNMD_Score > used_best_threshold, "NMD", "NonNMD"))
  # mall <- mall %>% left_join(sortedTFNMD.SUPPA, by = "AS_events") # 2024.09.20 Add AS.SUPPA
mallname2 <- paste(opt$Output, "FinalUniqNMDflag.NMD_score.txt", sep = ".")
  write.table(mall, file = mallname2, sep = "\t", col.names = NA)

#### 

cat("All done!")
