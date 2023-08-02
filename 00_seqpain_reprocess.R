#Do you need to process your Salmon quant files to ready for DESeq2?

library(tidyverse)
library(tximport)


#This method is a for loop that iterates over the list of files
#digests them one -by-one
#and makes a list

#first, list all the quant files

# dir <- "/home/shg618/taxolpathy/"
dir <- getwd()
data_path <- "salmon_quants/"
folder_list <- grep("salmon_quant", list.files(paste(dir, data_path, sep = "")), value = TRUE)
rm(file_list)
rm(temp_file_list)
for(folder in folder_list){
  temp_file_list <- data.frame(`folder` = folder,
                               filename = grep("quant.sf", list.files(paste(dir, data_path, folder, sep = "")), value = TRUE))
  if(exists("file_list")){
    file_list <- rbind(file_list, temp_file_list)
  }
  if(!exists("file_list")){
    file_list <- temp_file_list
  }
  rm(temp_file_list)
}
rm(folder)

#make an import manifest
#this was specific to the filenames I have. The data were sequenced over 3 separate runs
#with technical/sequencing replicates between two of the runs
coldata <- data.frame(folder = file_list$folder,
                      files = paste0(dir, data_path, file_list$folder, "/", file_list$filename)) %>%
  separate(folder, into = c("Run", "sample_number"), sep = "_S", extra = "merge", remove = TRUE) %>%
  dplyr::relocate(files, .after = sample_number) %>%
  mutate(sample_number = paste0("S", sample_number)) %>%
  separate(sample_number, into = c("sample_number"), sep = "_", extra = "drop", remove = FALSE) %>%
  mutate(seq_name = paste0(Run, "_", sample_number)) %>%
  mutate(seq_name = gsub("FC06905", "p6", seq_name)) %>%
  mutate(seq_name = gsub("FC07208", "p71", seq_name)) %>%
  mutate(seq_name = gsub("FC07215", "p72", seq_name)) %>%
  mutate(Run = factor(Run),
         seq_name = factor(seq_name))

rm(file_list)

#import the quant files
for(i in 1:length(levels(coldata$Run))){
  namevar <- levels(coldata$Run)[i]
  df <- coldata %>%
    dplyr::filter(Run %in% namevar) %>%
    droplevels(.)
  txi_list <- tximport::tximport(df$files, type = "salmon",
                                 countsFromAbundance = "no",
                                 txIn = TRUE,
                                 tx2gene = gene_tx_EnsDb[, c("tx_id_version", "tx_id", "gene_id")],
                                 ignoreTxVersion = FALSE)
  assign(paste0("txi_list_", namevar), txi_list, envir = .GlobalEnv)
  rm(namevar)
  rm(txi_list)
  rm(df)
}



counts_combined <- full_join(
  (txi_list_FC06905$counts %>% 
     round() %>% 
     data.frame %>%
     setNames(., (coldata %>%
                    dplyr::filter(Run == "FC06905") %>%
                    droplevels(.) %>%
                    dplyr::select(seq_name) %>%
                    simplify() %>%
                    as.character())) %>%
     rownames_to_column(., var = "tx_id")),
  (txi_list_FC07208$counts %>% 
     round() %>% 
     data.frame %>%
     setNames(., (coldata %>%
                    dplyr::filter(Run == "FC07208") %>%
                    droplevels(.) %>%
                    dplyr::select(seq_name) %>%
                    simplify() %>%
                    as.character())) %>%
     rownames_to_column(., var = "tx_id"))) %>%
  full_join(., 
            (txi_list_FC07215$counts %>% 
               round() %>% 
               data.frame %>%
               setNames(., (coldata %>%
                              dplyr::filter(Run == "FC07215") %>%
                              droplevels(.) %>%
                              dplyr::select(seq_name) %>%
                              simplify() %>%
                              as.character())) %>%
               rownames_to_column(., var = "tx_id"))
  ) 


counts_combined <- counts_combined %>%
  left_join(., gene_tx_EnsDb[, c("tx_id", "gene_id", "gene_name", "description", "gene_biotype")]) %>%
  relocate(any_of(c("gene_name", "description", "gene_biotype")), .after = gene_id) %>%
  dplyr::rename(., gene_symbol = gene_name)


rm(txi_list_FC06905)
rm(txi_list_FC07208)
rm(txi_list_FC07215)

write_tsv(counts_combined, file = "counts_combined.tsv.gz", quote = "none", col_names = TRUE)
