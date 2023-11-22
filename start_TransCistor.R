rm(list=ls())

library(rlang)
library(dplyr)
library(rlist)
library(data.table)
library(readxl)
library(ggrepel)
library(tidyverse)
library(hash)
setwd('./')
data_path <- paste0('RegulationFiles/')

metadata_file <- 'metadata.txt'
metadata <- read.table(metadata_file, header = TRUE, sep = "\t",na.strings=c("","NA"))

output_directory <- "Results/"
output_name <- 'TCResults_final.RData'

source( "transcistor.R")
enrichr_threshold <- 0

TADs = c()
Windows = c()
Distance = c()
Pvalues = c()

for(i in 1:nrow(metadata)){
    print(paste0(i, "_of_", (nrow(metadata))))
    file_name <- metadata[i,]$File_Name
    print(paste0("current_file : ", file_name))
    current_file <- read.table(paste0(data_path, file_name), header= FALSE, sep = '\t', comment.char = '')
    
    
   
    name <- metadata[i,]$SYMBOL
    chr <- metadata[i,]$CHR
    tss <- metadata[i,]$TSS
    strand <- metadata[i,]$strand
    file.type <- metadata[i,]$FILE_TYPE
    species <- metadata[i,]$SPECIES
  
    #TransCistor function call
    results <- TransCistor(input.file = current_file, 
                           id.type = file.type, 
                           species = species, 
                           TAD='All',
                           lncRNA.name = name,
                           lncRNA.chr = chr, 
                           lncRNA.tss = tss, 
                           lncRNA.strand = strand,
                           filename=file_name
                           )
    TADs <- rbind(TADs, results[[1]])
    Windows <- rbind(Windows, results[[2]])
    Distance <- rbind(Distance, results[[3]])
    Pvalues <- rbind(Pvalues, results[[4]])
    

}
save.image(paste0(output_directory, output_name))

dy=dplyr::filter(TADs, Proximal>0)
min(dy$Sizes)


