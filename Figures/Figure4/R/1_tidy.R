library(tidyverse)
library(dplyr)
library(cowplot)
library(paletteer)
library(viridis)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure4")
RNAdata <- c("~/Box/ZamanianLab/SeqLibraries/Mapping/")

########################################################
####### RAW/TPM Counts for rep 1
########################################################

dir <- c("O002_BmAF_RNAt/181220_BCCV00ANXX/")
counts_dir <- paste0(RNAdata,dir,"counts/")
sample_list <- list.files(path = counts_dir, full.names = FALSE, recursive = TRUE)

#pull in gtf file of gene_ids and lengths
gtf <- read.csv(paste0(RNAdata,"auxillary/Bma.geneset.csv"), header=F) 
colnames(gtf) <- c("chr","start","end","strand","gene_id")
gtf <- gtf %>% mutate(gene_len = abs(end-start)) %>% select(gene_id,gene_len)

#create blank count table using first sample file to pull gene ids
counts <- read.csv(paste0(counts_dir,sample_list[1]), header=F, sep="\t") %>% dplyr::slice(-(1:4)) %>% dplyr::select(1)
colnames(counts) <- c("gene_id")

#populate raw and tpm counts table
counts.raw <- counts
counts.tpm <- counts
for (i in sample_list){
  print(i) 
  counts.i <- read.csv(paste0(counts_dir,i), header=F, sep="\t") %>% dplyr::slice(-(1:4))
  colnames(counts.i) <- c("gene_id","count")
  
  #add tpm column
  counts.i <- left_join(counts.i, gtf, by = "gene_id") %>%
    mutate(rpk = count / (1000*gene_len))
  rpk_sum <- sum(counts.i$rpk, na.rm=TRUE)
  counts.i <- counts.i %>% 
    mutate(tpm = 1000000*rpk / rpk_sum) %>%
    select(gene_id,count,tpm)
  
  #add raw count and tpm data to dfs
  i <- str_replace(i, ".ReadsPerGene.tab","")
  
  counts.raw.i <- counts.i %>% select(gene_id,count)
  colnames(counts.raw.i) <- c("gene_id",i)
  counts.raw <- left_join(counts.raw,counts.raw.i, by = "gene_id")
  
  counts.tpm.i <- counts.i %>% select(gene_id,tpm)
  colnames(counts.tpm.i) <- c("gene_id",i)
  counts.tpm <- left_join(counts.tpm,counts.tpm.i, by = "gene_id")
  
}

# switch to long and add meta data
counts.raw.1 <- counts.raw %>%
  pivot_longer(2:ncol(counts.raw), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression") %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = "")) %>%
  mutate(rep = "R1")

counts.tpm.1 <- counts.tpm %>%
  pivot_longer(2:ncol(counts.tpm), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression") %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = "")) %>%
  mutate(rep = "R1")


########################################################
####### RAW/TPM Counts for reps 2 / 3
########################################################

dir <- c("Q002/200217_AHNHN3DMXX/")
counts_dir <- paste0(RNAdata,dir,"counts/")
sample_list <- list.files(path = counts_dir, full.names = FALSE, recursive = TRUE)

#pull in gtf file of gene_ids and lengths
gtf <- read.csv(paste0(RNAdata,"auxillary/Bma.geneset.csv"), header=F) 
colnames(gtf) <- c("chr","start","end","strand","gene_id")
gtf <- gtf %>% mutate(gene_len = abs(end-start)) %>% select(gene_id,gene_len)

#create blank count table using first sample file to pull gene ids
counts <- read.csv(paste0(counts_dir,sample_list[1]), header=F, sep="\t") %>% dplyr::slice(-(1:4)) %>% select(1)
colnames(counts) <- c("gene_id")

#populate raw and tpm counts table
counts.raw <- counts
counts.tpm <- counts
for (i in sample_list){
  print(i) 
  counts.i <- read.csv(paste0(counts_dir,i), header=F, sep="\t") %>% dplyr::slice(-(1:4))
  colnames(counts.i) <- c("gene_id","count")
  
  #add tpm column
  counts.i <- left_join(counts.i, gtf, by = "gene_id") %>%
    mutate(rpk = count / (1000*gene_len))
  rpk_sum <- sum(counts.i$rpk)
  counts.i <- counts.i %>% 
    mutate(tpm = 1000000*rpk / rpk_sum) %>%
    select(gene_id,count,tpm)
  
  #add raw count and tpm data to dfs
  i <- str_replace(i, ".ReadsPerGene.tab","")
  
  counts.raw.i <- counts.i %>% select(gene_id,count)
  colnames(counts.raw.i) <- c("gene_id",i)
  counts.raw <- left_join(counts.raw,counts.raw.i, by = "gene_id")
  
  counts.tpm.i <- counts.i %>% select(gene_id,tpm)
  colnames(counts.tpm.i) <- c("gene_id",i)
  counts.tpm <- left_join(counts.tpm,counts.tpm.i, by = "gene_id")
  
}

# switch to long and add meta data
counts.raw.23 <- counts.raw %>%
  pivot_longer(2:ncol(counts.raw), 
               names_to = c("rep","sample"), 
               names_pattern = "(.*?)_(.*?)_.*",
               values_to = "expression") %>%
  mutate(rep = gsub(x = rep, pattern = "P1", replacement = "R2")) %>%
  mutate(rep = gsub(x = rep, pattern = "P2", replacement = "R3")) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))

counts.tpm.23 <- counts.tpm %>%
  pivot_longer(2:ncol(counts.tpm), 
               names_to = c("rep","sample"), 
               names_pattern = "(.*?)_(.*?)_.*",
               values_to = "expression") %>%
  mutate(rep = gsub(x = rep, pattern = "P1", replacement = "R2")) %>%
  mutate(rep = gsub(x = rep, pattern = "P2", replacement = "R3")) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))


##########################################
### COMBINE Reps and export ##########
##########################################

counts.raw <- rbind(counts.raw.1, counts.raw.23)
counts.raw$sample <- factor(as.numeric(counts.raw$sample))
counts.raw$rep <- factor(counts.raw$rep)
saveRDS(counts.raw, "data/counts.raw.rds")

counts.tpm <- rbind(counts.tpm.1, counts.tpm.23)
counts.tpm$sample <- factor(as.numeric(counts.tpm$sample))
counts.tpm$rep <- factor(counts.tpm$rep)
saveRDS(counts.tpm, "data/counts.tpm.rds")




