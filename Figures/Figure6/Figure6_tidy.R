library(tidyverse)
library("ZamanianLabThemes")

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure6")

########################################################
####### (Intestine vs non-Intestine, Pharynx vs non-Pharynx) Seq Run 1
########################################################

RNAdata <- c("~/Box/ZamanianLab/SeqLibraries/Mapping/")

dir <- c("Q008/201218_BHMLMFDSXY/")
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
counts.raw <- counts.raw %>%
  pivot_longer(2:ncol(counts.raw), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

counts.tpm <- counts.tpm %>%
  pivot_longer(2:ncol(counts.tpm), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

# export
counts.raw$sample <- factor(counts.raw$sample)
saveRDS(counts.raw, "data/counts-ic-lcm.raw.rds")

counts.tpm$sample <- factor(counts.tpm$sample)
saveRDS(counts.tpm, "data/counts-ic-lcm.tpm.rds")


########################################################
####### (Intestine vs non-Intestine, Pharynx vs non-Pharynx) Seq Run 2
########################################################

dir <- c("Q008/210112_BHV7H5DMXX/")
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
counts.raw <- counts.raw %>%
  pivot_longer(2:ncol(counts.raw), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

counts.tpm <- counts.tpm %>%
  pivot_longer(2:ncol(counts.tpm), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

# export
counts.raw$sample <- factor(counts.raw$sample)
saveRDS(counts.raw, "data/counts-ic-lcm_b.raw.rds")

counts.tpm$sample <- factor(counts.tpm$sample)
saveRDS(counts.tpm, "data/counts-ic-lcm_b.tpm.rds")


########################################################
####### (Pharynx vs non-Pharynx) Additional Run
########################################################

dir <- c("Q009/210312_BHVN5VDSXY/")
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
counts.raw <- counts.raw %>%
  pivot_longer(2:ncol(counts.raw), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

counts.tpm <- counts.tpm %>%
  pivot_longer(2:ncol(counts.tpm), 
               names_to = c("sample"), 
               names_pattern = ".*?_(.*?)_.*",
               values_to = "expression")

# export
counts.raw$sample <- factor(counts.raw$sample)
saveRDS(counts.raw, "data/lcm_b.raw.rds")

counts.tpm$sample <- factor(counts.tpm$sample)
saveRDS(counts.tpm, "data/lcm_b.tpm.rds")







