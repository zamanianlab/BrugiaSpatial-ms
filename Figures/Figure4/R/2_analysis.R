library(tidyverse)
library(cowplot)
library(ggdendro)
library(scales)
library(DESeq2)
library(pheatmap)
library(wesanderson)
library(dtw)
library(viridis)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure4")

theme_zlab_white = function(base_size = 12, base_family = "Helvetica") {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
      axis.line = element_blank(),
      #axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = base_size*0.8, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.4),  
      axis.title.x = element_text(size = base_size, color = "black", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "black"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "black"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "white", color  =  NA),  
      panel.border = element_rect(fill = NA, color = "black"),  
      #panel.border = element_blank(),
      #panel.grid.major = element_line(color = "black"),  
      panel.grid.major = element_blank(), 
      #panel.grid.minor = element_line(color = "black"),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify faceting options
      strip.background = element_rect(fill = "white", color = "white"),  
      strip.text.x = element_text(size = base_size*0.8, color = "black"),  
      strip.text.y = element_text(size = base_size*0.8, color = "black",angle = -90),  
      # Specify plot options
      plot.background = element_rect(colour = NA, fill = NA),  
      plot.title = element_text(size = base_size*1.2, color = "black"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
}

##########################################
### OVERALL CLUSTERING ANALYSIS ##########
##########################################

# read in gene counts (raw and tpm)
counts.raw <- readRDS("counts/counts.raw.rds") 
counts.tpm <- readRDS("counts/counts.tpm.rds") 

counts <- counts.tpm

# TPM threshold 
counts.thr <- counts %>%
  group_by(gene_id) %>%
  summarise(sum = sum(expression)) %>%
  filter(sum > 20)
counts <- counts %>% filter(gene_id %in% unique(counts.thr$gene_id))


matrix_heatmap <- function (df) {
  df <- df %>%
    mutate(rep_sample = paste0(rep,"_",sample)) %>%
    dplyr::select(gene_id, rep_sample, expression) %>%
    pivot_wider(names_from = rep_sample, values_from = expression) %>%
    column_to_rownames(var = "gene_id")
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(log2(df.m+1)),center=TRUE,scale=TRUE)) #or log2(df.m+1)
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate gene distances and use hclust to cluster samples based on dist
geneDists <- dist(counts, method = "euclidean")
gclust.dist <- hclust(geneDists, method="ward.D2")

#get list order of clustered genes from hclust output
ord <- gclust.dist$order

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

heatmap <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_viridis() +
  ylab("Cryosection position") + xlab("B. malayi Genes") +
  theme_zlab_white() +
  #scale_y_reverse() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.title.y = element_text(size = 34),
    #axis.title.x = element_text(size = 34),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.x = element_blank(),
    #axis.text.x = element_text(angle = 45, hjust = 1, size = 6), + #element_text(angle = 45, hjust = 1, size = 6)
    legend.title=element_blank()) +
  facet_grid(rep ~ ., scale = "free", space="free")
ggsave("cluster.pdf", heatmap, width = 12, height = 10, units = "in")







##########################################
### Analyze / cluster specific genes #########
##########################################

# read in gene counts (TPM)
counts.r1 <- readRDS("data/counts-r1.tpm.rds") %>% mutate(rep = "R1")
counts.r2 <- readRDS("data/counts-r2.tpm.rds") %>% mutate(rep = "R2")
counts.r3 <- readRDS("data/counts-r3.tpm.rds") %>% mutate(rep = "R3")
counts.r4 <- readRDS("data/counts-r4.tpm.rds") %>% mutate(rep = "R4")

counts <- rbind(counts.r1,counts.r2,counts.r3,counts.r4) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))
counts$sample <- as.factor(as.numeric(counts$sample))
counts <- counts %>% filter(rep != "R2")

# TPM threshold 
counts.thr <- counts %>%
  group_by(gene_id) %>%
  summarise(sum = sum(expression)) %>%
  filter(sum > 20)
counts <- counts %>% filter(gene_id %in% unique(counts.thr$gene_id))

#pull in genes of interest
gene_list <- read.csv("lgic_list.tsv", header=F, sep = "") 
gene_list <- gene_list$V3 

#filter for gene list
counts <- counts %>%
  filter(gene_id %in% gene_list)

matrix_heatmap <- function (df) {
  df <- df %>%
    mutate(rep_sample = paste0(rep,"_",sample)) %>%
    dplyr::select(gene_id, rep_sample, expression) %>%
    pivot_wider(names_from = rep_sample, values_from = expression) %>%
    column_to_rownames(var = "gene_id")
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(log2(df.m+1)),center=TRUE,scale=TRUE)) #or log2(df.m+1)
  return(df.m)
}
counts <- matrix_heatmap(counts)

# correlation matrix on transposed counts (columns = genes)
counts.t <- t(counts)

#install.packages("Hmisc")
library(Hmisc)
cor <- rcorr(counts.t, type = c("pearson","spearman"))
cor.r <- cor$r
write.csv(cor.r, "plots/lgic_corr.csv")

#plot correlation
cor.df <- as.data.frame(cor.r) %>%
  rownames_to_column("gene_id") %>%
  gather(gene_id2,corr,2:17)
  
heatmap <- ggplot(cor.df, aes(gene_id, gene_id2)) +
  geom_tile(aes(fill = corr)) +
  scale_fill_viridis() +
  ylab("") + xlab("") + theme(axis.ticks = element_blank()) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.title.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 22)) +
  #theme(axis.ticks.y = element_blank()) +
  #theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_text(angle = 0, hjust = 0.5, size = 11)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 11)) +
  theme(legend.title=element_blank()) 
ggsave("plots/corr.pdf", heatmap, width = 12, height = 10, units = "in")







#calculate gene distances and use hclust to cluster samples based on dist
geneDists <- dist(counts, method = "euclidean")
gclust.dist <- hclust(geneDists, method="ward.D2")

#get list order of clustered genes from hclust output
ord <- gclust.dist$order

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "sample", values_to = "expression")
counts$sample <- factor(counts$sample)

heatmap <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_viridis() +
  ylab("Cryosection position") + xlab("LGIC Genes") + theme(axis.ticks = element_blank()) +
  theme_bw() +
  theme(
    plot.background = element_blank()
    ,panel.grid.major = element_blank()
    ,panel.grid.minor = element_blank()
    ,panel.border = element_blank()
  ) +
  theme(axis.title.y = element_text(size = 22)) +
  theme(axis.title.x = element_text(size = 22)) +
  #theme(axis.ticks.y = element_blank()) +
  #theme(axis.text.y = element_blank()) +
  theme(axis.ticks.x = element_blank()) +
  theme(axis.text.x = element_blank()) +
  theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 16)) +
  theme(legend.title=element_blank()) + 
  facet_grid(rep ~ .)
ggsave("plots/cluster2.pdf", heatmap, width = 20, height = 12, units = "in")





##########################################
### Plot specific genes #########
##########################################

# read in gene counts (TPM)
counts.r1 <- readRDS("data/counts-r1.tpm.rds") %>% mutate(rep = "R1")
counts.r2 <- readRDS("data/counts-r2.tpm.rds") %>% mutate(rep = "R2")
counts.r3 <- readRDS("data/counts-r3.tpm.rds") %>% mutate(rep = "R3")
counts.r4 <- readRDS("data/counts-r4.tpm.rds") %>% mutate(rep = "R4")

counts <- rbind(counts.r1,counts.r2,counts.r3,counts.r4) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))
counts$sample <- as.factor(as.numeric(counts$sample))
counts <- counts %>% filter(rep != "R2")

#pull in genes of interest
gene_list <- read.csv("lgic_list.tsv", header=F, sep = "") 
gene_list <- gene_list$V3 

#filter for gene list
counts <- counts %>%
  filter(gene_id %in% gene_list)

gene_plot <- ggplot(counts, aes(sample, log2(expression+1))) +
  geom_line(aes(colour = gene_id, group = gene_id)) +
  scale_fill_viridis() +
  ylab("Gene Expression log2(TPM+1)") + xlab("Cryosection Position") + theme(axis.ticks = element_blank()) +
  theme_bw() +
  theme(axis.title.y = element_text(size = 14)) +
  theme(axis.title.x = element_text(size = 14)) +
  #theme(axis.ticks.y = element_blank()) +
  #theme(axis.text.y = element_blank()) +
  #theme(axis.ticks.x = element_blank()) +
  #theme(axis.text.x = element_blank()) +
  #theme(axis.text.x = element_text(angle = 45, hjust = 0.5, size = 16)) +
  theme(legend.title=element_blank()) +
  facet_grid(rep ~ .)
ggsave("plots/lgic_plot.pdf", gene_plot, width = 10, height = 7, units = "in")


### STOPPED HERE


# plot dendogram
dend <- ggdendrogram(hc, rotate = FALSE, size = 40) + 
  theme(axis.text.x = element_blank()) +
  theme(axis.text.y = element_blank()) 
ggsave("temp_dend.pdf", dend, width = 39, height = 10, units = "in")
ggsave("temp_dend_Henr.pdf", dend, width = 20, height = 3, units = "in")

#### cutting up the tree and plotting subtrees
group <- cutree(hc, k = 20)
groupdf <- as.data.frame(group)
groupdf$gene_id <- rownames(groupdf)
groupdf <- groupdf %>% dplyr::rename(Cluster = group)
rownames(groupdf) <- seq_along(1:nrow(groupdf))

#plot normalized expression values (non-scaled and non-log transformed)
expdf <- pd
rownames(expdf) <-  seq_along(1:nrow(expdf))
expdf <- expdf %>% select(-id)
mergedf <- merge(expdf,groupdf, by="gene_id")

mergedf.plot <- mergedf #%>%
#filter(Cluster < 200)

grplots <- ggplot(mergedf.plot, aes(sample, value, group = gene_id, color = as.factor(Cluster))) +
  geom_line(size = 0.8, alpha =0.3) +
  #geom_ribbon(aes(ymin = value-0.075, ymax = value+0.075), alpha =0.2) +
  stat_summary(aes(group=Cluster), fun.data=mean_cl_normal, geom="smooth") +
  stat_summary(aes(group=Cluster), fun.y=mean, geom="line", linetype="dashed", colour="black", size = 0.7) +
  ylab("Normalized Expression") + xlab("") + 
  theme(axis.ticks = element_blank()) +
  theme(legend.position = "none") + theme_bw() + theme(legend.position = "none") +
  theme(
    strip.background = element_blank(),
    strip.text.x = element_blank()
  ) +
  facet_grid(Cluster~., scales ="free")
#grplots
ggsave("temp_grplots.pdf", grplots, width = 25, height = 49.5, units = "in")





##########################################
###### DESEQ-2 Sample Clustering from raw counts
##########################################

# read in gene counts (raw
counts.r1 <- readRDS("data/counts-r1.raw.rds") %>% mutate(rep = "R1")
counts.r2 <- readRDS("data/counts-r2.raw.rds") %>% mutate(rep = "R2")
counts.r3 <- readRDS("data/counts-r3.raw.rds") %>% mutate(rep = "R3")
counts.r4 <- readRDS("data/counts-r4.raw.rds") %>% mutate(rep = "R4")

counts <- rbind(counts.r1,counts.r2,counts.r3,counts.r4) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))
counts$sample <- as.factor(as.numeric(counts$sample))
counts <- counts %>% filter(rep != "R2")

#prep 
gene_count <- counts %>% column_to_rownames(var = "gene_id")
gene_count <- as.matrix(gene_count)

# read in raw gene counts (rep 1)
gene_count <- as.matrix(read.csv("~/Box/ZamanianLab/SeqLibraries/Mapping/Q003/counts/gene_count_matrix.csv", row.names = "gene_id"))
colnames_g <- gsub(x = colnames(gene_count), pattern = "_S\\d+_L\\d+", replacement = "")
colnames(gene_count) <- colnames_g
colnames_g <- gsub(x = colnames(gene_count), pattern = "P\\d+_S", replacement = "R1_S")
colnames(gene_count) <- colnames_g
gene_count.1 <- as.data.frame(gene_count) %>%
  rownames_to_column(var = "gene_id")

# read in raw gene counts (rep 2)
gene_count <- as.matrix(read.csv("~/Box/ZamanianLab/SeqLibraries/Mapping/Q004/counts/gene_count_matrix.csv", row.names = "gene_id"))
colnames_g <- gsub(x = colnames(gene_count), pattern = "_S\\d+_L\\d+", replacement = "")
colnames(gene_count) <- colnames_g
colnames_g <- gsub(x = colnames(gene_count), pattern = "P\\d+_S", replacement = "R2_S")
colnames(gene_count) <- colnames_g
gene_count.2 <- as.data.frame(gene_count) %>%
  rownames_to_column(var = "gene_id")

# read in raw gene counts (rep 3)
gene_count <- as.matrix(read.csv("~/Box/ZamanianLab/SeqLibraries/Mapping/Q005/counts/gene_count_matrix.csv", row.names = "gene_id"))
colnames_g <- gsub(x = colnames(gene_count), pattern = "_S\\d+_L\\d+", replacement = "")
colnames(gene_count) <- colnames_g
colnames_g <- gsub(x = colnames(gene_count), pattern = "P\\d+_S", replacement = "R3_S")
colnames(gene_count) <- colnames_g
gene_count.3 <- as.data.frame(gene_count) %>%
  rownames_to_column(var = "gene_id")

# merge counts across 3 reps and convertto matrix
gene_count <- left_join(gene_count.1,gene_count.2, by = "gene_id")
gene_count <- left_join(gene_count,gene_count.3, by = "gene_id") %>%
  column_to_rownames(var = "gene_id")
gene_count <- as.matrix(gene_count)

# make sample file (rownames = colnames(gene_count))
colData <- data.frame("sample_id"=colnames(gene_count)) %>%
  mutate("rep_sample" = sample_id) %>%
  separate(rep_sample, c("rep", "sample"), sep = "_", remove = TRUE) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = ""))
colData$sample <- factor(as.numeric(colData$sample))
rownames(colData) = colData$sample_id
colData <- colData %>% select(-sample_id)

#use vst (variance stabilizing transform) for sample clustering
gene_count <- gene_count + 1 #add pseudocount to allow de-seq to take geom mean
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = colData, design = ~ rep)
vst = vst(dds, blind=FALSE)
vsd=assay(vst)
pheatmap(cor(vsd), cluster_rows=T, cluster_cols=T)

Rdist <- dist(t(assay(vst)))
RdistMatrix <- as.matrix(Rdist)

Rheatmap <- pheatmap(RdistMatrix,cluster_rows=F, cluster_cols=F)
ggsave("heatmapd.pdf", Rheatmap, width = 20, height = 20, units = c("in"))

Rheatmap_cluster <- pheatmap(RdistMatrix,cluster_rows=T, cluster_cols=T)
ggsave("heatmap_clusterd.pdf", Rheatmap_cluster, width = 20, height = 20, units = c("in"))



########################
#### Plot gene subsets
########################

rna_counts <- readRDS("~/Box/ZamanianLab/Data/RNAseq/Bmalayi_AF_RNAt/expression_tidy/BmAF_RNAt.rds")

#read in annotated receptors and plot
annot <- read.csv2("~/Box/ZamanianLab/Data/Genetics/orthologs.csv", sep = ",")
annot <- annot %>%
  select(1:7) %>%
  dplyr::rename(Gene_ID = Bmalayi_Gene_ID) %>%
  filter(Type == "LGIC")
#filter(Type =="GPCR", Subtype == "Aminergic") 
gene_list <- unique(annot$Gene_ID)
gene_list <- c("WBGene00223839", "WBGene00228311", "WBGene00222703", "WBGene00221971")

# read in gene counts (TPM)
rna_counts.plot <- rna_counts %>%
  filter(metric == "TPM_g") %>%
  dplyr::select(gene_id, proj, sample, expression) %>%
  filter(gene_id %in% gene_list)

### Plot
plot <- ggplot(rna_counts.plot) +
  aes(x = sample, y = log2(expression), group = proj, colour = proj) + 
  theme_bw() +
  geom_point(size=1, alpha = 0.5) +
  geom_line() +
  theme(axis.title.x=element_blank()) +
  facet_grid(gene_id~., scales ="free")
plot

ggsave("temp_geneplots.pdf", plot, width = 20, height = 49.5, units = "in")

#transcript_list <- c("Bm2001")
gene_list <- c("WBGene00224213") #ov-16
gene_list <- c("WBGene00225800") #ASTMH lead
gene_list <- c("WBGene00224994","WBGene00229959")  #bma-btub-1,bma-btub-2
gene_list <- c("WBGene00227604","WBGene00221972") #tax-4, osm-9
gene_list <- c("WBGene00221972")
gene_list <- c("WBGene00221971") #avr-14
gene_list <- c("WBGene00225593") #vha-1 homolog
gene_list <- c("WBGene00225597") #Bma-gpa-7
gene_list <- c("WBGene00223643","WBGene00269103") #dop-1,ser-2

gene_list <- c("WBGene00227131","WBGene00225421") #mif-1, cpi-1, avr-14
gene_list <- c("WBGene00227865","WBGene00228788", "WBGene00226900") #ace-1,2,3
gene_list <- c("WBGene00223839", "WBGene00228311", "WBGene00222703", "WBGene00221971") #glc-2,3,4,avr-14
gene_list <- c("WBGene00227085") #bma-lad2









