library(tidyverse)
library(cowplot)
library(ggrepel)
library("ZamanianLabThemes")
library(paletteer)
library(DESeq2)
library(wesanderson)
library(viridis)
library(ggdendro)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
#BiocManager::install("ggtree")
#BiocManager::install("treeio")
library(ggtree)
library(treeio)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure4")
theme_zlab_white = function(base_size = 12, base_family = "Helvetica") {
  
  theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme( #top right bottom left
      # Specify axis options
      axis.line = element_blank(),
      #axis.line = element_line(color = "black"),
      axis.text.x = element_text(size = base_size*0.8, color = "black", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "black", lineheight = 0.9),  
      axis.ticks = element_line(color = "black", size  =  0.4),  
      axis.title.x = element_text(size = base_size, color = "black", margin = margin(5, 0, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "black", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "white"),  
      legend.key = element_rect(color = "white",  fill = "white"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "black"),  
      legend.title = element_text(size = base_size*0.8, hjust = 0, color = "black"),  #face = "bold",
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


####################################
## FIGURE 4
####################################


############
## Fig4.ill: illustration of RNAt
############

library(magick)
library(pdftools)
library(grConvert)
library(grImport2)

RNA.illustration <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/RNA_Tomography/Illustration_RNAt.pdf", density = 600)
Fig4.ill <- ggdraw() + 
  draw_image(RNA.illustration, scale = 1) 


############
## Read in run QC
############

run_stats <- readRDS("data/run_stats.rds") %>%
  filter(metric == "unique_perc")

#stats
r1.read.avg <- run_stats 
r1.read.avg$sample <- as.numeric(r1.read.avg$sample)
r1.read.avg <- r1.read.avg %>% filter(rep == "R1", sample > 11)
r1.read.m <- mean(r1.read.avg$value) #70%
r1.read.sd <- sd(r1.read.avg$value) #14%

############
## Robust genes and slices: Fig4.qc (cumulative gene, venn overlap) + stats (head-enriched overlap)
############

# read in raw gene counts
counts <- readRDS("data/counts.raw.rds")

# filter for rep and id robustly expressed genes and slices
quality_filter <- function (rep.filter) {
  counts <- readRDS("data/counts.raw.rds") %>% filter(rep == rep.filter) %>% arrange(sample)
  counts$sample <- factor(as.numeric(counts$sample))
  counts <- counts %>%
    dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
    dplyr::select(gene_id, rep_sample, expression) %>%
    pivot_wider(names_from = rep_sample, values_from = expression) %>%
    column_to_rownames(var = "gene_id")
  
  #gene filter 1 (require x reads across all samples)
  keep <- rowSums(counts) > 20
  counts <- counts[keep,]
  #gene filter 2 (require x reads in each of at least y samples)
  keep <- rowSums(counts > 10) >= 1
  counts <- counts[keep,]
  #surviving genes
  gene_list <- unique(rownames(counts)) 
  
  #slice filter (require x reads in each of at least y genes) - transpose
  counts.t <- t(counts)
  keep <- rowSums(counts.t > 10) >= 100
  slice_list <- colnames(counts[,keep])

  #export
  filtered <- list("gene_list" = gene_list, "slice_list" = slice_list)
  return(filtered)
}
gene_list_r1 <- quality_filter("R1")$gene_list #8870
gene_list_r2 <- quality_filter("R2")$gene_list #6253
gene_list_r3 <- quality_filter("R3")$gene_list #8138
gene_list_r13 <- intersect(gene_list_r1,gene_list_r3) #7676
gene_list_r123 <- intersect(gene_list_r13,gene_list_r2) #5810

#STAT: overlap with head-enriched transcripts from Figure1.R
#head-enriched total from Figure1.R (degs.head.list): 2406
#temp <- intersect(gene_list_r1,degs.head.list) #2375/2406
#temp <- intersect(gene_list_r123,degs.head.list) #2209/2406

slice_list_r1 <- quality_filter("R1")$slice_list #r1: slice 1-11
slice_list_r2 <- quality_filter("R2")$slice_list #r2: slice 1-8;15;21;41-43
slice_list_r3 <- quality_filter("R3")$slice_list #r3: slice 1-7;11;15 
slice_list_r12 <- append(slice_list_r1,slice_list_r2)
slice_list_r123 <- append(slice_list_r12, slice_list_r3)

#create venn diagram
#if (!require(devtools)) install.packages("devtools")
#devtools::install_github("gaospecial/ggVennDiagram")
library("ggVennDiagram")

x <- list(R1 = gene_list_r1, R2 = gene_list_r2, R3 = gene_list_r3)

venn.plot <- ggVennDiagram(x[1:3], label_alpha = 0, label_size=3, category.names = c("R1","R2","R3")) +
  ggplot2::scale_fill_gradient(low="white",high = "yellow") +
  ggplot2::scale_colour_manual(values = c("#E69F00", "#56B4E9" ,"#e600cb")) +
  theme(legend.position = "none")
venn.plot
Fig4.venn <- venn.plot

## Cumulative Gene Figure

# read in raw gene counts
counts <- readRDS("data/counts.raw.rds") %>%
  group_by(sample,rep,gene_id) %>%
  distinct(gene_id, .keep_all = TRUE) %>%
  ungroup()
counts$sample <- as.integer(counts$sample)
counts <- counts %>% 
  mutate(expressed = ifelse(expression >10, 1,0)) %>%
  filter(expressed != 0) 

#populate blank df with cumulative count (>0)
counts.cumul.final <- data.frame(sample=NA, rep=NA, cumul_count=NA)[numeric(0), ]
for (i in 1:48){
  counts.cumul.i <- counts %>%
    filter(sample <= i) %>%
    group_by(rep,gene_id) %>%
    distinct(rep,gene_id, .keep_all = TRUE) %>%
    group_by(rep) %>%
    summarise(cumul_count = sum(expressed)) %>%
    mutate(sample = i)
  counts.cumul.final <- rbind(counts.cumul.final,counts.cumul.i)
}

#plot
plot.gcount <- ggplot(counts.cumul.final, aes(x = sample, y = cumul_count, group = rep))+
  theme_zlab_white() + 
  geom_line(aes(colour = rep), size = 0.65, alpha=0.75) +
  ylab(expression(Cummulative~genes~detected)) + xlab("Cryosection position") +
  geom_hline(yintercept = 11866, linetype = "dashed", color = "black", 2) +
  geom_hline(yintercept = 8870, linetype = "dashed", color = "#E69F00", 2) +
  geom_hline(yintercept = 6253, linetype = "dashed", color = "#56B4E9", 2) +
  geom_hline(yintercept = 8138, linetype = "dashed", color = "#e600cb", 2) +
  annotate("text", x = 5, y = 12226, label = "11,896", color = "black", size = 3.5) + 
  annotate("text", x = 5, y = 9230, label = "8,900", color = "#E69F00", size = 3.5) + 
  annotate("text", x = 5, y = 6603, label = "6,283", color = "#56B4E9", size = 3.5) + 
  annotate("text", x = 5, y = 8498, label = "8,168", color = "#e600cb", size = 3.5) + 
  scale_fill_manual(values=c("#E69F00", "#56B4E9" ,"#e600cb"), 
                    name=" ",
                    labels=c("R1", "R2", "R3")) +
  scale_colour_manual(values=c("#E69F00", "#56B4E9", "#e600cb"), 
                    name=" ",
                    labels=c("R1", "R2", "R3")) 
plot.gcount
Fig4.gcount <- plot.gcount

#Fig 4 QC
Fig4.qc <- plot_grid(plot.gcount,venn.plot,ncol = 1,rel_heights = c(1,0.9))


############
## Supplemental Figure (cluster.samples): sample clustering (all reps)
############

# read in gene counts, filter for robust genes, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate z-score matrix
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(log2(df.m+1)),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate sample distances and cluster
sampleDists <- dist(t(counts), method = "euclidean") 
sclust.dist <- hclust(sampleDists, method="ward.D2")
ord <- sclust.dist$order #if you want to arrange samples by distance

#convert sample distance matrix to df and pivot long for plotting
sampleDists <- as.data.frame(as.matrix(sampleDists)) %>%
  tibble::rownames_to_column(var = "rep_sample_a") %>%
  tidyr::pivot_longer(!rep_sample_a,
               names_to = "rep_sample_b",
               values_to = "dist") %>%
  tidyr::separate(rep_sample_a, c("rep_a","sample_a")) %>%
  tidyr::separate(rep_sample_b, c("rep_b","sample_b"))
sampleDists$rep_a <- factor(sampleDists$rep_a, levels = c("R1","R2","R3"))
sampleDists$sample_a <- factor(as.integer(sampleDists$sample_a))
sampleDists$rep_b <- factor(sampleDists$rep_b, levels = c("R1","R2","R3"))
sampleDists$sample_b <- factor(as.integer(sampleDists$sample_b))

#plot sample distances
pal <- wes_palette("Zissou1", 10, type = "continuous")
cluster.ht <- ggplot(sampleDists, aes(x = sample_a, y = sample_b, fill = dist)) +
  theme_zlab_white() + xlab('') + ylab('') + labs(fill = "Distance") +
  geom_tile() + 
  theme(axis.text.x = element_text(angle=90)) +
  scale_fill_gradientn(colours = pal) +
  facet_grid(rep_a ~ rep_b)
cluster.ht

#Supplemental Figure (sample distances across reps)
ggsave("cluster.samples.pdf", cluster.ht, width = 24, height = 20, units = "in")


############
## Gene clustering (rep 1): Fig4.r1 (Rep 1 heatmap + dendogram + mapping rate) 
############

# read in gene counts, filter for robust genes, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate z-score matrix
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate sample distances and cluster
geneDists <- dist(counts, method = "euclidean") 
gclust.dist <- hclust(geneDists, method="ward.D2")
ord <- gclust.dist$order #if you want to arrange genes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

# plot gene heatmap (genes clustered by distance)
pal <- wes_palette("Zissou1", 10, type = "continuous")
heatmap.r1 <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colours = pal) +
  ylab("Cryosection position") + xlab(expression(paste(italic("Brugia malayi"), " genes"))) +
  guides(fill=guide_colourbar(title="z-score")) +
  scale_y_discrete(limits = rev(levels(counts$sample))) +
  theme_zlab_white() +
  annotate("rect", xmin = 5, xmax = 5805, ymin = 0.5, ymax = 37.5,
           alpha = 0, color = "#000066") +
  theme(
    #panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 8, color = "black", lineheight = 0.9),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.direction = "vertical", 
    legend.position = "left",
    plot.margin = unit(c(0.5, 0, 0.5, 0.5), "cm"))
#ggsave("cluster.pdf", heatmap.r1, width = 14, height = 7, units = "in")

# generate dendogram
dend <- ggdendrogram(gclust.dist, rotate = FALSE, size = 5) + 
  theme(axis.text.x = element_blank(),axis.text.y = element_blank())

#prep qc data and mark slices that should be filtered
run_stats.r1 <- run_stats %>% 
  filter(rep == "R1") %>%
  mutate(rep_sample = paste0(rep,"_",sample)) %>%
  mutate(bad_slice = ifelse(rep_sample %in% slice_list_r1,"no","yes"))

qc.r1 <- ggplot(run_stats.r1, aes(sample,100*value)) +
  geom_line(aes(group=rep, colour = bad_slice),size=0.75, show.legend = F, alpha = 0.75) +
  geom_point(aes(group=rep, colour = bad_slice),size=0.15, show.legend = F, alpha = 0.5) +
  scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("Mapping %") + xlab("") +
  coord_flip() + 
  scale_x_discrete(limits = rev(levels(counts$sample))) +
  scale_y_continuous(breaks=c(0,100),limits = c(0, 100)) +
  geom_hline(yintercept = 50, colour = "gray", linetype = "dashed") + 
  theme_zlab_white() +
  theme(
    panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    axis.text.x = element_text(size = 10, color = "black", lineheight = 0.9),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "top",
    plot.margin = unit(c(0.5, 0.75, 0.5, -0.75), "cm"))
qc.r1

# create Fig4b (R1 heatmap and qc) - fix plot.margin (top right bottom left)
Fig4.r1.bottom <- plot_grid(heatmap.r1,qc.r1, nrow = 1, labels = c('',''), align = "h", rel_widths = c(2,0.12), scale = 1) 
Fig4.r1 <- plot_grid(dend + theme(plot.margin = unit(c(0.25, 0.6, -1.05, 2.5), "cm")), Fig4.r1.bottom,ncol = 1, labels = c('',''), rel_heights = c(0.075,1), scale = 1) 
#ggsave("Fig4.r1.pdf", Fig4.r1, width =14, height = 7, units = "in")


############
## Fig4.grplots (cluster gene patterns R1 post-filtering of slices using log2(TPM))
############ 

# read in gene counts, filter for robust genes and slices, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r1) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate z-score matrix
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate sample distances and cluster
geneDists <- dist(counts, method = "euclidean") 
gclust.dist <- hclust(geneDists, method="ward.D2")
ord <- gclust.dist$order #if you want to arrange genes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

# cut up tree into cluster
cluster <- cutree(gclust.dist, k = 20)
cluster <- as.data.frame(cluster) %>%
  rownames_to_column("gene_id")

# join cluster id with normalized counts
counts <- left_join(counts, cluster, by = "gene_id")

Fig4.grplots <- ggplot(counts, aes(sample, expression, group = gene_id)) +
  geom_line(size = 0.3, alpha =0.1, colour = "#709CA2FF") +
  theme_zlab_white() +
  geom_hline(yintercept=0, color = "black") +
  stat_summary(aes(group=cluster), fun.data=mean_cl_normal, geom="smooth", colour = "#1D1F71FF", size = 1) +
  ylab("z-score (TPM)") + xlab("Cryosection position (anterior to posterior)") + 
  theme(
    strip.background = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0.1, "cm"),
    axis.text.x = element_blank(),
    strip.text.x = element_blank(),
    legend.position = "none"
  ) +
  facet_wrap(cluster ~ ., ncol = 4, scales ="free", labeller = labeller(cluster),strip.position = "left") 
#ggsave("Fig4.grplots.pdf", Fig4.grplots, width = 10, height = 8, units = "in")

#save supplemental file of clustered genes (cluster df)
counts <- counts %>% select(gene_id, cluster) %>%
  distinct(gene_id, .keep_all = TRUE)
write.csv(counts,"RNAtomography_clustering.csv", row.names = FALSE)




############
## Fig4.r2 - Gene clustering (rep 2) to make Fig4.r23
############ 

# read in gene counts, filter for robust genes, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R2") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate z-score matrix
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate sample distances and cluster
geneDists <- dist(counts, method = "euclidean") 
gclust.dist <- hclust(geneDists, method="ward.D2")
ord <- gclust.dist$order #if you want to arrange genes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

# plot gene heatmap (genes clustered by distance)
pal <- wes_palette("Zissou1", 10, type = "continuous")
heatmap.r2 <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colours = pal) +
  ylab("R2") + xlab("") +
  guides(fill=guide_colourbar(title="z-score")) +
  scale_y_discrete(limits = rev(levels(counts$sample))) +
  theme_zlab_white() +
  theme(
    #panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.direction = "vertical", 
    legend.position = "none",
    plot.margin = unit(c(0.5, 0, -0.5, 0.5), "cm"))

#prep qc data and mark slices that should be filtered
run_stats.r2 <- run_stats %>% 
  filter(rep == "R2") %>%
  mutate(rep_sample = paste0(rep,"_",sample)) %>%
  mutate(bad_slice = ifelse(rep_sample %in% slice_list_r2,"no","yes"))

qc.r2 <- ggplot(run_stats.r2, aes(sample,100*value)) +
  geom_line(aes(group=rep, colour = bad_slice),size=0.5, show.legend = F, alpha = 0.75) +
  geom_point(aes(group=rep, colour = bad_slice),size=0.1, show.legend = F, alpha = 0.5) +
  scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("") + xlab("") +
  coord_flip() + 
  scale_x_discrete(limits = rev(levels(counts$sample))) +
  scale_y_continuous(breaks=c(0,100),limits = c(0, 100)) +
  geom_hline(yintercept = 50, colour = "gray", linetype = "dashed") + 
  theme_zlab_white() +
  theme(
    panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, -0.5, -0.75), "cm"))
qc.r2

#Fig4.r2
Fig4.r2 <- plot_grid(heatmap.r2,qc.r2, nrow = 1, labels = c('',''), align = "h", rel_widths = c(2,0.2), scale = 1) 


############
## Fig4.r3 - Gene clustering (rep 3) to make Fig4.r23
############ 

# read in gene counts, filter for robust genes, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R3") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate z-score matrix
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

#calculate sample distances and cluster
geneDists <- dist(counts, method = "euclidean") 
gclust.dist <- hclust(geneDists, method="ward.D2")
ord <- gclust.dist$order #if you want to arrange genes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

# plot gene heatmap (genes clustered by distance)
pal <- wes_palette("Zissou1", 10, type = "continuous")
heatmap.r3 <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colours = pal) +
  ylab("R3") + xlab("") +
  guides(fill=guide_colourbar(title="z-score")) +
  scale_y_discrete(limits = rev(levels(counts$sample))) +
  theme_zlab_white() +
  theme(
    #panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.direction = "vertical", 
    legend.position = "none",
    plot.margin = unit(c(0.5, 0, -0.5, 0.5), "cm"))

#prep qc data and mark slices that should be filtered
run_stats.r3 <- run_stats %>% 
  filter(rep == "R3") %>%
  mutate(rep_sample = paste0(rep,"_",sample)) %>%
  mutate(bad_slice = ifelse(rep_sample %in% slice_list_r3,"no","yes"))

qc.r3 <- ggplot(run_stats.r3, aes(sample,100*value)) +
  geom_line(aes(group=rep, colour = bad_slice),size=0.5, show.legend = F, alpha = 0.75) +
  geom_point(aes(group=rep, colour = bad_slice),size=0.1, show.legend = F, alpha = 0.5) +
  scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("") + xlab("") +
  coord_flip() + 
  scale_x_discrete(limits = rev(levels(counts$sample))) +
  scale_y_continuous(breaks=c(0,100),limits = c(0, 100)) +
  geom_hline(yintercept = 50, colour = "gray", linetype = "dashed") + 
  theme_zlab_white() +
  theme(
    panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0.5, 0.5, -0.5, -0.75), "cm"))
qc.r3

# Fig4.r3
Fig4.r3 <- plot_grid(heatmap.r3,qc.r3, nrow = 1, labels = c('',''), align = "h", rel_widths = c(2,0.2), scale = 1) 

# Fig4.r23: combine Fig4.r2 / r3 (top right bottom left)
Fig4.r23 <- plot_grid(Fig4.r2 + theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")), 
                      Fig4.r3 + theme(plot.margin = unit(c(-0.5, 0, -0.5, 0), "cm")), 
                      nrow = 2, labels = c('',''), align = "h", rel_heights = c(1,1), scale = 1) 


#########################
##### Figure 4 export ##########
#########################


Fig4.T <- plot_grid(Fig4.ill, Fig4.r1, nrow = 1, labels = c('A','B'), align = "h", rel_widths = c(0.25,1), scale = 1) 
Fig4.M <- plot_grid(Fig4.r23, Fig4.gcount, Fig4.venn , nrow = 1, labels = c('C','D','E'), align = "h", rel_widths = c(0.4,0.35,0.25), scale = 1) 
Fig4.B <- plot_grid(Fig4.grplots, nrow = 1, labels = c('F'), align = "h", scale = 1) 


Fig4 <- plot_grid(Fig4.T, Fig4.M, Fig4.B, ncol = 1, align = "v", rel_heights = c(1,0.6,0.7), scale = 0.99) 
ggsave("Fig4.pdf", Fig4, width = 16, height = 15, units = "in")







####################################
## FIGURE 5
####################################


############
## Fig5.ill: horizontal RNAt
############

library(magick)
library(pdftools)
library(grConvert)
library(grImport2)

RNA.illustration <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/RNA_Tomography/Illustration_RNAt_horizontal.pdf", density = 600)
Fig5.ill <- ggdraw() + 
  draw_image(RNA.illustration, scale = 1) 


##########################################
### Load R1 counts data, filter for core (R1,2,3) gene list, show highly expressed genes across combined cryosections
##########################################

# read in gene counts (TPM), filter for robust slices and genes
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>% 
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r1) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  separate(rep_sample, c("rep","sample"))
counts$sample <- as.integer(counts$sample)

# load in and join gene metadata
Bma.TM <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bma.TM.csv",
                   header = FALSE, sep = ",") %>%
  dplyr::rename("protein_id" = V1, "gene_id" = V2, "gene_name" = V3, "TM_count" = V4) %>% select(-TM_count) %>%
  distinct(gene_id, .keep_all=TRUE) %>%
  mutate(protein_id = str_remove(protein_id, "a")) %>%
  mutate(gene_name = str_replace(gene_name, "Bma-", "")) %>%
  mutate(gene_name = ifelse(is.na(gene_name),protein_id,gene_name)) %>%
  mutate(gene_name = ifelse(gene_name == "Bm4233","val-1",
                            ifelse(gene_name == "Bm2001","nlp-77",
                                   ifelse(gene_name == "Bm14007","rpl-30",
                                          ifelse(gene_name == "Bm2632","mlc-1/2",gene_name)))))
counts <- left_join(counts, Bma.TM, by = "gene_id") %>%
  dplyr::filter(!is.na(protein_id))


#sample range
min_sample <- as.integer(min(counts$sample))
max_sample <- as.integer(max(counts$sample))

#assigning to sets of 5 and finding highest genes that are distinctive (<5 groups)
counts.5 <- counts %>%
  mutate(gr=cut(sample, breaks= seq(min_sample-1, max_sample, by = 5))) %>%
  group_by(gr,gene_id) %>%
  mutate(mean.expression = mean(expression)) %>%
  distinct(gr,gene_id,.keep_all=TRUE) %>% select(-sample,-expression)

counts.5$gr.n <- as.integer(counts.5$gr)
gr.last <- as.integer(max(counts.5$gr.n, na.rm = TRUE)) + 1
counts.5 <- counts.5 %>%
  mutate(gr.n = ifelse(is.na(gr.n),gr.last,gr.n)) 
  
counts.5 <- counts.5 %>%
  group_by(gr.n) %>%
  slice_max(order_by = mean.expression, n = 100) %>% #grab top 100 genes for each group
  ungroup() %>% 
  group_by(gene_id) %>%
  mutate(n = n()) %>%
  filter(n < 5) %>% #keep if gene is found in fewer than 5 other groups
  ungroup() %>%
  group_by(gr.n) %>%
  slice_max(order_by = mean.expression, n = 10) %>% #grab top 10 per group
  ungroup() %>%
  mutate(group_size = "five")

#assigning to sets of 10 and finding highest genes that are distinctive (<4 groups)
counts.10 <- counts %>%
  mutate(gr=cut(sample, breaks= seq(min_sample-1, max_sample, by = 10))) %>%
  group_by(gr,gene_id) %>%
  mutate(mean.expression = mean(expression)) %>%
  distinct(gr,gene_id,.keep_all=TRUE) %>% select(-sample,-expression)

counts.10$gr.n <- as.integer(counts.10$gr)
gr.last <- as.integer(max(counts.10$gr.n, na.rm = TRUE)) + 1
counts.10 <- counts.10 %>%
  mutate(gr.n = ifelse(is.na(gr.n),gr.last,gr.n)) 

counts.10 <- counts.10 %>%
  group_by(gr.n) %>%
  slice_max(order_by = mean.expression, n = 100) %>% #grab top 100 genes for each group
  ungroup() %>% 
  group_by(gene_id) %>%
  mutate(n = n()) %>%
  filter(n < 4) %>% #keep if gene is found in fewer than 4 other groups
  ungroup() %>%
  group_by(gr.n) %>%
  slice_max(order_by = mean.expression, n = 20) %>% #grab top 20 per group
  ungroup() %>%
  mutate(group_size = "ten")

#top genes across entire tomography run
counts.all <- counts %>%
  mutate(gr = NA) %>%
  group_by(gene_id) %>%
  mutate(mean.expression = mean(expression)) %>%
  distinct(gene_id,.keep_all=TRUE) %>% select(-sample,-expression) %>%
  ungroup() %>%
  slice_max(order_by = mean.expression, n = 20) %>% #grab top 20
  ungroup() %>%
  mutate(n = 1, gr.n = 1, group_size = "all") 

#join all tilings
counts.tiles <- rbind(counts.5,counts.10, counts.all)
counts.tiles$group_size <- factor(counts.tiles$group_size, levels = c("five","ten","all"))
counts.tiles <- counts.tiles %>%
  group_by(group_size) %>%
  mutate(width.scale = 10/max(gr.n)) %>% ungroup()

#plot tiling of grouped sections with highly-expressed markers
tile.plot <- ggplot(counts.tiles) + #
  geom_rect(data = filter(counts.tiles,group_size == "five"), 
            mapping=aes(xmin=(gr.n*5)-4, xmax=(gr.n*5)+0.5, ymin=0, ymax=0.5),fill = "#01295F", color="white", alpha=0.1) +
  geom_rect(data = filter(counts.tiles,group_size == "ten"), 
            mapping=aes(xmin=(gr.n*10)-9, xmax=(gr.n*10)+0.5, ymin=0, ymax=0.5), fill = "#437F97", color="white", alpha=0.1) +
  geom_rect(data = filter(counts.tiles,group_size == "all"), 
            mapping=aes(xmin=1, xmax=40.5, ymin=0, ymax=0.5), fill = "#849324", color="white", alpha=0.1) +
  #geom_text(color = "white") + # add white text in the middle
  #scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  theme_void() + 's '
  theme(strip.background = element_blank(),  
        strip.text = element_text(color = "transparent")) +
  facet_grid (group_size ~ .)# remove any axis markings
tile.plot



library(ggtext)
#plot tiling of grouped sections with highly-expressed markers
tile.plot <- ggplot(counts.tiles) + #
  geom_textbox(data = filter(counts.tiles,group_size == "five"), aes(x=((gr.n*5)-4), y = 0, label = gr.n),
               width = grid::unit(0.12, "npc"), fill = "#01295F", color="white", alpha=0.5, hjust = 0, vjust = 1, halign = 0.5) +
  geom_textbox(data = filter(counts.tiles,group_size == "ten"), aes(x=((gr.n*5)-4), y = 0, label = gr.n),
               width = grid::unit(0.12, "npc"), fill = "#01295F", color="white", alpha=0.5, hjust = 0, vjust = 1, halign = 0.5) +
  
  #geom_richtext(data = filter(counts.tiles,group_size == "five"), aes(x=((gr.n*5)-4)/10, y = 0, label = gene_id),
  #              angle = 0, fill = "#01295F", color="white", alpha=0.1) +
  geom_richtext(data = filter(counts.tiles,group_size == "ten"), aes(x=((gr.n*10)-9)/10, y = 0, label = gene_id),
                angle = 0, fill = "#01295F", color="white", alpha=0.1) +
  #geom_text(color = "white") + # add white text in the middle
  #scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  theme_void() + 
  xlim(0,42) +
  theme(strip.background = element_blank(),  
        strip.text = element_text(color = "transparent")) +
  facet_grid (group_size ~ .)# remove any axis markings
tile.plot






#########################
##### Figure 5 export
#########################

Fig5.L <- plot_grid(Fig5.ill, NA, tile.plot, NA, NA, ncol = 1, labels = c('A','','','','B'), rel_heights = c(0.25,0.01,1,0.01,0.6), scale = 1) 
Fig5.L
ggsave("Fig5L.pdf", Fig5.L, width = 14, height = 10, units = "in")


Fig5.R <- NA

Fig5 <- plot_grid(Fig5.L, NA, Fig5.R, nrow = 1, rel_widths = c(1,0.01,0.25), scale = 0.99) 
ggsave("Fig5.pdf", Fig5, width = 14, height = 10, units = "in")








tile.plot.10 <- ggplot(counts.10,aes(gr.n, y = 0, fill = gr.n, label = gene_name)) +
  geom_tile(width = .9, height = .9) + # make square tiles
  geom_text(color = "white") + # add white text in the middle
  scale_fill_identity(guide = "none") + # color the tiles with the colors in the data frame
  coord_fixed() + # make sure tiles are square
  theme_void() # remove any axis markings


##########################################
### Load Antigens, Proteome, Drug Targets and join into gene_list df
##########################################

antigens <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Antigens.csv",
                     header = TRUE, sep = ",") %>%
  dplyr::slice(1:36) %>% select(gene_id,gene_name,Vaccine) %>% mutate(antigen = 1)

drug_targets <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Drug_Targets.csv",
                         header = TRUE, sep = ",") %>%
  mutate(gene_name = ifelse(is.na(name),transcript_id,name)) %>%
  select(gene_id,gene_name,type,subtype)%>% mutate(drug_target = 1)

gene_list <- full_join(antigens,drug_targets,by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name.x),gene_name.y,gene_name.x)) %>%
  select(-`gene_name.x`,-`gene_name.y`) %>%
  mutate(gene_name = str_replace(gene_name, "Bma-", "")) %>%
  filter(is.na(subtype) | subtype != "Chemosensory")

gene_list.antigens <- gene_list %>% filter(antigen == 1)
gene_list.antigens <- unique(gene_list.antigens$gene_id)

gene_list.drug <- gene_list %>% filter(drug_target == 1)
gene_list.drug <- unique(gene_list.drug$gene_id)


##########################################
### Load R1 counts data, filter for gene_list and cluster (pasted from grplots)
##########################################

# read in gene counts, filter for robust slices and genes of interest, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>% 
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r1) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate matrix (no scale normalization)
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  df.m <- log2(df.m+1) #df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)




#calculate sample distances and cluster
geneDists <- dist(counts, method = "euclidean") 
gclust.dist <- hclust(geneDists, method="ward.D2")
ord <- gclust.dist$order #if you want to arrange genes by distance

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]))
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$rep <- factor(counts$rep, levels = c("R1","R2","R3"))
counts$sample <- factor(as.integer(counts$sample))

# cut up tree into cluster
cluster <- cutree(gclust.dist, k = 10)
cluster <- as.data.frame(cluster) %>%
  rownames_to_column("gene_id")

# join cluster id with normalized counts
counts <- left_join(counts, cluster, by = "gene_id")

grplots <- ggplot(counts, aes(sample, expression, group = gene_id)) +
  geom_line(aes(colour = gene_id), size = 0.3, alpha =0.4) +
  theme_zlab_white() +
  geom_hline(yintercept=0, color = "black") +
  stat_summary(aes(group=cluster), fun.data=mean_cl_normal, geom="smooth", colour = "#1D1F71FF", size = 0.5) +
  ylab("z-score (TPM)") + xlab("Cryosection position (anterior to posterior)") + 
  theme(
    strip.background = element_blank(),
    #axis.ticks.y = element_blank(),
    #axis.text.y = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.ticks.length.x = unit(0.1, "cm"),
    #axis.text.x = element_blank(),
    #strip.text.x = element_blank(),
    legend.position = "right"
  ) +
  facet_wrap(cluster ~ ., ncol = 4, scales ="free", labeller = labeller(cluster),strip.position = "left") 
#ggsave("Fig4.grplots.pdf", grplots, width = 10, height = 8, units = "in")


#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(!gene_id, names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$sample <- factor(as.integer(counts$sample))

#join with gene_list
counts <- left_join(counts,gene_list, by = "gene_id")





########################
## Fig5.drug (known and putative drug target patterns) / mark if head-enriched
########################

# reinstalling all of the packages from the source code
#remotes::install_github("YuLab-SMU/tidytree", force = TRUE)
#remotes::install_github("YuLab-SMU/ggtree", force = TRUE)
#remotes::install_github("tidyverse/dplyr@v1.0.5", force = TRUE)

#loading the packages
library(ggtree)
library(ape)
library(dplyr)
library(stringr)

#load in Bma and Cel ids
Bma.id <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bma.Proteins.csv",
                      header = FALSE, sep = ",") %>% 
  dplyr::rename("protein_id"=V1, "gene_id"=V2, "gene_name"=V3) %>%
  group_by(gene_id,gene_name) %>%
  distinct(gene_id, .keep_all=TRUE)
Bma.protein.list <- unique(Bma.id$protein_id)

Cel.id <- read.csv("/Users/mzamanian/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Cel.Proteins.csv",
                   header = FALSE, sep = ",") %>% 
  dplyr::rename("protein_id"=V1, "gene_id"=V2, "gene_name"=V3) %>%
  group_by(gene_id,gene_name) %>%
  distinct(gene_id, .keep_all=TRUE)
Cel.protein.list <- unique(Cel.id$protein_id)

# read in iqtree file
lgic.phylo <- read.tree("~/Desktop/phylo/tree/LGIC_trim_final.aln.treefile")

#convert tree (phylo object) to tibble and add species column
x <- as_tibble(lgic.phylo) %>% 
  mutate(species = ifelse(label %in% Bma.protein.list,"Bma","Cel"))
x$species[153:302] <- NA

#add useful columns for Bma and Cel
Bma.id <- Bma.id %>% dplyr::rename("label" = protein_id)
Cel.id <- Cel.id %>% dplyr::rename("label" = gene_id)
x <- left_join(x, Bma.id, by = "label") 
x <- left_join(x, Cel.id, by = "label") 
x <- x %>%
  mutate(gene_id = ifelse(species == "Cel",label,gene_id)) %>%
  mutate(protein_id = ifelse(species == "Bma", label, protein_id)) %>%
  mutate(gene_name = ifelse(species == "Bma", gene_name.x, NA)) %>%
  mutate(gene_name = ifelse(species == "Cel", gene_name.y, gene_name)) %>%
  select(-gene_name.x,-gene_name.y)
x <- x %>%
  mutate(new_label = ifelse(is.na(gene_name),protein_id,gene_name)) 
  
#replace label with new label (gene name and protein id as backup)
lgic.phylo = rename_taxa(lgic.phylo, x, label, new_label)
#convert tibble back to tree (phylo object)
#lgic.phylo <- as.phylo(x)

### get coexpression data
#get all brugia gene ids
lgics <- x %>% filter(species == "Bma")
lgic.list <- unique(lgics$gene_id)

# read in gene counts, filter for robust slices and genes of interest, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% lgic.list) %>% 
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r1) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  pivot_wider(names_from = rep_sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# generate matrix (no scale normalization)
matrix_heatmap <- function (df) {
  df.m <- data.matrix(df, rownames.force = TRUE)
  ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
  df.m <- df.m[!ind,]
  #df.m <- log2(df.m+1) 
  df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
  return(df.m)
}
counts <- matrix_heatmap(counts)

# correlation matrix on transposed counts (columns = genes)
library(Hmisc)
counts.t <- t(counts)
cor <- rcorr(counts.t, type = c("pearson","spearman"))
cor.r <- as.data.frame(cor$r) %>%
  rownames_to_column(var = "gene_id") %>%
  pivot_longer(cols=2:51, names_to = "gene_id_2", values_to = "corr") %>%
  transmute(from = pmin(gene_id, gene_id_2), to = pmax(gene_id, gene_id_2), corr) %>%
  distinct() 

#STUCK HERE
cor.r <- left_join(cor.r,lgics,by="gene_id") %>%
  select(-gene_id) %>%
  dplyr::rename("from" = gene_name, "gene_id" =y)
cor.r <- left_join(cor.r,lgics,by="gene_id") %>%
  select(-gene_id) %>%
  dplyr::rename("to" = gene_name)
cor.r <- cor.r %>% filter(from != to) %>%
  distinct(corr,from,to)

  


#data on LGIC connections
dat <- data.frame(from=c("WBGene00221971", "Y71D11A.5", "K11G12.7"), 
                  to=c("WBGene00228955", "F55D10.5", "F55D10.5"), 
                  h=c(1, 1, 0.1), 
                  type=c("t1", "t2", "t3"), 
                  s=c(2, 1, 2))

dat <- data.frame(from=cor.r$from, 
                  to=cor.r$to, 
                  h=rep(c(1),length(cor.r$from)),
                  type=rep(c("t1"), length(cor.r$from)), 
                  s=rep(c(2), length(cor.r$from)))
dat <- dat %>% slice_head(n=1)



#LGIC phylogenetic tree 
lgic.tree <- ggtree(lgic.phylo, layout="circular", branch.length="none") + 
  geom_treescale() + theme_tree() + 
  geom_tiplab(size=2.5) +
  geom_hilight(mapping=aes(subset = node %in% c(38, 48, 58, 36),
                           node = node,
                           fill = as.factor(node))) +
  #geom_taxalink(data=cor.r, mapping=aes(taxa1=x, taxa2=y), 
       #         alpha=1, offset=0.05, ncp=5) +
 geom_taxalink(data=dat, mapping=aes(taxa1=from, taxa2=to), 
               alpha=1, offset=0.05, ncp=5)
lgic.tree





#DO ALL THE SAME FOR C ELEGANS (gene names...)
#COLLAPSE isoforms (a/b/c truncate as new column in x...etc)







Fig5.drug <- ggplot(data=filter(counts,drug_target == 1 ), aes(sample,log2(expression+1))) +
  geom_line(aes(group=gene_id, colour = type),size=0.5, show.legend = F, alpha = 0.5) +
  #geom_point(aes(group=rep, colour = bad_slice),size=0.1, show.legend = F, alpha = 0.5) +
  #scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("log2(TPM + 1)") + xlab("Cryosection position") +
  theme_zlab_white() +
  theme(
    #panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    legend.position = "none") +
  facet_grid(type ~ .)
Fig5.drug


############
## Fig5.antigen (antigen-vaccine targets) / mark if head-enriched
############ 

Fig5.antigen <- ggplot(data=filter(counts,antigen == 1 ), aes(sample,log2(expression+1))) +
  geom_line(aes(group=gene_name, colour = gene_name),size=0.5, show.legend = F, alpha = 0.5) +
  #geom_point(aes(group=rep, colour = bad_slice),size=0.1, show.legend = F, alpha = 0.5) +
  #scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("log2(TPM + 1)") + xlab("Cryosection position") +
  theme_zlab_white() +
  theme(
    #panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    legend.position = "none")
Fig5.antigen


Fig5.antigen <- ggplot(data=filter(counts, antigen == 1 ), aes(sample,expression)) +
  geom_line(aes(group=gene_id, colour = Vaccine),size=0.5, show.legend = F, alpha = 0.5) +
  #geom_point(aes(group=rep, colour = bad_slice),size=0.1, show.legend = F, alpha = 0.5) +
  #scale_colour_manual(values=c("no"="#000066","yes"="#CC0033")) +
  ylab("z-score") + xlab("Cryosection position") +
  theme_zlab_white() +
  theme(
    #panel.border = element_rect(color = "gray", fill = NULL, size = NULL, linetype = NULL),
    legend.position = "none") +
  facet_grid(Vaccine ~ .)
Fig5.antigen



############
## Fig5.lgic (subunit correlation - all reps)
############ 








#########################
##### Figure 5 export ##########
#########################




















#### JUNK BELOW




### OVERALL SAMPLE DISTANCE / CLUSTERING ANALYSIS USING DESEQ ##########

# read in gene counts (raw)
counts <- readRDS("data/counts.raw.rds") %>%
  filter(rep == "R1") %>% select(-rep) %>%
  arrange(sample)
counts$sample <- factor(as.numeric(counts$sample))

# prep for matrix conversion
counts <- counts %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# make metadata file
sample_list <- colnames(counts)
samples <- data.frame("sample" = sample_list, stringsAsFactors = FALSE)
saveRDS(samples, "data/samples-rnat.rds")
samples <- readRDS("data/samples-rnat.rds")

# convert to matrix and create deseq2 object
gene_count <- as.matrix(counts)
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ 1)

#filter 1 (require x reads total across all samples)
keep <- rowSums(counts(dds)) > 20
dds <- dds[keep,]

#filter 2 (require x reads in each of at least y samples)
keep <- rowSums(counts(dds) >= 20) >= 1
dds <- dds[keep,]
nrow(dds) #11718 > 8491

#pick vst transform
vsd <- vst(dds, blind = FALSE)

temp <- t(assay(vsd))


#calculate sample distances and use hclust to cluster samples based on dist
sampleDists <- dist(t(assay(vsd)), method = "euclidean") 
sampleDists <- dist(t(counts), method = "euclidean") 
clust.dist <- hclust(sampleDists, method="ward.D2")

#get list order of clustered samples from hclust output
ord <- clust.dist$order

#convert original distance matrix to df, re-order, and tidy for plotting
sampleDists.df <- as.data.frame(as.matrix(sampleDists)) %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = 2:49,
               names_to = "sample_2",
               values_to = "dist")
sampleDists.df$sample <- factor(as.integer(sampleDists.df$sample))
sampleDists.df$sample_2 <- factor(as.integer(sampleDists.df$sample_2))

#plot heatmap cluster
pal <- wes_palette("Zissou1", 25, type = "continuous")
cluster.ht <- ggplot(sampleDists.df, aes(x = sample, y = sample_2, fill = dist)) +
  theme_zlab() + xlab('') + ylab('') + labs(fill = "") +
  geom_tile() +
  scale_fill_gradientn(colours = pal)
cluster.ht


########## TRANSCRIPT CLUSTERING ##########

# read in counts (vsd normalized) 
counts <- as.data.frame(assay(vsd)) 

# read in counts (tpm)
counts <- readRDS("data/counts.tpm.rds") %>%
  filter(rep == "R1") %>% arrange(sample)

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


# plot heatmap
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
heatmap <- ggplot(counts, aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_viridis() +
  #scale_fill_gradientn(colours = pal, "Z-Score",
 #                      guide = guide_colorbar(
  #                       direction = "horizontal",
 #                        title.position = "left",
 #                        label.position = "bottom"
 #                      )) +
 # scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
 # scale_x_discrete(labels = c("Body 1" = "","Body 2" = "Body","Body 3" = "",
                          #    "Head 1" = "","Head 2" = "Head","Head 3" = ""),
  #                 position = "top") +
  #geom_vline(xintercept = 3.5, linetype="dashed", color = "black", size=0.5) +
  ylab("Cryosection position") + xlab("B. malayi Genes") +
  theme_zlab_white() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(face = "italic", size = 8, angle =0, hjust=1),
    strip.background = element_rect(fill="white", color = "white"),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing=unit(-0.75, "lines"),
    legend.position = "right") 
ggsave("cluster.pdf", heatmap, width = 14, height = 7, units = "in")










