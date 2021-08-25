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
  #dplyr::filter(rep_sample %in% slice_list_r123) %>% #optional (cluster all vs QC'd)
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
#cluster.ht

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
## Fig4.grplots (cluster gene patterns R1 post-filtering of slices using scale-normalized TPM)
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

# add gene metadata
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

# group and order by single-cryo uptick vs diffuse uptick
cluster.markers <- counts %>%
  group_by(cluster) %>%
  slice_max(order_by = expression, n = 3) %>%
  select(cluster,gene_name,sample) %>% 
  dplyr::rename(max.slice = "sample") %>%
  dplyr::summarise(marker=paste(gene_name, collapse="  "), max.slice=max.slice) %>%
  distinct()
cluster.order <- cluster.markers %>% select(cluster, max.slice)
# discrete: 3,12,19,1,13,4,11,16,17,8,6,2 / diffuse: 9,10,14,20,15,18,5,7

# add marker information and order
counts <- left_join (counts, cluster.markers, by = "cluster")
counts$cluster <- factor(counts$cluster, levels = c("3","12","19","1","13","4","11","16","17","8","6","2",
                                                    "9","10","14","20","15","18","5","7"))
discrete <- c("3","12","19","1","13","4","11","16","17","8","6","2")
counts <- counts %>%
  mutate(type = ifelse(cluster %in% discrete, "discrete", "diffuse"))

counts.labels <- counts %>%
  filter(type == "discrete") %>%
  select(cluster,marker) %>% distinct()

# plot clusters with gene marker labels
Fig4.grplots <- ggplot(counts, aes(sample, expression)) +
  geom_line(size = 0.3, alpha =0.1, aes(colour = type, group = gene_id)) + # colour = "#709CA2FF"
  scale_color_manual(values=c(discrete="#78BC61",diffuse="#FFA500"))+
  theme_zlab_white() +
  geom_hline(yintercept=0, color = "black") +
  stat_summary(aes(group=cluster), fun.data=mean_cl_normal, geom="smooth", colour = "#1D1F71FF", size = 1) +
  geom_text(data = counts.labels, aes(x = 14, y = 5, label = marker, fontface = "italic"), 
           hjust = 0, size = 3, color = "black") +
  ylab("z-score") + xlab("Cryosection position (anterior to posterior)") + 
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
  facet_wrap(cluster ~ ., ncol = 4, scales ="free", labeller = labeller(cluster), strip.position = "left") 
ggsave("Fig4.grplots.pdf", Fig4.grplots, width = 12, height = 4, units = "in")

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
## Fig5.ill: horizontal RNAt illustration
############

library(magick)
library(pdftools)
library(grConvert)
library(grImport2)

RNA.illustration <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/RNA_Tomography/Illustration_RNAt_horizontal.pdf", density = 600)
Fig5.ill <- ggdraw() + 
  draw_image(RNA.illustration, scale = 1) 


##########################################
### Plot top expressed genes per 5 sections and antigen / drug target distributions
##########################################

# read in gene counts (TPM), filter for robust slices and genes
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::filter(gene_id %in% gene_list_r123) %>% 
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r123) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  separate(rep_sample, c("rep","sample"))
counts$sample <- as.integer(counts$sample)

#update cryosection indexing based on qc (R1 = 12, R2 = 9, R3 = 8)
counts <- counts %>%
  mutate(sample = ifelse(rep == "R1", sample - 11, ifelse(rep =="R2", sample - 8, ifelse(rep == "R3", sample - 7, sample))))

#sample range
min_sample <- as.integer(min(counts$sample))
max_sample <- as.integer(max(counts$sample))

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
                                          ifelse(gene_name == "Bm10669", "cpi-2",
                                              ifelse(gene_name == "Bm2632","mlc-1/2",gene_name))))))
counts <- left_join(counts, Bma.TM, by = "gene_id") %>%
  dplyr::filter(!is.na(protein_id))

# plot: coalescing in groups of 5 cryosections and finding highest expressed genes for each region
counts.gr <- counts %>%
  mutate(gr=cut(sample, breaks= seq(min_sample-1, max_sample, by = 5))) %>%
  group_by(gr,gene_id) %>%
  mutate(mean.expression = mean(expression)) %>%
  distinct(gr,gene_id,.keep_all=TRUE) %>% select(-sample,-expression)

counts.gr$gr.n <- as.integer(counts.gr$gr)
gr.last <- as.integer(max(counts.gr$gr.n, na.rm = TRUE)) + 1
counts.gr <- counts.gr %>%
  mutate(gr.n = ifelse(is.na(gr.n),gr.last,gr.n)) 

counts.gr <- counts.gr %>%
  group_by(gr.n) %>%
  slice_max(order_by = mean.expression, n = 10) %>% 
  mutate(rank=as.numeric(row_number())) %>%
  ungroup()

tile.plot <- ggplot(counts.gr) + 
  geom_rect(mapping=aes(xmin=gr.n-0.4, xmax=gr.n+0.4, ymin=0, ymax=10+1), fill = "gray95", color="white", size = 0, alpha = 0.2) +
  geom_vline(aes(xintercept = gr.n-0.5), linetype="dotted", color="#01295F", alpha=0.25) +
  geom_vline(aes(xintercept = 8.5), linetype="dotted", color="#01295F", alpha=0.25) +
  xlim(0.5,8.5) +
  #geom_line(aes(x=gr.n,y = rank, group = gene_name, colour = gene_name), alpha = 0.2) +
  #geom_text(aes(x = gr.n, y = rank, label = gene_name, colour = gene_name), size = 3, fontface = "italic") +
  geom_text(aes(x = gr.n, y = rank, label = gene_name), colour = "grey25", size = 3, fontface = "italic") +
  theme_zlab_white() +
  scale_y_reverse() +
  theme_void() +
  theme(
    strip.background =element_rect(fill="#01295F"),
    strip.text = element_text(colour = 'white'),
    panel.background = element_rect(fill = 'white', colour = 'white'),
    plot.background = element_rect(fill = 'white', colour = 'white'),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "none")
tile.plot


##########################################
### Plot drug target distributions
##########################################

# read in gene counts (TPM), filter for robust slices and genes
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r123) %>%
  dplyr::select(gene_id, rep_sample, expression) %>%
  separate(rep_sample, c("rep","sample"))
counts$sample <- as.integer(counts$sample)

#update cryosection indexing based on qc (R1 = 12, R2 = 9, R3 = 8)
counts <- counts %>%
  mutate(sample = ifelse(rep == "R1", sample - 11, ifelse(rep =="R2", sample - 8, ifelse(rep == "R3", sample - 7, sample))))

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
                                          ifelse(gene_name == "Bm10669", "cpi-2",
                                                 ifelse(gene_name == "Bm2632","mlc-1/2",gene_name))))))
counts <- left_join(counts, Bma.TM, by = "gene_id") %>%
  dplyr::filter(!is.na(protein_id))

# gpcrs, trps, lgics and merge with counts
drugtargets.1 <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Drug_Targets.csv",
                         header = TRUE, sep = ",") %>% filter(type %in% c("GPCR","TRP")) %>% select(gene_id, type,subtype)
drugtargets.2 <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bm.LGICs.txt",
                          header = FALSE, sep = ",") %>% dplyr::rename(protein_id = "V1") %>% mutate(type = "LGIC")
counts.dr <- left_join(counts, drugtargets.1, by="gene_id")
counts.dr <- left_join(counts.dr, drugtargets.2, by="protein_id")
counts.dr <- counts.dr %>% filter(!is.na(type.x) | !is.na(type.y))
counts.dr <- counts.dr %>%
  mutate(type = ifelse(is.na(type.x), type.y , type.x)) %>% select(-type.x,-type.y)

# mark which ones are highest in ES pore, nerve ring, and vulva regions
counts.dr.max <- counts.dr %>% 
  group_by(gene_name) %>% #identify slice with highest TPM
  slice_max(order_by = expression, n = 1) %>%
  ungroup() %>%
  select(gene_id,type,sample,expression) %>%
  dplyr::rename(sample.max = "sample", expression.max = "expression") %>%
  group_by(type) %>% # order (rank) the genes in each type by location of sample max (max slice)
  arrange(sample.max) %>%
  mutate(rank=as.numeric(row_number())) %>%
  ungroup() %>%
  filter(expression.max > 10) %>% select(-type) # remove if you don't hit at least 10 tpm

counts.dr <- left_join(counts.dr, counts.dr.max, by = "gene_id") %>%
  filter(!is.na(sample.max))

counts.dr <- counts.dr %>%
  mutate(enr.loc = ifelse(sample.max >= 160/20 - 1 & sample.max <= 260/20 + 1, "nring",
                      ifelse(sample.max >= 380/20 - 1 & sample.max <= 480/20 + 1, "secp",
                            ifelse(sample.max >= 556/20 -1 & sample.max <= 720/20 + 1, "vulva","diffuse"))))


# order and factor for plotting (GPCRs: 10, TRPs: 10, LGICs: )
counts.dr$type <- factor(counts.dr$type, levels = c("GPCR","LGIC","TRP")) 

gpcr.list <- filter(counts.dr, type == "GPCR")
gpcr.list <- sort(unique(gpcr.list$gene_name))
trp.list <- filter(counts.dr, type == "TRP")
trp.list <- sort(unique(trp.list$gene_name))
lgic.list <- filter(counts.dr, type == "LGIC")
lgic.list <- sort(unique(lgic.list$gene_name))
receptor.list <- c(gpcr.list,trp.list,lgic.list)

counts.dr$gene_name <- factor(counts.dr$gene_name, levels = c(receptor.list))

# gpcrs
gpcr.plot <- ggplot(data=filter(counts.dr,type %in% c("GPCR")), aes(sample,expression, group = gene_name)) +
  theme_zlab_white() +
  annotate("rect", xmin=160/20, xmax=260/20, ymin=-Inf, ymax=Inf, fill = "#547DBF", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=380/20, xmax=480/20, ymin=-Inf, ymax=Inf, fill = "#3FFFBB", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=556/20, xmax=720/20, ymin=-Inf, ymax=Inf, fill = "#E574FF", color="white", size = 0, alpha = 0.3) +
  geom_line(aes(colour = enr.loc), size = 0.5, alpha = 0.6) +
  scale_color_manual("", values=c(nring="#0000FF",secp="#006837",vulva="#662D91",diffuse="gray54")) +
  scale_y_discrete(breaks = NULL) +
  scale_x_discrete(position = "top") +
  geom_text(data=filter(counts.dr,type %in% c("GPCR") & sample == "1"), aes(x = 25, y = 0.8*expression.max, label = gene_name, fontface = "italic"), 
            hjust = 0, size = 3, color = "black") +
  ylab('')+ xlab('GPCRs') + 
  theme(legend.position = "none",
        axis.title.x = element_text(size =10, face = "bold"),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = 8),
        strip.text.x = element_blank()) +
  facet_wrap(gene_name ~ ., ncol = 5, strip.position = "top", scales ="free_y") 
#gpcr.plot

# trps
trp.plot <- ggplot(data=filter(counts.dr,type %in% c("TRP")), aes(sample,expression, group = gene_name)) +
  theme_zlab_white() +
  annotate("rect", xmin=160/20, xmax=260/20, ymin=-Inf, ymax=Inf, fill = "#547DBF", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=380/20, xmax=480/20, ymin=-Inf, ymax=Inf, fill = "#3FFFBB", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=556/20, xmax=720/20, ymin=-Inf, ymax=Inf, fill = "#E574FF", color="white", size = 0, alpha = 0.3) +
  geom_line(aes(colour = enr.loc), size = 0.5, alpha = 0.6) +
  scale_color_manual("", values=c(nring="#0000FF",secp="#006837",vulva="#662D91",diffuse="gray54")) +
  scale_y_discrete(breaks = NULL) +
  scale_x_discrete(position = "top") +
  geom_text(data=filter(counts.dr,type %in% c("TRP") & sample == "1"), aes(x = 25, y = 0.8*expression.max, label = gene_name, fontface = "italic"), 
            hjust = 0, size = 3, color = "black") +
  ylab('')+ xlab('TRPs') +
  theme(legend.position = "none",
        axis.title.x = element_text(size =10, face = "bold"),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = 8),
        strip.text.x = element_blank()) +
  facet_wrap(gene_name ~ ., ncol = 5, strip.position = "top", scales ="free_y") 
#trp.plot

# lgics
lgic.plot <- ggplot(data=filter(counts.dr,type %in% c("LGIC")), aes(sample,expression, group = gene_name)) +
  theme_zlab_white() +
  annotate("rect", xmin=160/20, xmax=260/20, ymin=-Inf, ymax=Inf, fill = "#547DBF", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=380/20, xmax=480/20, ymin=-Inf, ymax=Inf, fill = "#3FFFBB", color="white", size = 0, alpha = 0.3) +
  annotate("rect", xmin=556/20, xmax=720/20, ymin=-Inf, ymax=Inf, fill = "#E574FF", color="white", size = 0, alpha = 0.3) +
  geom_line(aes(colour = enr.loc), size = 0.5, alpha = 0.6) +
  scale_color_manual("", values=c(nring="#0000FF",secp="#006837",vulva="#662D91",diffuse="gray54")) +
  scale_y_discrete(breaks = NULL) +
  scale_x_discrete(position = "top") +
  geom_text(data=filter(counts.dr,type %in% c("LGIC") & sample == "1"), aes(x = 25, y = 0.8*expression.max, label = gene_name, fontface = "italic"), 
            hjust = 0, size = 3, color = "black") +
  ylab('')+ xlab('LGICs') +
  theme(legend.position = "none",
        axis.title.x = element_text(size =10, face = "bold"),
        axis.title.y = element_text(size = 9),
        axis.text.y = element_text(size = 6),
        axis.ticks.x = element_blank(),
        axis.text.x = element_blank(),
        #axis.text.x = element_text(size = 8),
        strip.text.x = element_blank()) +
  facet_wrap(gene_name ~ ., ncol = 5, strip.position = "top", scales ="free_y") 
#lgic.plot

#strip.text.y = element_text(size = 10, face = "italic")) + 
#facet_wrap(gene_name ~ ., ncol = 5, labeller = labeller(gene_name), strip.position = "none", scales ="free_y") 




#########################
##### Alternatively, use scale-normalized data for the drug target plots (similar results)
#########################

# read in gene counts, filter for robust genes, prep for matrix conversion
counts <- readRDS("data/counts.tpm.rds") %>% arrange(sample) %>%
  dplyr::filter(rep == "R1") %>%
  dplyr::mutate(rep_sample = paste0(rep,"_",sample)) %>%
  dplyr::filter(rep_sample %in% slice_list_r123) %>%
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

#convert count matrix to df, re-order based on clustering, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")
counts <- counts %>%
  pivot_longer(2:ncol(counts), names_to = "rep_sample", values_to = "expression") %>%
  separate(rep_sample, c("rep","sample"))
counts$sample <- as.integer(counts$sample)

#update cryosection indexing based on qc (R1 = 12, R2 = 9, R3 = 8)
counts <- counts %>%
  mutate(sample = ifelse(rep == "R1", sample - 11, ifelse(rep =="R2", sample - 8, ifelse(rep == "R3", sample - 7, sample))))









########################
## Fig5.phylo / LGIC correlations 
########################

### CORRELATION MATRIX
#read lgics in from tree
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

# write and read LGIC correlations
saveRDS(cor.r, file = "data/lgic.cor.rds")

# don't run above, just load the RDS
setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure4")
cor.r <- readRDS(file = "data/lgic.cor.rds")

### TREE
# reinstalling all of the packages from the source code
# remotes::install_github("YuLab-SMU/tidytree", force = TRUE)
# remotes::install_github("YuLab-SMU/ggtree", force = TRUE)
# devtools::install_github("GuangchuangYu/treeio")
# remotes::install_github("tidyverse/dplyr@v1.0.5", force = TRUE)
library(ggtree)
library(tidytree)
library(treeio)
library(ape)
library(dplyr)
library(stringr)

# load in Bma and Cel ids
Bma.id <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bma.Proteins.csv",
                      header = FALSE, sep = ",") %>% 
  rename("protein_id" = V1, "gene_id" = V2, "gene_name" = V3) %>%
  group_by(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)
Bma.protein.list <- unique(Bma.id$protein_id)

Cel.id <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Cel.Proteins.csv",
                   header = FALSE, sep = ",") %>% 
  rename("protein_id" = V1, "gene_id" = V2, "gene_name" = V3) %>%
  group_by(gene_id,gene_name) %>%
  distinct(gene_id, .keep_all=TRUE)
Cel.gene.list <- unique(Cel.id$gene_id)
Cel.protein.list <- unique(Cel.id$protein_id)

all_id <- bind_rows(Bma.id, Cel.id) %>% 
  mutate(label = case_when(
    str_detect(protein_id, 'Bm') == TRUE ~ protein_id,
    TRUE ~ gene_id)
  )

# read in iqtree file
lgic.phylo <- read.iqtree("data/LGIC_trim_final.aln.treefile") 

# convert tree (phylo object) to tibble (relabel) and generate d1 with other data
d1 <- as_tibble(lgic.phylo) %>% 
  mutate(species = case_when(
    label %in% Bma.protein.list ~ 'Bma',
    label %in% Cel.gene.list ~ 'Cel'
  )) %>% 
  left_join(., all_id) %>% 
  mutate(tiplab = case_when(
           is.na(gene_name) == TRUE ~ protein_id,
           TRUE ~ gene_name
         ))


#reroot
lgic.phylo <- ape::root(lgic.phylo, node = 214)

(node_labels <- ggtree(lgic.phylo, layout = "circular", branch.length = "none") %<+% d1 +
    geom_text2(aes(subset = !isTip, label = node), size = 2, hjust = -.3) +
    geom_text2(aes(subset = !isTip, label = UFboot, color = UFboot), size = 2, hjust = -.3, vjust = 3) +
    geom_tiplab(aes(label = tiplab)) +
    theme_tree2() +
    NULL)

# prepare correlation data
cor_tree <- cor.r %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('from' = 'gene_id')) %>% 
  rename(from_gene_id = tiplab) %>% 
  left_join(., select(d1, gene_id, tiplab), by = c('to' = 'gene_id')) %>%
  rename(to_gene_id = tiplab) %>% 
  select(from = from_gene_id, to = to_gene_id, corr)

# subset tree (glc)
glc.phylo <- tree_subset(lgic.phylo, node = 161, levels_back = 0)  #163 originally
glc.tibble <- as_tibble(glc.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:10) %>%
  mutate(newlab = str_remove(tiplab, "Bma-"))

glc_cor <- cor_tree %>% 
  filter(from %in% glc.tibble$tiplab & to %in% glc.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>%
  mutate(color = cut(corr, 9, labels = LETTERS[1:9]))

glc.phylo@phylo$tip.label <- glc.tibble$tiplab
glc.phylo@data$SH_aLRT <- 0.1

(glc.tree <- ggtree(glc.phylo, layout = 'inward_circular', branch.length = 'SH_aLRT', xlim = 4) %<+% glc.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = glc_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8) +
    geom_tiplab(aes(label = newlab),
                size = 3,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = 0,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty") +
    NULL)

save_plot('glc_tree.pdf', glc.tree, base_height = 12)

# subset tree (acc)
acc.phylo <- tree_subset(lgic.phylo, node = 191, levels_back = 0) 
acc.tibble <- as_tibble(acc.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:14) %>%
  mutate(newlab = str_remove(tiplab, "Bma-"))

acc_cor <- cor_tree %>% 
  filter(from %in% acc.tibble$tiplab & to %in% acc.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>% 
  mutate(color = cut(corr, 10, labels = LETTERS[1:10]))

acc.phylo@phylo$tip.label <- acc.tibble$tiplab

(acc.tree <- ggtree(acc.phylo, layout = 'inward_circular', xlim = 4) %<+% acc.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = acc_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8) +
    geom_tiplab(aes(label = newlab),
                size = 3,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = -0.,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty") +
    NULL)

save_plot('acc_tree.pdf', plot_grid(acc.tree), base_height = 12)

# subset tree (nachr)

nachr.phylo <- tree_subset(lgic.phylo, node = 241, levels_back = 0) 
nachr.tibble <- as_tibble(nachr.phylo) %>% 
  left_join(., d1, by = 'label') %>% 
  select(!contains('.y')) %>% 
  rename(node = node.x) %>% 
  dplyr::slice(1:30) %>%
  mutate(newlab = str_remove(tiplab, "Bma-"))

nachr_cor <- cor_tree %>% 
  filter(from %in% nachr.tibble$tiplab & to %in% nachr.tibble$tiplab) %>%
  filter(from != to) %>% 
  filter(corr > 0) %>% 
  # remove isoform links
  separate(from, into = c('from', 'from_isoform'), sep = '\\.') %>% 
  separate(to, into = c('to', 'to_isoform'), sep = '\\.') %>% 
  filter(from != to) %>% 
  mutate(to_isoform = case_when(
    is.na(to_isoform) == TRUE ~ '',
    TRUE ~ to_isoform
  )) %>% 
  mutate(from_isoform = case_when(
    is.na(from_isoform) == TRUE ~ '',
    TRUE ~ from_isoform
  )) %>% 
  mutate(
    from = str_c(from, from_isoform, sep = '.'),
    to = str_c(to, to_isoform, sep = '.')
  ) %>% 
  mutate(
    from = str_remove(from, '\\.$'), 
    to = str_remove(to, '\\.$')
  ) %>% 
  select(from, to, corr) %>% 
  mutate(color = cut(corr, 10, labels = LETTERS[1:10])) 

nachr.phylo@phylo$tip.label <- nachr.tibble$tiplab

(nachr.tree <- ggtree(nachr.phylo, layout = 'inward_circular', xlim = 5) %<+% nachr.tibble +
    geom_tippoint(aes(fill = species),
                  shape = 21,
                  size = 2) +
    geom_taxalink(data = nachr_cor, mapping = aes(taxa1 = from, taxa2 = to,color = color),
                  size = 1,
                  ncp = 2,
                  offset = 1.1,
                  outward = FALSE,
                  alpha = 0.8
    ) +
    geom_tiplab(aes(label = newlab),
                size = 2.5,
                align = TRUE,
                linesize = 0,
                linetype = 0,
                offset = -1,
                hjust = 0,
                fontface = 'italic') +
    scale_color_brewer(palette = 'Reds') +
    scale_fill_manual(values = c('black', 'white')) +
    theme(legend.position = "empty") +
    NULL)

save_plot('nachr_tree.pdf', plot_grid(nachr.tree), base_height = 3)

fig5c <- plot_grid(acc.tree, glc.tree, nachr.tree, nrow = 3)

ggsave('Fig5C.pdf', fig5c, height = 9.5, width = 3)


# supp fig 2
(all_corr <- cor_tree %>% 
    filter(corr > 0) %>% 
    mutate(class = case_when(
      from %in% acc.tibble$gene_name & to %in% acc.tibble$gene_name ~ 'ACC',
      from %in% glc.tibble$gene_name & to %in% glc.tibble$gene_name ~ 'GLC',
      from %in% nachr.tibble$gene_name & to %in% nachr.tibble$gene_name ~ 'N-AchR',
    )) %>% 
    # filter(!is.na(class)) %>% 
    ggplot() +
    geom_tile(aes(x = from, y = to, fill = corr)) +
    scale_fill_viridis() +
    labs(fill = 'Correlation') +
    # facet_wrap(facets = vars(class)) +
    theme(
      axis.title = element_blank(),
      axis.text = element_text(face = 'italic'),
      axis.text.x = element_text(angle = 45, hjust = 1),
      panel.background = element_rect(fill = 'black'),
      panel.grid = element_blank(),
      legend.position = 'bottom'
    ) +
    NULL)

save_plot('SupplementaryFigure2_correlation.pdf', all_corr, base_width = 9.5, base_height = 9.5)


#########################
##### Figure 5 export
#########################

Fig5a.1 <- tile.plot 
Fig5a.2 <- Fig5.ill 
Fig5a <- plot_grid(Fig5a.1 + theme(plot.margin = unit(c(0.35,3.45,-0.1,-0.35), "cm")), 
                   Fig5a.2 + theme(plot.margin = unit(c(0,0,0,0.1), "cm")), ncol = 1, labels = c('',''), rel_heights = c(1,1), scale = 1)

Fig5b <- plot_grid(gpcr.plot + theme(plot.margin = unit(c(0.3,0,0,0), "cm")), 
                   NA,
                   trp.plot + theme(plot.margin = unit(c(0.3,0,0,0), "cm")),
                   NA,
                   lgic.plot + theme(plot.margin = unit(c(0.3,0,0,0), "cm")), 
                   ncol = 1, labels = c('','',''), rel_heights = c(1,0.02,1,0.02,3.5), scale = 1) 

Fig5c <- plot_grid(acc.tree + theme(plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm")), 
                   glc.tree + theme(plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm")), 
                   nachr.tree + theme(plot.margin = unit(c(-.5,-.5,-.5,-.5), "cm")), 
                   labels = c(' ACC',' GluCl',' nAChR'), ncol = 1, scale = 1, 
                   label_size = 10.5, label_fontface = "plain", label_x = 0.5)

Fig5L <- plot_grid(Fig5a, NA, Fig5b, ncol= 1, labels = c('A','','B'), rel_heights = c(1.5,0.02,3), scale = 0.99)
Fig5R <- plot_grid(Fig5c, ncol= 1, labels = c('C'), scale = 1)

Fig5 <- plot_grid(Fig5L,NA,Fig5R, nrow = 1, labels = c('','',''), rel_widths = c(1, 0, 0.4), scale = 0.99)

ggsave("Fig5.pdf", Fig5 , width = 13, height = 10.5, units = "in")






