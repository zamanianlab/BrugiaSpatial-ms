library(tidyverse)
library(cowplot)
library(ggrepel)
library("ZamanianLabThemes")
library(paletteer)
library(DESeq2)
library(hexbin)
library(wesanderson)

setwd("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure6")
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


############
## Fig6.ill: illustration of lcm / intestine dissections
############

library(magick)
library(pdftools)
library(grImport2)

RNA.illustration <- image_read_pdf("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/LCM_Pharynx/Illustration_IntestinePharynx.pdf", density = 300)
Fig6.ill <- ggdraw() + 
  draw_image(RNA.illustration, scale = 1) 


##########################################
### Intestine / non-Intestine analysis (DESeq/Volcano Plot)
##########################################

# Bma AF Intestine/non-Intestine (S07-S12) sample list 
# (S07/S08 = intestine 1 / non_intestine 1 ...etc)
sample_filter <- c('S07','S08','S09','S10','S11','S12')

# read in counts (raw) and filter for samples
counts_a.raw <- readRDS("data/counts-ic-lcm.raw.rds") %>%
  filter(sample %in% sample_filter)
counts_b.raw <- readRDS("data/counts-ic-lcm_b.raw.rds") %>%
  filter(sample %in% sample_filter) 

counts.raw <- left_join(counts_a.raw,counts_b.raw, by = c("gene_id","sample"))
counts.raw[is.na(counts.raw)] = 0
counts.raw <- counts.raw %>% 
  mutate(expression = expression.x + expression.y) %>%
  select(gene_id,sample,expression)

# convert to matrix
counts.raw <- counts.raw %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# make metadata file
sample_list <- colnames(counts.raw)
samples <- data.frame("sample_id" = sample_list, stringsAsFactors = FALSE)
samples <- samples %>%
  mutate(cond = rep(c("Intestine", "N_Intestine"), times = 3)) %>%
  mutate(worm = c("1","1","2","2","3","3"))
samples$cond <- as.factor(samples$cond)
samples$worm <- as.factor(samples$worm)
saveRDS(samples, "deg/samples-ic.rds")

# read in sample metadata df > samples
samples <- readRDS("deg/samples-ic.rds")

##########################################
## Transform and Cluster Samples
##########################################

gene_count <- as.matrix(counts.raw)
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond + worm)

#filter 1 (require x reads total across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

#filter 2 (require x reads in each of at least y samples)
keep <- rowSums(counts(dds) >= 20) >= 3
dds <- dds[keep,]
nrow(dds) #11718 > 6075

#pick vst transform
vsd <- vst(dds, blind = FALSE)

#calculate sample distances and use hclust to cluster samples based on dist
sampleDists <- dist(t(assay(vsd)), method = "euclidean") #chose vsd over rld
clust.dist <- hclust(sampleDists, method="ward.D2")

#get list order of clustered samples from hclust output
ord <- clust.dist$order

#convert original distance matrix to df, re-order, and tidy for plotting
sampleDists.df <- as.data.frame(as.matrix(sampleDists)) %>%
  rownames_to_column(var = "sample") 
sampleDists.df$sample <- factor(sampleDists.df$sample, 
                                   levels = c(sampleDists.df$sample[ord]),
                                   labels = c("Non-Intestine 3", "Non_Intestine 1", "Non_Intestine 2", "Intestine 2", "Intestine 1", "Intestine 3"))
sampleDists.df <- sampleDists.df %>%
  pivot_longer(cols = 2:ncol(sampleDists.df),
               names_to = "sample_2",
               values_to = "dist")
sampleDists.df$sample_2 <- factor(sampleDists.df$sample_2, 
                                     levels = c(sampleDists.df$sample_2[ord]), 
                                     labels = c("Non-Intestine 3", "Non_Intestine 1", "Non_Intestine 2", "Intestine 2", "Intestine 1", "Intestine 3"))

#plot heatmap cluster
pal <- wes_palette("Zissou1", 25, type = "continuous")
cluster.ic <- ggplot(sampleDists.df, aes(x = sample, y = sample_2, fill = dist)) +
  theme_zlab() + xlab('') + ylab('') + labs(fill = "") +
  geom_tile() +
  scale_fill_gradientn(colours = pal)
cluster.ic

##########################################
## PCA analysis
##########################################

#Run PCA on vst transformed data and return output for ggplot
pca <- plotPCA(vsd, intgroup = c("cond"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

pca$cond <- factor(pca$cond, levels = c("Intestine","N_Intestine"), labels = c("Intestine", "Non-intestine"))
#plot pca
pca.ic <- ggplot(pca, aes(x = PC1, y = PC2, shape = cond)) +
  geom_point(size =2.5, alpha = 0.75) +
  theme_zlab_white() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank())+
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  #scale_color_manual(values = c("#904494","#dfc765")) +
  theme(legend.position= "top",
        legend.direction = "horizontal",
        legend.title = element_blank(),
        legend.margin=margin(-10,0,0,0),
        legend.justification="center")
pca.ic

##########################################
# Identify DEGs
##########################################

gene_count <- as.matrix(counts.raw)

#Read in sample metadata df > samples
samples <- readRDS("deg/samples-ic.rds")

# run deseq to compare baselines (ctl) across strains
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond + worm)
keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep,]
dds <- DESeq(dds) #resultsNames(dds)
res <- results(dds, contrast = c("cond", "Intestine", "N_Intestine")) #summary(res)
res.df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>% 
  arrange(padj) %>%
  mutate(sig = ifelse((abs(log2FoldChange) > 1) & padj < 0.01, "yes", "no"))
saveRDS(res.df, "deg/ic.res.rds")

#Supplementary file 
write.csv(res.df,"IC_DEGs.csv", row.names = FALSE)

#plotCounts(dds, gene="WBGene00220555", intgroup="cond")

# get list of DEGs
degs <- res.df %>% filter(sig == "yes") #2759
degs.int <- degs %>% filter(log2FoldChange >= 1) #1077
degs.nint <- degs %>% filter(log2FoldChange <= 1) #1682
degs.int.list <- unique(degs.int$gene_id) 


##########################################
# Annotated Volcano Plot
##########################################

# Prep TM(+) gene list
Bma.TM <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bma.TM.csv",
                      header = FALSE, sep = ",") %>%
  dplyr::rename("protein_id" = V1, "gene_id" = V2, "gene_name" = V3, "TM_count" = V4) %>%
  distinct(gene_id, .keep_all=TRUE)

# Join TM data and add columns for plotting
res.df <- left_join(res.df,Bma.TM, by="gene_id")
res.df <- res.df %>%
  dplyr::mutate(gene_name = str_replace(gene_name, "Bma-", "")) %>%
  dplyr::mutate(TM = ifelse(TM_count > 0,"yes","no"))

# Join intestinal proteom data
prot.int.t20 <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bm.int.proteome.t20.txt",
                   header = FALSE, sep = ",")
prot.int.t20 <- prot.int.t20$V1
res.df <- res.df %>% mutate(prot.int.t20 = ifelse(gene_id %in% prot.int.t20, "yes", "no"))

prot.int.only <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bm.int.proteome.only.txt",
                         header = FALSE, sep = ",")
prot.int.only <- prot.int.only$V1
res.df <- res.df %>% mutate(prot.int.only = ifelse(gene_id %in% prot.int.only, "yes", "no"))


# Add in TPM counts
sample_filter <- c('S07','S08','S09','S10','S11','S12')
counts_a.tpm <- readRDS("data/counts-ic-lcm.tpm.rds") %>%
  filter(sample %in% sample_filter)
counts_b.tpm <- readRDS("data/counts-ic-lcm_b.tpm.rds") %>%
  filter(sample %in% sample_filter) 
counts.tpm <- left_join(counts_a.tpm,counts_b.tpm, by = c("gene_id","sample"))
counts.tpm[is.na(counts.tpm)] = 0
counts.tpm <- counts.tpm %>% 
  mutate(expression = expression.x + expression.y) %>%
  select(gene_id,sample,expression) %>%
  pivot_wider(names_from = sample, values_from = expression)
res.df.tpm <- left_join(res.df, counts.tpm, by = "gene_id") %>%
  mutate(tm.int = ifelse(TM == "yes" & sig == "yes" & log2FoldChange >= 1, "yes","no")) %>%
  mutate(mean.int.tpm = ((S07+S09+S11)/3),mean.nint.tpm = ((S08+S10+S12)/3)) %>%
  mutate(hid.ant = ifelse(sig == "yes" & TM == "yes" & mean.int.tpm > 100 & (mean.int.tpm/mean.nint.tpm) >10, "yes", "no"))
#merge gene_name/protein_id
res.df.tpm <- res.df.tpm %>%
  mutate(gene_name = ifelse(is.na(gene_name),protein_id,gene_name))

# export Supplementary Table
write.csv(res.df.tpm,"IC_DEGs_TPM.csv", row.names = FALSE)

# Volcano plot
volcano.ic <- ggplot(data=res.df.tpm, aes(x=log2FoldChange, y=-log10(padj))) +
  theme_zlab_white() +
  geom_vline(xintercept=c(0), col="gray") +
  geom_vline(xintercept=c(-1, 1), linetype='dashed', col="black", size = 0.25) +
  geom_hline(yintercept=2, linetype='dashed', col="black", size = 0.25) +
  geom_point(data=subset(res.df.tpm, tm.int == "yes"), shape=16, colour ="#D1A86F", alpha=0.9, size=1.75) + 
  geom_point(data=subset(res.df.tpm, tm.int == "no"), shape=16, colour ="black", alpha=0.25, size=1) + 
  geom_point(data=subset(res.df.tpm, hid.ant == "yes"), shape=16, colour ="#6600CC", alpha=1, size=2.25) + 
  #geom_text_repel(data=subset(res.df.tpm, hid.ant == "yes"),aes(label=gene_name, fontface=3),
  #                 segment.alpha=0.4,segment.size=0.4, segment.color = "#6600CC",
  #                 nudge_y=0,nudge_x=0, colour ="black", max.iter = 100000, max.overlaps = 20 , force = 3, size=3) +
  annotate(geom="text", x=4, y= 120, label="Intestine-enriched",color="black", size = 4, hjust = 0) + 
  ylim(0,120) + xlim(-10,10) +
  xlab(expression(log["2"]~'(Fold Change)')) + ylab(expression(-log["10"]~italic(P)))
volcano.ic

#Stats
# all degs 
int.degs <- res.df.tpm %>% filter(sig == "yes")  #2759
int.degs <- int.degs$gene_id
# number of intestinal-enriched TM proteins
tm.int <- res.df.tpm %>% filter(tm.int == "yes") #489
tm.int <- tm.int$gene_id
# number of prioritized intestinal-enriched
hid.ant <- res.df.tpm %>% filter(hid.ant == "yes") #64
hid.ant <- hid.ant$gene_id
#how many of t20 proteome hits (prot.int.t20) were re-identified
all.int <- res.df.tpm %>% filter(sig == "yes" & log2FoldChange >1)
all.int <- all.int$gene_id #1077
t20.capture <- intersect(all.int, prot.int.t20) #7 of 17
#how many of intestine-only proteome hits (prot.int.only) were re-identified
only.capture <- intersect(all.int, prot.int.only) #1 of 4




##########################################
### Heatmap of prioritized hidden vaccine candidates with overlay of high in pharynx
##########################################

# read in counts (vsd normalized) 
vsd_degs <- as.data.frame(assay(vsd)) 
vsd_degs <- vsd_degs %>%
  rownames_to_column(var = "gene_id") %>%
  #filter(gene_id %in% hid.ant) %>%  # option: filter for prioritized genes (hidden antigens)
  #filter(gene_id %in% int.degs) %>% # option: filter for all DEGs
  column_to_rownames(var = "gene_id")
counts <- vsd_degs
colnames(counts) <- c("Intestine_1", "NIntestine_1", "Intestine_2", "NIntestine_2", "Intestine_3", "NIntestine_3")

# read in sample metadata df > samples
samples <- readRDS("deg/samples-ic.rds")

#join counts and sample info, widen, declare and normalize matrix
df.m <- data.matrix(counts, rownames.force = TRUE)
ind <- apply(df.m, 1, var) == 0  #remove genes with no variance 
df.m <- df.m[!ind,]
df.m <- t(scale(t(df.m),center=TRUE,scale=TRUE))
counts <- df.m

#calculate gene distances and use hclust to cluster samples based on dist
geneDists <- dist(counts, method = "euclidean")
gclust.dist <- hclust(geneDists, method="ward.D2")

#get list order of clustered genes from hclust output
ord <- gclust.dist$order
ord_labels <- gclust.dist$labels

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")

#add annotations, factor based on clustering order
counts <- left_join(counts,Bma.TM, by="gene_id")
counts <- counts %>% 
  mutate(gene_name = ifelse(is.na(gene_name),protein_id,gene_name)) %>%
  mutate(gene_name = str_replace(gene_name, "Bma-", ""))
counts$gene_id <- factor(counts$gene_id, levels = counts$gene_id[ord], labels = counts$gene_name[ord])

#long form for plotting
counts <- counts %>%
  pivot_longer(2:7, names_to = "sample", values_to = "expression") %>%
  tidyr::separate(sample, into = "group", sep = "_", remove = FALSE)
counts$sample <- factor(counts$sample, levels = c("Intestine_1", "Intestine_2", "Intestine_3","NIntestine_1", "NIntestine_2", "NIntestine_3"), 
                        labels = c("Intestine 1", "Intestine 2", "Intestine 3", "Non-intestine 1", "Non-intestine 2", "Non-intestine 3"))
counts$group <- factor(counts$group, levels = c("Intestine", "NIntestine"), 
                       labels = c("Intestine", "Non-intestine"))

# number of prioritized intestinal-enriched
hid.ant <- res.df.tpm %>% filter(hid.ant == "yes") #64
hid.ant <- hid.ant$gene_name

# plot heatmap
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
heatmap <- ggplot(filter(counts, gene_id %in% hid.ant), aes(gene_id, sample)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colours = pal, "z-score",
                       guide = guide_colorbar(
                         direction = "horizontal",
                         title.position = "left",
                         label.position = "bottom"
                       )) +
  #scale_x_discrete(guide = guide_axis(n.dodge = 1)) +
  scale_y_discrete(labels = c("Intestine 1" = "","Intestine 2" = "Intestine","Intestine 3" = "",
                              "Non-intestine 1" = "","Non-intestine 2" = "Non-intestine", "Non-intestine 3" = ""),
                   position = "left") +
  geom_hline(yintercept = 3.5, linetype="dashed", color = "black", size=0.5) +
  xlab("") + ylab("") + 
  theme_zlab_white() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.y = element_text(size = 12, angle = 90),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x = element_text(face = "italic", angle = 45, size = 9, hjust =1, vjust = 1),
    strip.background = element_rect(fill="white", color = "white"),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing=unit(-0.75, "lines"),
    legend.position = "bottom",
    legend.key.width=unit(1,"cm"),
    legend.margin=margin(-15,0,0,0)) 
heatmap



##########################################
### Pharynx / non-Pharynx analysis (DESeq/Volcano Plot)
##########################################

# Bma AF Pharynx/non-Pharynx (S05-S06) sample list 
sample_filter <- c('S05','S06')

# read in counts (raw) and filter and combine for samples from first two runs
counts_a.raw <- readRDS("data/counts-ic-lcm.raw.rds") %>%
  filter(sample %in% sample_filter)
counts_b.raw <- readRDS("data/counts-ic-lcm_b.raw.rds") %>%
  filter(sample %in% sample_filter) 

counts.raw <- left_join(counts_a.raw,counts_b.raw, by = c("gene_id","sample"))
counts.raw[is.na(counts.raw)] = 0
counts.raw <- counts.raw %>% 
  mutate(expression = expression.x + expression.y) %>%
  select(gene_id,sample,expression)

# read in counts (raw) and filter for sample (S1-S4) from second run (4 samples)
#S1_b/S2_b = non-pharynx/pharynx ... 
counts_add.raw <- readRDS("data/lcm_b.raw.rds") %>%
  mutate(sample = paste0(sample,"_b"))

#join all the samples
counts.raw <- rbind (counts.raw, counts_add.raw)

# convert to matrix
counts.raw <- counts.raw %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# make metadata file
sample_list <- colnames(counts.raw)
samples <- data.frame("sample_id" = sample_list, stringsAsFactors = FALSE)
samples <- samples %>%
  mutate(cond = c("Pharynx", "Nonpharynx","Nonpharynx","Pharynx", "Nonpharynx", "Pharynx")) %>%
  mutate(worm = c("1","1","2","2","3","3")) 
samples$cond <- as.factor(samples$cond)
samples$worm <- as.factor(samples$worm)
saveRDS(samples, "deg/samples-lcm.rds")

# read in sample metadata df > samples
samples <- readRDS("deg/samples-lcm.rds")


##########################################
## Transform and Cluster Samples / Supplementary Figure
##########################################

gene_count <- as.matrix(counts.raw)
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond + worm)

#filter 1 (require x reads total across all samples)
keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]

#filter 2 (require x reads in each of at least y samples)
keep <- rowSums(counts(dds) >= 10) >= 3
dds <- dds[keep,]
nrow(dds) #11718 > 2450 (1649)

#pick vst transform
vsd <- vst(dds, blind = FALSE)

#calculate sample distances and use hclust to cluster samples based on dist
sampleDists <- dist(t(assay(vsd)), method = "euclidean") #chose vsd over rld
clust.dist <- hclust(sampleDists, method="ward.D2")

#get list order of clustered samples from hclust output
ord <- clust.dist$order

#convert original distance matrix to df, re-order, and tidy for plotting
sampleDists.df <- as.data.frame(as.matrix(sampleDists)) %>%
  rownames_to_column(var = "sample") 
sampleDists.df$sample <- factor(sampleDists.df$sample, 
                                levels = c(sampleDists.df$sample[ord]),
                                labels = c("Non-Pharynx 1", "Non-Pharynx 2", "Non-Pharynx 3", "Pharynx 1", "Pharynx 2", "Pharynx 3"))
sampleDists.df <- sampleDists.df %>%
  pivot_longer(cols = 2:ncol(sampleDists.df),
               names_to = "sample_2",
               values_to = "dist")
sampleDists.df$sample_2 <- factor(sampleDists.df$sample_2, 
                                  levels = c(sampleDists.df$sample_2[ord]),
                                  labels = c("Non-Pharynx 1", "Non-Pharynx 2", "Non-Pharynx 3", "Pharynx 1", "Pharynx 2", "Pharynx 3"))

#plot heatmap cluster
pal <- wes_palette("Zissou1", 25, type = "continuous")
cluster.pn <- ggplot(sampleDists.df, aes(x = sample, y = sample_2, fill = dist)) +
  theme_zlab_white() + xlab('') + ylab('') +
  geom_tile() +
  scale_y_discrete(labels = c("Non-Pharynx 1" = "","Non-Pharynx 2" = "Non-Pharynx","Non-Pharynx 3" = "",
                              "Pharynx 1" = "","Pharynx 2" = "Pharynx","Pharynx 3" = ""),
                   position = "left") +
  scale_x_discrete(labels = c("Non-Pharynx 1" = "","Non-Pharynx 2" = "Non-Pharynx","Non-Pharynx 3" = "",
                              "Pharynx 1" = "","Pharynx 2" = "Pharynx","Pharynx 3" = ""),
                   position = "bottom") +
  scale_fill_gradientn(colours = pal, "Distance",
                       guide = guide_colorbar(
                         direction = "horizontal",
                         title.position = "left",
                         label.position = "bottom")) +
  geom_vline(xintercept = 3.5, linetype="dashed", color = "black", size=0.5) +
  geom_hline(yintercept = 3.5, linetype="dashed", color = "black", size=0.5) +
  theme(axis.ticks = element_blank(),
        axis.text.y = element_text(angle=90,hjust = 0.5),
        axis.text.x = element_text(angle=0,hjust = 0.5),
        legend.key.width=unit(1,"cm"),
        legend.position = "top")
cluster.pn
ggsave("pn_cluster.pdf", cluster.pn,width = 4, height = 4.5)



##########################################
# Identify DEGs - IGNORE
##########################################

gene_count <- as.matrix(counts.raw)

#Read in sample metadata df > samples
samples <- readRDS("deg/samples-lcm.rds")

# run deseq to compare baselines (ctl) across strains
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond + worm)
keep <- rowSums(counts(dds)) >= 30 
dds <- dds[keep,]
dds <- DESeq(dds) #resultsNames(dds)
res <- results(dds, contrast = c("cond", "Pharynx", "Nonpharynx")) #summary(res)
res.df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>% 
  arrange(padj) %>%
  mutate(sig = ifelse((abs(log2FoldChange) > 1) & padj < 0.01, "yes", "no"))
saveRDS(res.df, "deg/lcm.res.rds")


##########################################
# TPM prioritization and Venn Diagram
##########################################

# Bma AF Pharynx/non-Pharynx (S05-S06) sample list 
sample_filter <- c('S05','S06')

# read in counts (raw) and filter and combine for samples from first two runs
counts_a.tpm <- readRDS("data/counts-ic-lcm.tpm.rds") %>%
  filter(sample %in% sample_filter)
counts_b.tpm <- readRDS("data/counts-ic-lcm_b.tpm.rds") %>%
  filter(sample %in% sample_filter) 

counts.tpm <- left_join(counts_a.tpm,counts_b.tpm, by = c("gene_id","sample"))
counts.tpm[is.na(counts.tpm)] = 0
counts.tpm <- counts.tpm %>% 
  mutate(expression = expression.x + expression.y) %>%
  select(gene_id,sample,expression)

# read in counts (raw) and filter for sample (S1-S4) from second run (4 samples)
#S1_b/S2_b = non-pharynx/pharynx ... 
counts_add.tpm <- readRDS("data/lcm_b.tpm.rds") %>%
  mutate(sample = paste0(sample,"_b"))

#join all the samples
counts.tpm <- rbind (counts.tpm, counts_add.tpm)

#populate with gene id info and tpm totals/thresholds
counts.tpm <- left_join(counts.tpm, Bma.TM, by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name),protein_id,gene_name))
counts.tpm <- counts.tpm %>%
  pivot_wider(names_from = sample, values_from = expression)
colnames(counts.tpm) <- c("gene_id","protein_id","gene_name","TM","Pharynx_1","NPharynx_1","NPharynx_2", "Pharynx_2","NPharynx_3","Pharynx_3")

counts.tpm <- counts.tpm %>%
  mutate(mean.pharynx.tpm = ((Pharynx_1 + Pharynx_2 + Pharynx_3)/3),
         mean.npharynx.tpm = ((NPharynx_1 + NPharynx_2 + NPharynx_3)/3)) %>%
  mutate(high.pharynx = ifelse(TM > 0 & Pharynx_1 > 100 & Pharynx_2 > 100 & Pharynx_3 > 100, "yes", "no")) %>%
  mutate(ph.enrch = (mean.pharynx.tpm+1)/(mean.npharynx.tpm+1)) %>%
  mutate(prioritize = ifelse(high.pharynx == "yes" & ph.enrch > 5, "yes","no"))

#Supplementary file 
write.csv(counts.tpm,"PN_expression.csv", row.names = FALSE)


# Venn diagram 
library("ggVennDiagram")

#list of top-50 pharynx transmembrane genes
ph.t50 <- counts.tpm %>% 
  filter(high.pharynx == "yes") %>% 
  dplyr::top_n(50, mean.pharynx.tpm) 
ph.t50 <- ph.t50$gene_id

#list of prioritized pharynx transmembrane genes
ph.ha <- counts.tpm %>% 
  filter(prioritize == "yes")
ph.ha <- ph.ha$gene_id

# list of prioritized intestinal-enriched
hid.ant <- res.df.tpm %>% filter(hid.ant == "yes") #64
hid.ant <- hid.ant$gene_id

# overlap
intersect(hid.ant,ph.ha) #"WBGene00226064" hpo-8 "WBGene00227683" pcp-1

# venn plot
x <- list(R1 = ph.t50, R2 = ph.ha, R3 = hid.ant)

venn.plot <- ggVennDiagram(x[1:3], label_alpha = 0, set_size = 3.5, label_size=3.5,  category.names = c("Pharynx High","                  Pharynx Restricted","Intestinal HA                  ")) +
  ggplot2::scale_fill_gradient(low="white",high = "orange") +
  ggplot2::scale_colour_manual(values = c("#88A80D", "#07610D" ,"#D1A86F")) +
  theme(legend.position = "none") +
  geom_segment(aes(x = 1.5, y = 2, xend = -3.5, yend = 6)) +
  ggplot2::annotate(geom="text", x=-6, y= 7, label=c("italic(Bma-hpo-8)"), parse = TRUE, color="black", size = 3.25, hjust = 0) +
  ggplot2::annotate(geom="text", x=-6, y= 6.5, label=c("italic(Bma-pcp-1)"), parse = TRUE, color="black", size = 3.25, hjust = 0) 
venn.plot
Fig6.venn <- venn.plot





##########################################
### Figure 6
##########################################
library(Cairo)
theme_set(theme_cowplot(font_family = "Helvetica") + 
            theme(text = element_text(colour = "black")))

Fig6a <- Fig6.ill
Fig6b <- volcano.ic + annotation_custom(ggplotGrob(pca.ic), xmin = -11, xmax = -1, ymin = 70, ymax = 127)
Fig6T <- plot_grid(Fig6a, NULL, Fig6b, labels = c('A','','B'), nrow = 1, rel_widths = c(1,0.01,1.35), scale = 1) 

Fig6c <- heatmap
Fig6d <- Fig6.venn
Fig6B <- plot_grid(Fig6c, NULL, Fig6d, labels = c('C','','D'), nrow = 1, rel_widths = c(1,0.01, 0.45), scale = 1)

Fig6 <- plot_grid(Fig6T, NULL, Fig6B, labels = c('','',''), ncol = 1, rel_heights = c(1,0.01,0.65), scale = 0.99) 


ggsave("Fig6.pdf", Fig6, width = 14, height = 10, units = "in")







