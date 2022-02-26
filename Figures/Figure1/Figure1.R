library(tidyverse)
library(cowplot)
library(ggrepel)
library("ZamanianLabThemes")
library(paletteer)
library(DESeq2)
library(hexbin)
library(wesanderson)

setwd("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure1")

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

############
## FIG 1A Illustration of Head/Body dissection
############

library(magick)
library(pdftools)
#devtools::install_github("sjp/grConvert")
library(grConvert)
library(grImport2)

illust.hb <- image_read_pdf("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/HeadBody/Illustration_HeadBody.pdf", density = 600)
Fig1a <- ggdraw() + 
  draw_image(illust.hb, scale = 1) 


############
## FIG 1B/C Head/Body RNAseq analysis (B = volcano + pca; C = annotated heat map)
############

########################################################
####### READ IN AND SAVE RAW/TPM counts (Head vs Tail)
########################################################

RNAdata <- c("~/Library/CloudStorage/Box-Box/ZamanianLab/SeqLibraries/Mapping/")

dir <- c("O001_BmAF_HT/190114_BCCV49ANXX/")
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
saveRDS(counts.raw, "counts/ht.raw.rds")

counts.tpm$sample <- factor(counts.tpm$sample)
saveRDS(counts.tpm, "counts/ht.tpm.rds")


##########################################
### Prep raw counts for DESEQ2 analysis of AF H/T 
##########################################

# read in counts (raw) and filter for samples
counts.raw <- readRDS("counts/ht.raw.rds") 

# convert to matrix
counts.raw <- counts.raw %>%
  pivot_wider(names_from = sample, values_from = expression) %>%
  column_to_rownames(var = "gene_id")

# make metadata file
sample_list <- colnames(counts.raw)
samples <- data.frame("sample" = sample_list, stringsAsFactors = FALSE)
samples <- samples %>%
  mutate(cond = rep(c("Head", "Body"), times = 3)) %>%
  mutate(worm = c("1","1","2","2","3","3"))
samples$cond <- as.factor(samples$cond)
samples$worm <- as.factor(samples$worm)
saveRDS(samples, "deg/samples-ht.rds")

# read in sample metadata df > samples
samples <- readRDS("deg/samples-ht.rds")


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
nrow(dds) #11718 > 10098

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
                                   labels = c("Body 1", "Body 2", "Body 3", "Head 1", "Head 2", "Head 3"))
sampleDists.df <- sampleDists.df %>%
  pivot_longer(cols = 2:ncol(sampleDists.df),
               names_to = "sample_2",
               values_to = "dist")
sampleDists.df$sample_2 <- factor(sampleDists.df$sample_2, 
                                     levels = c(sampleDists.df$sample_2[ord]), 
                                     labels = c("Body 1", "Body 2", "Body 3", "Head 1", "Head 2", "Head 3"))

#plot heatmap cluster
pal <- wes_palette("Zissou1", 25, type = "continuous")
cluster.ht <- ggplot(sampleDists.df, aes(x = sample, y = sample_2, fill = dist)) +
  theme_zlab() + xlab('') + ylab('') + labs(fill = "") +
  geom_tile() +
  scale_fill_gradientn(colours = pal)
cluster.ht


##########################################
## PCA analysis
##########################################

#Run PCA on vst transformed data and return output for ggplot
pca <- plotPCA(vsd, intgroup = c("cond"), returnData = TRUE)
percentVar <- round(100 * attr(pca, "percentVar"))

#plot pca
pca.ht <- ggplot(pca, aes(x = PC1, y = PC2, shape = cond)) +
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
pca.ht


##########################################
# Identify DEGs
##########################################

gene_count <- as.matrix(counts.raw)

#Read in sample metadata df > samples
samples <- readRDS("deg/samples-ht.rds")

# run deseq to compare baselines (ctl) across strains
dds <- DESeqDataSetFromMatrix(countData = gene_count, colData = samples, design = ~ cond + worm)
keep <- rowSums(counts(dds)) >= 30
dds <- dds[keep,]
dds <- DESeq(dds) #resultsNames(dds)
res <- results(dds, contrast = c("cond", "Head", "Body")) #summary(res)
res.df <- as.data.frame(res) %>%
  rownames_to_column(var = "gene_id") %>% 
  arrange(padj) %>%
  mutate(sig = ifelse((abs(log2FoldChange) > 1) & padj < 0.01, "yes", "no"))
saveRDS(res.df, "deg/ht.res.rds")

# integrate gene names
Bma.id <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Bma.Proteins.csv",
                   header = FALSE, sep = ",") %>% 
  dplyr::rename("protein_id" = V1, "gene_id" = V2, "gene_name" = V3) %>%
  group_by(gene_id, gene_name) %>%
  distinct(gene_id, .keep_all = TRUE)
res.df <- left_join(res.df,Bma.id,by="gene_id")

#Supplementary file 
write.csv(res.df,"HeadBody_DEGs.csv", row.names = FALSE)

# get list of DEGs
degs <- res.df %>% filter(sig == "yes")
degs.head <- degs %>% filter(log2FoldChange >= 1) #2406
degs.body <- degs %>% filter(log2FoldChange <= 1) #2876
degs <- unique(degs$gene_id) 
degs.head.list <- unique(degs.head$gene_id) #use in Figure4.R


##########################################
# Annotated Volcano Plot
##########################################

res.df <- readRDS("deg/ht.res.rds")

# Prep proteome gene lists
proteome1 <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Secretome_Proteins_Bennuru_All.csv",
                      header = TRUE, sep = ",") %>%
  select(gene_id,gene_name) 

proteome2 <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Secretome_Proteins_Bennuru_AF.csv",
                      header = TRUE, sep = ",") %>%
  select(gene_id,gene_name) 

proteome3 <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Secretome_Proteins_Moreno_MF_AF_AM.csv",
                      header = TRUE, sep = ",") %>%
  select(gene_id,gene_name)

proteome <- rbind(proteome1,proteome2,proteome3) %>%
  distinct(gene_id,gene_name) %>% mutate(proteome = 1)

# Prep antigen gene list
antigens <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Antigens.csv",
                      header = TRUE, sep = ",") %>%
  dplyr::slice(1:36) %>% select(1:8) %>% mutate(antigen = 1)

# Join proteome/antigen lists
res.df <- left_join(res.df,proteome, by="gene_id")
res.df <- left_join(res.df,antigens, by="gene_id")
res.df <- res.df %>%
  mutate(gene_name = ifelse(is.na(gene_name.x),gene_name.y,gene_name.x)) %>%
  select(-`gene_name.x`,-`gene_name.y`) %>%
  mutate(gene_name = str_replace(gene_name, "Bma-", ""))

#Stats (how many genes in each quadrant)
#head-enriched proteome and antigens
he.prot <- res.df %>% filter(proteome == "1" & log2FoldChange >= 1) #141
he.antigen <- res.df %>% filter(antigen == "1" & log2FoldChange >= 1) #25

#body-enriched proteome and antigens
be.prot <- res.df %>% filter(proteome == "1" & log2FoldChange <= -1) #103
be.antigen <- res.df %>% filter(antigen == "1" & log2FoldChange <= -1) #4

# Volcano plot
volcano.ht <- ggplot(data=res.df, aes(x=log2FoldChange, y=-log10(padj))) +
  theme_zlab_white() +
  geom_vline(xintercept=c(0), col="gray") +
  geom_vline(xintercept=c(-1, 1), linetype='dashed', col="black", size = 0.25) +
  geom_hline(yintercept=2, linetype='dashed', col="black", size = 0.25) +
  geom_point(data=subset(res.df, sig == "no"), shape=16, colour = "grey", alpha = 0.025, size = 0.75) +
  geom_point(data=subset(res.df, sig == "yes"), shape=16, alpha = 0.1, size = 0.75) +
  geom_point(data=subset(res.df, proteome == "1" & sig == "yes"), shape=16, colour ="lightseagreen", alpha=0.5, size=1) + 
  geom_point(data=subset(res.df, proteome == "1" & sig == "no"), shape=16, colour ="lightseagreen", alpha=0.05, size=0.75) + 
  geom_point(data=subset(res.df, antigen == "1" & sig == "yes"), shape=16, colour ="#6600CC", alpha=0.9, size=2) + 
  geom_point(data=subset(res.df, antigen == "1" & sig == "no"), shape=16, colour ="#6600CC", alpha=0.1, size=1) + 
  geom_text_repel(data=subset(res.df, antigen == "1" & sig == "yes"),aes(label=gene_name, fontface=3),
                   segment.alpha=0.4,segment.size=0.4, segment.color = "#6600CC",
                   nudge_y=0,nudge_x=0, colour ="black", max.iter = 100000, max.overlaps = 20 , force = 3, size=3) +
  annotate(geom="text", x=7, y= 70, label="Head-enriched",color="black", size = 4, hjust = 0) + 
  annotate(geom="text", x=-10, y= 70, label="Body-enriched",color="black", size = 4, hjust = 0) + 
  #annotate(geom="text", x=-10, y= 65, label="103",color="lightseagreen", size = 4, hjust = 0) + 
  annotate(geom="text", x=8, y= 65, label="141 (of 244)",color="lightseagreen", size = 4, hjust = 0) + 
  #annotate(geom="text", x=-10, y= 62, label="4",color="#6600CC", size = 4, hjust = 0) + 
  annotate(geom="text", x=8, y= 62, label="25 (of 29)",color="#6600CC", size = 4, hjust = 0) + 
  ylim(0,70) + xlim(-11,11) +
  xlab(expression(log["2"]~'(Fold Change)')) + ylab(expression(-log["10"]~italic(P)))
volcano.ht



##########################################
### Load Antigens, Proteome, Drug Targets and join into gene_list df
##########################################

antigens <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Antigens.csv",
                     header = TRUE, sep = ",") %>%
  dplyr::slice(1:36) %>% select(gene_id,gene_name,Vaccine) %>% mutate(antigen = 1)

#current drug targets
drug_targets <- read.csv("~/Library/CloudStorage/Box-Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/Gene_Lists/Drug_Targets.csv",
                         header = TRUE, sep = ",") %>%
  filter(existing == "yes") %>%
  mutate(gene_name = ifelse(is.na(name),transcript_id,name)) %>%
  select(gene_id,gene_name,type,subtype)%>% mutate(drug_target = 1)

gene_list <- full_join(antigens,drug_targets,by = "gene_id") %>%
  mutate(gene_name = ifelse(is.na(gene_name.x),gene_name.y,gene_name.x)) %>%
  select(-`gene_name.x`,-`gene_name.y`) %>%
  mutate(gene_name = str_replace(gene_name, "Bma-", ""))

#add adjusted p-value from DE analysis above
padj <- res.df %>% select(gene_id,padj)
gene_list <- left_join(gene_list, padj, by="gene_id") %>%
  mutate(siglevel = ifelse(padj < 0.001, '***',
                            ifelse(padj < 0.01, '**',
                              ifelse(padj < 0.05, '*',''))))


##########################################
### 1B. Heatmap of genes of interest (vaccine and drug targets) - VSD
##########################################

# read in counts (vsd normalized) 
vsd_degs <- as.data.frame(assay(vsd)) 
vsd_degs <- vsd_degs %>%
  rownames_to_column(var = "gene_id") %>%
  filter(gene_id %in% unique(gene_list$gene_id)) %>%  # option: filter for gene list
  column_to_rownames(var = "gene_id")
counts <- vsd_degs
colnames(counts) <- c("Head_1", "Body_1", "Head_2", "Body_2", "Head_3", "Body_3")

# read in sample metadata df > samples
samples <- readRDS("deg/samples-ht.rds")

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

#convert genecount matrix to df, re-order, and tidy (long form + metadata) for plotting
counts <- as.data.frame(counts) %>%
  rownames_to_column(var = "gene_id")

#add annotations (switch labels to gene_name)
counts <- left_join(counts, gene_list, by = "gene_id")
counts$gene_id <- factor(counts$gene_id, levels = c(counts$gene_id[ord]), labels = c(counts$gene_name[ord]))

#long form for plotting
counts <- counts %>%
  pivot_longer(2:7, names_to = "sample", values_to = "expression") %>%
  separate(sample, into = "group", sep = "_", remove = FALSE)
counts$sample <- factor(counts$sample, levels = c("Body_1", "Body_2", "Body_3", "Head_1", "Head_2", "Head_3"), 
                        labels = c("Body 1", "Body 2", "Body 3", "Head 1", "Head 2", "Head 3"))
counts$group <- factor(counts$group, levels = c("Body", "Head"), 
                       labels = c("Body", "Head"))
counts$Vaccine <- factor(counts$Vaccine)
counts$antigen <- factor(counts$antigen)
counts$drug_target <- factor(counts$drug_target)
counts$type <- factor(counts$type)
counts$subtype <- factor(counts$subtype)

# plot heatmap
library(wesanderson)
pal <- wes_palette("Zissou1", 100, type = "continuous")
heatmap <- ggplot(counts, aes(sample, gene_id)) +
  geom_tile(aes(fill = expression)) +
  scale_fill_gradientn(colours = pal, "Z-Score",
                       guide = guide_colorbar(
                         direction = "horizontal",
                         title.position = "left",
                         label.position = "bottom"
                       )) +
  scale_y_discrete(guide = guide_axis(n.dodge = 1)) +
  scale_x_discrete(labels = c("Body 1" = "","Body 2" = "Body","Body 3" = "",
                              "Head 1" = "","Head 2" = "Head","Head 3" = ""),
                   position = "top") +
  geom_vline(xintercept = 3.5, linetype="dashed", color = "black", size=0.5) +
  xlab("") + ylab("") + 
  theme_zlab_white() +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.text.x = element_text(size = 12, angle = 0),
    axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(face = "italic", size = 9, angle =0, hjust=1),
    strip.background = element_rect(fill="white", color = "white"),
    strip.text.x = element_text(size = 10),
    strip.text.y = element_text(size = 10),
    plot.margin = unit(c(0, 0, 0, 0), "cm"),
    panel.spacing=unit(-0.75, "lines"),
    legend.position = "bottom",
    legend.margin=margin(-5,0,0,0)) 
heatmap

heatmap_ann1 <- ggplot(counts) +
  geom_tile(aes(x = "Antigen", y = gene_id, fill = antigen)) +
  geom_tile(aes(x = "Vaccine", y = gene_id, fill = Vaccine)) +
  scale_fill_manual(values=c("black"),na.value="white") +
  scale_x_discrete(limits = c("Antigen","Vaccine")) +
  xlab("") + ylab("") + 
  theme_zlab_white() +
  guides(x=guide_axis(angle = 45)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
heatmap_ann1

library(ggnewscale)
pal_type <- c("lightsteelblue","seagreen2","firebrick4")
heatmap_ann2 <- ggplot(counts) +
  geom_tile(aes(x = "Type", y = gene_id, fill = as.character(factor(type)))) +
  scale_fill_manual(breaks = levels(counts$type), labels=c("\u03B2-tubulin","BK","LGIC"), values = pal_type, na.value="white") +
  scale_x_discrete(labels=c("Anthelmintic")) +
  xlab("") + ylab("") +
  theme_zlab_white() +
  guides(x=guide_axis(angle = 45)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    #axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    #axis.ticks.x = element_blank(),
    axis.ticks.y = element_blank(),
    strip.background = element_rect(fill="white", color = "white"),
    legend.title = element_blank(),
    legend.text=element_text(size=10),
    legend.position = "right", legend.box = 'vertical', 
    legend.justification = 'left', 
    legend.box.just = 'center',
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  ) 
heatmap_ann2

heatmap_sig <- ggplot(counts,aes(x = "p-value", y = gene_id)) +
  geom_tile(fill = "white") +
  scale_fill_manual(values=c("white"),na.value="white") +
  scale_x_discrete(limits = c("p-value")) +
  geom_text(aes(label=siglevel), size = 2.5) +
  xlab("") + ylab("") + 
  theme_zlab_white() +
  guides(x=guide_axis(angle = 45)) +
  theme(
    plot.background = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    legend.title = element_blank(),
    legend.position = "none",
    plot.margin = unit(c(0, 0, 0, 0), "cm")
  )
heatmap_sig


##########################################
### Figure 1
##########################################
#Fig1b <- volcano.ht + annotation_custom(ggplotGrob(pca.ht), xmin = -12, xmax = -2, ymin = 45, ymax = 70)

theme_set(theme_cowplot(font_family = "Helvetica") + 
            theme(text = element_text(colour = "black")))

library(Cairo)
Fig1b <- pca.ht
Fig1c <- volcano.ht 
Fig1d <- plot_grid(heatmap, NULL, heatmap_sig, NULL,heatmap_ann1, NULL, heatmap_ann2, nrow = 1, 
                   rel_widths = c(1.25, -0.05, 0.15, 0 ,0.25,-.05,0.52), axis = "tb", align = "h")

Fig1ab <- plot_grid(Fig1a,Fig1b, nrow =1, labels = c('A','B'), rel_widths = c(1.3,1))
FigL <- plot_grid(Fig1ab, Fig1c, nrow =2, labels = c('','C'), rel_heights = c(0.6,1))

Fig1 <- plot_grid(FigL, NULL, Fig1d, labels = c('','','D'), nrow = 1, 
                  rel_widths = c(0.9,0.01,1), scale = 0.99) +
  theme(plot.background = element_rect(fill = "white"))

#Fig1
#ggsave("Fig1.pdf", Fig1, width = 14.5, height = 8.25, units = "in")
ggsave("Fig1.png", Fig1, width = 14.5, height = 8.25, units = "in")
ggsave("Fig1.tiff", Fig1, width = 14.5, height = 8.25, units = "in")




##########################################
### Supp - GO enrichment of DEGs
##########################################

library(topGO)
library(conflicted)
library(Rgraphviz)
conflict_prefer("select", "dplyr")
conflict_prefer("filter", "dplyr")
conflict_prefer("rename", "dplyr")

### WormBase ParaSite species
library(biomaRt)
mart <- useMart("parasite_mart", dataset = "wbps_gene", host = "https://parasite.wormbase.org", port = 443)

# brugia
brugia_go <- getBM(mart = mart, 
                   filters = c("species_id_1010"),
                   value = list("brmalaprjna10729"),
                   attributes = c("wbps_gene_id", "wbps_transcript_id", "go_accession", "go_name_1006", "go_definition_1006", "go_linkage_type", "go_namespace_1003")) %>%
  janitor::clean_names()

# Note: we use gene_id for WBP species
brugia_go_out <- dplyr::select(brugia_go, gene_id = wbps_gene_id, go_id = go_accession) %>%
  group_by(gene_id) %>%
  distinct() %>%
  filter(go_id != "") %>%
  # mutate(transcript_id = str_remove(transcript_id, '\\.[0-9]*$')) %>%
  summarise(go_ids = list(go_id)) %>%
  mutate(go_ids = str_c(go_ids)) %>%
  mutate(go_ids = str_remove_all(go_ids, "c\\(")) %>%
  mutate(go_ids = str_remove_all(go_ids, '\\"')) %>%
  mutate(go_ids = str_remove_all(go_ids, '\\)'))

write.table(brugia_go_out, 'go/brugia_gene_go.txt', 
            sep = "\t",
            row.names = FALSE,
            col.names = FALSE,
            quote = FALSE)

## GO enrichment
# load in the mappings (this must be a file, can't be a df)
gene_go <- topGO::readMappings('go/brugia_gene_go.txt')

# read the df for later perusal
genes_go <- read_tsv('go/brugia_gene_go.txt', col_names = c('gene_id', 'go_ids')) %>%
  separate_rows(go_ids, sep = ', ') %>%
  rename(go_id = go_ids)

# get the list of possible gene_ids
gene_ids <- names(gene_go)

# read in your list of genes/transcripts of interest (upregulated)
interest_genes <- readRDS("deg/ht.res.rds") %>%
  filter(sig == "yes" & log2FoldChange < -1) %>% select(gene_id) %>% distinct(gene_id)
interest_go <- left_join(interest_genes, genes_go, by = c('gene_id')) 

go_summary <- group_by(interest_go, go_id) %>%
  summarize(n = n()) %>%
  filter(!is.na(go_id))

# the final data.frame needs to have one column with all the transcript_ids
# and a second column denoting whether or not it is a transcript of interest

final_genes <- distinct(select(genes_go, gene_id)) %>%
  mutate(interest = case_when(
    gene_id %in% interest_genes$gene_id ~ 1,
    TRUE ~ 0
  )) 

# get the interest column as a factor
final_genes_tg <- as.factor(final_genes$interest)

# convert to a factor with names
names(final_genes_tg) <- final_genes$gene_id

# create the topGOdata objects
# MF == molecular function, BP == biological process, CC == cellular component
go_data_mf <- new("topGOdata", ontology = "MF", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_bp <- new("topGOdata", ontology = "BP", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)
go_data_cc <- new("topGOdata", ontology = "CC", allGenes = final_genes_tg, annot = annFUN.gene2GO, gene2GO = gene_go)

# create the statistical test
fisher_test <- new("classicCount", testStatistic = GOFisherTest, name = "Fisher test")

all_go_data <- tibble(top_go_object = c(go_data_mf, go_data_bp, go_data_cc)) %>% # make a tibble of all 3 topGOobjects
  mutate(test = map(top_go_object, getSigGroups, fisher_test)) %>% # run the fisher test on each topGOobject
  mutate(result = map2(top_go_object, test, GenTable, ranksOf = "classic", topNodes = 10)) %>% # extract significant GO IDs to a df
  mutate(result = map2(top_go_object, result, ~mutate(.y, class = .x@ontology))) # add the GO class as a column for each nested df

plot_data <- select(all_go_data, result) %>%
  unnest(cols = c(result)) %>%
  janitor::clean_names() %>%
  rename(result = result1) %>%
  mutate(class = case_when(
    class == 'BP' ~ 'Biological Process',
    class == 'MF' ~ 'Molecular Function',
    class == 'CC' ~ 'Cellular Component'
  ))

go_plot <- ggplot(plot_data) + 
    geom_point(aes(x = term, y = -log10(as.numeric(result))), color = "black", size = 5) +
    facet_grid(rows = vars(class), scales = "free") +
    labs(y = "-log10(p-value)", x = "GO Term") + 
    coord_flip() +
    theme_minimal(base_size = 16, base_family = "Helvetica") +
    theme(legend.position = "none",
          axis.text.y = element_text(size = 8))
go_plot



