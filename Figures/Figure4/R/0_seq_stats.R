library(tidyverse)
library(cowplot)
library(gdata) #for read.xls
library(wesanderson)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure4")
RNAdata <- c("~/Box/ZamanianLab/SeqLibraries/Mapping/")

# run bash one-liners to aggregate STAR mapping stats for each seq run

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

#############
#load in STAR / qc stats for rep 1 - 4
#############

#rep 1
dir <- c("O002_BmAF_RNAt/181220_BCCV00ANXX/")
star_dir <- paste0(RNAdata,dir,"star/")
star_total <- read.csv(paste0(star_dir,"STAR_total_reads.csv"), header=F)
star_unique <- read.csv(paste0(star_dir,"STAR_unique_mapped.csv"), header=F)
star_multi <- read.csv(paste0(star_dir,"STAR_multi_mapped.csv"), header=F)
star_stats <- left_join(star_total,star_unique, by = "V1") %>% left_join(.,star_multi, by ="V1")
colnames(star_stats) <- c("sample","total_reads","unique_reads","multi_reads")
star_r1 <- star_stats %>% mutate(rep = "R1") %>%
  separate(sample, c("rep", "sample"), sep = "_", remove = TRUE) %>%
  mutate(rep = "R1") %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = "")) 


#rep 2/3
dir <- c("Q002/200217_AHNHN3DMXX/")
star_dir <- paste0(RNAdata,dir,"star/")
star_total <- read.csv(paste0(star_dir,"STAR_total_reads.csv"), header=F)
star_unique <- read.csv(paste0(star_dir,"STAR_unique_mapped.csv"), header=F)
star_multi <- read.csv(paste0(star_dir,"STAR_multi_mapped.csv"), header=F)
star_stats <- left_join(star_total,star_unique, by = "V1") %>% left_join(.,star_multi, by ="V1")
colnames(star_stats) <- c("sample","total_reads","unique_reads","multi_reads")
star_r23 <- star_stats %>% mutate(rep = "R1") %>%
  separate(sample, c("rep", "sample"), sep = "_", remove = TRUE) %>%
  mutate(rep = ifelse(rep == "P1", "R2", "R3")) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = "")) 


#combine reps
star_stats <- rbind(star_r1,star_r23) %>%
  mutate(unique_perc = unique_reads/total_reads, mapped_perc = (unique_reads+multi_reads)/total_reads) %>%
  gather(metric,value,total_reads:multi_reads,unique_perc:mapped_perc) 
star_stats$sample <- factor(as.integer(star_stats$sample))


#############
#load in library qc stats for rep 2/3 (don't have for r1)
#############

#rep 1 (don't have normal excel file)
#load in total RNA yield (bioanalyzer, 400-9500bp range, pg/uL) for just P1
qc_r1 <- read.csv("~/Box/ZamanianLab/SeqLibraries/Submission/O002_BmAF_RNAt/QC/RNAyield_stats.csv", header = FALSE, sep = ",")
colnames(qc_r1) <- c("sample","rna")
qc_r1 <- qc_r1 %>% 
  arrange(sample) %>%
  mutate(rep = "R1")
qc_r1$sample <- as.numeric(qc_r1$sample)
qc_r1 <- qc_r1 %>% 
  mutate(sample = ifelse(sample > 10, sample -1, sample))
qc_r1$sample <- as.factor(as.numeric(qc_r1$sample))

#rep 2/3
qc_r23 <- read.xls("~/Box/ZamanianLab/SeqLibraries/Submission/Q002_20191223/Q002.xlsx", sheet = 1, header = TRUE, stringsAsFactors=FALSE) %>%
  select(3,18:21)
colnames(qc_r23) <- c("id","pcr_1","qc_1","pcr_2","qc_2") 
qc_r23 <- qc_r23 %>%
  separate(id, c("rep", "sample"), sep = "_", remove = TRUE) %>%
  mutate(rep = ifelse(rep == "P1", "R2", "R3")) %>%
  mutate(sample = gsub(x = sample, pattern = "S", replacement = "")) 
qc_r23$qc_1 <- as.numeric(qc_r23$qc_1)
qc_r23$sample <- as.factor(as.numeric(qc_r23$sample))

#combine reps
qc_r1 <- qc_r1 %>%
  gather(metric,value,rna)
qc_r23 <- qc_r23 %>%
  gather(metric,value,pcr_1:qc_2)

library_qc <- rbind(qc_r1,qc_r23)
library_qc$sample <- factor(as.integer(library_qc$sample))




##########
# Combine QC data and plot
##########

#combine
run_stats <- rbind(star_stats,library_qc)

#order and rename metrics
run_stats$metric <- factor(run_stats$metric, levels=c("rna","pcr_1","qc_1","pcr_2","qc_2","total_reads","unique_reads","multi_reads","unique_perc","mapped_perc"))

#export for figure generation
saveRDS(run_stats, "data/run_stats.rds")





#plot total reads and frac
plot_stats <- ggplot(run_stats, aes(x = sample, y = value, group = rep()))+
  theme_zlab_white() +
  geom_line(aes(group = rep, colour = rep), size = 0.5, alpha=0.75) +
  #geom_line(size = 0.5, alpha=0.75, linetype = "dashed", colour = "darkorange4") +
  geom_point(aes(group = rep, colour = rep, fill = rep), shape=21, size=2.5, alpha=0.75) +
  #geom_point(shape=21, size=3, alpha=0.75, colour = "darkorange4", fill = "darkorange2") +
  #ylab(expression(Total~Reads~(10^{"6"}))) +
  ylab('') +
  xlab("Slice Position") +
  theme(axis.text.x = element_text(size=5),
        axis.text.y = element_text(face="bold", size=8),
        axis.title.x  = element_text(face="bold", size=12, vjust = 0.1), 
        axis.title.y  = element_text(face="bold", angle=90, size=12)) +
  theme(axis.title.x=element_blank()) +
  facet_grid(metric ~ rep, scales = "free_y")
plot_stats

ggsave("seq_stats.pdf", plot_stats, width = 14, height = 10, units = c("in"))


