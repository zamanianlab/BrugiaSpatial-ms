library(tidyverse)
library(cowplot)
library(ggrepel)
library("ZamanianLabThemes")
library(paletteer)
library(DESeq2)
library(hexbin)
library(wesanderson)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure3")

theme_zlab_black = function(base_size = 12, base_family = "Helvetica") {
  
    theme_bw(base_size = base_size, base_family = base_family) %+replace%
    theme(
      # Specify axis options
        #axis.line = element_blank(),
      axis.line = element_line(color = "white"),
      axis.text.x = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.text.y = element_text(size = base_size*0.8, color = "white", lineheight = 0.9),  
      axis.ticks = element_line(color = "white", size  =  0.4),  
      axis.title.x = element_text(size = base_size, color = "white", margin = margin(0, 10, 0, 0)),  
      axis.title.y = element_text(size = base_size, color = "white", angle = 90, margin = margin(0, 10, 0, 0)),  
      axis.ticks.length = unit(0.3, "lines"),   
      # Specify legend options
      legend.background = element_rect(color = NA, fill = "black"),  
      legend.key = element_rect(color = "white",  fill = "black"),  
      legend.key.size = unit(1.2, "lines"),  
      legend.key.height = NULL,  
      legend.key.width = NULL,      
      legend.text = element_text(size = base_size*0.8, color = "white"),  
      legend.title = element_text(size = base_size*0.8, face = "bold", hjust = 0, color = "white"),  
      legend.position = "right",  
      legend.text.align = NULL,  
      legend.title.align = NULL,  
      legend.direction = "vertical",  
      legend.box = NULL, 
      # Specify panel options
      panel.background = element_rect(fill = "black", color  =  NA),  
        #panel.border = element_rect(fill = NA, color = "white"),  
      panel.border = element_blank(),
        #panel.grid.major = element_line(color = "grey35"),  
      panel.grid.major = element_blank(), 
        #panel.grid.minor = element_line(color = "grey20"),  
      panel.grid.minor = element_blank(),  
      panel.spacing = unit(0.5, "lines"),   
      # Specify faceting options
      strip.background = element_rect(fill = "grey30", color = "grey10"),  
      strip.text.x = element_text(size = base_size*0.8, color = "white"),  
      strip.text.y = element_text(size = base_size*0.8, color = "white",angle = -90),  
      # Specify plot options
      plot.background = element_rect(color = "black", fill = "black"),  
      plot.title = element_text(size = base_size*1.2, color = "white"),  
      plot.margin = unit(rep(1, 4), "lines")
      
    )
  
}

############
## FIG 3A 3D cryosection render + panel
############

library(magick)
library(pdftools)
library(grConvert)
library(grImport2)

RNAt.3d <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/RNA_Tomography/FluorescentSectionsLocations.pdf", density = 600)
Fig3a <- ggdraw() + 
  draw_image(RNAt.3d, scale = 1) 


############
## FIG 3B/C Nuclei density plots from cryosection and LS replicates
############

nuclei.cryo <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/RNAseq/RNA_Tomography/FluorescentSlices_NucleiCounts.csv", header=T) 
# adult based on vulva location, worm 3 vulva = 35, worm 4 = 36, worm 6 = 39, worm 7 = 60
nuclei.cryo <- nuclei.cryo %>%
  mutate(Slice_adj = ifelse(Worm == 3, Slice-35,
                            ifelse(Worm == 4, Slice - 36,
                                   ifelse(Worm == 6, Slice - 39,
                                          ifelse(Worm == 7, Slice -60, Slice)))))

nuclei.ls <- read.csv("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Microscopy/LightSheet/Lightsheet_Nuclei_Counts.csv", header=T) 
nuclei.ls <- nuclei.lcm %>% select(1:5) %>%
  mutate(intercept = ifelse(Sex == "Female", 600, NA))
nuclei.ls$Worm <- factor(nuclei.lcm$Worm)
nuclei.ls$Sex <- factor(nuclei.lcm$Sex)

#plot nuclei distribution data for cryo
nuclei.cryo.plot <- ggplot(nuclei.cryo, aes(x = Slice_adj, y = Nuclei)) +
  geom_vline(xintercept = 0, color = "#E574FF") +
  #annotate(geom="text", x=-1, y= 30, label="Vulva",color="#E574FF", size = 4) + 
  geom_line(aes(group=Worm), colour = "grey", size =0.5, alpha = 0.5, linetype = "dotted") +
  geom_boxplot(aes(group=Slice_adj), colour = "white", fill = "grey", size =0.5, alpha = 0.75) +
  geom_smooth(se = FALSE, color = "red",span = 0.5, size = 0.5, linetype = "solid") +
  scale_x_reverse() + coord_flip() +
  xlim(25,-30) + ylim(0,40) +
  theme_zlab_black() +
  xlab("Cryosection position relative to vulva (20 µm units)") + ylab("Nuclei per cryosection") +
  theme(legend.position= "none") 
nuclei.cryo.plot


#devtools::install_github("easystats/see")
#install.packages("ggdist")
library(see)
library(ggdist)
#plot nuclei distribution data for lcm
nuclei.ls.plot <- ggplot(nuclei.ls, aes(x = Worm, y = Distance)) +
  geom_jitter(data = filter(nuclei.ls,Worm %in% c("0516_E4","0329_E6")), colour = "white", size = 0.75, alpha = 0.8, height = 0, width = 0.1) +
  stat_histinterval(data = filter(nuclei.ls,Worm %in% c("0516_E4","0329_E6")), alpha = 1, slab_type = "histogram", 
                    breaks = 30, slab_color = "white", outline_bars = FALSE, n=100, width = 0.25, justification = -0.75, interval_size = 0, interval_alpha = 0) +
  theme_zlab_black() +
  scale_y_reverse() +
  ylim(600,0) +
  geom_hline(aes(yintercept=intercept), color = "#E574FF") + 
  theme(axis.text.x = element_blank(),  
        axis.ticks.x = element_blank()) +
        #axis.title.x = element_blank() +
  ylab("Distance from anterior (µm)") + xlab("Lightsheet nuclei density") +
  theme(legend.position= "none")  +
  facet_grid(. ~ Sex , scale = "free", space="free")
nuclei.ls.plot


# nuclei density stats (cryo)
nuclei.cryo.stats <- nuclei.cryo %>% 
  filter(Slice_adj < 0) %>%
  group_by(Worm) %>%
  summarise(nuclei.mean = mean(Nuclei))
mean.cryo = mean(nuclei.cryo.stats$nuclei.mean) #10.15
sd.cryo = sd(nuclei.cryo.stats$nuclei.mean) #1.46

# nuclei density stats (ls)

#female stats
nuclei.ls.stats <- nuclei.ls %>% 
  filter(Sex == "Female") %>%
  group_by(Worm) %>%
  summarise(nuclei.total = n(), dist.total = max(Distance)) %>%
  mutate(nuclei.dens = 20*(nuclei.total / dist.total))
# mean nuclei density for female (9.70)

#male stats
nuclei.ls.stats <- nuclei.ls %>% 
  filter(Sex == "Male") %>%
  group_by(Worm) %>%
  summarise(nuclei.total = n(), dist.total = max(Distance)) %>%
  mutate(nuclei.dens = 20*(nuclei.total / dist.total))
mean.ls.male = mean(nuclei.ls.stats$nuclei.dens) #9.23
sd.ls.male = sd(nuclei.ls.stats$nuclei.dens) #1.86



##########################################
### Figure 3
##########################################

Fig3a <- Fig3a 
Fig3b <- nuclei.cryo.plot
Fig3c <- nuclei.ls.plot

theme_set(theme_cowplot(font_family = "Helvetica") + 
            theme(text = element_text(colour = "white")))

Fig3 <- plot_grid(Fig3a, NULL, Fig3b, NULL, Fig3c, nrow = 1, labels = c('A','','B','','C'), rel_widths = c(0.6,0.04,1,0.04,1), scale = 0.99)

ggsave("Fig3.pdf", Fig3, width =14, height = 8, units = "in", bg = "black")



