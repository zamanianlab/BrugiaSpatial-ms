library(tidyverse)
library(dplyr)
library(cowplot)
library(LaCroixColoR)

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure2")

theme_zlab <- function(base_size = 16, base_family = "helvetica") {
  (ggthemes::theme_foundation(base_size = base_size, base_family = base_family) +
     ggplot2::theme_bw(base_size = 16, base_family = "Helvetica") +
     ggplot2::theme(
       plot.title = ggplot2::element_text(size = 11, colour = "gray29", face = "bold", hjust = 0.5),
       # axes
       axis.text.x = ggplot2::element_text(size = 10, angle = 45, hjust = 1),
       axis.text.y = ggplot2::element_text(size = 10),
       axis.title.x = ggplot2::element_text(size = 10),
       axis.title.y = ggplot2::element_text(angle = 90, size = 10),
       axis.ticks = ggplot2::element_line(size = 0.25),
       #axis.line = ggplot2::element_line(size = 0.75, colour = "black"),
       
       # facets
       strip.text = ggplot2::element_text(face = "bold", size = 11),
       strip.background = element_rect(fill="white"),
       
       # panels
       panel.grid.minor = ggplot2::element_blank(),
       panel.grid.major = ggplot2::element_blank(),
       panel.background = ggplot2::element_blank(),
       
       # legend
       legend.title = element_text(size = 8),
       legend.text = element_text(size = 8)
       #legend.position = "none"
     ))
}

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
      #panel.border = element_rect(fill = NA, color = "black"),  
      panel.border = element_blank(),
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
## FIG 2A Clade III Comparative Anatomy 
############

setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Microscopy/")

anat <- read.csv("Comparative_Anatomy.csv", header = TRUE, sep = ",")
anat$Stage <- factor(anat$Stage, levels = c("mf","L2/3","AF"), 
                     labels = c("mf/L1","L2/L3","Adult Female")) 
anat <- anat %>% dplyr::select(Node,Species,Stage,Feature,Start=Start_R,End=End_R,Order) %>%
  filter(Species %in% c('Dracunculus globacephalus','Dracunculus lutrae','Parascaris equorum','Oxyspirura mansoni','Anisakis simplex','Toxocara canis','Gongylonema pulchrum') == FALSE)

anat$Species <- fct_reorder(anat$Species, anat$Node)
anat$Species <- fct_rev(anat$Species)
anat$Order <- factor(anat$Order, levels = c("Spirurida","Ascaridida","Oxyurida"))
anat$Order <- factor(anat$Order, levels = c("Spirurida","Ascaridida","Oxyurida"))

#plot counts
anat.plot <- ggplot(data=anat, aes(y=Species)) +
  theme_zlab_white() +
  theme(axis.text.y = element_text(face = "italic", hjust =1),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        legend.title = element_blank(),
        panel.spacing = unit(1, "lines"),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75)) +
  scale_color_manual(values=lacroix_palette("PassionFruit",type = "continuous", n=6)) +
  geom_point(data = filter(anat, Stage %in% c("Adult Female")),aes(x = (Start+End)/2, y = Species, color = Feature), size = 5, alpha = 0.75) +
  geom_point(data = filter(anat, Stage %in% c("mf/L1","L2/L3")),aes(x = 500*(Start), y = Species, color = Feature), size = 4.5, alpha = 0.75, show.legend = FALSE) +
  #geom_segment(data = filter(anat, Stage %in% c("mf/L1","L2/L3")),
  #             aes(x = Start*400, y = Species, xend = End*400, yend = Species,color = Feature), size = 5, alpha = 0.75, show.legend = FALSE) +
  ylab("") + xlab("Anterior (L) to Posterior (R)") +
  facet_grid(Order ~ Stage,  scales = "free", space="free")
anat.plot


############
## FIG 2B LIGHT-SHEET IMAGE ES PORE
############

library(magick)
library(pdftools)
library(grConvert)
library(grImport2)

LS.pulse <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Microscopy/Lightsheet/LightSheet_Pulse2.pdf", density = 600)
Fig2b <- ggdraw() + 
  draw_image(LS.pulse, scale = 1) 


############
## FIG 2C/D SBF-SEM / TEM
############

SBFSEM <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Microscopy/Bmalayi_AF_SBFSEM/Figure2C_SBF_SEM.pdf", density = 600)
Fig2c <- ggdraw() + 
  draw_image(SBFSEM, scale = 1) 

TEM <- image_read_pdf("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Microscopy/Bmalayi_AF_SBFSEM/Figure2D_TEM.pdf", density = 600)
Fig2d <- ggdraw() + 
  draw_image(TEM, scale = 1) 



##########################################
### Figure 2
##########################################
setwd("~/Box/ZamanianLab/Data/Airs_Experiments/Bm_spatial_trans-ms/Figures/Figure2/")

Fig2a <- anat.plot

Fig2top <- plot_grid(Fig2a + theme(plot.margin = unit(c(0.1, 0.1, 0.3, -0.5), "cm")), NULL, Fig2b, nrow = 1, labels = c('A','','B'), hjust = 0, label_x = 0.01, rel_widths = c(1.2,0.02,0.5), scale = 1)

theme_set(theme_cowplot(font_family = "Helvetica") + 
            theme(text = element_text(colour = "white")))
Fig2bot <- plot_grid(Fig2c, NULL, Fig2d, nrow = 1, labels = c('C','','D'), hjust = -2.5, label_x = 0.01, rel_widths = c(0.69,0.01,1), scale = 1)

Fig2 <- plot_grid(Fig2top, Fig2bot, nrow = 2, labels = c('',''), rel_heights = c(1,1.1), scale = 0.99)

theme_set(theme_cowplot(font_family = "Helvetica") + 
            theme(text = element_text(colour = "black")))
ggsave("Figure2.pdf", Fig2, width =12, height = 13.5, units = "in")



