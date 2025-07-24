##### Libraries #####
library("dplyr")
library("ggplot2")
library("ggpubr")

##### Barplot ######
##COI##
#occurrence
Phyla_Cobble <- AllSp_Cobble_COI[,c(2,3)]
Phyla_Cobble$nrMOTUs <- 1
Phyla_Cobble <- aggregate(. ~ Phylum, data = Phyla_Cobble, sum)
Phyla_Cobble$Method <- "Cobble"
Phyla_water <- AllSp_water_COI[,c(2,3)]
Phyla_water$nrMOTUs <- 1
Phyla_water <- aggregate(. ~ Phylum, data = Phyla_water, sum)
Phyla_water$Method <- "water"
Phyla_DAM <- AllSp_DAM_COI[,c(2,3)]
Phyla_DAM$nrMOTUs <- 1
Phyla_DAM <- aggregate(. ~ Phylum, data = Phyla_DAM, sum)
Phyla_DAM$Method <- "DAM"
Phyla_BAM <- AllSp_BAM_COI[,c(2,3)]
Phyla_BAM$nrMOTUs <- 1
Phyla_BAM <- aggregate(. ~ Phylum, data = Phyla_BAM, sum)
Phyla_BAM$Method <- "BAM"

Phyla_All <- rbind(Phyla_BAM,Phyla_Cobble,Phyla_DAM,Phyla_water)
colnames(Phyla_All) <- c("Phyla", "Reads", "MOTUs", "Methods")

Phyla_All_percent <- Phyla_All[,c(1,2,3)]
Phyla_All_percent <- aggregate(. ~ Phyla, data = Phyla_All_percent, sum)
Sum_MOTUs <- sum(Phyla_All_percent$MOTUs)
Sum_reads <- sum(Phyla_All_percent$Reads)
Phyla_All_percent$Reads <- Phyla_All_percent$Reads/Sum_reads
Phyla_All_percent$MOTUs <- Phyla_All_percent$MOTUs/Sum_MOTUs

palette_Phyla = c("#bcbd22", "#aa3c9d",  "darkblue", "darkcyan", "#ffc300", "turquoise", "#8c564b",   "#ff7f0e","#f7b6d2", "blue",
                  "#dbdb8d", "red", "#ffbb89",  "darkgreen", "tomato", "deepskyblue",  "darkgoldenrod","darkred", "#393b79",  
                  "darkviolet", "plum", "#98df8a", "#e5716b",  "#7f7f7f", "#538136", "lightblue", "orange", "forestgreen", 
                  "yellow", "#2d5da1",  "#f5d55c", "#17becf", "#c70039","#9467bd", "lawngreen", "#e0a18e")

Barplot_Phyla_MOTUs <- ggplot(Phyla_All[,-2], aes(fill=Phyla, y= MOTUs, x=Methods)) + 
  scale_fill_manual(values = palette_Phyla) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(colour = "black", size = 10, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 10), 
        panel.background = element_blank(), 
        legend.position = "none")

Barplot_Phyla_abund <- ggplot(Phyla_All[,-3], aes(fill=Phyla, y= Reads, x=Methods)) + 
  scale_fill_manual(values = palette_Phyla) +
  geom_bar(position="fill", stat="identity") +
  theme(axis.text.x = element_text(colour = "black", size = 10, vjust = 0.3, hjust = 1), 
        axis.text.y = element_text(colour = "black",  size = 10), 
        panel.background = element_blank(), 
        legend.position = "none")

Barplot_Phyla_all = ggarrange(Barplot_Phyla_MOTUs, Barplot_Phyla_abund, 
                              nrow = 1, ncol = 2, common.legend = TRUE, legend="right")


