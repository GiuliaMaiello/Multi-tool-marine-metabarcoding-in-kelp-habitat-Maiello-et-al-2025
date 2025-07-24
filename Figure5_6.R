##### Libraries #####
library("dplyr")
library("ggforce")
library("ggnewscale")
library("ggpubr")
library("ggplot2")
library("ggrepel")
library("gridExtra")
library("indicspecies")
library("plyr")
library("reshape2")
library("stringr")
library("vegan")
library("tibble")

##### Bubble Plot #####
## 12S ##
CHNMS_data_bubble <- CHNMS_data[-which(CHNMS_data$Source == 'FRAM'),]
CHNMS_data_bubble_NoFRAM <- CHNMS_data_bubble[,-c(1,2)]
CHNMS_data_bubble_NoFRAM <- CHNMS_data_bubble_NoFRAM[,-which(apply(CHNMS_data_bubble_NoFRAM,2,sum) == 0)]
CHNMS_data_bubble <- cbind(CHNMS_data_bubble[,1:2],CHNMS_data_bubble_NoFRAM)
CHNMS_data_bubble = melt(CHNMS_data_bubble, id = c("Site", "Source"))
CHNMS_data_bubble$value <- sqrt(CHNMS_data_bubble$value)
CHNMS_data_bubble$Site <- factor(CHNMS_data_bubble$Site,levels=unique(CHNMS_data_bubble$Site))
colnames(CHNMS_data_bubble) <- c("Site", "Source", "Species", "Reads")
CHNMS_data_bubble$Species<- as.character(CHNMS_data_bubble$Species)
CHNMS_data_bubble <- CHNMS_data_bubble[order(CHNMS_data_bubble$Species),]
CHNMS_data_bubble$Class <- as.character(CHNMS_data_bubble$Class)
CHNMS_data_bubble <- CHNMS_data_bubble[order(CHNMS_data_bubble$Class),]

Bubble_12S <- ggplot(CHNMS_data_bubble, aes(x = Site, y = ordered(Species, 
                                                                 levels=sort(unique(Species), decreasing = F)))) + 
  scale_y_discrete(limits = rev) +
  geom_point(aes(size = reads, fill = Source), alpha = 0.8, shape = 21) + 
  scale_fill_manual(values = c("darkorange2", "mediumpurple2", "deepskyblue2"), guide = guide_legend(override.aes = list(size=5))) +
  scale_size_continuous(limits = c(1, 400), range = c(1,10), breaks = c(1,100,200,300)) +  
  labs( x= "", y = "", size = expression(sqrt("Nr. of reads")), fill = "")  + 
  facet_grid(Class ~ Region, scales = "free", space = "free") +
  theme_light() +
  theme(strip.background.y = element_rect(fill= "lightgrey"),
        strip.background.x = element_rect(fill= "lightgrey"),
        strip.text.x= element_text(colour = "#000000"),
        panel.grid.major = element_line(colour = "white"),
        panel.grid.minor = element_line(colour = "white"), 
        axis.ticks = element_line(colour = "#000000"),
        strip.text.y = element_text(colour = "#000000", angle = 360),
        axis.text.x = element_text(colour = "#000000", angle = 90, hjust = 1),
        axis.text.y = element_text(colour = "#000000", face = "italic"),
        text = element_text(size = 10),
        legend.direction = "vertical", 
        legend.box = "vertical") 

##### nMDS #####
## 12S ##
## dataset preparation
CHNMSSp_NMDS <- melt(CHNMS_data, id.vars = c("Site", "Source"))
colnames(CHNMSSp_NMDS) <- c("Site", "Source", "Species", "Reads")
CHNMSSp_NMDS$Reads <- sqrt(CHNMSSp_NMDS$Reads)
CHNMSSp_NMDS$Site <- substr(CHNMSSp_NMDS$Site, 1, 6)
CHNMSSp_NMDS <- CHNMSSp_NMDS[-which(CHNMSSp_NMDS$Source == "FRAM"),]
CHNMS_data_nMDS <- dcast(CHNMSSp_NMDS, Site + Source ~ Species, value.var = "Reads", fun.aggregate = mean)

CHNMS_data_nMDS_NoFRAM <- CHNMS_data_nMDS[,-c(1,2)]
CHNMS_data_nMDS_NoFRAM <- CHNMS_data_nMDS_NoFRAM[,-which(apply(CHNMS_data_nMDS_NoFRAM,2,sum) == 0)]
CHNMS_data_nMDS <- cbind(CHNMS_data_nMDS[,1:2],CHNMS_data_nMDS_NoFRAM)
values <- c("North", "South")
repetitions <- c(15, 6)
CHNMS_data_nMDS$GeoLocation <- rep(values, repetitions)

## nMDS calucation and plotting 
NMDS_CHNMS = metaMDS(CHNMS_data_nMDS[,-c(1:2, 46)], distance="bray", k = 3, trymax = 200, autotransform = F)

NMDSspecies = data.frame(Species=rownames(NMDS_CHNMS$species),
                         x=as.numeric(NMDS_CHNMS$species[,1]),
                         y=as.numeric(NMDS_CHNMS$species[,2]))

NMDSstaz = data.frame(Site=CHNMS_data_nMDS$Site,  Source = CHNMS_data_nMDS$Source, Group = CHNMS_data_nMDS$GeoLocation,
                      x=as.numeric(NMDS_CHNMS$points[,1]), 
                      y=as.numeric(NMDS_CHNMS$points[,2]))

NMDSstaz$Source = as.character(NMDSstaz$Source)
NMDSstaz$NrSpecies = rowSums(CHNMS_data_nMDS[,-c(1,2, 46)] > 0)

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls_Group <- ddply(NMDSstaz, .(Group), find_hull)
hulls_Group$Group = factor(hulls_Group$Group,
                           levels = sort(unique(hulls_Group$Group)))


NMDS_CHNMS_plot =  ggplot(data = NMDSstaz, aes(x=x, y=y)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(data = NMDSstaz, aes(x=x, y=y, fill = Source, size = NrSpecies), shape = 21) + 
  scale_fill_manual(values = c("darkorange2", "mediumpurple2", "deepskyblue2")) +
  geom_polygon(data = hulls_Group, aes(x=x, y=y, group=Group, colour = Group), 
               size = 1.2, linetype = 2, fill = NA, alpha = 0.6) +
  scale_color_manual(values = c("darkgoldenrod1", "forestgreen")) + 
  geom_text_repel(data = NMDSstaz, aes(x=x, y=y, label=Site), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = -0.7, y = 0.9, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_CHNMS$stress,2)))

Xs = CHNMS_data_nMDS[,c(1,2)]
Xs$Location <- NMDSstaz$Group
Ys = CHNMS_data_nMDS[,-c(1:2,46)]
permanova = adonis2(Ys ~ Xs$Source + Xs$Location, permutations = 9999, method = "bray")
permanova_env = adonis2(Ys ~ Xs$Location, permutations = 9999, method = "bray")

### Indicator Species Analysis
IndicSp <- CHNMS_data_nMDS[,c(1,3:46)]

NorthSp <- IndicSp[c(1:15),]
SpRichenessNorth <- apply(NorthSp[,-c(1,45)],2,sum)
SpRichenessNorth <- which(SpRichenessNorth > 0)

SouthSp <- IndicSp[c(16:21),]
SpRichenessSouth <- apply(SouthSp[,-c(1,45)],2,sum)
SpRichenessSouth <- which(SpRichenessSouth > 0)

abund = IndicSp[,2:44]
GeoLocation <- IndicSp$GeoLocation

inv = multipatt(abund, IndicSp$GeoLocation, control = how(nperm=9999))
summary(inv)

## COI all##
## dataset preparation
CHNMSSp_NMDS_COI <- melt(CHNMS_COIdata, id.vars = c("Site", "Source"))
colnames(CHNMSSp_NMDS_COI) <- c("Site", "Source", "Species", "Reads")
CHNMSSp_NMDS_COI$Reads <- sqrt(CHNMSSp_NMDS_COI$Reads)
CHNMSSp_NMDS_COI$Reads = 1*(CHNMSSp_NMDS_COI$Reads>0) 

CHNMS_data_nMDS_COI <- dcast(CHNMSSp_NMDS_COI, Site + Source ~ Species, value.var = "Reads", fun.aggregate = mean)

#jaccard
CHNMSSp_NMDS_COI <- melt(CHNMS_data_nMDS_COI, id.vars = c("Site", "Source"))
colnames(CHNMSSp_NMDS_COI) <- c("Site", "Source", "Species", "Reads")
CHNMSSp_NMDS_COI$Reads <- 1*(CHNMSSp_NMDS_COI$Reads>0)
CHNMS_data_nMDS_COI <- dcast(CHNMSSp_NMDS_COI, Site + Source ~ Species, value.var = "Reads", fun.aggregate = sum)

values <- c("North", "South")
repetitions <- c(21, 7)
CHNMS_data_nMDS_COI$GeoLocation <- rep(values, repetitions)

## nMDS calucation and plotting 
NMDS_CHNMS_COI = metaMDS(CHNMS_data_nMDS_COI[,-c(1:2, 458)], distance="jaccard", k = 3, trymax = 200, autotransform = F)

NMDSspecies = data.frame(Species=rownames(NMDS_CHNMS_COI$species),
                         x=as.numeric(NMDS_CHNMS_COI$species[,1]),
                         y=as.numeric(NMDS_CHNMS_COI$species[,2]))

NMDSstaz = data.frame(Site=CHNMS_data_nMDS_COI$Site,  Source = CHNMS_data_nMDS_COI$Source, Group = CHNMS_data_nMDS_COI$GeoLocation,
                      x=as.numeric(NMDS_CHNMS_COI$points[,1]), 
                      y=as.numeric(NMDS_CHNMS_COI$points[,2]))

NMDSstaz$Source = as.character(NMDSstaz$Source)
NMDSstaz$NrSpecies = rowSums(CHNMS_data_nMDS_COI[,-c(1,2, 458)] > 0)

NMDS_CHNMS_plot_COI =  ggplot(data = NMDSstaz, aes(x=x, y=y, xmin -3.75, xmax = 2, ymin= -1, ymax=1)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(data = NMDSstaz, aes(x=x, y=y, fill = Source, size = NrSpecies), shape = 21) + 
  scale_fill_manual(values = c("orange", "mediumseagreen","mediumpurple2", "deepskyblue2")) +
  geom_text_repel(data = NMDSstaz, aes(x=x, y=y, label=Site), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 1.5, y = 1.1, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_CHNMS_COI$stress,2)))

Xs = CHNMS_data_nMDS_COI[,c(1,2)]
Ys = CHNMS_data_nMDS_COI[,-c(1:2,458)]
permanova = adonis2(Ys ~  Xs$Source, permutations = 9999, method = "jaccard")

## COI eDNA##
## dataset preparation
CHNMSSp_NMDS_eDNA <- CHNMSSp_NMDS_COI[-which(CHNMSSp_NMDS_COI$Source == 'Cobble'),]
CHNMSSp_NMDS_eDNA <- CHNMSSp_NMDS_eDNA[-which(CHNMSSp_NMDS_eDNA$Source == 'DAM'),]
CHNMS_data_nMDS_eDNA <- dcast(CHNMSSp_NMDS_eDNA, Site + Source  ~ Species, value.var = "Reads", fun.aggregate = mean)

CHNMS_data_nMDS_species <- CHNMS_data_nMDS_eDNA[,-c(1,2)]
CHNMS_data_nMDS_species = CHNMS_data_nMDS_species[, colSums(CHNMS_data_nMDS_species) > 0]
CHNMS_data_nMDS_eDNA <- cbind(CHNMS_data_nMDS_eDNA[,c(1,2)], CHNMS_data_nMDS_species)
values <- c("North", "South")
repetitions <- c(9, 7)
CHNMS_data_nMDS_eDNA$GeoLocation <- rep(values, repetitions)

## nMDS calucation and plotting 
NMDS_CHNMS_COI = metaMDS(CHNMS_data_nMDS_eDNA[,-c(1,2, 175)], distance="jaccard", k = 3, trymax = 200, autotransform = F)

NMDSspecies_COI = data.frame(Species=rownames(NMDS_CHNMS_COI$species),
                             x=as.numeric(NMDS_CHNMS_COI$species[,1]),
                             y=as.numeric(NMDS_CHNMS_COI$species[,2]))

NMDSstaz_COI = data.frame(Site=CHNMS_data_nMDS_eDNA$Site, Source =CHNMS_data_nMDS_eDNA$Source, Group = CHNMS_data_nMDS_eDNA$GeoLocation,
                          x=as.numeric(NMDS_CHNMS_COI$points[,1]), 
                          y=as.numeric(NMDS_CHNMS_COI$points[,2]))

NMDSstaz_COI$Group= as.character(NMDSstaz_COI$Group)
NMDSstaz_COI$NrSpecies = rowSums(CHNMS_data_nMDS_eDNA[,-c(1,2,175)] > 0)

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls_Group <- ddply(NMDSstaz_COI, .(Group), find_hull)
hulls_Group$Group = factor(hulls_Group$Group,
                           levels = sort(unique(hulls_Group$Group)))

NMDS_CHNMS_plot_COI_eDNA =  ggplot(data = NMDSstaz_COI, aes(x=x, y=y, xmin = -0.7, xmax = 0.7, ymin = -0.85, ymax = 0.85))+
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_point(data = NMDSstaz_COI, aes(x=x, y=y, fill = Source, size = NrSpecies), shape = 21) + 
  scale_fill_manual(values = c("darkorange2", "deepskyblue2")) +
  geom_polygon(data = hulls_Group, aes(x=x, y=y, group=Group, colour = Group), 
               size = 1.2, linetype = 2, fill = NA, alpha = 0.6) +
  scale_color_manual(values = c("darkgoldenrod1", "forestgreen")) + 
  geom_text_repel(data = NMDSstaz_COI, aes(x=x, y=y, label=Site), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 0.5, y = 0.7, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_CHNMS_COI$stress,2)))

Xs = CHNMS_data_nMDS_eDNA[,c(1,2)]
Xs$Location <- NMDSstaz_COI$Group
Ys = CHNMS_data_nMDS_eDNA[,-c(1:2,175)]
permanova = adonis2(Ys ~  Xs$Source, permutations = 9999, method = "jaccard")
permanova_env = adonis2(Ys ~ Xs$Location, permutations = 9999, method = "jaccard")

### Indicator Species Analysis
IndicSp_eDNA <- CHNMS_data_nMDS_eDNA[,c(1,3:175)]
IndicSp_eDNA$GeoLocation <- CHNMS_data_nMDS_eDNA$GeoLocation

abund = IndicSp_eDNA[,2:173]

inv = multipatt(abund, IndicSp_eDNA$GeoLocation, control = how(nperm=9999))
summary(inv)

## COI only water##
## dataset preparation
CHNMSSp_NMDS_water <- CHNMSSp_NMDS_COI[which(CHNMSSp_NMDS_COI$Source == 'water'),]
CHNMS_data_nMDS_water <- dcast(CHNMSSp_NMDS_water, Site + Source ~ Species, value.var = "Reads", fun.aggregate = mean)

CHNMS_data_nMDS_species <- CHNMS_data_nMDS_water[,-c(1,2)]
CHNMS_data_nMDS_species = CHNMS_data_nMDS_species[, colSums(CHNMS_data_nMDS_species) > 0]
CHNMS_data_nMDS_water <- cbind(CHNMS_data_nMDS_water[,c(1,2)], CHNMS_data_nMDS_species)
values <- c("North", "South")
repetitions <- c(6, 3)
CHNMS_data_nMDS_water$GeoLocation <- rep(values, repetitions)

## nMDS calucation and plotting 
NMDS_CHNMS_COI = metaMDS(CHNMS_data_nMDS_water[,-c(1,2, 119)], distance="jaccard", k = 3, trymax = 200, autotransform = F)

NMDSspecies_COI = data.frame(Species=rownames(NMDS_CHNMS_COI$species),
                             x=as.numeric(NMDS_CHNMS_COI$species[,1]),
                             y=as.numeric(NMDS_CHNMS_COI$species[,2]))

NMDSstaz_COI = data.frame(Site=CHNMS_data_nMDS_water$Site, Source = CHNMS_data_nMDS_water$Source, Group = CHNMS_data_nMDS_water$GeoLocation,
                          x=as.numeric(NMDS_CHNMS_COI$points[,1]), 
                          y=as.numeric(NMDS_CHNMS_COI$points[,2]))

NMDSstaz_COI$Group= as.character(NMDSstaz_COI$Group)
NMDSstaz_COI$NrSpecies = rowSums(CHNMS_data_nMDS_water[,-c(1,2,119)] > 0)

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls_Group <- ddply(NMDSstaz_COI, .(Group), find_hull)
hulls_Group$Group = factor(hulls_Group$Group,
                           levels = sort(unique(hulls_Group$Group)))

NMDS_CHNMS_plot_COI_water =  ggplot(data = NMDSstaz_COI, aes(x=x, y=y, xmin = -1.1, xmax=1.1, ymax = 0.5, ymin = -0.5)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_polygon(data = hulls_Group, aes(x=x, y=y, group=Group, colour = Group), 
               size = 1, linetype = 2, fill = NA, alpha = 0.8) +
  scale_color_manual(values = c("goldenrod1", "forestgreen")) + 
  geom_point(data = NMDSstaz_COI, aes(x=x, y=y, fill = Source, size = NrSpecies), shape = 21) + 
  scale_fill_manual(values = c("deepskyblue2")) +
  geom_text_repel(data = NMDSstaz_COI, aes(x=x, y=y, label=Site), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 0.7, y = 0.4, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_CHNMS_COI$stress,2)))

Xs = CHNMS_data_nMDS_water[,c(1,2)]
Xs$Location <- NMDSstaz_COI$Group
Ys = CHNMS_data_nMDS_water[,-c(1:2,119)]
permanova = adonis2(Ys ~ Xs$Site, permutations = 9999, method = "jaccard")
permanova_env = adonis2(Ys ~ Xs$Location, permutations = 9999, method = "jaccard")

### Indicator Species Analysis
IndicSp_water <- CHNMS_data_nMDS_water[,c(1,3:118)]
IndicSp_water$GeoLocation <- CHNMS_data_nMDS_water$GeoLocation

abund = IndicSp_water[,2:117]

inv = multipatt(abund, IndicSp_water$GeoLocation, control = how(nperm=9999))
summary(inv)

## COI only BAM##
## dataset preparation
CHNMSSp_NMDS_BAM <- CHNMSSp_NMDS_COI[which(CHNMSSp_NMDS_COI$Source == 'BAM'),]
CHNMS_data_nMDS_BAM <- dcast(CHNMSSp_NMDS_BAM, Site + Source ~ Species, value.var = "Reads", fun.aggregate = mean)

CHNMS_data_nMDS_species <- CHNMS_data_nMDS_BAM[,-c(1,2)]
CHNMS_data_nMDS_species = CHNMS_data_nMDS_species[, colSums(CHNMS_data_nMDS_species) > 0]
CHNMS_data_nMDS_BAM <- cbind(CHNMS_data_nMDS_BAM[,c(1,2)], CHNMS_data_nMDS_species)
values <- c("North", "South")
repetitions <- c(3, 4)
CHNMS_data_nMDS_BAM$GeoLocation <- rep(values, repetitions)

## nMDS calucation and plotting 
NMDS_CHNMS_COI = metaMDS(CHNMS_data_nMDS_BAM[,-c(1,2, 148)], distance="jaccard", k = 3, trymax = 200, autotransform = F)

NMDSspecies_COI = data.frame(Species=rownames(NMDS_CHNMS_COI$species),
                             x=as.numeric(NMDS_CHNMS_COI$species[,1]),
                             y=as.numeric(NMDS_CHNMS_COI$species[,2]))

NMDSstaz_COI = data.frame(Site=CHNMS_data_nMDS_BAM$Site, Source = CHNMS_data_nMDS_BAM$Source, Group = CHNMS_data_nMDS_BAM$GeoLocation,
                          x=as.numeric(NMDS_CHNMS_COI$points[,1]), 
                          y=as.numeric(NMDS_CHNMS_COI$points[,2]))

NMDSstaz_COI$Group= as.character(NMDSstaz_COI$Group)
NMDSstaz_COI$NrSpecies = rowSums(CHNMS_data_nMDS_BAM[,-c(1,2,148)] > 0)

find_hull <- function(df) df[chull(df$x, df$y), ]
hulls_Group <- ddply(NMDSstaz_COI, .(Group), find_hull)
hulls_Group$Group = factor(hulls_Group$Group,
                           levels = sort(unique(hulls_Group$Group)))


NMDS_CHNMS_plot_COI_BAM =  ggplot(data = NMDSstaz_COI, aes(x=x, y=y, xmin = -0.6, xmax = 0.6, ymin = -0.4, ymax = 0.4)) +
  geom_hline(yintercept = 0, color = "darkgrey", size=0.5) + 
  geom_vline(xintercept = 0, color = "darkgrey", size=0.5) + 
  xlab("NMDS1") + ylab("NMDS2") +
  geom_polygon(data = hulls_Group, aes(x=x, y=y, group=Group, colour = Group), 
               size = 1, linetype = 2, fill = NA, alpha = 0.8) +
  scale_color_manual(values = c("goldenrod1", "forestgreen")) + 
  geom_point(data = NMDSstaz_COI, aes(x=x, y=y, fill = Source, size = NrSpecies), shape = 21) + 
  scale_fill_manual(values = c("darkorange2")) +
  geom_text_repel(data = NMDSstaz_COI, aes(x=x, y=y, label=Site), size = 3.5, force = 10,
                  fontface = "bold", inherit.aes = F, max.overlaps = 20) + 
  theme(legend.position = "right") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  annotate("text", x = 0.2, y = 0.4, size = 4, fontface = "bold",
           label = paste0("Stress: ",round(NMDS_CHNMS_COI$stress,2)))

Xs = CHNMS_data_nMDS_BAM[,c(1,2)]
Xs$Location <- NMDSstaz_COI$Group
Ys = CHNMS_data_nMDS_BAM[,-c(1:2,148)]
permanova = adonis2(Ys ~ Xs$Site, permutations = 9999, method = "jaccard")
permanova_env = adonis2(Ys ~ Xs$Location, permutations = 9999, method = "jaccard")

### Indicator Species Analysis
IndicSp_BAM <- CHNMS_data_nMDS_BAM[,c(1,3:147)]
IndicSp_BAM$GeoLocation <- CHNMS_data_nMDS_BAM$GeoLocation

abund = IndicSp_BAM[,2:146]

inv = multipatt(abund, IndicSp_BAM$GeoLocation, control = how(nperm=9999))
summary(inv)