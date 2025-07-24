##### Libraries #####
library("dplyr")
library("ggplot2")
library("ggpubr")
library("iNEXT")
library("reshape2") 
library("tidyr")

##### Split data by sampling type#####
## 12S ##
BAM_12S <- CHNMS_data[which(CHNMS_data$Source == 'BAM'),]
DAM_12S <- CHNMS_data[which(CHNMS_data$Source == 'DAM'),]
FRAM_12S <- CHNMS_data[which(CHNMS_data$Source == 'FRAM'),]
water_12S <- CHNMS_data[which(CHNMS_data$Source == 'water'),]

BAM_12S_Allsp <- BAM_12S[,3:54]
BAM_12S_Allsp <- BAM_12S_Allsp[,- which(apply(BAM_12S_Allsp, 2, sum) == 0)]
BAM_12S <- cbind(BAM_12S[,c(1,2)], BAM_12S_Allsp)

DAM_12S_Allsp <- DAM_12S[,3:54]
DAM_12S_Allsp <- DAM_12S_Allsp[,- which(apply(DAM_12S_Allsp, 2, sum) == 0)]
DAM_12S <- cbind(DAM_12S[,c(1,2)], DAM_12S_Allsp)

FRAM_12S_Allsp <- FRAM_12S[,3:54]
FRAM_12S_Allsp <- FRAM_12S_Allsp[,- which(apply(FRAM_12S_Allsp, 2, sum) == 0)]
FRAM_12S <- cbind(FRAM_12S[,c(1,2)], FRAM_12S_Allsp)

water_12S_Allsp <- water_12S[,3:54]
water_12S_Allsp <- water_12S_Allsp[,- which(apply(water_12S_Allsp, 2, sum) == 0)]
water_12S <- cbind(water_12S[,c(1,2)], water_12S_Allsp)

## COI ##
BAM_COI <- CHNMS_COIdata[which(CHNMS_COIdata$Source == 'BAM'),]
DAM_COI <- CHNMS_COIdata[which(CHNMS_COIdata$Source == 'DAM'),]
water_COI <- CHNMS_COIdata[which(CHNMS_COIdata$Source == 'water'),]
Cobble_COI <- CHNMS_COIdata[which(CHNMS_COIdata$Source == 'Cobble'),]

BAM_COI_Allsp <- BAM_COI[,3:457]
BAM_COI_Allsp <- BAM_COI_Allsp[,- which(apply(BAM_COI_Allsp, 2, sum) == 0)]
BAM_COI <- cbind(BAM_COI[,c(1,2)], BAM_COI_Allsp)
AllSp_BAM_COI <- as.data.frame(colnames(BAM_COI)[-c(1,2)])
colnames(AllSp_BAM_COI) <- "scientific_name"
COI_Sp_Phyla <- rbind(COIeDNA[,c(2,6)], COICobble[,c(2,6)])
AllSp_BAM_COI <- left_join(AllSp_BAM_COI, COI_Sp_Phyla)
AllSp_BAM_COI <- AllSp_BAM_COI[!duplicated(AllSp_BAM_COI[,c('scientific_name')]),]
AllSp_BAM_COI$reads <- apply(BAM_COI[,-c(1,2)],2,sum)

DAM_COI_Allsp <- DAM_COI[,3:457]
DAM_COI_Allsp <- DAM_COI_Allsp[,- which(apply(DAM_COI_Allsp, 2, sum) == 0)]
DAM_COI <- cbind(DAM_COI[,c(1,2)], DAM_COI_Allsp)
AllSp_DAM_COI <- as.data.frame(colnames(DAM_COI)[-c(1,2)])
colnames(AllSp_DAM_COI) <- "scientific_name"
AllSp_DAM_COI <- left_join(AllSp_DAM_COI, COI_Sp_Phyla)
AllSp_DAM_COI <- AllSp_DAM_COI[!duplicated(AllSp_DAM_COI[,c('scientific_name')]),]
AllSp_DAM_COI$reads <- apply(DAM_COI[,-c(1,2)],2,sum)

water_COI_Allsp <- water_COI[,3:457]
water_COI_Allsp <- water_COI_Allsp[,- which(apply(water_COI_Allsp, 2, sum) == 0)]
water_COI <- cbind(water_COI[,c(1,2)], water_COI_Allsp)
AllSp_water_COI <- as.data.frame(colnames(water_COI)[-c(1,2)])
colnames(AllSp_water_COI) <- "scientific_name"
AllSp_water_COI <- left_join(AllSp_water_COI, COI_Sp_Phyla)
AllSp_water_COI <- AllSp_water_COI[!duplicated(AllSp_water_COI[,c('scientific_name')]),]
AllSp_water_COI$reads <- apply(water_COI[,-c(1,2)],2,sum)

Cobble_COI_Allsp <- Cobble_COI[,3:457]
Cobble_COI_Allsp <- Cobble_COI_Allsp[,- which(apply(Cobble_COI_Allsp, 2, sum) == 0)]
Cobble_COI <- cbind(Cobble_COI[,c(1,2)], Cobble_COI_Allsp)
AllSp_Cobble_COI <- as.data.frame(colnames(Cobble_COI)[-c(1,2)])
colnames(AllSp_Cobble_COI) <- "scientific_name"
AllSp_Cobble_COI <- left_join(AllSp_Cobble_COI, COI_Sp_Phyla)
AllSp_Cobble_COI <- AllSp_Cobble_COI[!duplicated(AllSp_Cobble_COI[,c('scientific_name')]),]
AllSp_Cobble_COI$reads <- apply(Cobble_COI[,-c(1,2)],2,sum)

## COI all MOTUS##
BAM_COI_allMOTUS <- CHNMS_COIdata_allMOTUS[which(CHNMS_COIdata_allMOTUS$Source == 'BAM'),]
DAM_COI_allMOTUS <- CHNMS_COIdata_allMOTUS[which(CHNMS_COIdata_allMOTUS$Source == 'DAM'),]
water_COI_allMOTUS <- CHNMS_COIdata_allMOTUS[which(CHNMS_COIdata_allMOTUS$Source == 'water'),]
Cobble_COI_allMOTUS <- CHNMS_COIdata_allMOTUS[which(CHNMS_COIdata_allMOTUS$Source == 'Cobble'),]

BAM_COI_Allsp_allMOTUS <- BAM_COI_allMOTUS[,3:4904]
BAM_COI_Allsp_allMOTUS <- BAM_COI_Allsp_allMOTUS[,- which(apply(BAM_COI_Allsp_allMOTUS, 2, sum) == 0)]
BAM_COI_allMOTUS <- cbind(BAM_COI_allMOTUS[,c(1,2)], BAM_COI_Allsp_allMOTUS)
AllSp_BAM_COI_allMOTUS <- as.data.frame(colnames(BAM_COI_allMOTUS)[-c(1,2)])
colnames(AllSp_BAM_COI_allMOTUS) <- "id"
COI_Sp_Phyla_allMOTUS <- rbind(COIeDNA_allMOTUS[,c(1,6)], COICobble_allMOTUS[,c(1,6)])
AllSp_BAM_COI_allMOTUS <- left_join(AllSp_BAM_COI_allMOTUS, COI_Sp_Phyla_allMOTUS)
AllSp_BAM_COI_allMOTUS$reads <- apply(BAM_COI_allMOTUS[,-c(1,2)],2,sum)

DAM_COI_Allsp_allMOTUS <- DAM_COI_allMOTUS[,3:4904]
DAM_COI_Allsp_allMOTUS <- DAM_COI_Allsp_allMOTUS[,- which(apply(DAM_COI_Allsp_allMOTUS, 2, sum) == 0)]
DAM_COI_allMOTUS <- cbind(DAM_COI_allMOTUS[,c(1,2)], DAM_COI_Allsp_allMOTUS)
AllSp_DAM_COI_allMOTUS <- as.data.frame(colnames(DAM_COI_allMOTUS)[-c(1,2)])
colnames(AllSp_DAM_COI_allMOTUS) <- "id"
AllSp_DAM_COI_allMOTUS <- left_join(AllSp_DAM_COI_allMOTUS, COI_Sp_Phyla_allMOTUS)
AllSp_DAM_COI_allMOTUS$reads <- apply(DAM_COI_allMOTUS[,-c(1,2)],2,sum)

water_COI_Allsp_allMOTUS <- water_COI_allMOTUS[,3:4904]
water_COI_Allsp_allMOTUS <- water_COI_Allsp_allMOTUS[,- which(apply(water_COI_Allsp_allMOTUS, 2, sum) == 0)]
water_COI_allMOTUS <- cbind(water_COI_allMOTUS[,c(1,2)], water_COI_Allsp_allMOTUS)
AllSp_water_COI_allMOTUS <- as.data.frame(colnames(water_COI_allMOTUS)[-c(1,2)])
colnames(AllSp_water_COI_allMOTUS) <- "id"
AllSp_water_COI_allMOTUS <- left_join(AllSp_water_COI_allMOTUS, COI_Sp_Phyla_allMOTUS)
AllSp_water_COI_allMOTUS$reads <- apply(water_COI_allMOTUS[,-c(1,2)],2,sum)

Cobble_COI_Allsp_allMOTUS <- Cobble_COI_allMOTUS[,3:4904]
Cobble_COI_Allsp_allMOTUS <- Cobble_COI_Allsp_allMOTUS[,- which(apply(Cobble_COI_Allsp_allMOTUS, 2, sum) == 0)]
Cobble_COI_allMOTUS <- cbind(Cobble_COI_allMOTUS[,c(1,2)], Cobble_COI_Allsp_allMOTUS)
AllSp_Cobble_COI_allMOTUS <- as.data.frame(colnames(Cobble_COI_allMOTUS)[-c(1,2)])
colnames(AllSp_Cobble_COI_allMOTUS) <- "id"
AllSp_Cobble_COI_allMOTUS <- left_join(AllSp_Cobble_COI_allMOTUS, COI_Sp_Phyla_allMOTUS)
AllSp_Cobble_COI_allMOTUS$reads <- apply(Cobble_COI_allMOTUS[,-c(1,2)],2,sum)

##### BoxPlots #####
nr_sp_BAM_12S <- data.frame(BAM_12S$Site, BAM_12S$Source, rowSums(BAM_12S[,-c(1,2)] > 0), apply(BAM_12S[,-c(1,2)],1,sum))
colnames(nr_sp_BAM_12S) <- c("site", "group","nr_sp_12S", "nr_reads_12S")
nr_sp_DAM_12S <- data.frame(DAM_12S$Site, DAM_12S$Source, rowSums(DAM_12S[,-c(1,2)] > 0), apply(DAM_12S[,-c(1,2)],1,sum))
colnames(nr_sp_DAM_12S) <- c("site", "group","nr_sp_12S", "nr_reads_12S")
nr_sp_FRAM_12S <- data.frame(FRAM_12S$Site, FRAM_12S$Source, rowSums(FRAM_12S[,-c(1,2)] > 0), apply(FRAM_12S[,-c(1,2)],1,sum))
colnames(nr_sp_FRAM_12S) <- c("site", "group","nr_sp_12S", "nr_reads_12S")
nr_sp_water_12S <- data.frame(water_12S$Site, water_12S$Source, rowSums(water_12S[,-c(1,2)] > 0), apply(water_12S[,-c(1,2)],1,sum))
colnames(nr_sp_water_12S) <- c("site", "group","nr_sp_12S", "nr_reads_12S")

nr_sp_12S <- rbind(nr_sp_BAM_12S, nr_sp_DAM_12S, nr_sp_FRAM_12S, nr_sp_water_12S)

nr_sp_BAM_COI <- data.frame(BAM_COI$Site, BAM_COI$Source, rowSums(BAM_COI[,-c(1,2)] > 0), apply(BAM_COI[,-c(1,2)],1,sum))
colnames(nr_sp_BAM_COI) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_DAM_COI <- data.frame(DAM_COI$Site, DAM_COI$Source, rowSums(DAM_COI[,-c(1,2)] > 0), apply(DAM_COI[,-c(1,2)],1,sum))
colnames(nr_sp_DAM_COI) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_water_COI <- data.frame(water_COI$Site, water_COI$Source, rowSums(water_COI[,-c(1,2)] > 0), apply(water_COI[,-c(1,2)],1,sum))
colnames(nr_sp_water_COI) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_Cobble_COI <- data.frame(Cobble_COI$Site, Cobble_COI$Source, rowSums(Cobble_COI[,-c(1,2)] > 0), apply(Cobble_COI[,-c(1,2)],1,sum))
colnames(nr_sp_Cobble_COI) <- c("site", "group","nr_sp_COI", "nr_reads_COI")

nr_sp_COI <- rbind(nr_sp_BAM_COI, nr_sp_Cobble_COI, nr_sp_DAM_COI, nr_sp_water_COI)

nr_sp_BAM_COI_allMOTUS <- data.frame(BAM_COI_allMOTUS$Site, BAM_COI_allMOTUS$Source, rowSums(BAM_COI_allMOTUS[,-c(1,2)] > 0), apply(BAM_COI_allMOTUS[,-c(1,2)],1,sum))
colnames(nr_sp_BAM_COI_allMOTUS) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_DAM_COI_allMOTUS <- data.frame(DAM_COI_allMOTUS$Site, DAM_COI_allMOTUS$Source, rowSums(DAM_COI_allMOTUS[,-c(1,2)] > 0), apply(DAM_COI_allMOTUS[,-c(1,2)],1,sum))
colnames(nr_sp_DAM_COI_allMOTUS) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_water_COI_allMOTUS <- data.frame(water_COI_allMOTUS$Site, water_COI_allMOTUS$Source, rowSums(water_COI_allMOTUS[,-c(1,2)] > 0), apply(water_COI_allMOTUS[,-c(1,2)],1,sum))
colnames(nr_sp_water_COI_allMOTUS) <- c("site", "group","nr_sp_COI", "nr_reads_COI")
nr_sp_Cobble_COI_allMOTUS <- data.frame(Cobble_COI_allMOTUS$Site, Cobble_COI_allMOTUS$Source, rowSums(Cobble_COI_allMOTUS[,-c(1,2)] > 0), apply(Cobble_COI_allMOTUS[,-c(1,2)],1,sum))
colnames(nr_sp_Cobble_COI_allMOTUS) <- c("site", "group","nr_sp_COI","nr_reads_COI")

nr_sp_COI_allMOTUS <- rbind(nr_sp_BAM_COI_allMOTUS, nr_sp_Cobble_COI_allMOTUS, nr_sp_DAM_COI_allMOTUS, nr_sp_water_COI_allMOTUS)

color_12S = c("darkorange2", "mediumpurple2", "red3", "deepskyblue2")
color_COI = c("darkorange2", "mediumseagreen",  "mediumpurple2", "deepskyblue2")

boxplot_nr_sp_12S <- ggboxplot(nr_sp_12S[,-1], x = "group", y = "nr_sp_12S", 
                               color = "black", fill = "group", alpha = 0.7, palette = color_12S, size = 0.5, add = "jitter",
                               ylab = "Number of species", xlab = "Methods")

kruskal.test(nr_sp_12S ~ group, data = nr_sp_12S)
pairwise.wilcox.test(nr_sp_12S$nr_sp_12S, nr_sp_12S$group,
                     p.adjust.method = "BH")

wilcox.test(nr_sp_12S$nr_sp_12S[7:14], nr_sp_12S$nr_sp_12S[23:39], alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE, p.adjust.method = "BH")

boxplot_nr_sp_COI <- ggboxplot(nr_sp_COI[,-1], x = "group", y = "nr_sp_COI", 
                               color = "black", fill = "group", alpha = 0.7, palette = color_COI, size = 0.5, add = "jitter",
                               #order = c("Ethanol", "Silica", "Slush"),
                               ylab = "Number of species", xlab = "Methods")

kruskal.test(nr_sp_COI ~ group, data = nr_sp_COI)
pairwise.wilcox.test(nr_sp_COI$nr_sp_COI, nr_sp_COI$group,
                     p.adjust.method = "BH")
wilcox.test(nr_sp_COI$nr_sp_COI[1:44], nr_sp_COI$nr_sp_COI[44:59], alternative = "two.sided", paired = FALSE, exact = FALSE, correct = TRUE, p.adjust.method = "BH")

boxplot_nr_sp_COI_allMOTUS <- ggboxplot(nr_sp_COI_allMOTUS, x = "group", y = "nr_sp_COI", 
                                        color = "black", fill = "group", alpha = 0.7, palette = color_COI, size = 0.5, add = "jitter",
                                        #order = c("Ethanol", "Silica", "Slush"),
                                        ylab = "Number of MOTUs", xlab = "Methods")

kruskal.test(nr_sp_COI ~ group, data = nr_sp_COI_allMOTUS)
pairwise.wilcox.test(nr_sp_COI_allMOTUS$nr_sp_COI, nr_sp_COI_allMOTUS$group,
                     p.adjust.method = "BH")

##### Species accumulation curves #####
# Function to convert a column to binary based on a threshold
binary_converter <- function(column, threshold) {
  binary_data <- ifelse(column >= threshold, 1, 0)
  return(binary_data)
}
# Apply the binary_converter function to all columns
AccCurve_12S <- CHNMS_data[,c(3:54)]
thresholds <- 1
AccCurve_12S <- data.frame(lapply(AccCurve_12S, binary_converter, threshold = thresholds))
AccCurve_12S <- data.frame(CHNMS_data[,c(1,2)], AccCurve_12S)

## Accumulation curves 12S
AccBAM <- AccCurve_12S[which(AccCurve_12S$Source == 'BAM'),-2]
AccBAM <- as.data.frame(t(AccBAM))
colnames(AccBAM) <- AccBAM[1,]
AccBAM <- AccBAM[-1,]
AccBAM$species <- rownames(AccBAM)
rownames(AccBAM) <- NULL

AccDAM <- AccCurve_12S[which(AccCurve_12S$Source == 'DAM'),-2]
AccDAM <- as.data.frame(t(AccDAM))
colnames(AccDAM) <- AccDAM[1,]
AccDAM <- AccDAM[-1,]
AccDAM$species <- rownames(AccDAM)
rownames(AccDAM) <- NULL

AccFRAM <- AccCurve_12S[which(AccCurve_12S$Source == 'FRAM'),-2]
AccFRAM <- as.data.frame(t(AccFRAM))
colnames(AccFRAM) <- AccFRAM[1,]
AccFRAM <- AccFRAM[-1,]
AccFRAM$species <- rownames(AccFRAM)
rownames(AccFRAM) <- NULL

Accwater <- AccCurve_12S[which(AccCurve_12S$Source == 'water'),-2]
Accwater <- as.data.frame(t(Accwater))
colnames(Accwater) <- Accwater[1,]
Accwater <- Accwater[-1,]
Accwater$species <- rownames(Accwater)
rownames(Accwater) <- NULL

BAM <- as.matrix(apply(AccBAM[,-7],2,as.integer))
row.names(BAM) <- AccBAM[,7]
DAM <- as.matrix(apply(AccDAM[,-9],2,as.integer))
row.names(DAM) <- AccDAM[,9]
FRAM <- as.matrix(apply(AccFRAM[,-5],2,as.integer))
row.names(FRAM) <- AccFRAM[,5]
water <- as.matrix(apply(Accwater[,-18],2,as.integer))
row.names(water) <- Accwater[,18]

AccCurve12S = list(BAM = BAM, DAM = DAM, FRAM = FRAM, water = water)

Acc12S_next <- iNEXT(AccCurve12S, datatype="incidence_raw", endpoint=50)
Acc12S <- ggiNEXT(Acc12S_next ) +
  scale_fill_manual(values=color_12S) +
  scale_color_manual(values = color_12S) +
  xlab("Number of sampling sites") + ylab("Species richness") +
  theme_bw()
Acc12S <- Acc12S + theme(legend.position = "bottom")

## Accumulation curves COI
AccCurve_COI <- CHNMS_COIdata[,c(3:457)]
AccCurve_COI <- data.frame(lapply(AccCurve_COI, binary_converter, threshold = thresholds))
AccCurve_COI <- data.frame(CHNMS_COIdata[,c(1,2)], AccCurve_COI)

AccBAM_COI <- AccCurve_COI[which(AccCurve_COI$Source == 'BAM'),-2]
AccBAM_COI <- as.data.frame(t(AccBAM_COI))
colnames(AccBAM_COI) <- AccBAM_COI[1,]
AccBAM_COI <- AccBAM_COI[-1,]
AccBAM_COI$species <- rownames(AccBAM_COI)
rownames(AccBAM_COI) <- NULL

AccDAM_COI <- AccCurve_COI[which(AccCurve_COI$Source == 'DAM'),-2]
AccDAM_COI <- as.data.frame(t(AccDAM_COI))
colnames(AccDAM_COI) <- AccDAM_COI[1,]
AccDAM_COI <- AccDAM_COI[-1,]
AccDAM_COI$species <- rownames(AccDAM_COI)
rownames(AccDAM_COI) <- NULL

Accwater_COI <- AccCurve_COI[which(AccCurve_COI$Source == 'water'),-2]
Accwater_COI <- as.data.frame(t(Accwater_COI))
colnames(Accwater_COI) <- Accwater_COI[1,]
Accwater_COI <- Accwater_COI[-1,]
Accwater_COI$species <- rownames(Accwater_COI)
rownames(Accwater_COI) <- NULL

AccCobble_COI <- AccCurve_COI[which(AccCurve_COI$Source == 'Cobble'),-2]
AccCobble_COI <- as.data.frame(t(AccCobble_COI))
colnames(AccCobble_COI) <- AccCobble_COI[1,]
AccCobble_COI <- AccCobble_COI[-1,]
AccCobble_COI$species <- rownames(AccCobble_COI)
rownames(AccCobble_COI) <- NULL

BAM_COI <- as.matrix(apply(AccBAM_COI[,-9],2,as.integer))
row.names(BAM_COI) <- AccBAM_COI[,9]
DAM_COI <- as.matrix(apply(AccDAM_COI[,-10],2,as.integer))
row.names(DAM_COI) <- AccDAM_COI[,10]
water_COI <- as.matrix(apply(Accwater_COI[,-16],2,as.integer))
row.names(water_COI) <- Accwater_COI[,16]
Cobble_COI <- as.matrix(apply(AccCobble_COI[,-28],2,as.integer))
row.names(Cobble_COI) <- AccCobble_COI[,28]

AccCurveCOI = list(BAM = BAM_COI, DAM = DAM_COI, water = water_COI, Cobble = Cobble_COI)

AccCOI_next <- iNEXT(AccCurveCOI, datatype="incidence_raw", endpoint=200)
AccCOI <- ggiNEXT(AccCOI_next) +
  scale_fill_manual(values=color_COI) +
  scale_color_manual(values = color_COI) +
  xlab("Number of sampling sites") + ylab("Species richness") +
  theme_bw()
AccCOI <- AccCOI + theme(legend.position = "bottom")

## Accumulation curves COI allMOTUS
AccCurve_COI_allMOTUS <- CHNMS_COIdata_allMOTUS[,c(3:4904)]
AccCurve_COI_allMOTUS <- data.frame(lapply(AccCurve_COI_allMOTUS, binary_converter, threshold = thresholds))
AccCurve_COI_allMOTUS <- data.frame(CHNMS_COIdata_allMOTUS[,c(1,2)], AccCurve_COI_allMOTUS)

AccBAM_COI_allMOTUS <- AccCurve_COI_allMOTUS[which(AccCurve_COI_allMOTUS$Source == 'BAM'),-2]
AccBAM_COI_allMOTUS <- as.data.frame(t(AccBAM_COI_allMOTUS))
colnames(AccBAM_COI_allMOTUS) <- AccBAM_COI_allMOTUS[1,]
AccBAM_COI_allMOTUS <- AccBAM_COI_allMOTUS[-1,]
AccBAM_COI_allMOTUS$species <- rownames(AccBAM_COI_allMOTUS)
rownames(AccBAM_COI_allMOTUS) <- NULL

AccDAM_COI_allMOTUS <- AccCurve_COI_allMOTUS[which(AccCurve_COI_allMOTUS$Source == 'DAM'),-2]
AccDAM_COI_allMOTUS <- as.data.frame(t(AccDAM_COI_allMOTUS))
colnames(AccDAM_COI_allMOTUS) <- AccDAM_COI_allMOTUS[1,]
AccDAM_COI_allMOTUS <- AccDAM_COI_allMOTUS[-1,]
AccDAM_COI_allMOTUS$species <- rownames(AccDAM_COI_allMOTUS)
rownames(AccDAM_COI_allMOTUS) <- NULL

Accwater_COI_allMOTUS <- AccCurve_COI_allMOTUS[which(AccCurve_COI_allMOTUS$Source == 'water'),-2]
Accwater_COI_allMOTUS <- as.data.frame(t(Accwater_COI_allMOTUS))
colnames(Accwater_COI_allMOTUS) <- Accwater_COI_allMOTUS[1,]
Accwater_COI_allMOTUS <- Accwater_COI_allMOTUS[-1,]
Accwater_COI_allMOTUS$species <- rownames(Accwater_COI_allMOTUS)
rownames(Accwater_COI_allMOTUS) <- NULL

AccCobble_COI_allMOTUS <- AccCurve_COI_allMOTUS[which(AccCurve_COI_allMOTUS$Source == 'Cobble'),-2]
AccCobble_COI_allMOTUS <- as.data.frame(t(AccCobble_COI_allMOTUS))
colnames(AccCobble_COI_allMOTUS) <- AccCobble_COI_allMOTUS[1,]
AccCobble_COI_allMOTUS <- AccCobble_COI_allMOTUS[-1,]
AccCobble_COI_allMOTUS$species <- rownames(AccCobble_COI_allMOTUS)
rownames(AccCobble_COI_allMOTUS) <- NULL

BAM_COI_allMOTUS <- as.matrix(apply(AccBAM_COI_allMOTUS[,-10],2,as.integer))
row.names(BAM_COI_allMOTUS) <- AccBAM_COI_allMOTUS[,10]
DAM_COI_allMOTUS <- as.matrix(apply(AccDAM_COI_allMOTUS[,-10],2,as.integer))
row.names(DAM_COI_allMOTUS) <- AccDAM_COI_allMOTUS[,10]
water_COI_allMOTUS <- as.matrix(apply(Accwater_COI_allMOTUS[,-18],2,as.integer))
row.names(water_COI_allMOTUS) <- Accwater_COI_allMOTUS[,18]
Cobble_COI_allMOTUS <- as.matrix(apply(AccCobble_COI_allMOTUS[,-29],2,as.integer))
row.names(Cobble_COI_allMOTUS) <- AccCobble_COI_allMOTUS[,29]

AccCurveCOI_allMOTUS = list(BAM = BAM_COI_allMOTUS, DAM = DAM_COI_allMOTUS, water = water_COI_allMOTUS, Cobble = Cobble_COI_allMOTUS)

AccCOI_next_allMOTUS <- iNEXT(AccCurveCOI_allMOTUS, datatype="incidence_raw", endpoint=200)
AccCOI_allMOTUS <- ggiNEXT(AccCOI_next_allMOTUS) +
  scale_fill_manual(values=color_COI) +
  scale_color_manual(values = color_COI) +
  xlab("Number of sampling sites") + ylab("MOTUs richness") +
  theme_bw()
AccCOI_allMOTUS <- AccCOI_allMOTUS + theme(legend.position = "bottom")
