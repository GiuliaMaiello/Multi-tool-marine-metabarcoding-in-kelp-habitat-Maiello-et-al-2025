##### Libraries #####
library("dplyr")
library("VennDiagram")

##### Venn Diagrams ####
## Comparison 12S_DAM/12S_water
SpeciesDAM_12S <- colnames(DAM_12S_Allsp)
Specieswater_12S <- colnames(water_12S_Allsp)
SpeciesBAM_12S <- colnames(BAM_12S_Allsp)
SpeciesFRAM_12S <- colnames(FRAM_12S_Allsp)

Venn_water_DAM <- draw.pairwise.venn(area1 = length(SpeciesDAM_12S), 
                                     area2 = length(Specieswater_12S), 
                                     cross.area = length(intersect(SpeciesDAM_12S, Specieswater_12S)), 
                                     category = c("DAM", "water"),
                                     col=c("mediumpurple2", "deepskyblue2"), lwd=4, font.face=rep("bold",7),
                                     fill = c(alpha("mediumpurple2",0.3), alpha("deepskyblue2",0.3)),
                                     fontfamily = rep("Arial", 3), cat.fontfamily = "Arial")

OnlyWater_12S = setdiff(unique(Specieswater_12S), intersect(Specieswater_12S, SpeciesDAM_12S))
OnlyDAM_12S = setdiff(unique(SpeciesDAM_12S), intersect(Specieswater_12S, SpeciesDAM_12S))
DAM_Water_12S = intersect(Specieswater_12S, SpeciesDAM_12S)

## Comparison BAM/DAM/water 12S
Venn_12S <- draw.triple.venn(area1 = length(SpeciesDAM_12S), 
                             area2 = length(Specieswater_12S), 
                             area3 = length(SpeciesBAM_12S),
                             n12 = length(intersect(SpeciesDAM_12S, Specieswater_12S)),
                             n23 = length(intersect(Specieswater_12S, SpeciesBAM_12S)),
                             n13 = length(intersect(SpeciesDAM_12S, SpeciesBAM_12S)), 
                             n123 = length(intersect(SpeciesDAM_12S, intersect(Specieswater_12S, SpeciesBAM_12S))),
                             euler.d = TRUE, scaled = TRUE,
                             category = c("DAM", "water", "BAM"),
                             col=c("mediumpurple2", "deepskyblue2", "darkorange2"), lwd=3, font.face=rep("bold",7),
                             fill = c(alpha("mediumpurple2",0.3), alpha("deepskyblue2",0.3), alpha("darkorange2",0.3)),
                             fontfamily = rep("Arial", 7), cat.fontfamily = "Arial"
)


## Comparison COI_eDNA/COI_Cobble
eDNA_COI <- CHNMS_COIdata[-which(CHNMS_COIdata$Source == 'Cobble'),]

eDNA_COI_Allsp <- eDNA_COI[,3:457]
eDNA_COI_Allsp <- eDNA_COI_Allsp[,- which(apply(eDNA_COI_Allsp, 2, sum) == 0)]
eDNA_COI <- cbind(eDNA_COI[,c(1,2)], eDNA_COI_Allsp)
AllSp_eDNA_COI <- as.data.frame(colnames(eDNA_COI)[-c(1,2)])
colnames(AllSp_eDNA_COI) <- "scientific_name"
COI_Sp_Phyla <- rbind(COIeDNA[,c(2,6)], COICobble[,c(2,6)])
COI_Sp_Phyla_unique <- COI_Sp_Phyla[!duplicated(COI_Sp_Phyla[,c('scientific_name')]),]
AllSp_eDNA_COI <- left_join(AllSp_eDNA_COI, COI_Sp_Phyla)
AllSp_eDNA_COI <- AllSp_eDNA_COI[!duplicated(AllSp_eDNA_COI[,c('scientific_name')]),]
AllSp_eDNA_COI$reads <- apply(eDNA_COI[,-c(1,2)],2,sum)


SpeciesCobble_COI <- unique(AllSp_Cobble_COI$scientific_name)
SpecieseDNA_COI <- unique(AllSp_eDNA_COI$scientific_name)

Venn_COI_Cob_eDNA <- draw.pairwise.venn(area1 = length(SpecieseDNA_COI), 
                                        area2 = length(SpeciesCobble_COI), 
                                        cross.area = length(intersect(SpecieseDNA_COI, SpeciesCobble_COI)), 
                                        col=c("tomato4", "mediumseagreen"), lwd=3, font.face=rep("bold",7),
                                        fill = c(alpha("tomato4",0.4), alpha("mediumseagreen",0.4)),
                                        fontfamily = rep("Arial", 3), cat.fontfamily = "Arial")

OnlyeDNA_COI = setdiff(unique(SpecieseDNA_COI), intersect(SpecieseDNA_COI, SpeciesCobble_COI))
OnlyCobble_COI = setdiff(unique(SpeciesCobble_COI), intersect(SpecieseDNA_COI, SpeciesCobble_COI))
eDNA_Cobble_COI = intersect(SpecieseDNA_COI, SpeciesCobble_COI)

Abund_eDNACOI <- data.frame(SpecieseDNA_COI, AllSp_eDNA_COI$reads)
Abund_CobCOI <- data.frame(SpeciesCobble_COI, AllSp_Cobble_COI$reads)

Abund_OnlyeDNA <- Abund_eDNACOI[which(Abund_eDNACOI$SpecieseDNA_COI %in% OnlyeDNA_COI),]
Abund_OnlyCobble <- Abund_CobCOI[which(Abund_CobCOI$SpeciesCobble_COI %in% OnlyCobble_COI),]

TotReadseDNA <- sum(Abund_eDNACOI$AllSp_eDNA_COI.reads)
TotReadsCobble <- sum(Abund_CobCOI$AllSp_Cobble_COI.reads)
TotReadsCOI <- TotReadsCobble + TotReadseDNA
PercOnlyeDNA <- sum(Abund_OnlyeDNA$AllSp_eDNA_COI.reads)/TotReadsCOI
PercOnlyCob <- sum(Abund_OnlyCobble$AllSp_Cobble_COI.reads)/TotReadsCOI

## Comparison COI_BAM/COI_Cobble/COI_DAM/COI_water
SpeciesBAM_COI <- unique(AllSp_BAM_COI$scientific_name)
SpeciesCobble_COI <- unique(AllSp_Cobble_COI$scientific_name)
SpeciesDAM_COI <- unique(AllSp_DAM_COI$scientific_name)
Specieswater_COI <- unique(AllSp_water_COI$scientific_name)

VennALL_COI <- draw.quad.venn(area1 = length(unique(SpeciesBAM_COI)), 
                              area2 = length(unique(SpeciesCobble_COI)), 
                              area3 = length(unique(SpeciesDAM_COI)),
                              area4 = length(unique(Specieswater_COI)),
                              n12 = length(intersect(unique(SpeciesBAM_COI), unique(SpeciesCobble_COI))),
                              n13 = length(intersect(unique(SpeciesBAM_COI), unique(SpeciesDAM_COI))), 
                              n14 = length(intersect(unique(SpeciesBAM_COI), unique(Specieswater_COI))),
                              n23 = length(intersect(unique(SpeciesCobble_COI), unique(SpeciesDAM_COI))),
                              n24 = length(intersect(unique(SpeciesCobble_COI), unique(Specieswater_COI))),
                              n34 = length(intersect(unique(SpeciesDAM_COI), unique(Specieswater_COI))),
                              n123 = length(intersect(unique(SpeciesBAM_COI), intersect(unique(SpeciesCobble_COI), unique(SpeciesDAM_COI)))),
                              n124 = length(intersect(unique(SpeciesBAM_COI), intersect(unique(SpeciesCobble_COI), unique(Specieswater_COI)))),
                              n134 = length(intersect(unique(SpeciesBAM_COI), intersect(unique(SpeciesDAM_COI), unique(Specieswater_COI)))),
                              n234 = length(intersect(unique(SpeciesCobble_COI), intersect(unique(SpeciesDAM_COI), unique(Specieswater_COI)))),
                              n1234 = length(intersect(intersect(unique(SpeciesBAM_COI), unique(SpeciesCobble_COI)), intersect(unique(SpeciesDAM_COI), unique(Specieswater_COI)))),
                              euler.d = TRUE, scaled = TRUE,
                              category = c("BAM", "Cobble", "DAM", "water"),
                              col = c("darkorange2", "mediumseagreen",  "mediumpurple2", "deepskyblue2"), lwd=3, font.face=rep("bold",7),
                              cat.fontfamily = "Arial")
