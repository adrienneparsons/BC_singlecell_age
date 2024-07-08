## Code for rankNet - glm - bubbleplots ER

# Load the required libraries
library(CellChat)
library(patchwork)
library(VennDiagram)
library(RColorBrewer)
library(ComplexHeatmap)
library(xlsx)
library(tidyr)

# Load the cellchat objects saved in Single-cell_breast_cancer > Github > BC_singlecell_age > Data
# Change directory!!
cellchatER.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatER.young_computeCommunProb_psize_TRUE.Rdata")
cellchatER.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatER.old_computeCommunProb_psize_TRUE.Rdata")

object.list <- list(ER.young = cellchatER.young, ER.old = cellchatER.old)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# liftCellchat() necessary for ER
lapply(object.list, function(x) length(levels(x@idents)))
object.list[[1]] <- liftCellChat(object.list[[1]], levels(object.list[[2]]@idents))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

## RankNet ER
## Selected cell types: iCAF, myCAF, lumA, CD4, CD8, macrophages, PVL cycling, ACKR1

## iCAF
# iCAF-iCAF
icaf_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-iCAF"))
icaf_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_icaf_data <- icaf_icaf_data$signaling.contribution

# iCAF-myCAF
icaf_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-myCAF"))
icaf_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_mycaf_data <- icaf_mycaf_data$signaling.contribution

# iCAF-lumA
icaf_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-Cancer lumA"))
icaf_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_lumA_data <- icaf_lumA_data$signaling.contribution

# iCAF-T cells CD4
icaf_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-CD4"))
icaf_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_CD4_data <- icaf_CD4_data$signaling.contribution

# iCAF-T cells CD8
icaf_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-CD8"))
icaf_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_CD8_data <- icaf_CD8_data$signaling.contribution

# iCAF-macrophages
icaf_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-Macrophages"))
icaf_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_macrophages_data <- icaf_macrophages_data$signaling.contribution

# iCAF-PVL cycling
icaf_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-PVL cycling"))
icaf_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_PVLc_data <- icaf_PVLc_data$signaling.contribution

# iCAF-Endothelial ACKR1
icaf_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-ACKR1"))
icaf_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_ACKR1_data <- icaf_ACKR1_data$signaling.contribution

## myCAF
# myCAF-iCAF
mycaf_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-iCAF"))
mycaf_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_icaf_data <- mycaf_icaf_data$signaling.contribution

# myCAF-myCAF
mycaf_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-myCAF"))
mycaf_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_mycaf_data <- mycaf_mycaf_data$signaling.contribution

# myCAF-lumA
mycaf_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-Cancer lumA"))
mycaf_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_lumA_data <- mycaf_lumA_data$signaling.contribution

# myCAF-T cells CD4
mycaf_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-CD4"))
mycaf_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_CD4_data <- mycaf_CD4_data$signaling.contribution

# myCAF-T cells CD8
mycaf_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-CD8"))
mycaf_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_CD8_data <- mycaf_CD8_data$signaling.contribution

# myCAF-macrophages
mycaf_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-Macrophages"))
mycaf_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_macrophages_data <- mycaf_macrophages_data$signaling.contribution

# myCAF-PVL cycling
mycaf_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-PVL cycling"))
mycaf_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_PVLc_data <- mycaf_PVLc_data$signaling.contribution

# myCAF-Endothelial ACKR1
mycaf_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-ACKR1"))
mycaf_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_ACKR1_data <- mycaf_ACKR1_data$signaling.contribution

## lumA
# lumA-iCAF
lumA_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-iCAF"))
lumA_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_icaf_data <- lumA_icaf_data$signaling.contribution

# lumA-myCAF
lumA_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-myCAF"))
lumA_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_mycaf_data <- lumA_mycaf_data$signaling.contribution

# lumA-lumA
lumA_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-Cancer lumA"))
lumA_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_lumA_data <- lumA_lumA_data$signaling.contribution

# lumA-T cells CD4
lumA_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-CD4"))
lumA_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_CD4_data <- lumA_CD4_data$signaling.contribution

# lumA-T cells CD8
lumA_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-CD8"))
lumA_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_CD8_data <- lumA_CD8_data$signaling.contribution

# lumA-macrophages
lumA_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-Macrophages"))
lumA_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_macrophages_data <- lumA_macrophages_data$signaling.contribution

# lumA-PVL cycling
lumA_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-PVL cycling"))
lumA_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_PVLc_data <- lumA_PVLc_data$signaling.contribution

# lumA-Endothelial ACKR1
lumA_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer lumA-ACKR1"))
lumA_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
lumA_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(8), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
lumA_ACKR1_data <- lumA_ACKR1_data$signaling.contribution

## T cells CD4
# T cells CD4-iCAF
CD4_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-iCAF"))
CD4_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_icaf_data <- CD4_icaf_data$signaling.contribution

# T cells CD4-myCAF
CD4_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-myCAF"))
CD4_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_mycaf_data <- CD4_mycaf_data$signaling.contribution

# T cells CD4-lumA
CD4_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-Cancer lumA"))
CD4_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_lumA_data <- CD4_lumA_data$signaling.contribution

# T cells CD4-T cells CD4
CD4_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-CD4"))
CD4_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_CD4_data <- CD4_CD4_data$signaling.contribution

# T cells CD4-T cells CD8
CD4_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-CD8"))
CD4_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_CD8_data <- CD4_CD8_data$signaling.contribution

# T cells CD4-macrophages
CD4_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-Macrophages"))
CD4_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_macrophages_data <- CD4_macrophages_data$signaling.contribution

# T cells CD4-PVL cycling -- no inferred communications
# CD4_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-PVL cycling"))
# CD4_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
# CD4_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
# CD4_PVLc_data <- CD4_PVLc_data$signaling.contribution

# T cells CD4-Endothelial ACKR1
CD4_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-ACKR1"))
CD4_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_ACKR1_data <- CD4_ACKR1_data$signaling.contribution

## T cells CD8
# T cells CD8-iCAF
CD8_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-iCAF"))
CD8_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_icaf_data <- CD8_icaf_data$signaling.contribution

# T cells CD8-myCAF
CD8_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-myCAF"))
CD8_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_mycaf_data <- CD8_mycaf_data$signaling.contribution

# T cells CD8-lumA
CD8_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-Cancer lumA"))
CD8_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_lumA_data <- CD8_lumA_data$signaling.contribution

# T cells CD8-T cells CD4
CD8_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-CD4"))
CD8_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_CD4_data <- CD8_CD4_data$signaling.contribution

# T cells CD8-T cells CD8
CD8_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-CD8"))
CD8_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_CD8_data <- CD8_CD8_data$signaling.contribution

# T cells CD8-macrophages
CD8_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-Macrophages"))
CD8_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_macrophages_data <- CD8_macrophages_data$signaling.contribution

# T cells CD8-PVL cycling
CD8_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-PVL cycling"))
CD8_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_PVLc_data <- CD8_PVLc_data$signaling.contribution

# T cells CD8-Endothelial ACKR1
CD8_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-ACKR1"))
CD8_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_ACKR1_data <- CD8_ACKR1_data$signaling.contribution

## Macrophages
# Macrophages-iCAF
macrophages_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-iCAF"))
macrophages_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_icaf_data <- macrophages_icaf_data$signaling.contribution

# Macrophages-myCAF
macrophages_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-myCAF"))
macrophages_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_mycaf_data <- macrophages_mycaf_data$signaling.contribution

# Macrophages-lumA
macrophages_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-Cancer lumA"))
macrophages_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_lumA_data <- macrophages_lumA_data$signaling.contribution

# Macrophages-T cells CD4
macrophages_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-CD4"))
macrophages_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_CD4_data <- macrophages_CD4_data$signaling.contribution

# Macrophages-T cells CD8
macrophages_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-CD8"))
macrophages_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_CD8_data <- macrophages_CD8_data$signaling.contribution

# Macrophages-macrophages
macrophages_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-Macrophages"))
macrophages_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_macrophages_data <- macrophages_macrophages_data$signaling.contribution

# Macrophages-PVL cycling
macrophages_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-PVL cycling"))
macrophages_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_PVLc_data <- macrophages_PVLc_data$signaling.contribution

# Macrophages-Endothelial ACKR1
macrophages_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-ACKR1"))
macrophages_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_ACKR1_data <- macrophages_ACKR1_data$signaling.contribution

## PVL cycling
# PVL cycling-iCAF
PVLc_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-iCAF"))
PVLc_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_icaf_data <- PVLc_icaf_data$signaling.contribution

# PVL cycling-myCAF
PVLc_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-myCAF"))
PVLc_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_mycaf_data <- PVLc_mycaf_data$signaling.contribution

# PVL cycling-lumA
PVLc_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-Cancer lumA"))
PVLc_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_lumA_data <- PVLc_lumA_data$signaling.contribution

# PVL cycling-T cells CD4
PVLc_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-CD4"))
PVLc_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_CD4_data <- PVLc_CD4_data$signaling.contribution

# PVL cycling-T cells CD8
PVLc_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-CD8"))
PVLc_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_CD8_data <- PVLc_CD8_data$signaling.contribution

# PVL cycling-macrophages
PVLc_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-Macrophages"))
PVLc_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_macrophages_data <- PVLc_macrophages_data$signaling.contribution

# PVL cycling-PVL cycling
PVLc_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-PVL cycling"))
PVLc_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_PVLc_data <- PVLc_PVLc_data$signaling.contribution

# PVL cycling-Endothelial ACKR1
PVLc_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for PVL cycling-ACKR1"))
PVLc_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
PVLc_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(10), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
PVLc_ACKR1_data <- PVLc_ACKR1_data$signaling.contribution

## Endothelial ACKR1
# Endothelial ACKR1-iCAF
ACKR1_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-iCAF"))
ACKR1_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_icaf_data <- ACKR1_icaf_data$signaling.contribution

# Endothelial ACKR1-myCAF
ACKR1_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-myCAF"))
ACKR1_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_mycaf_data <- ACKR1_mycaf_data$signaling.contribution

# Endothelial ACKR1-lumA
ACKR1_lumA <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(8), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-Cancer lumA"))
ACKR1_lumA + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_lumA_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(8), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_lumA_data <- ACKR1_lumA_data$signaling.contribution

# Endothelial ACKR1-T cells CD4
ACKR1_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-CD4"))
ACKR1_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_CD4_data <- ACKR1_CD4_data$signaling.contribution

# Endothelial ACKR1-T cells CD8
ACKR1_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-CD8"))
ACKR1_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_CD8_data <- ACKR1_CD8_data$signaling.contribution

# Endothelial ACKR1-macrophages
ACKR1_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-Macrophages"))
ACKR1_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_macrophages_data <- ACKR1_macrophages_data$signaling.contribution

# Endothelial ACKR1-PVL cycling
ACKR1_PVLc <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(10), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-PVL cycling"))
ACKR1_PVLc + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_PVLc_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(10), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_PVLc_data <- ACKR1_PVLc_data$signaling.contribution

# Endothelial ACKR1-Endothelial ACKR1
ACKR1_ACKR1 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(14), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for ACKR1-ACKR1"))
ACKR1_ACKR1 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
ACKR1_ACKR1_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(14), targets.use = c(14), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
ACKR1_ACKR1_data <- ACKR1_ACKR1_data$signaling.contribution


## Signaling pathways ER
pathwaysER <- c(unique(icaf_icaf_data$name),unique(icaf_mycaf_data$name),
                unique(icaf_lumA_data$name),unique(icaf_CD4_data$name),
                unique(icaf_CD8_data$name),unique(icaf_macrophages_data$name),
                unique(icaf_PVLc_data$name),unique(icaf_ACKR1_data$name),
                unique(mycaf_icaf_data$name),unique(mycaf_mycaf_data$name),
                unique(mycaf_lumA_data$name),unique(mycaf_CD4_data$name),
                unique(mycaf_CD8_data$name),unique(mycaf_macrophages_data$name),
                unique(mycaf_PVLc_data$name),unique(mycaf_ACKR1_data$name),
                unique(lumA_icaf_data$name),unique(lumA_mycaf_data$name),
                unique(lumA_lumA_data$name),unique(lumA_CD4_data$name),
                unique(lumA_CD8_data$name),unique(lumA_macrophages_data$name),
                unique(lumA_PVLc_data$name),unique(lumA_ACKR1_data$name),
                unique(CD4_icaf_data$name),unique(CD4_mycaf_data$name),
                unique(CD4_lumA_data$name),unique(CD4_CD4_data$name),
                unique(CD4_CD8_data$name),unique(CD4_macrophages_data$name),
                unique(CD4_ACKR1_data$name),unique(CD8_icaf_data$name),
                unique(CD8_mycaf_data$name),unique(CD8_lumA_data$name),
                unique(CD8_CD4_data$name),unique(CD8_CD8_data$name),
                unique(CD8_macrophages_data$name),unique(CD8_PVLc_data$name),
                unique(CD8_ACKR1_data$name),unique(macrophages_icaf_data$name),
                unique(macrophages_mycaf_data$name),unique(macrophages_lumA_data$name),
                unique(macrophages_CD4_data$name),unique(macrophages_CD8_data$name),
                unique(macrophages_macrophages_data$name),unique(macrophages_PVLc_data$name),
                unique(macrophages_ACKR1_data$name),unique(PVLc_icaf_data$name),
                unique(PVLc_mycaf_data$name),unique(PVLc_lumA_data$name),
                unique(PVLc_CD4_data$name),unique(PVLc_CD8_data$name),
                unique(PVLc_macrophages_data$name),unique(PVLc_PVLc_data$name),
                unique(PVLc_ACKR1_data$name),unique(ACKR1_icaf_data$name),
                unique(ACKR1_mycaf_data$name),unique(ACKR1_lumA_data$name),
                unique(ACKR1_CD4_data$name),unique(ACKR1_CD8_data$name),
                unique(ACKR1_macrophages_data$name),unique(ACKR1_PVLc_data$name),
                unique(ACKR1_ACKR1_data$name))

unique(pathwaysER)
length(unique(pathwaysER)) #102
table_pathwaysER <- table(pathwaysER)[order(table(pathwaysER), decreasing=T)]
names_pathwaysER <- names(table_pathwaysER[table_pathwaysER != 0])


## Bubble plots - we have to create them to obtain the first table, but the definitive bubble plots from fig 6 are coming later
# bubble plot
# sources.use & targets.use are all the previously specified cells (change as needed)
netVisual_bubble(cellchat, signaling = unique(pathwaysER), sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29), comparison = c(1, 2), angle.x = 45)
# bubble plot associated data (return.data = TRUE)
bubble_data_all <- netVisual_bubble(cellchat, signaling = unique(pathwaysER), sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29), comparison = c(1, 2), angle.x = 45, return.data = TRUE)
# LRpair <- netAnalysis_contribution(object.list[[1]], signaling = unique(pathwaysER), sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29), return.data = TRUE)

# Keep the ones that have p<0.01 (pval == 3)
bubble_data_all <- bubble_data_all$communication[bubble_data_all$communication$pval == 3,]

# Create the row when the probability is 0
bubble_data_all <- bubble_data_all %>% arrange(interaction_name, desc(group.names))

# Create a row when the probability is 0; (also in script Function_for_ratio.R)
fun <- function(data){
  i <- 1
  final_data <- data[1,]
  final_data <- final_data[-1,]
  
  while(i < nrow(data)){
    print(i)
    group.name <- data$group.names[i]
    int.name <- data$interaction_name[i]
    oy <- data$dataset[i]
    
    n <- i+1
    
    if(data$group.names[n] == group.name & data$interaction_name[n] == int.name & data$dataset[n] != oy){
      final_data <- rbind(final_data, data[c(i, n),])
      i <- n+1
    } 
    
    else if((data$group.names[n] == group.name & data$interaction_name[n] == int.name & data$dataset[n] != oy) == F){
      print(group.name)
      print(int.name)
      print(oy)
      
      if(oy == "ER.old"){
        print("adding young")
        final_data <- rbind(final_data, data[i,])
        final_data <- rbind(final_data, data[i,])
        final_data$dataset[nrow(final_data)] <- "ER.young"
        final_data$prob[nrow(final_data)] <- 0
        i <- i + 1
      }
      
      if(oy == "ER.young"){
        print("adding old")
        final_data <- rbind(final_data, data[i,])
        final_data <- rbind(final_data, data[i,])
        final_data$dataset[nrow(final_data)] <- "ER.old"
        final_data$prob[nrow(final_data)] <- 0
        i <- i + 1
      }
    }
    
  }
  
  return(final_data)
  
  
}

bubble_data_all_2 <- fun(bubble_data_all) 

## Logistic regression
bubble_data_all_3 <- bubble_data_all_2 %>% arrange(interaction_name, desc(group.names), dataset)
# Table with all the communication probabilities for each L-R pair in each cell-cell interaction (ER young & old)
write.xlsx(bubble_data_all_3, "Bubbledata_ER.xlsx", sheetName = "Sheet1")

df <- data.frame(pathway_name = character(0),
                 cell_interaction = character(0),
                 ER.young = numeric(0),
                 ER.old = numeric(0))
for (i in 1:nrow(bubble_data_all_3)){
  sp <- bubble_data_all_3$pathway_name[i]
  int <- bubble_data_all_3$group.names[i]
  prob_ER_young <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "ER.young" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  prob_ER_old <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "ER.old" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  df <- rbind(df, data.frame(pathway_name = sp, cell_interaction = int, ER.young = prob_ER_young, ER.old = prob_ER_old))
}
df <- unique(df)

dfyoung <- cbind(df[,1:3],rep(c("Young"),859))
colnames(dfyoung) <- c("pathway_name","cell_interaction","avg_prob","dataset")
dfold <- cbind(df[,c(1:2,4)],rep(c("Old"),859))
colnames(dfold) <- c("pathway_name","cell_interaction","avg_prob","dataset")
df2 <- rbind(dfyoung,dfold)

df3 <- pivot_wider(df2, names_from = pathway_name, values_from = avg_prob)
df3[is.na(df3)] <- 0
# Table including the cell-cell interactions in rows and the signaling pathways in columns
write.xlsx(df3, "Bubbledata_ER_glm.xlsx", sheetName = "Sheet1")

df3$dataset <- factor(df3$dataset, levels=c("Young","Old"), labels = c(0,1))
# summary(glm(dataset ~ ., data = df3[,2:104], family = "binomial"))
colnames(df3) <- make.names(colnames(df3))

# Univariate logistic regression model
res <- lapply(3:104, function(i) summary(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

res_confint <- lapply(3:104, function(i) confint(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

# Statistically significant signaling pathways (p<0.05)
sign_res <- lapply(res, function(x) {
  ifelse(x$coefficients[, "Pr(>|z|)"] < 0.05, TRUE, FALSE)
})

for (i in seq_along(sign_res)) {
  if (any(sign_res[[i]])) {
    cat("Model", i, "significant coefficients:\n")
    print(sign_res[[i]])
    print(res[[i]]$coefficients[,"Estimate"][2])
    cat("\n")
  }
}

# Supplementary table
suppltable <- data.frame(pathway_name = character(0),
                         estimate = numeric(0),
                         pvalue = numeric(0))
for (x in res) {
  suppltable <- as.data.frame(rbind(suppltable,c(rownames(x$coefficients)[2],x$coefficients[2, "Estimate"], x$coefficients[2, "Pr(>|z|)"])))
}
colnames(suppltable) <- c("pathway_name","estimate","pvalue")

ci95 <- data.frame()
for (x in res_confint) {
  ci95 <- as.data.frame(rbind(ci95,c(rownames(x)[2],x[2,])))
}
colnames(ci95) <- c("pathway_name","inf CI 95%","sup CI 95%")

suppltable_values <- merge(suppltable, ci95, by = "pathway_name")
suppltable_values <- suppltable_values[,c(1:2,4:5,3)]

table_pathwaysER <- as.data.frame(table_pathwaysER)
table_pathwaysER$pathwaysER <- make.names(table_pathwaysER$pathwaysER)
suppltable_values_def <- merge(suppltable_values,table_pathwaysER,by.x="pathway_name",by.y="pathwaysER")
suppltable_values_def$included_fig5 <- ifelse(as.numeric(suppltable_values_def$pvalue) < 0.05 & suppltable_values_def$Freq >= 15, TRUE, FALSE)

# Table with the signaling pathways and their results from the logistic regression model (estimate, ci95 and p-value) and the rankNet frequency
write.xlsx(suppltable_values_def,"Bubbledata_ER_glm_res.xlsx",sheetName = "Sheet1")

# Signaling pathways obtained from logistic regression
sp_from_lr <- c("PERIOSTIN","FN1","THBS","CD99","NOTCH","GAP","ADGRE")
# Signaling pathways obtained from logistic regression & rankNet frequency >=15
sp_from_lr_ranknet <- c("FN1","THBS","CD99","GAP","ADGRE")

# bubble plot
# sources.use & targets.use are all the previously specified cells (change as needed)
# sp_from_lr_ranknet includes the signaling pathways obtained from logistic regression & rankNet frequency >=15
netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29), comparison = c(1, 2), angle.x = 45, thresh = 0.01)
# bubble plot associated data (return.data = TRUE)
databubble <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29), comparison = c(1, 2), angle.x = 45, return.data = TRUE)

# Increased signaling in older group
gg1 <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in older group", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg1

# Increased signaling in younger group
gg2 <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,8,10,14,19,28,29), targets.use = c(3,4,8,10,14,19,28,29),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in younger group", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg2
