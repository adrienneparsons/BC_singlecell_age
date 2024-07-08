## Code for rankNet - glm - bubbleplots TNBC

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
cellchatTNBC.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatTN.young_computeCommunProb_psize_TRUE.Rdata")
cellchatTNBC.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatTN.old_computeCommunProb_psize_TRUE.Rdata")

object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

## RankNet TNBC
## Selected cell types: iCAF, myCAF, basal, CD4, CD8, macrophages, monocytes

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

# iCAF-basal
icaf_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-Cancer basal"))
icaf_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_basal_data <- icaf_basal_data$signaling.contribution

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

# iCAF-monocytes
icaf_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for iCAF-Monocytes"))
icaf_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
icaf_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(3), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
icaf_monocytes_data <- icaf_monocytes_data$signaling.contribution

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

# myCAF-basal
mycaf_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-Cancer basal"))
mycaf_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_basal_data <- mycaf_basal_data$signaling.contribution

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

# myCAF-monocytes
mycaf_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for myCAF-Monocytes"))
mycaf_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
mycaf_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(4), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
mycaf_monocytes_data <- mycaf_monocytes_data$signaling.contribution

## basal
# basal-iCAF
basal_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-iCAF"))
basal_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_icaf_data <- basal_icaf_data$signaling.contribution

# basal-myCAF
basal_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-myCAF"))
basal_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_mycaf_data <- basal_mycaf_data$signaling.contribution

# basal-basal
basal_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-Cancer basal"))
basal_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_basal_data <- basal_basal_data$signaling.contribution

# basal-T cells CD4
basal_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-CD4"))
basal_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_CD4_data <- basal_CD4_data$signaling.contribution

# basal-T cells CD8
basal_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-CD8"))
basal_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_CD8_data <- basal_CD8_data$signaling.contribution

# basal-macrophages
basal_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-Macrophages"))
basal_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_macrophages_data <- basal_macrophages_data$signaling.contribution

# basal-monocytes
basal_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Cancer basal-Monocytes"))
basal_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
basal_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(5), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
basal_monocytes_data <- basal_monocytes_data$signaling.contribution

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

# T cells CD4-basal
CD4_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-Cancer basal"))
CD4_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_basal_data <- CD4_basal_data$signaling.contribution

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

# T cells CD4-monocytes
CD4_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD4-Monocytes"))
CD4_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD4_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(28), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD4_monocytes_data <- CD4_monocytes_data$signaling.contribution

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

# T cells CD8-basal
CD8_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-Cancer basal"))
CD8_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_basal_data <- CD8_basal_data$signaling.contribution

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

# T cells CD8-monocytes
CD8_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for CD8-Monocytes"))
CD8_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
CD8_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(29), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
CD8_monocytes_data <- CD8_monocytes_data$signaling.contribution

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

# Macrophages-basal
macrophages_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-Cancer basal"))
macrophages_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_basal_data <- macrophages_basal_data$signaling.contribution

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

# Macrophages-monocytes
macrophages_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Macrophages-Monocytes"))
macrophages_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
macrophages_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(19), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
macrophages_monocytes_data <- macrophages_monocytes_data$signaling.contribution

## Monocytes
# Monocytes-iCAF
monocytes_icaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(3), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-iCAF"))
monocytes_icaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_icaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(3), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_icaf_data <- monocytes_icaf_data$signaling.contribution

# Monocytes-myCAF
monocytes_mycaf <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(4), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-myCAF"))
monocytes_mycaf + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_mycaf_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(4), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_mycaf_data <- monocytes_mycaf_data$signaling.contribution

# Monocytes-basal
monocytes_basal <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(5), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-Cancer basal"))
monocytes_basal + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_basal_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(5), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_basal_data <- monocytes_basal_data$signaling.contribution

# Monocytes-T cells CD4
monocytes_CD4 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(28), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-CD4"))
monocytes_CD4 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_CD4_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(28), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_CD4_data <- monocytes_CD4_data$signaling.contribution

# Monocytes-T cells CD8
monocytes_CD8 <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(29), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-CD8"))
monocytes_CD8 + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_CD8_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(29), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_CD8_data <- monocytes_CD8_data$signaling.contribution

# Monocytes-macrophages
monocytes_macrophages <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(19), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-Macrophages"))
monocytes_macrophages + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_macrophages_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(19), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_macrophages_data <- monocytes_macrophages_data$signaling.contribution

# Monocytes-monocytes
monocytes_monocytes <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(21), do.stat = FALSE, show.raw = FALSE, title = paste0("Signaling pathways for Monocytes-Monocytes"))
monocytes_monocytes + theme(axis.text.y = element_text(size=8), title = element_text(size=10))
monocytes_monocytes_data <- rankNet(cellchat, measure = "weight", mode = "comparison", stacked = F, sources.use = c(21), targets.use = c(21), do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
monocytes_monocytes_data <- monocytes_monocytes_data$signaling.contribution

## Signaling pathways TNBC
pathwaysTNBC <- c(unique(icaf_icaf_data$name),unique(icaf_mycaf_data$name),
                  unique(icaf_basal_data$name),unique(icaf_CD4_data$name),
                  unique(icaf_CD8_data$name),unique(icaf_macrophages_data$name),
                  unique(icaf_monocytes_data$name),unique(mycaf_icaf_data$name),
                  unique(mycaf_mycaf_data$name),unique(mycaf_basal_data$name),
                  unique(mycaf_CD4_data$name),unique(mycaf_CD8_data$name),
                  unique(mycaf_macrophages_data$name),unique(mycaf_monocytes_data$name),
                  unique(basal_icaf_data$name),unique(basal_mycaf_data$name),
                  unique(basal_basal_data$name),unique(basal_CD4_data$name),
                  unique(basal_CD8_data$name),unique(basal_macrophages_data$name),
                  unique(basal_monocytes_data$name),unique(CD4_icaf_data$name),
                  unique(CD4_mycaf_data$name),unique(CD4_basal_data$name),
                  unique(CD4_CD4_data$name),unique(CD4_CD8_data$name),
                  unique(CD4_macrophages_data$name),unique(CD4_monocytes_data$name),
                  unique(CD8_icaf_data$name),unique(CD8_mycaf_data$name),
                  unique(CD8_basal_data$name),unique(CD8_CD4_data$name),
                  unique(CD8_CD8_data$name),unique(CD8_macrophages_data$name),
                  unique(CD8_monocytes_data$name),unique(macrophages_icaf_data$name),
                  unique(macrophages_mycaf_data$name),unique(macrophages_basal_data$name),
                  unique(macrophages_CD4_data$name),unique(macrophages_CD8_data$name),
                  unique(macrophages_macrophages_data$name),unique(macrophages_monocytes_data$name),
                  unique(monocytes_icaf_data$name),unique(monocytes_mycaf_data$name),
                  unique(monocytes_basal_data$name),unique(monocytes_CD4_data$name),
                  unique(monocytes_CD8_data$name),unique(monocytes_macrophages_data$name),
                  unique(monocytes_monocytes_data$name))

unique(pathwaysTNBC)
length(unique(pathwaysTNBC)) #71
table_pathwaysTNBC <- table(pathwaysTNBC)[order(table(pathwaysTNBC), decreasing=T)]
names_pathwaysTNBC <- names(table_pathwaysTNBC[table_pathwaysTNBC != 0])

## Bubble plots - we have to create them to obtain the first table, but the definitive bubble plots from fig 5 are coming later
# bubble plot
# sources.use & targets.use are all the previously specified cells (change as needed)
netVisual_bubble(cellchat, signaling = unique(pathwaysTNBC), sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29), comparison = c(1, 2), angle.x = 45)
# bubble plot associated data (return.data = TRUE)
bubble_data_all <- netVisual_bubble(cellchat, signaling = unique(pathwaysTNBC), sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29), comparison = c(1, 2), angle.x = 45, return.data = TRUE)
# LRpair <- netAnalysis_contribution(object.list[[1]], signaling = unique(pathwaysTNBC), sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29), return.data = TRUE)

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
      
      if(oy == "TNBC.old"){
        print("adding young")
        final_data <- rbind(final_data, data[i,])
        final_data <- rbind(final_data, data[i,])
        final_data$dataset[nrow(final_data)] <- "TNBC.young"
        final_data$prob[nrow(final_data)] <- 0
        i <- i + 1
      }
      
      if(oy == "TNBC.young"){
        print("adding old")
        final_data <- rbind(final_data, data[i,])
        final_data <- rbind(final_data, data[i,])
        final_data$dataset[nrow(final_data)] <- "TNBC.old"
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
# Table with all the communication probabilities for each L-R pair in each cell-cell interaction (TNBC young & old)
write.xlsx(bubble_data_all_3, "Bubbledata_TNBC.xlsx", sheetName = "Sheet1")

df <- data.frame(pathway_name = character(0),
                 cell_interaction = character(0),
                 TNBC.young = numeric(0),
                 TNBC.old = numeric(0))
for (i in 1:nrow(bubble_data_all_3)){
  sp <- bubble_data_all_3$pathway_name[i]
  int <- bubble_data_all_3$group.names[i]
  prob_TNBC_young <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "TNBC.young" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  prob_TNBC_old <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "TNBC.old" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  df <- rbind(df, data.frame(pathway_name = sp, cell_interaction = int, TNBC.young = prob_TNBC_young, TNBC.old = prob_TNBC_old))
}
df <- unique(df)

dfyoung <- cbind(df[,1:3],rep(c("Young"),646))
colnames(dfyoung) <- c("pathway_name","cell_interaction","avg_prob","dataset")
dfold <- cbind(df[,c(1:2,4)],rep(c("Old"),646))
colnames(dfold) <- c("pathway_name","cell_interaction","avg_prob","dataset")
df2 <- rbind(dfyoung,dfold)

df3 <- pivot_wider(df2, names_from = pathway_name, values_from = avg_prob)
df3[is.na(df3)] <- 0
# Table including the cell-cell interactions in rows and the signaling pathways in columns
write.xlsx(df3, "Bubbledata_TNBC_glm.xlsx", sheetName = "Sheet1")

df3$dataset <- factor(df3$dataset, levels=c("Young","Old"), labels = c(0,1))
# summary(glm(dataset ~ ., data = df3[,2:73], family = "binomial"))
colnames(df3) <- make.names(colnames(df3))

# Univariate logistic regression model
res <- lapply(3:73, function(i) summary(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

res_confint <- lapply(3:73, function(i) confint(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

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

table_pathwaysTNBC <- as.data.frame(table_pathwaysTNBC)
table_pathwaysTNBC$pathwaysTNBC <- make.names(table_pathwaysTNBC$pathwaysTNBC)
suppltable_values_def <- merge(suppltable_values,table_pathwaysTNBC,by.x="pathway_name",by.y="pathwaysTNBC")
suppltable_values_def$included_fig5 <- ifelse(as.numeric(suppltable_values_def$pvalue) < 0.05 & suppltable_values_def$Freq >= 15, TRUE, FALSE)

# Table with the signaling pathways and their results from the logistic regression model (estimate, ci95 and p-value) and the rankNet frequency
write.xlsx(suppltable_values_def,"Bubbledata_TNBC_glm_res.xlsx",sheetName = "Sheet1")

# Signaling pathways obtained from logistic regression
sp_from_lr <- c("GALECTIN","CypA","PLAU","FN1","THBS","Cholesterol","APP","ICAM","MHC-I","MPZ","VCAM","ADGRE")
# Signaling pathways obtained from logistic regression & rankNet frequency >=15
sp_from_lr_ranknet <- c("GALECTIN","CypA","PLAU","FN1","THBS","APP","MHC-I","MPZ","ADGRE")

# bubble plot
# sources.use & targets.use are all the previously specified cells (change as needed)
# sp_from_lr_ranknet includes the signaling pathways obtained from logistic regression & rankNet frequency >=15
netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29), comparison = c(1, 2), angle.x = 45, thresh = 0.01)
# bubble plot associated data (return.data = TRUE)
databubble <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29), comparison = c(1, 2), angle.x = 45, return.data = TRUE)

# Increased signaling in older group
gg1 <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29),  comparison = c(1, 2), max.dataset = 2, title.name = "Increased signaling in older group", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg1

# Increased signaling in younger group
gg2 <- netVisual_bubble(cellchat, signaling = sp_from_lr_ranknet, sources.use = c(3,4,5,19,21,28,29), targets.use = c(3,4,5,19,21,28,29),  comparison = c(1, 2), max.dataset = 1, title.name = "Increased signaling in younger group", angle.x = 45, remove.isolate = T, thresh = 0.01)
gg2