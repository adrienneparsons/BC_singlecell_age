# Figure 4: Comparison analysis of multiple datasets using CellChat
###########################################################
# Esther Sauras Colon, 2024-06-11

# Load the required libraries
library(CellChat)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)

## MAJOR CELL TYPES
# Load the cellchat objects previously generated
# Create an object list
# DO A SEPARATE ANALYSIS FOR TNBC AND ER

#--------------------------------------------------------------------------------------------------------------------------------
# TNBC
# Change directory!!
cellchatTNBC.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/cellchatTN.young_computeCommunProb_psize_TRUE_celltype_major.Rdata")
cellchatTNBC.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/cellchatTN.old_computeCommunProb_psize_TRUE_celltype_major.Rdata")
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# ER
# Change directory!!
cellchatER.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/cellchatER.young_computeCommunProb_psize_TRUE_celltype_major.Rdata")
cellchatER.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/cellchatER.old_computeCommunProb_psize_TRUE_celltype_major.Rdata")
object.list <- list(ER.young = cellchatER.young, ER.old = cellchatER.old)
#--------------------------------------------------------------------------------------------------------------------------------

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Interaction strength among different cell populations for both groups (young and old)
weight.max <- getMaxWeight(object.list, attribute = c("idents","weight"))
par(mfrow = c(1,2), xpd=TRUE)
for (i in 1:length(object.list)) {
  netVisual_circle(object.list[[i]]@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Strength of interactions - ", names(object.list)[i]))
}

# Circle plot with the differential interaction strength among different cell populations
par(mfrow = c(1,1), xpd=TRUE)
netVisual_diffInteraction(cellchat, weight.scale = T, measure = "weight")

## MINOR CELL TYPES
# Load the cellchat objects previously generated

# Create an object list
# DO A SEPARATE ANALYSIS FOR TNBC AND ER

#--------------------------------------------------------------------------------------------------------------------------------
# TNBC
# Change directory!!
cellchatTNBC.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatTN.young_computeCommunProb_psize_TRUE.Rdata")
cellchatTNBC.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatTN.old_computeCommunProb_psize_TRUE.Rdata")
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# ER
# Change directory!!
cellchatER.young <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatER.young_computeCommunProb_psize_TRUE.Rdata")
cellchatER.old <- readRDS("C:/Users/18458505Q/Desktop/HTVC Esther/BWH/RData/cellchat_files/cellchatER.old_computeCommunProb_psize_TRUE.Rdata")
object.list <- list(ER.young = cellchatER.young, ER.old = cellchatER.old)
#--------------------------------------------------------------------------------------------------------------------------------

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# liftCellchat() if necessary 
lapply(object.list, function(x) length(levels(x@idents)))
# Remove the following code comment for ER (ER young has no B cells Naive, so we need to lift up)
# object.list[[1]] <- liftCellChat(object.list[[1]], levels(object.list[[2]]@idents))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap with the differential interaction strength among different cell populations
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg2

# Compare the major sources and targets in 2D space
# Compute the network centrality scores
object.list <- lapply(object.list, function(x) {netAnalysis_computeCentrality(x, slot.name = "netP")})
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

num.link <- sapply(object.list, function(x) {rowSums(x@net$count) + colSums(x@net$count)-diag(x@net$count)})
weight.MinMax <- c(min(num.link), max(num.link)) # control the dot size in the different datasets
gg <- list()
for (i in 1:length(object.list)) {
  gg[[i]] <- netAnalysis_signalingRole_scatter(object.list[[i]], title = names(object.list)[i], weight.MinMax = weight.MinMax, label.size=4) + theme(aspect.ratio = 1) + scale_y_continuous(limits = c(0,0.25)) + scale_x_continuous(limits = c(0,0.15))
}
# The x and y axes limits change
# TNBC: scale_y_continuous(limits = c(0,0.25)) + scale_x_continuous(limits = c(0,0.15))
# ER: scale_y_continuous(limits = c(0,0.125)) + scale_x_continuous(limits = c(0,0.15))
patchwork::wrap_plots(plots = gg)
