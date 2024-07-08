## All CellChat data: ER/ER.old/ER.young/TN/TN.old/TN.young

# Load the required libraries
library(CellChat)
library(patchwork)

# Load the Seurat objects
seu_def_er <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defER_subset.Rdata")
seu_def_er_old <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defERold_subset.Rdata")
seu_def_er_young <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defERyoung_subset.Rdata")
seu_def_tnbc <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defTNBC_subset.Rdata")
seu_def_tnbc_old <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defTNBCold_subset.Rdata")
seu_def_tnbc_young <- readRDS("C:/Users/esthe/Desktop/PhD/BWH/RData/seurat_CCA_integrated_defTNBCyoung_subset.Rdata")

# Create CellChat objects
# Define group.by 
# celltype_minor or celltype_major
cellgroup = "celltype_minor"
cellchatER <- createCellChat(object = seu_def_er, group.by = cellgroup, assay = "RNA")
cellchatER.old <- createCellChat(object = seu_def_er_old, group.by = cellgroup, assay = "RNA")
cellchatER.young <- createCellChat(object = seu_def_er_young, group.by = cellgroup, assay = "RNA")
cellchatTN <- createCellChat(object = seu_def_tnbc, group.by = cellgroup, assay = "RNA")
cellchatTN.old <- createCellChat(object = seu_def_tnbc_old, group.by = cellgroup, assay = "RNA")
cellchatTN.young <- createCellChat(object = seu_def_tnbc_young, group.by = cellgroup, assay = "RNA")

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 

# Set the used database in the object
cellchatER@DB <- CellChatDB
cellchatER.old@DB <- CellChatDB
cellchatER.young@DB <- CellChatDB
cellchatTN@DB <- CellChatDB
cellchatTN.old@DB <- CellChatDB
cellchatTN.young@DB <- CellChatDB

# Create a list with all the cellchat objects
cc.list <- list(ER = cellchatER, ER.old = cellchatER.old,
                ER.young = cellchatER.young, TN = cellchatTN,
                TN.young = cellchatTN.young, TN.old = cellchatTN.old)

# Subset the expression data of signaling genes for saving computation cost
cc.list <- lapply(cc.list, function(x) subsetData(x))

future::plan("multisession", workers = 4) # do parallel

# Part II: Inference of cell-cell communication network
# Compute the communication probability and infer cellular communication network
cc.list <- lapply(cc.list, function(x) identifyOverExpressedGenes(x))
cc.list <- lapply(cc.list, function(x) identifyOverExpressedInteractions(x))
# For celltype_minor:
# The number of highly variable ligand-receptor pairs used for signaling inference is 2498 
# The number of highly variable ligand-receptor pairs used for signaling inference is 2243
# The number of highly variable ligand-receptor pairs used for signaling inference is 2451
# The number of highly variable ligand-receptor pairs used for signaling inference is 2459
# The number of highly variable ligand-receptor pairs used for signaling inference is 2362
# The number of highly variable ligand-receptor pairs used for signaling inference is 2275 

# For celltype_major:
# The number of highly variable ligand-receptor pairs used for signaling inference is 2448
# The number of highly variable ligand-receptor pairs used for signaling inference is 2210
# The number of highly variable ligand-receptor pairs used for signaling inference is 2384
# The number of highly variable ligand-receptor pairs used for signaling inference is 2421
# The number of highly variable ligand-receptor pairs used for signaling inference is 2324
# The number of highly variable ligand-receptor pairs used for signaling inference is 2211

# Project gene expression data onto PPI (Optional)
# cellchat <- projectData(cellchat, PPI.human)

# Define population size 
# TRUE for unsorted single-cell transcriptomics or FALSE for sorted-enriched single cells
psize = TRUE
options(future.globals.maxSize = 600 * 1024^2) # default: 500 * 1024^2 = 500 MiB
cc.list <- lapply(cc.list, function(x) computeCommunProb(x, population.size = psize))

lapply(names(cc.list), function(x) saveRDS(cc.list[[x]], paste0("C:/Users/esthe/Desktop/PhD/BWH/cellchat",x,"_computeCommunProb_psize_",psize,"_",cellgroup,".Rdata")))
