# Cellchat object creation from Seurat objects
###########################################################
# Esther Sauras Colon & Adrienne Parsons, 2025-02-28

# Prerequisite packages
library(dplyr)
library(Seurat)
library(patchwork)
library(ggplot2)
library(xlsx)
library(stringr)
library(cowplot)
library(CellChat)
library(patchwork)
library(VennDiagram)
library(RColorBrewer)
library(ComplexHeatmap)

# options(Seurat.object.assay.version = "v5")
# Load the SWARBRICK dataset, available for download through the Broad Single-Cell portal at https://singlecell.broadinstitute.org/single_cell/study/SCP1039
# Also available at GEO Accession # GSE176078

data_path <- "<YOUR DATA PATH>"
options(Seurat.object.assay.version = "v5")

# Function to recursively read scRNAseq data using read10x and create seurat objects including metadata
datasheet <- "metadata.csv"
read_scRNAseq_data <- function(directory_path) {
  # Get a list of the folders in the directory
  file_list <- list.files(path = directory_path, recursive = FALSE, full.names = TRUE)
  # Create a list to store Seurat objects
  seurat_list <- list()
  # Loop through each file
  for (file_path in file_list) {
    # File with the associated metadata
    metadatafiles <- read.csv(file.path(file_path, datasheet))
    # Read the 10x Genomics data using Read10X
    seurat_read <- Read10X(data.dir = file_path, gene.column=1)
    # Create the seurat object
    seurat_obj <- CreateSeuratObject(counts = seurat_read, project = "swb", 
                                     meta.data = metadatafiles,
                                     min.cells = 3, min.features = 200)
    # Add the Seurat object to the list
    seurat_list[[file_path]] <- list(seurat_read, seurat_obj)
  }
  # Return the list of Seurat objects
  return(seurat_list)
}

# Specify the path to the 'Data_Swarbrick' folder
data_folder_path <- paste0(data_path,"/Data_Swarbrick")
# Call the function to read scRNAseq data recursively
scRNAseq_data <- read_scRNAseq_data(data_folder_path)

matrix_all <- sapply(scRNAseq_data,function(x) x[1])
seurat_object_all <- sapply(scRNAseq_data,function(x) x[2])

# Merge all the seurat objects in one
seu_all <- merge(x = seurat_object_all[[1]], y = seurat_object_all[-1])
Layers(seu_all) # Same as Layers(seu_all[["RNA]])

# Select ER+ & TNBC (21 patients = 80753 samples)
seu_er_tnbc <- subset(seu_all, subset = subtype == "ER+" | subtype == "TNBC")
# Change ER+ to ER
seu_er_tnbc$subtype <- replace(seu_er_tnbc$subtype,seu_er_tnbc$subtype=="ER+","ER")

# Subset ER
seu_er <- subset(seu_er_tnbc, subset = subtype == "ER")

# Subset TNBC
seu_tnbc <- subset(seu_er_tnbc, subset = subtype == "TNBC")

# QC (nFeature & percent.mt)
summary(seu_er_tnbc$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max.
# 201     870    1421    1822    2420    9683
summary(seu_er$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201     795    1532    1835    2464    9683 
summary(seu_tnbc$nFeature_RNA)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 201.0   905.8  1332.0  1810.1  2380.0  8312.0 
df <- data.frame(seu_er_tnbc$nFeature_RNA,seu_er_tnbc$celltype_minor)
df1 <- data.frame(seu_er$nFeature_RNA,seu_er$celltype_minor)
df2 <- data.frame(seu_tnbc$nFeature_RNA,seu_tnbc$celltype_minor)
ggplot(df, aes(x=seu_er_tnbc.celltype_minor,
               y=seu_er_tnbc.nFeature_RNA,
               fill=seu_er_tnbc.celltype_minor)) + 
  geom_boxplot() + labs(title="Plot of nFeature_RNA per cell type",x="Cell type",y="nFeature_RNA",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))
ggplot(df1, aes(x=seu_er.celltype_minor,
               y=seu_er.nFeature_RNA,
               fill=seu_er.celltype_minor)) + 
  geom_boxplot() + labs(title="Plot of nFeature_RNA per cell type ER",x="Cell type",y="nFeature_RNA",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))
ggplot(df2, aes(x=seu_tnbc.celltype_minor,
                y=seu_tnbc.nFeature_RNA,
                fill=seu_tnbc.celltype_minor)) + 
  geom_boxplot() + labs(title="Plot of nFeature_RNA per cell type TNBC",x="Cell type",y="nFeature_RNA",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))

summary(seu_er_tnbc$percent.mito)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   3.285   5.458   6.324   8.527  19.997
summary(seu_er$percent.mito)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   4.047   6.574   7.227   9.624  19.989
summary(seu_tnbc$percent.mito)
# Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
# 0.000   2.756   4.620   5.511   7.207  19.997

dff <- data.frame(seu_er_tnbc$percent.mito,seu_er_tnbc$celltype_minor)
dff1 <- data.frame(seu_er$percent.mito,seu_er$celltype_minor)
dff2 <- data.frame(seu_tnbc$percent.mito,seu_tnbc$celltype_minor)
ggplot(dff, aes(x=seu_er_tnbc.celltype_minor,
               y=seu_er_tnbc.percent.mito,
               fill=seu_er_tnbc.celltype_minor)) + 
  geom_violin() + labs(title="Plot of percent.mito per cell type",x="Cell type",y="percent.mito",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))
ggplot(dff1, aes(x=seu_er.celltype_minor,
                y=seu_er.percent.mito,
                fill=seu_er.celltype_minor)) + 
  geom_violin() + labs(title="Plot of percent.mito per cell type ER",x="Cell type",y="percent.mito",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))
ggplot(dff2, aes(x=seu_tnbc.celltype_minor,
                 y=seu_tnbc.percent.mito,
                 fill=seu_tnbc.celltype_minor)) + 
  geom_violin() + labs(title="Plot of percent.mito per cell type TNBC",x="Cell type",y="percent.mito",fill="Cell type") +
  theme_classic() + 
  theme(axis.text.x = element_text(angle = 90))
seu_def <- subset(seu_er_tnbc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 20)
# seu_def (21 patients = 77873 samples)

# Perform analysis without integration
# Run standard analysis workflow
seu_def <- NormalizeData(seu_def)
seu_def <- FindVariableFeatures(seu_def)
seu_def <- ScaleData(seu_def)
seu_def <- RunPCA(seu_def)

seu_def <- FindNeighbors(seu_def, dims = 1:30, reduction = "pca")
seu_def <- FindClusters(seu_def, resolution = 2, cluster.name = "unintegrated_clusters")

seu_def <- RunUMAP(seu_def, dims = 1:30, reduction = "pca", reduction.name = "umap.unintegrated")
DimPlot(seu_def, reduction = "umap.unintegrated", group.by = c("orig.ident", "seurat_clusters"))

# Perform integration
seu_def <- IntegrateLayers(object = seu_def, method = CCAIntegration, 
                           orig.reduction = "pca", 
                           new.reduction = "integrated.cca", verbose = FALSE)

# Re-join layers after integration
seu_def[["RNA"]] <- JoinLayers(seu_def[["RNA"]])

seu_def <- FindNeighbors(seu_def, reduction = "integrated.cca", dims = 1:30)
seu_def <- FindClusters(seu_def, resolution = 1)

seu_def <- RunUMAP(seu_def, dims = 1:30, reduction = "integrated.cca")

# Visualization
DimPlot(seu_def, reduction = "umap", group.by = c("orig.ident"))
DimPlot(seu_def, reduction = "umap", group.by = c("celltype_major"))
DimPlot(seu_def, reduction = "umap", group.by = c("celltype_minor"))
DimPlot(seu_def, reduction = "umap", group.by = c("celltype_subset"))

# Add patients' age to the object
agedb <- read.xlsx2(paste0(data_path,"/wu,swarbrick2021 - Supplement.xlsx"), sheetIndex = 1, startRow = 4)
agedb <- agedb[agedb$Subtype.by.IHC == "ER+"|agedb$Subtype.by.IHC == "TNBC",]
agedb <- agedb[,c(1,3)]
agedb$Case.ID <- str_sub(agedb$Case.ID,1,4)
colnames(agedb) <- c("ID","age")

seu_def[["ID"]] <- str_sub(seu_def$orig.ident,4,7)
seu_def[["age"]] <- merge(seu_def@meta.data,agedb,by="ID")
seu_def[["ID"]] <- NULL

# Save seurat objects
saveRDS(seu_def,
        paste0(data_path,"/seurat_CCA_integrated_def_subset.Rdata"))

# Subset ER+ & TNBC in two objects
seu_def_er <- subset(seu_def, subset = subtype == "ER")
seu_def_tnbc <- subset(seu_def, subset = subtype == "TNBC")

# Subset old & young within each subtype
summary(as.numeric(seu_def_er$age)) # median 55
summary(as.numeric(seu_def_tnbc$age)) # median 54

groupageER <- c()
for (i in 1:length(seu_def_er$age)){
  if(seu_def_er$age[i] <= 55){
    groupageER <- c(groupageER,"young")
  } else {
    groupageER <- c(groupageER, "old")
  }
}
seu_def_er$groupage <- groupageER

groupageTN <- c()
for (i in 1:length(seu_def_tnbc$age)){
  if(seu_def_tnbc$age[i] <= 55){
    groupageTN <- c(groupageTN,"young")
  } else {
    groupageTN <- c(groupageTN, "old")
  }
}
seu_def_tnbc$groupage <- groupageTN

saveRDS(seu_def_er,
        paste0(data_path,"/seurat_CCA_integrated_defER_subset.Rdata"))
saveRDS(seu_def_tnbc,
        paste0(data_path,"/seurat_CCA_integrated_defTNBC_subset.Rdata"))

seu_def_tnbc_old <- subset(seu_def_tnbc, subset = groupage == "old")
seu_def_tnbc_young <- subset(seu_def_tnbc, subset = groupage == "young")
saveRDS(seu_def_tnbc_old,
        paste0(data_path,"/seurat_CCA_integrated_defTNBCold_subset.Rdata"))
saveRDS(seu_def_tnbc_young,
        paste0(data_path,"/seurat_CCA_integrated_defTNBCyoung_subset.Rdata"))

seu_def_er_old <- subset(seu_def_er, subset = groupage == "old")
seu_def_er_young <- subset(seu_def_er, subset = groupage == "young")
saveRDS(seu_def_er_old,
        paste0(data_path,"/seurat_CCA_integrated_defERold_subset.Rdata"))
saveRDS(seu_def_er_young,
        paste0(data_path,"/seurat_CCA_integrated_defERyoung_subset.Rdata"))

## All CellChat data: ER/ER.old/ER.young/TN/TN.old/TN.young

# Load the Seurat objects
seu_def_er <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defER_subset.Rdata"))
seu_def_er_old <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defERold_subset.Rdata"))
seu_def_er_young <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defERyoung_subset.Rdata"))
seu_def_tnbc <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defTNBC_subset.Rdata"))
seu_def_tnbc_old <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defTNBCold_subset.Rdata"))
seu_def_tnbc_young <- readRDS(paste0(data_path,"/seurat_CCA_integrated_defTNBCyoung_subset.Rdata"))

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

lapply(names(cc.list), function(x) saveRDS(cc.list[[x]], paste0(data_path,"/",x,"_computeCommunProb_psize_",psize,"_",cellgroup,".Rdata")))

# ___________________________________________________________________________________________
# Creating single-patient CellChat objects

# options(Seurat.object.assay.version = "v5")
# Load the SWARBRICK dataset, available for download through the Broad Single-Cell portal at https://singlecell.broadinstitute.org/single_cell/study/SCP1039
# Also available at GEO Accession # GSE176078
data_path <- "<YOUR DATA PATH>"
options(Seurat.object.assay.version = "v5")

# Function to recursively read scRNAseq data using read10x and create seurat objects including metadata
datasheet <- "metadata.csv"
read_scRNAseq_data <- function(directory_path) {
  # Get a list of the folders in the directory
  file_list <- list.files(path = directory_path, recursive = FALSE, full.names = TRUE)
  # Create a list to store Seurat objects
  seurat_list <- list()
  # Loop through each file
  for (file_path in file_list) {
    # File with the associated metadata
    metadatafiles <- read.csv(file.path(file_path, datasheet))
    # Read the 10x Genomics data using Read10X
    seurat_read <- Read10X(data.dir = file_path, gene.column=1)
    # Create the seurat object
    seurat_obj <- CreateSeuratObject(counts = seurat_read, project = "swb", 
                                     meta.data = metadatafiles,
                                     min.cells = 3, min.features = 200)
    # Add the Seurat object to the list
    seurat_list[[file_path]] <- list(seurat_read, seurat_obj)
  }
  # Return the list of Seurat objects
  return(seurat_list)
}

# Specify the path to the 'Data_Swarbrick' folder
data_folder_path <- data_path
# Call the function to read scRNAseq data recursively
scRNAseq_data <- seurat_list

matrix_all <- sapply(scRNAseq_data,function(x) x[1])
seurat_object_all <- sapply(scRNAseq_data,function(x) x[2])

# Merge all the seurat objects in one
seu_all <- merge(x = seurat_object_all[[1]], y = seurat_object_all[-1])
Layers(seu_all) # Same as Layers(seu_all[["RNA]])

# Select ER+ & TNBC (21 patients = 80753 samples)
seu_er_tnbc <- subset(seu_all, subset = subtype == "ER+" | subtype == "TNBC")
# Change ER+ to ER
seu_er_tnbc$subtype <- replace(seu_er_tnbc$subtype,seu_er_tnbc$subtype=="ER+","ER")

# QC selection
seu_def <- subset(seu_er_tnbc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 20)

# Add patients' age to the object
agedb <- read.xlsx2(paste0(directory_path,"/wu,swarbrick2021 - Supplement.xlsx"), sheetIndex = 1, startRow = 4)
agedb <- agedb[agedb$Subtype.by.IHC == "ER+"|agedb$Subtype.by.IHC == "TNBC",]
agedb <- agedb[,c(1,3)]
agedb$Case.ID <- str_sub(agedb$Case.ID,1,4)
colnames(agedb) <- c("ID","age")

seu_def[["ID"]] <- str_sub(seu_def$orig.ident,4,7)
seu_def[["age"]] <- merge(seu_def@meta.data,agedb,by="ID")
seu_def[["ID"]] <- NULL

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 

# scRNAseq preprocessing
seu_def <- NormalizeData(seu_def)
seu_def <- FindVariableFeatures(seu_def)
seu_def <- ScaleData(seu_def)
seu_def <- RunPCA(seu_def)
seu_def <- FindNeighbors(seu_def, dims = 1:30, reduction = "pca")
seu_def <- FindClusters(seu_def, resolution = 2, cluster.name = "unintegrated_clusters")

# Running in parallel
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 600 * 1024^2)

# For each patient:
for(ident in unique(seu_def$orig.ident)){
  print(ident)
  
  # Subset to that patient
  seu <- subset(seu_def, subset = orig.ident == ident)
  
  # For the celltype_major and celltype_minor annotations:
  for(cellgroup in c("celltype_major", "celltype_minor")){
    print(cellgroup)
    
    # Make a cellchat object of the single patient for that annotation group
    cellchat <- createCellChat(object = seu, group.by = cellgroup, assay = "RNA")
    cellchat@DB <- CellChatDB
    
    # CellChat object preprocessing
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, population.size = T)
    
    # Save
    saveRDS(cellchat, paste0(ident,"_computeCommunProb_psize_T_" ,cellgroup,".Rdata"))
    
  }
}# options(Seurat.object.assay.version = "v5")
# Load the SWARBRICK dataset, available for download through the Broad Single-Cell portal at https://singlecell.broadinstitute.org/single_cell/study/SCP1039
# Also available at GEO Accession # GSE176078

data_path <- "<YOUR DATA PATH>"
options(Seurat.object.assay.version = "v5")

# Function to recursively read scRNAseq data using read10x and create seurat objects including metadata
datasheet <- "metadata.csv"
read_scRNAseq_data <- function(directory_path) {
  # Get a list of the folders in the directory
  file_list <- list.files(path = directory_path, recursive = FALSE, full.names = TRUE)
  # Create a list to store Seurat objects
  seurat_list <- list()
  # Loop through each file
  for (file_path in file_list) {
    # File with the associated metadata
    metadatafiles <- read.csv(file.path(file_path, datasheet))
    # Read the 10x Genomics data using Read10X
    seurat_read <- Read10X(data.dir = file_path, gene.column=1)
    # Create the seurat object
    seurat_obj <- CreateSeuratObject(counts = seurat_read, project = "swb", 
                                     meta.data = metadatafiles,
                                     min.cells = 3, min.features = 200)
    # Add the Seurat object to the list
    seurat_list[[file_path]] <- list(seurat_read, seurat_obj)
  }
  # Return the list of Seurat objects
  return(seurat_list)
}

# Specify the path to the 'Data_Swarbrick' folder
data_folder_path <- data_path
# Call the function to read scRNAseq data recursively
scRNAseq_data <- seurat_list

matrix_all <- sapply(scRNAseq_data,function(x) x[1])
seurat_object_all <- sapply(scRNAseq_data,function(x) x[2])

# Merge all the seurat objects in one
seu_all <- merge(x = seurat_object_all[[1]], y = seurat_object_all[-1])
Layers(seu_all) # Same as Layers(seu_all[["RNA]])

# Select ER+ & TNBC (21 patients = 80753 samples)
seu_er_tnbc <- subset(seu_all, subset = subtype == "ER+" | subtype == "TNBC")
# Change ER+ to ER
seu_er_tnbc$subtype <- replace(seu_er_tnbc$subtype,seu_er_tnbc$subtype=="ER+","ER")

# QC selection
seu_def <- subset(seu_er_tnbc, subset = nFeature_RNA > 200 & nFeature_RNA < 5000 & percent.mito < 20)

# Add patients' age to the object
agedb <- read.xlsx2(paste0(directory_path,"/wu,swarbrick2021 - Supplement.xlsx"), sheetIndex = 1, startRow = 4)
agedb <- agedb[agedb$Subtype.by.IHC == "ER+"|agedb$Subtype.by.IHC == "TNBC",]
agedb <- agedb[,c(1,3)]
agedb$Case.ID <- str_sub(agedb$Case.ID,1,4)
colnames(agedb) <- c("ID","age")

seu_def[["ID"]] <- str_sub(seu_def$orig.ident,4,7)
seu_def[["age"]] <- merge(seu_def@meta.data,agedb,by="ID")
seu_def[["ID"]] <- NULL

# Set the ligand-receptor interaction database
CellChatDB <- CellChatDB.human 

# scRNAseq preprocessing
seu_def <- NormalizeData(seu_def)
seu_def <- FindVariableFeatures(seu_def)
seu_def <- ScaleData(seu_def)
seu_def <- RunPCA(seu_def)
seu_def <- FindNeighbors(seu_def, dims = 1:30, reduction = "pca")
seu_def <- FindClusters(seu_def, resolution = 2, cluster.name = "unintegrated_clusters")

# Running in parallel
future::plan("multisession", workers = 4)
options(future.globals.maxSize = 600 * 1024^2)

# For each patient:
for(ident in unique(seu_def$orig.ident)){
  print(ident)
  
  # Subset to that patient
  seu <- subset(seu_def, subset = orig.ident == ident)
  
  # For the celltype_major and celltype_minor annotations:
  for(cellgroup in c("celltype_major", "celltype_minor")){
    print(cellgroup)
    
    # Make a cellchat object of the single patient for that annotation group
    cellchat <- createCellChat(object = seu, group.by = cellgroup, assay = "RNA")
    cellchat@DB <- CellChatDB
    
    # CellChat object preprocessing
    cellchat <- subsetData(cellchat)
    cellchat <- identifyOverExpressedGenes(cellchat)
    cellchat <- identifyOverExpressedInteractions(cellchat)
    cellchat <- computeCommunProb(cellchat, population.size = T)
    
    # Save
    saveRDS(cellchat, paste0(ident,"_computeCommunProb_psize_T_" ,cellgroup,".Rdata"))
    
  }
}







