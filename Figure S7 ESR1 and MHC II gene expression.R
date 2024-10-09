# Figure S7: ESR1 and MHC II gene analysis
###########################################################
# Adrienne Parsons, 2024-05-30

# Prerequisites -----------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(Seurat)
library(gdata)
library(readxl)
library(msigdbr)
library(dplyr)
library(gage)
library(fgsea)
library(ggplot2)
library(ggforce)

# Set up ------------------------------------------------------------------------------------------
setwd("<YOUR DATA PATH>")
#setwd("~/DropboxMGB/Projects/Single-cell_breast_cancer/AnalysisAdrienne") # for Peter
rm(list=ls())

# Function to split character string
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------

# Load supplemental table, downloaded from here: https://www.nature.com/articles/s41588-021-00911-1#Sec39
sup.tib <- read_excel("./wu,swarbrick2021 - Supplement.xlsx", skip = 3)

# Load expression and metadata, downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
folders <- list.files("../Data_Swarbrick", full.names = T)
folders <- folders[!grepl("xlsx", folders)]
data.tib <- tibble(Sample = cutf(folders, d = "/", f = 3),
                   Folder = folders)

# Generate empty list 
seu.ls <- vector(mode = "list", length = nrow(data.tib))

# Load Seurat objects
for (n in 1:nrow(data.tib)) {
  #n = 1
  print(data.tib$Sample[n])
  
  # Define folder
  folder <- data.tib$Folder[n]
  
  # Load data
  counts.mat <- Read10X(folder, gene.column = 1)
  metadata.df <- read.csv(paste0(folder, "/metadata.csv"), row.names = 1)
  
  # Create Seurat object
  seu <- CreateSeuratObject(counts = counts.mat, 
                            meta.data = metadata.df,
                            min.cells = 10, min.features = 100, names.field = 2, project = data.tib$Sample[n])
  
  # Add it to the list
  seu.ls[[n]] <- seu
}

# Remove objects from the loop to prevent any confusion
rm(counts.mat, metadata.df, seu)

# Merge
seu.all <- merge(x = seu.ls[[1]], y = seu.ls[-1])

# Generate a matrix of the donors' ages
age_matrix <- matrix(pull(sup.tib, 3))
rownames(age_matrix) <- gsub("-", "", paste0("CID", sup.tib$`Case ID`))
rownames(age_matrix) <- gsub("CID4290", "CID4290A", rownames(age_matrix))
rownames(age_matrix) <- gsub("CID4530", "CID4530N", rownames(age_matrix))
colnames(age_matrix) <- "age"

# Get the ER+ donors
seu.ER <- subset(seu.all, subset = subtype == "ER+")

#Annotate donor ages and age groups
seu.ER$Age <- NA
seu.ER$Age_group <- NA

for(donor in unique(seu.ER$orig.ident)){
  seu.ER$Age[seu.ER$orig.ident == donor] <- age_matrix[which(rownames(age_matrix) == donor), 1]
}

seu.ER$Age_group[seu.ER$Age > 55] <- "Older"
seu.ER$Age_group[seu.ER$Age <= 55] <- "Young"
seu.ER <- NormalizeData(seu.ER)


# Set the Idents to the celltype to plot by cell type
Idents(seu.ER) <- "celltype_minor"

setwd("<YOUR RESULTS PATH>")

# Get the ER+ data for ESR1 and format the cell types in a new data frame
esr1_data <- GetAssayData(seu.ER, slot = "data") %>% .["ESR1", ]

er_data <- data.frame("Celltype" = seu.ER$celltype_minor, "Age" = seu.ER$Age_group, "Data" = esr1_data)

er_data$Celltype <- gsub("CAFs MSC ", "", er_data$Celltype)
er_data$Celltype <- gsub("CAFs ", "", er_data$Celltype)
er_data$Celltype <- gsub(" SC", "", er_data$Celltype)
er_data$Celltype <- gsub("Endothelial ", "", er_data$Celltype)
er_data$Celltype <- gsub(" cells", "", er_data$Celltype)
er_data$Celltype <- gsub("-cells", "", er_data$Celltype)
er_data$Celltype <- gsub("_", " ", er_data$Celltype)

# Reorder the cell types to match other plots
er_data$Celltype <- factor(er_data$Celltype, levels=c("B Memory", "B Naive", "Plasmablasts",
                                                          "iCAF-like", "myCAF-like", "Cancer Basal",
                                                          "Cancer Her2", "Cancer LumA", "Cancer LumB", 
                                                          "Cancer Cycling",
                                                          "Luminal Progenitors", "Mature Luminal", "Myoepithelial",
                                                          "ACKR1", "CXCL12",
                                                          "Lymphatic LYVE1", "RGS5",
                                                          "DCs", "Macrophage", "Monocyte", "Cycling Myeloid",
                                                          "PVL Differentiated",
                                                          "PVL Immature", "Cycling PVL", "NK", "NKT", "T CD4+", "T CD8+",
                                                          "Cycling T"))

# Plot the ESR1 Sina plots
p <- ggplot(er_data, aes(x = Celltype, y = Data, colour = Age))+ geom_sina(scale = "width", size = .1, stroke = .05)+
  geom_violin(draw_quantiles = 0.5, scale = "width", fill = NA, linewidth = 0.2)+
  ylab("ESR1 Expression")+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6), aspect.ratio = .09)

ggsave("ESR1_sinaPlot.pdf", p, width = 6.5, height = 1.5)

# Repeat for MHC and TNBC

# First, subset the data to just the TNBC donors
seu.TNBC <- subset(seu.all, subset = subtype == "TNBC")

#Annotate donor ages and age groups
seu.TNBC$Age <- NA
seu.TNBC$Age_group <- NA

for(donor in unique(seu.TNBC$orig.ident)){
  seu.TNBC$Age[seu.TNBC$orig.ident == donor] <- age_matrix[which(rownames(age_matrix) == donor), 1]
}

seu.TNBC$Age_group[seu.TNBC$Age > 55] <- "Older"
seu.TNBC$Age_group[seu.TNBC$Age <= 55] <- "Young"
seu.TNBC <- NormalizeData(seu.TNBC)

# Set the idents to celltype_minor to group sina plots by cell type
Idents(seu.TNBC) <- "celltype_minor"


setwd("<YOUR RESULTS PATH>")

# For the MHC genes of interest that are present in the dataset:
for(feature in c("HLA-DRB1", "HLA-DRB5", "HLA-DRA", "HLA-DQA1", "HLA-DQB1", "HLA-DOB", 
                 "HLA-DOA", "HLA-DMA", "HLA-DMB",
                 "HLA-DPA1", "HLA-DPB1")){
  
  # Get the data for that gene
  genedata <- GetAssayData(seu.TNBC, slot = "data") %>% .[feature, ]
  
  # Make a new data frame with cell type, age group, and the gene's data
  TNBC_data <- data.frame("Celltype" = seu.TNBC$celltype_minor, "Age" = seu.TNBC$Age_group, "Data" = genedata)
  
  # Format the cell type names and re-order to match other plots
  TNBC_data$Celltype <- gsub("CAFs MSC ", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub("CAFs ", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub(" SC", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub("Endothelial ", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub(" cells", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub("-cells", "", TNBC_data$Celltype)
  TNBC_data$Celltype <- gsub("_", " ", TNBC_data$Celltype)
  
  TNBC_data$Celltype <- factor(TNBC_data$Celltype, levels=c("B Memory", "B Naive", "Plasmablasts",
                                                  "iCAF-like", "myCAF-like", "Cancer Basal",
                                                  "Cancer Her2", "Cancer LumA", "Cancer LumB", 
                                                  "Cancer Cycling",
                                                  "Luminal Progenitors", "Mature Luminal", "Myoepithelial",
                                                  "ACKR1", "CXCL12",
                                                  "Lymphatic LYVE1", "RGS5",
                                                  "DCs", "Macrophage", "Monocyte", "Cycling Myeloid",
                                                  "PVL Differentiated",
                                                  "PVL Immature", "Cycling PVL", "NK", "NKT", "T CD4+", "T CD8+",
                                                  "Cycling T"))
  
  # Generate the sina plot for that gene
  p <- ggplot(TNBC_data, aes(x = Celltype, y = Data, colour = Age))+ geom_sina(scale = "width", size = .1, stroke = .05)+
    geom_violin(draw_quantiles = 0.5, scale = "width", fill = NA, linewidth = 0.2)+
    ylab(paste0(feature, " Expression"))+
    theme_bw()+
    theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 1, size = 6), aspect.ratio = .09)
  
  ggsave(paste0(feature, "_MHCII_sinaPlot.pdf"), p, width = 6.5, height = 1.5)
}

# Statistics to determine which comparisons are significant
Idents(seu.ER) <- "Age_group"
Idents(seu.TNBC) <- "Age_group"

# Make a dataframe to populate with the significant MHC data
MHC <- data.frame("Celltype" = NA, "p_val" = NA, "avg_log2FC" = NA, "pct.1" = NA, "pct.2" = NA, "p_val_adj" = NA, "gene" = NA)
MHC <- MHC[-1,]

# For each cell type, determine if the MHC markers are significantly DE between young and old
for(type in unique(seu.TNBC$celltype_minor)){
  print(type)
  obj <- subset(seu.TNBC, subset = celltype_minor == type)
  DE <- FindMarkers(obj, features = c("HLA-DRB1", "HLA-DRB5", "HLA-DRA", "HLA-DQA1", "HLA-DQB1", "HLA-DOB", 
                                            "HLA-DOA", "HLA-DMA", "HLA-DMB",
                                            "HLA-DPA1", "HLA-DPB1"), ident.1 = "Older", ident.2 = "Young")
  # Add the cell type and statistics to the data frame
  DE$Celltype <- type
  DE$gene <- rownames(DE)
  
  DE <- DE[, c(7, 1:6)]
  
  MHC <- rbind(MHC, DE)
}

# Subset to just the statistically significant comparisons
MHC2 <- MHC[MHC$p_val_adj < 0.05,]


# Statistics to determine which ER+ cell types have significant differences in ESR1 expression between 
# young and old
ESR <- data.frame("Celltype" = NA, "p_val" = NA, "avg_log2FC" = NA, "pct.1" = NA, "pct.2" = NA, "p_val_adj" = NA, "gene" = NA)
ESR <- ESR[-1,]

# For each cell type, subsset the data to that cell type and do DE gene expression analysis fir ESR1
for(type in unique(seu.ER$celltype_minor)){
  print(type)
  obj <- subset(seu.ER, subset = celltype_minor == type)
  
  if("Older" %in% Idents(obj) & "Young" %in% Idents(obj)){
    DE <- FindMarkers(obj, features = c("ESR1"), ident.1 = "Older", ident.2 = "Young")
    
    # If ESR1 is expressed and statistics are available, add those statistics to the data frame
    if(nrow(DE) > 0){
      DE$Celltype <- type
      DE$gene <- rownames(DE)
      
      DE <- DE[, c(7, 1:6)]
      
      ESR <- rbind(ESR, DE)
  }
  
  
  }
  
}

# Subset to just the significant cell types
ESR2 <- ESR[ESR$p_val_adj < 0.05,]

# Write the significant cell types to a csv
write.csv(file = "MHC_DEgenes.csv", MHC2, row.names = F)
write.csv(file = "ESR1_DEgenes.csv", ESR2, row.names = F)

