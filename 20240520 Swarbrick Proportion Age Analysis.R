# Prerequisites -----------------------------------------------------------------------------------
library(tidyverse)
library(data.table)
library(Seurat)
library(gdata)
library(ggforce)
library(readxl)
library(msigdbr)
library(dplyr)
library(gage)
library(fgsea)
library(ggplot2)
library(ggpubr)

# Set up ------------------------------------------------------------------------------------------
setwd("/Users/addie/Dropbox (Personal)/Single-cell_breast_cancer/AnalysisAdrienne")
#setwd("~/DropboxMGB/Projects/Single-cell_breast_cancer/AnalysisAdrienne") # for Peter
rm(list=ls())

# Function to split character string
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------

# Load supplemental table
sup.tib <- read_excel("../Data_Swarbrick/wu,swarbrick2021 - Supplement.xlsx", skip = 3)

# Load expression and metadata
folders <- list.files("../Data_Swarbrick", full.names = T)
folders <- folders[!grepl("xlsx", folders)]
data.tib <- tibble(Sample = cutf(folders, d = "/", f = 3),
                   Folder = folders)

# Generate empty list 
seu.ls <- vector(mode = "list", length = nrow(data.tib))

# Load Seurat objects
for (n in 1:nrow(data.tib)) {
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

##############################################################################
# Generate a matrix of the donors' ages
age_matrix <- matrix(pull(sup.tib, 3))
rownames(age_matrix) <- gsub("-", "", paste0("CID", sup.tib$`Case ID`))
rownames(age_matrix) <- gsub("CID4290", "CID4290A", rownames(age_matrix))
rownames(age_matrix) <- gsub("CID4530", "CID4530N", rownames(age_matrix))
colnames(age_matrix) <- "age"

# TNBC Proportion analysis
# New Working Directory
setwd("/Users/addie/desktop/Proportions/TNBC")

# minor to the TNBC donors and add ages to metadata
seu.TNBC <- minor(seu.all, minor = subtype == "TNBC")
seu.TNBC@meta.data$Age <- NA

# Two donors have the same age of 52; make slightly different for indexing
age_matrix[5,1] <- 52.00001

# Add age metadata to Seurat object
for(donor in unique(seu.TNBC$orig.ident)){
  seu.TNBC@meta.data$Age[seu.TNBC@meta.data$orig.ident == donor] <- as.character(age_matrix[rownames(age_matrix) == donor,1])
}

# Find the maximum number of unique minor cell types for each major cell type
for(type in unique(seu.TNBC$celltype_major)){
  print(type)
  type_obj <- minor(seu.TNBC, minor = celltype_major == type)
  print(length(unique(type_obj$celltype_minor)))
}

# The maximum is 13; choose 5 high contrast colors for plotting stacked proportion bar charts
colors <- c("#9b5de5", "#ff9f1c", "#00bbf9", "#fee440", "#f15bb5", "#8bc34a", "#f44336", "#072ac8",
            "#b52dc8", "#dbf679", "#e76f51", "#a0c4ff", "#51696c")

# For each major cell type:
for(type in unique(seu.TNBC$celltype_major)){
  print(type)
  
  # minor the seurat object to just be the majojr cell type
  type_obj <- minor(seu.TNBC, minor = celltype_major == type)
  
  # Make a data frame of minor cell types within that major cell type, populate with the unique minor cell types
  # and donors' ages as column names
  minors_df <- data.frame(matrix(NA, nrow = length(unique(type_obj$celltype_minor)), ncol = length(unique(type_obj$Age))+1))
  colnames(minors_df) <- c("celltype", unique(type_obj$Age))
  minors_df$celltype <- unique(type_obj$celltype_minor)
  
  # For each minor cell type and for each age represented in the data frame
  for(celltype in minors_df$celltype){
    for(age in colnames(minors_df[,2:length(minors_df)])){
      # calculate the proportion of cells of that minor cell type as a proportion of the major cell type for a given donor
      n <- nrow(type_obj@meta.data[type_obj@meta.data$Age == age & type_obj@meta.data$celltype_minor == celltype,])/nrow(type_obj@meta.data[type_obj@meta.data$Age == age,])
      
      # Add that proportion to the dataframe for the correct cell type/age combo
      minors_df[minors_df$celltype == celltype, age] <- n
    }
  }
  
    # Prep for plotting by making data long
      minors_df_long <- pivot_longer(minors_df, cols = 2:ncol(minors_df))
      
      #plot the stacked bar chart
      p_minors <- ggplot(minors_df_long, aes(x = name, y = value, fill = celltype))+
        geom_bar(stat = "identity", position = "fill")+
        xlab("Age")+
        ylab("Proportion")+
        scale_fill_manual(values = colors[1:length(unique(minors_df_long$celltype))])+
        theme(aspect.ratio = 1)
      
      # Save
      ggsave(paste0(type, "_barchart.pdf"), p_minors, device = "pdf", height = 3, width = 4)
    
      # For each minor cell type make a scatter plot of proportion to age
      for(types in unique(minors_df_long$celltype)){
        p <- ggplot(minors_df_long[minors_df_long$celltype == types,], aes(x = as.numeric(name), y = as.numeric(value)))+
          geom_point()+
          theme(aspect.ratio = 1)+
          geom_smooth(method = "lm", se = T)+
          stat_cor()+
          xlab("Age")+
          ylab("Proportion")
        
        # Print the correlation test
        data <- minors_df_long[minors_df_long$celltype == types,]
        tryCatch({cor <- cor.test(as.numeric(data$name), y = as.numeric(data$value))
          print(types)
          print(cor)
        },
        error = function(e){
          print("")
        })

        # Save the scatter plot
        ggsave(paste0("TNBC_", types, "_cor_to_age.pdf"), p, device = "pdf", height = 2, width =2)
      }
      
      
      }
    

# Change working directory for ER+
setwd("/Users/addie/desktop/Proportions/ER")

# minor the data to just be the ER+ donors and add donorss ages to Seurat object metadata
seu.ER <- minor(seu.all, minor = subtype == "ER+")
seu.ER@meta.data$Age <- NA
for(donor in unique(seu.ER$orig.ident)){
  seu.ER@meta.data$Age[seu.ER@meta.data$orig.ident == donor] <- as.character(age_matrix[rownames(age_matrix) == donor,1])
}

# Check to ensure there are 11 unique ages represented
print(length(unique(seu.ER$Age)))


# Same as TNBC, check the maximum number of minor cell types within each major cell type
for(type in unique(seu.ER$celltype_major)){
  print(type)
  type_obj <- minor(seu.ER, minor = celltype_major == type)
  print(length(unique(type_obj$celltype_minor)))
}

# For each major cell type:
for(type in unique(seu.ER$celltype_major)){
  print(type)
  
  # minor the ER+ seurat object to be just that majojr cell type
  type_obj <- minor(seu.ER, minor = celltype_major == type)
  
  # Make a data frame of minor cell types within that major cell type, populate with the unique minor cell types
  # and donors' ages as column names
  minors_df <- data.frame(matrix(NA, nrow = length(unique(type_obj$celltype_minor)), ncol = length(unique(type_obj$Age))+1))
  colnames(minors_df) <- c("celltype", unique(type_obj$Age))
  minors_df$celltype <- unique(type_obj$celltype_minor)
  
  # For each minor cell type and for each age represented in the data frame
  for(celltype in minors_df$celltype){
    for(age in colnames(minors_df[,2:length(minors_df)])){
      # calculate the proportion of cells of that minor cell type as a proportion of the major cell type for a given donor
      n <- nrow(type_obj@meta.data[type_obj@meta.data$Age == age & type_obj@meta.data$celltype_minor == celltype,])/nrow(type_obj@meta.data[type_obj@meta.data$Age == age,])
      
      # Add that proportion to the dataframe for the correct cell type/age combo
      minors_df[minors_df$celltype == celltype, age] <- n
    }
  }
  
  # Prep for plotting by making data long
  minors_df_long <- pivot_longer(minors_df, cols = 2:ncol(minors_df))
  
  #plot the stacked bar chart
  p_minors <- ggplot(minors_df_long, aes(x = name, y = value, fill = celltype))+
    geom_bar(stat = "identity", position = "fill")+
    xlab("Age")+
    ylab("Proportion")+
    scale_fill_manual(values = colors[1:length(unique(minors_df_long$celltype))])+
    theme(aspect.ratio = 1)
  
  # Save
  ggsave(paste0(type, "_barchart.pdf"), p_minors, device = "pdf", height = 3, width = 4)
  
  # For each minor cell type make a scatter plot of proportion to age
  for(types in unique(minors_df_long$celltype)){
    p <- ggplot(minors_df_long[minors_df_long$celltype == types,], aes(x = as.numeric(name), y = as.numeric(value)))+
      geom_point()+
      theme(aspect.ratio = 1)+
      geom_smooth(method = "lm", se = T)+
      stat_cor()+
      xlab("Age")+
      ylab("Proportion")
    
    # Print the correlation test
    data <- minors_df_long[minors_df_long$celltype == types,]
    tryCatch({cor <- cor.test(as.numeric(data$name), y = as.numeric(data$value))
    print(types)
    print(cor)
    },
    error = function(e){
      print("")
    })
    
    # Save the scatter plot
    ggsave(paste0("ER_", types, "_cor_to_age.pdf"), p, device = "pdf", height = 2, width =2)
  }
  
  
}

