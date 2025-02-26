# ASPEN with generalized options for other uses
###########################################################
# Adrienne Parsons & Peter van Galen, 2025-02-26

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
# Ensure proper version of Seurat
options(Seurat.object.assay.version = "v4")

setwd("<YOUR DATA PATH>")
rm(list=ls())

# Function to split character string
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------

# Data should be a Seurat object with a metadata column that is your cell type annotations
# Merge any individual samples into one object
# The continuous variable of interest should be stored as a matrix with the rownames the sample IDs

seu <- "<YOUR SEURAT OBJECT>"
seu <- seu[["RNA3"]] <- as(object = seu[["RNA"]], Class = "Assay")

# Store important variables for generalizability here:
annots <- "<THE METADATA COLUMN NAME THAT INCLUDES CELL TYPE ANNOTATIONS>"
orig.ident <- "<THE METADATA COLUMN NAME THAT INCLUDES SAMPLE IDs>"
cont.mat <- "<THE MATRIX CONTAINING THE CONTINUOUS VARIABLE OF INTEREST>"
# Also decide on a descriptive string for some function arguments; denoted in script to enter as "<YOUR STRING>"

##############################################################################

# GSEA function (adapted from https://bioinformaticsbreakdown.com/how-to-gsea/)
GSEA <- function(gene_list, GO_file, pval) {
  set.seed(54321)
  library(dplyr)
  library(gage)
  library(fgsea)
  
  if ( any( duplicated(names(gene_list)) )  ) {
    warning("Duplicates in gene names")
    gene_list = gene_list[!duplicated(names(gene_list))]
  }
  if  ( !all( order(gene_list, decreasing = TRUE) == 1:length(gene_list)) ){
    warning("Gene list not sorted")
    gene_list = sort(gene_list, decreasing = TRUE)
  }
  
  # Run fGSEA for the given GO file (gmt)
  myGO = fgsea::gmtPathways(GO_file)
  
  
  # AT THIS STAGE, SUBSET myGO TO PATHWAYS OF INTEREST IF NOT INTERESTED IN ENTIRE GENE SET!
  
  print("running fgsea")
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  if(dim(fgRes)[1]!=0){
    
    # Filter FGSEA by using gage results. Must be significant and in same direction to keep 
    print("running gage")
    gaRes = gage::gage(gene_list, gsets=myGO, same.dir=TRUE, set.size =c(15,600))
    
    ups = as.data.frame(gaRes$greater) %>% 
      tibble::rownames_to_column("Pathway") %>% 
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")
    
    downs = as.data.frame(gaRes$less) %>% 
      tibble::rownames_to_column("Pathway") %>% 
      dplyr::filter(!is.na(p.geomean) & q.val < pval ) %>%
      dplyr::select("Pathway")
    
    print("finalizing")
    if(dim(rbind(ups,downs))[1]!=0){
      
      ## Define up / down pathways which are significant in both tests
      keepups = fgRes[fgRes$NES > 0 & !is.na(match(fgRes$pathway, ups$Pathway)), ]
      keepdowns = fgRes[fgRes$NES < 0 & !is.na(match(fgRes$pathway, downs$Pathway)), ]
      
      fgRes = fgRes[ !is.na(match(fgRes$pathway, 
                                  c( keepups$pathway, keepdowns$pathway))), ] %>% 
        arrange(desc(NES))
      fgRes$pathway = stringr::str_replace(fgRes$pathway, "GO_" , "")
      
      fgRes$Enrichment = ifelse(fgRes$NES > 0, "Up-regulated", "Down-regulated")
      filtRes = rbind(head(fgRes, n = 10),
                      tail(fgRes, n = 10 ))
      
      # Prep for plotting and plot
      upcols =  colorRampPalette(colors = c("red4", "red1", "lightpink"))( sum(filtRes$Enrichment == "Up-regulated"))
      downcols =  colorRampPalette(colors = c( "lightblue", "blue1", "blue4"))( sum(filtRes$Enrichment == "Down-regulated"))
      colos = c(upcols, downcols)
      names(colos) = 1:length(colos)
      filtRes$Index = as.factor(1:nrow(filtRes))
      
      g = ggplot(filtRes, aes(reorder(pathway, NES), NES)) +
        geom_col( aes(fill = Index )) +
        scale_fill_manual(values = colos ) +
        coord_flip() +
        labs(x="Pathway", y="Normalized Enrichment Score",
             title="GSEA - Biological Processes") + 
        theme_minimal()
      
      
      output = list("Results" = fgRes, "Plot" = g)
      return(output)
      # If there is an error due to lack of results, add 0 values for each pathway
    }else{
      resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 6))
      rownames(resultsnone) <- names(myGO)
      resultsnone[,1] <- names(myGO)
      return(list("Results"= resultsnone))
    }
    # If there is an error due to lack of results, add 0 values for each pathway
  }else{
    resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 6))
    rownames(resultsnone) <- names(myGO)
    resultsnone[,1] <- names(myGO)
    return(list("Results"= resultsnone))
  }
}

########################################################################################
# Get the Hallmark gene sets from the msigdbr package
h_gene_sets = msigdbr(species = "<YOUR SPECIES OF INTEREST>", category = "<YOUR GENE SET OF INTEREST>") 

# make a smaller dataframe of the data we need that will be subsetted out 
hallmark_gene_set_names <- unique(h_gene_sets$gs_name)
gene_set_df <- data.frame(h_gene_sets$gs_name, h_gene_sets$human_gene_symbol)

# AT THIS STAGE, SUBSET gene_set_df TO PATHWAYS OF INTEREST IF NECESSARY

# Function for signature scoring
signature_scoring <- function(gene_set_names, obj, gsea_data, subtype_str){
  # start by getting a list of genes for each hallmark pathway (gene_list)
  # Then create a new Seurat object with a column
  # of module scores added to the metadata (hallmark_MS)
  # Create a data frame called final_data for each pathway that will be populated
  # with the average module score for each donor/cell type combination later
  donors <- unique(obj@meta.data[,orig.ident])
  matrix.sub <- subset(cont.mat, subset = rownames(cont.mat) %in% donors)
  
  for( i in gene_set_names){
    print(i)
    gene_list <- list(subset(gene_set_df$h_gene_sets.human_gene_symbol, gene_set_df$h_gene_sets.gs_name == i))
    hallmark_MS <- AddModuleScore(obj, features = gene_list)
    celltype.ls <- unique(obj@meta.data[,annots])
    final_data <- as.data.frame(matrix(NA, nrow = length(celltype.ls), ncol = length(donors)))
    colnames(final_data) <- donors
    rownames(final_data) <- celltype.ls
    
    # This loop takes each donor and creates a data frame (cell_score_df) that is
    # the module score for each cell type for each donor
    # it then creates another data frame for each donor called all_score_df which has 
    # all of the cells and their corresponding module scores
    
    for(j in 1:length(donors)){
      cell_score_df <- data.frame(matrix(NA, nrow = length(celltype.ls), ncol = 1))
      rownames(cell_score_df) <- celltype.ls

      indices <- which(hallmark_MS@meta.data[,orig.ident] == j)


      each_donor <- hallmark_MS@meta.data[indices,]
      all_score_df <- data.frame(matrix(NA, nrow = nrow(each_donor), ncol = 2))
      all_score_df[,1] <- each_donor[,annots]
      all_score_df[,2] <- each_donor$Cluster1
      
      # Given a data frame of module scores for all cells for one donor
      # create a new data frame (data1) that is all of the module scores for one cell type
      # then take the mean module score for each cell of that cell type for that donor
      # Then populate a new data frame with the mean module score for each cell type
      
      
      for(k in 1:length(celltype.ls)){
        data1 <- subset(all_score_df, all_score_df[,1] == celltype.ls[k])
        means <- mean(data1[,2])
        if(celltype.ls[k] %in% rownames(cell_score_df)){
          cell_score_df[k,1] <- means
        }
      }
      
      # Going back to the final_data data frame:
      # for each column of final data, assign a donor and its mean module score for
      # each cell type
      # then, going back to the first loop, make a data frame for each hallmark pathway
      # the columns are the donors, the rows are the cell types, and the cells are
      # the mean module scores for each donor/cell type combination for that pathway
      
      final_data[j] <- assign(x = donors[j], cell_score_df)
    }
    assign(x = paste0(i,"_df"), final_data)
  }
  
  # create a new data frame to populate with the correlation values for each
  # pathway and cell type
  
  correlation_df <- as.data.frame(matrix(data = NA, nrow = length(hallmark_gene_set_names), ncol = length(celltype.ls)))
  rownames(correlation_df) <-gene_set_names
  colnames(correlation_df) <- celltype.ls
  
  
  # Create another data frame that is the p-values for each correlation 
  # each pathway and cell type)
  p_val_df <- as.data.frame(matrix(data = NA, nrow = length(hallmark_gene_set_names), ncol = length(celltype.ls)))
  rownames(p_val_df) <- gene_set_names
  colnames(p_val_df) <- celltype.ls
  
  # Populate those data frames!
  # First for loop calls each of the 50 data frames created from the last series of loops
  # The second for loop populates the correlation and p-value data frames
  # "use = 'complete.obs'" makes it so that only points with paired data (no NAs) are 
  # used in the correlation analysis.
  
  for(f in 1:length(gene_set_names)){
    df <- get(paste0(gene_set_names[f], "_df"))
    df[is.na(df)] <- 0
    
    for(g in 1:nrow(df)){
      tryCatch({
        cor_vector <- cor(as.numeric(df[g,]), as.numeric(matrix.sub), use="complete.obs")
        cor_result <- cor.test(as.numeric(df[g,]), as.numeric(matrix.sub), use="complete.obs")
        p_val <- cor_result$p.value
        p_val_df[f,g] <- p_val
        correlation_df[f,g] <- cor_vector
      },
      error=function(e){
        print(rownames(df)[g])
        print(paste0("not enough finite observations"))
      })
    }
  }
  
  
  adj_pval_df <- as.data.frame(matrix(data = NA, nrow = length(gene_set_names), ncol = length(celltype.ls)))
  rownames(adj_pval_df) <- gene_set_names
  colnames(adj_pval_df) <- celltype.ls
  
  for(v in 1:length(p_val_df)){
    cell_list <- as.list(p_val_df[,v])
    adj_pvals <- p.adjust(cell_list, method = "BH", n = length(cell_list))
    adj_pval_df[,v] <- adj_pvals
  }
  
  # Make the final data frame for the GSEA analysis
  data_wide5 <- t(gsea_data[-nrow(gsea_data),])
  
  # Make all the data "long" so it can be merged
  long_cor <- reshape2::melt(as.matrix(correlation_df), value.name = "Cor", na.rm = FALSE)
  long_gsea1 <- reshape2::melt(data_wide5, value.name = "enrich_score")
  long_p_val <- reshape2::melt(as.matrix(adj_pval_df), value.name = "p_value", na.rm = FALSE)
  
  # merge the data into one giant master data set
  data_master <- merge(long_cor, long_gsea1, by = c("Var1", "Var2"), all = TRUE)
  data_master <- merge(data_master, long_p_val, by = c("Var1", "Var2"), all = TRUE)
  data_master$abs <- abs(data_master$Cor)
  
  colnames(data_master) <- c("Pathway", "Celltype", "Cor", "Enrich_Score", "p_Value", "abs")
  # Trim down the row names so they aren't redundant
  data_master$Pathway <- sub("HALLMARK_", "", data_master$Pathway)
  data_master$Pathway <- gsub("_", " ", data_master$Pathway)
  
  # Change NAs up to this point to 0
  data_master[is.na(data_master)] <- 0
  data_master$Celltype <- as.character(data_master$Celltype)
  
  # If no correlation was run on a given cell type, make all pathway entries for that 
  # cell type NA (grey bubbles in plot)
  for(celltype in unique(data_master$Celltype)){
    if(is.numeric(get(paste0(celltype, subtype_str))) == F) {
      data_master$Enrich_Score[data_master$Celltype == celltype] <- NA
    }
  }
  
  # Format column names
  data_master$Celltype <- gsub("_", " ", data_master$Celltype)
  
  # Return the final data frame
  return(data_master)
}


# Function for running ASPEN analysis
subset_masterdata <- function(subtype_str){ # Note: the first subsetting step and subtype_str may not apply;
                                            # If this is the case, then have "subtype_str" be a descriptive string
                                            # for your data 
  # create a Seurat object with all of the subtype-specific donors
  # also create some variables that will be used later
  obj <- "<YOUR SEURAT OBJECT>" # Include a subsetting step if necessary here
  obj <- NormalizeData(obj)
  donors <- unique(obj@meta.data[, orig.ident])
  cell_type <- unique(obj@meta.data[, annots])
  matrix.sub <- subset(cont.mat, subset = rownames(cont.mat) %in% donors)
  
  
  #Iterate through the cell types and make a matrix of gene expression for each cell type
  for(i in 1:length(cell_type)){
    
    matrix2 <- matrix.sub[,1]
    
    print(i)
    print(cell_type[i])
    
    # Generate empty list of matrices
    expr.ls <- vector(mode = "list", length = length(unique(obj@meta.data[,orig.ident])))
    names(expr.ls) <- unique(obj@meta.data[,orig.ident])
    noCell <- FALSE
    
    # For each donor, get the assay data for the cell type being iterated over
    for (n in names(expr.ls)) {
      print(n)
      tryCatch({
        expr.ls[[n]] <- GetAssayData(subset(obj, cells = rownames(obj@meta.data)[obj@meta.data[[orig.ident]] == n &
                                                                                            obj@meta.data[[annots]] == cell_type[i]]), 
                                     slot = "data")
      },
      # Not all cell types are present in each donor; this accounts for that
      error=function(e){
        print(paste0(n," : no cells found"))
        noCell <- TRUE
      })
      if(noCell == TRUE){
        # Remove the donor spaces and expression spaces from the corresponding
        # lists if a given donor doesn't have the cell type
        matrix2[n] <<- NA
        expr.ls[[n]] <- NA
        noCell <- FALSE
      }}
    print("removing values from expr.ls")
    expr.ls <<- unlist(expr.ls)
    print(length(expr.ls))
    
    print("removing values from age_matrix2")
    matrix2 <<- na.omit(matrix2)
    print(length(matrix2))
    print("made expr.ls")
    
    # If fewer than half of the donors have the cell type, don't run the correlation
    if(length(unlist(expr.ls)) <= 0.5*length(donors)){
      print("less than half of donors represented")
      corresults <- NA
      assign(x = paste0(cell_type[i], subtype_str), corresults, envir = .GlobalEnv)
    } else{
      
      print("expr.ls formatted")
      # Take the row means of each element in the list
      means.ls <- lapply(unlist(expr.ls), rowMeans)
      #print("mean.ls made")
      
      # Generate a matrix with genes x samples
      mean_expr.mat <- do.call(cbind, means.ls)
      mean_expr.mat <- as.data.frame(mean_expr.mat)
      # Only use genes that are expressed in >6 samples (>50% of donors)
      
      matrix2 <- matrix2[colnames(mean_expr.mat)]
      
      # Make a matrix of correlation coefficients for each gene
      correlation_matrix <- t( cor(matrix(matrix2)[,1], t(mean_expr.mat))) 
      correlation_matrix[is.na(correlation_matrix)] <- 0
      
      # Order the genes by correlation coefficient
      sorted_correlation_matrix <- correlation_matrix[order(correlation_matrix[,1],decreasing=TRUE),,drop=F]
      sorted_correlation_matrix <- na.omit(sorted_correlation_matrix)
      colnames(sorted_correlation_matrix) <- "cor"
      
      corresults <- c(sorted_correlation_matrix[,1])
      names(corresults) <- row.names(sorted_correlation_matrix)
      
      # Remove any correlation coefficients of 0; these are arbitrarily ranked
      # and will affect outcome of GSEA if included
      corresults <- corresults[corresults != 0]
      assign(x = paste0(cell_type[i], subtype_str), corresults, envir = .GlobalEnv)}
    } 
  
  # Set the working directory and import the GSEA file for analysis
  # gmt file found in Data folder on GitHub or at: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
  setwd("<YOUR DATA PATH>")
  GO_file <- "<YOUR GMT FILE>"
  
  # Run GSEA
  results_df=data.frame(Cell="trash", Pathway="trash", NES="trash", padj = "trash")
  for(j in 1:length(cell_type)){
    temp2 <- list("Results" = NA, "Plot" = NA)
    if(is.numeric(get(paste0(cell_type[j], subtype_str)))){
        print(cell_type[j])
      
      tryCatch({
        temp2 <- GSEA(get(paste0(cell_type[j], subtype_str)), GO_file, pval=0.05)
      }, error = function(e){
        print("no significant pathways")
        temp2 <- list("Results" = NA, "Plot" = NA)
      })
        
      if(class(temp2$Results)=="data.frame"){
        x=cbind(cell_type[j], temp2$Results[,1], temp2$Results[,6], temp2$Results[,3])
        colnames(x)=c("Cell", "Pathway", "NES", "padj")
        results_df=rbind(results_df, x)} else{
          myGO = fgsea::gmtPathways(GO_file)
          resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 4))
          resultsnone[,2] <- names(myGO)
          resultsnone[,1] <- rep(cell_type[j], 50)
          resultsnone[,3] <- rep(0, 50)
          resultsnone[,4] <- rep("NS", 50)
          colnames(resultsnone) <- c("Cell", "Pathway", "NES", "padj")
          results_df=rbind(results_df, resultsnone)
        }
      
      
    } else{
      print("no vector found")
      myGO = fgsea::gmtPathways(GO_file)
      resultsnone <- data.frame(matrix(NA, nrow = length(myGO), ncol = 4))
      resultsnone[,2] <- names(myGO)
      resultsnone[,1] <- rep(cell_type[j], 50)
      resultsnone[3] <- rep(NA, 50)
      resultsnone[4] <- rep(NA, 50)
      colnames(resultsnone) <- c("Cell", "Pathway", "NES", "padj")
      results_df=rbind(results_df, resultsnone)}
    
  }
  
  
  # Prep for bubble plot
  results_df2=results_df
  results_df2=results_df2[-1,]
  results_df3 <- results_df2
  results_df3$padj <- NULL
  data_wide <- spread(results_df3, key = "Pathway", value = "NES", fill = NA)
  
  # Get p-values for data tables
  p_wide <- spread(results_df2, key = "Pathway", value = "padj", fill = NA)
  ps <- pivot_longer(p_wide, cols = 3:length(p_wide))
  
  assign(paste0(subtype_str, "_pvals"), na.omit(ps), envir = .GlobalEnv)
  
  # Index which columns have all NA values
  index <- rowSums(is.na(data_wide)) != ncol(data_wide)-1
  index2 <- which(index == "TRUE")
  
  # Make all of those values 0
  for(i in index2){
    if(rowSums(is.na(data_wide[1,])) != ncol(data_wide)-1){
      nas <- which(is.na(data_wide[i,]) == TRUE)
      data_wide[i, nas] <- 0
      
    }
  }
  
  # Prep data for plotting
  row.names(data_wide)=data_wide[,1]
  data_wide=data_wide[,-1]
  data_wide2=as.matrix(data_wide)
  data_wide3<- matrix(as.numeric(data_wide2),    # Convert to numeric matrix
                      ncol = ncol(data_wide2))
  row.names(data_wide3)=row.names(data_wide2)
  colnames(data_wide3)=colnames(data_wide2)
  data_wide3=rbind(data_wide3, colSums(abs(data_wide3)))
  data_wide4=data_wide3[,order(-data_wide3[nrow(data_wide3),])]
  
  # Arm 2 of ASPEN: Seurat Signature scoring
  data_master <- signature_scoring(hallmark_gene_set_names, obj, data_wide4, subtype_str)
  
  # Remove unneeded variables from the environment
  for(type in cell_type){
    rm(list = paste0(type, subtype_str))
  }
  return(data_master)
}

# Run the ASPEN function
data_master <- subset_masterdata("<YOUR STRING>")


# Remove unneeded variables from the environment
for(type in unique(seu@meta.data[[annots]])){
  rm(list = paste0(type, "<YOUR STRING>"))
}

# Get the maximum NES
maxboth <- max(abs(data_master$Enrich_Score), na.rm = T)

# Bubble plotting function
bubble_plotting <- function(data){
  
  pathway_df <- data
  
  # Plot everything!
  p <- ggplot(pathway_df, mapping = aes(x = Celltype,
                                        y= Pathway, size=abs)) +
    geom_point(alpha = 1, aes(fill = Enrich_Score), pch = 21
    ) +  scale_fill_gradientn(colours = c("blue","white","red"), limits = c(maxboth*-1, maxboth))+
    labs(
      size = "Magnitude GSEA Cor.",
      fill = "GSEA NES",
      x = "",
      y = "") +
    
    guides(colour = guide_legend(order = 1),
           size = guide_legend(order = 2)) +
    scale_shape_identity() +
    theme(axis.text.x = element_text(size = 6, angle = 45, hjust=1, vjust = 1)) +
    theme(axis.text.y = element_text(size = 6))+
    theme(legend.position = "right")+
    theme(legend.key.size = unit(0.3, 'cm'))+
    scale_size(range = c(0, 2))
  
  ggsave(paste0("ASPEN_results.pdf"), p, device = "pdf",
         height = 3, width = 12, units = "in") # Specify length and width that work
}

setwd("<YOUR RESULTS PATH>")

# Run the bubble plotting functions
bubble_plotting(data_master)

