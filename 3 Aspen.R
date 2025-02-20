# Figure 3: ASPEN analysis on Hallmark pathways
###########################################################
# Adrienne Parsons & Peter van Galen, 2024-09-10

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

# Load supplemental table, from manuscript supplementary info, found at link: https://www.nature.com/articles/s41588-021-00911-1#Sec39
sup.tib <- read_excel("./wu,swarbrick2021 - Supplement.xlsx", skip = 3)

# Load expression and metadata, downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
folders <- list.files(".", full.names = T)
folders <- folders[!grepl("xlsx", folders)]
data.tib <- tibble(Sample = cutf(folders, d = "/", f = 2),
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

# GSEA analysis

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
h_gene_sets = msigdbr(species = "human", category = "H")

# make a smaller dataframe of the data we need that will be subsetted out 
hallmark_gene_set_names <- unique(h_gene_sets$gs_name)
gene_set_df <- data.frame(h_gene_sets$gs_name, h_gene_sets$human_gene_symbol)

# Function for signature scoring
signature_scoring <- function(gene_set_names, obj, gsea_data, subtype_str){
  # start by getting a list of genes for each hallmark pathway (gene_list)
  # Then create a new Seurat object with a column
  # of module scores added to the metadata (hallmark_MS)
  # Create a data frame called final_data for each pathway that will be populated
  # with the average module score for each donor/cell type combination later
  donors <- unique(obj@meta.data$orig.ident[obj@meta.data$subtype == subtype_str])
  age_matrix.sub <- subset(age_matrix, subset = rownames(age_matrix) %in% donors)
  
  for( i in gene_set_names){
    print(i)
    gene_list <- list(subset(gene_set_df$h_gene_sets.human_gene_symbol, gene_set_df$h_gene_sets.gs_name == i))
    hallmark_MS <- AddModuleScore(obj, features = gene_list)
    celltype.ls <- unique(obj@meta.data$celltype_minor)
    final_data <- as.data.frame(matrix(NA, nrow = length(celltype.ls), ncol = length(donors)))
    colnames(final_data) <- donors
    rownames(final_data) <- celltype.ls
    
    # This loop takes each donor and creates a data frame (cell_score_df) that is
    # the module score for each cell type for each donor
    # it then creates another data frame for each donor called all_score_df which has 
    # all of the cells and their corresponding module scores
    
    for(j in 1:length(donors)){
      cell_score_df <- data.frame(matrix(NA, nrow = length(celltype.ls), ncol = 1))
      rownames(cell_score_df) <- unique(obj$celltype_minor)
      each_donor <- subset(hallmark_MS@meta.data, hallmark_MS@meta.data$orig.ident == donors[j])
      all_score_df <- data.frame(matrix(NA, nrow = nrow(each_donor), ncol = 2))
      all_score_df[,1] <- each_donor$celltype_minor
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
        cor_vector <- cor(as.numeric(df[g,]), as.numeric(age_matrix.sub), use="complete.obs")
        cor_result <- cor.test(as.numeric(df[g,]), as.numeric(age_matrix.sub), use="complete.obs")
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
subset_masterdata <- function(subtype_str){
  # create a Seurat object with all of the subtype-specific donors
  # also create some variables that will be used later
  obj <- subset(x = seu.all, subset = subtype == subtype_str)
  obj <- NormalizeData(obj)
  donors <- unique(seu.all@meta.data$orig.ident[seu.all@meta.data$subtype == subtype_str])
  cell_type <- unique(obj@meta.data$celltype_minor)
  age_matrix.sub <- subset(age_matrix, subset = rownames(age_matrix) %in% donors)
  
  
  #Iterate through the cell types and make a matrix of gene expression for each cell type
  for(i in 1:length(cell_type)){
    
    age_matrix2 <- age_matrix.sub[,1]
    
    print(i)
    print(cell_type[i])
    
    # Generate empty list of matrices
    expr.ls <- vector(mode = "list", length = length(unique(obj$orig.ident)))
    names(expr.ls) <- unique(obj$orig.ident)
    noCell <- FALSE
    
    # For each donor, get the assay data for the cell type being iterated over
    for (n in unique(obj$orig.ident)) {
      print(n)
      tryCatch({
        expr.ls[[n]] <- GetAssayData(subset(obj, orig.ident == n & celltype_minor == cell_type[i]
        ), slot = "data")
      },
      # Not all cell types are present in each donor; this accounts for that
      error=function(e){
        print(paste0(n," : no cells found"))
        noCell <- TRUE
      })
      if(noCell == TRUE){
        # Remove the donor spaces and expression spaces from the corresponding
        # lists if a given donor doesn't have the cell type
        age_matrix2[n] <<- NA
        expr.ls[[n]] <- NA
        noCell <- FALSE
      }}
    print("removing values from expr.ls")
    expr.ls <<- unlist(expr.ls)
    print(length(expr.ls))
    
    print("removing values from age_matrix2")
    age_matrix2 <<- na.omit(age_matrix2)
    print(length(age_matrix2))
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
      
      age_matrix2 <- age_matrix2[colnames(mean_expr.mat)]
      
      # Make a matrix of correlation coefficients for each gene
      correlation_matrix <- t( cor(matrix(age_matrix2)[,1], t(mean_expr.mat))) 
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
  GO_file <- "h.all.v7.4.symbols.gmt"
  
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

# Run the ASPEN functions for TNBC and ER+
#options(future.globals.maxSize = 8000 * 1024^2)

data_master_TNBC <- subset_masterdata("TNBC")
data_master_ER <- subset_masterdata("ER+")

# Remove unneeded variables from the environment
for(type in unique(seu.all$celltype_minor)){
  rm(list = paste0(type, "ER+"))
  rm(list = paste0(type, "TNBC"))
}

# Get the maximum NES for both TNBC and ER+
maxboth <- max(c(abs(data_master_TNBC$Enrich_Score), 
                 abs(data_master_ER$Enrich_Score)), na.rm = T)

# Generate an emmpty data frame with the TNBC master data column names
data_master_TNBC2 <- data_master_TNBC[1,]
data_master_TNBC2 <- data_master_TNBC2[-1,]

# Subset the data to just be complete cases
for(pathway in unique(data_master_TNBC$Pathway)){
  data1 <- data_master_TNBC[data_master_TNBC$Pathway == pathway,]
  data2 <- data_master_ER[data_master_ER$Pathway == pathway,]
  if(sum(na.omit(data2$Enrich_Score), na.omit(data1$Enrich_Score)) != 0){
    data_master_TNBC2 <- rbind(data_master_TNBC2, data2)
  }
}

# Repeat for ER+
data_master_ER2 <- data_master_ER[1,]
data_master_ER2 <- data_master_ER2[-1,]

# Subset the data to just be complete cases
for(pathway in unique(data_master_ER$Pathway)){
  data1 <- data_master_TNBC[data_master_TNBC$Pathway == pathway,]
  data2 <- data_master_ER[data_master_ER$Pathway == pathway,]
  if(sum(na.omit(data2$Enrich_Score), na.omit(data1$Enrich_Score)) != 0){
    data_master_ER2 <- rbind(data_master_ER2, data2)
  }
}

# Prep for data tables
data_master_TNBC3 <- data_master_TNBC[,c("Pathway", "Celltype", "Cor", "Enrich_Score")]
data_master_TNBC3$Padj <- "> 0.05"
TNBC_pvals$name <- gsub("HALLMARK_", "", TNBC_pvals$name)
TNBC_pvals$name <- gsub("_", " ", TNBC_pvals$name)
TNBC_pvals$Cell <- gsub("_", " ", TNBC_pvals$Cell)
for(i in 1:nrow(TNBC_pvals)){
  data_master_TNBC3$Padj[data_master_TNBC3$Pathway == TNBC_pvals$name[i] & 
                           data_master_TNBC3$Celltype == TNBC_pvals$Cell[i]] <- TNBC_pvals$value[i]
}

data_master_ER3 <- data_master_ER[,c("Pathway", "Celltype", "Cor", "Enrich_Score")]
data_master_ER3$Padj <- "> 0.05"
`ER+_pvals`$name <- gsub("HALLMARK_", "", `ER+_pvals`$name)
`ER+_pvals`$name <- gsub("_", " ", `ER+_pvals`$name)
`ER+_pvals`$Cell <- gsub("_", " ", `ER+_pvals`$Cell)
for(i in 1:nrow(`ER+_pvals`)){
  data_master_ER3$Padj[data_master_ER3$Pathway == `ER+_pvals`$name[i] & 
                           data_master_ER3$Celltype == `ER+_pvals`$Cell[i] &
                           data_master_ER3$Enrich_Score == `ER+_pvals`$NES[i]] <- `ER+_pvals`$value[i]
}

setwd("<YOUR RESULTS PATH>")
write.csv(data_master_TNBC3, file = "TNBC_ASPENResults.csv", row.names = F, quote = F)
write.csv(data_master_ER3, file = "ER_ASPENResults.csv", row.names = F, quote = F)

# Bubble plotting function
bubble_plotting <- function(data, subtype_str, max){
  
  # Define groupings of Hallmark Pathways (all 50)
  Miscellaneous <- c("XENOBIOTIC METABOLISM", "WNT BETA CATENIN SIGNALING", "SPERMATOGENESIS",
                     "PANCREAS BETA CELLS", "NOTCH SIGNALING", "MYOGENESIS",
                      "HEME METABOLISM", "HEDGEHOG SIGNALING", 
                     "COAGULATION", "CHOLESTEROL HOMEOSTASIS", "BILE ACID METABOLISM",
                     "ANGIOGENESIS")
  Cell_stress <- c("UV RESPONSE UP", "UV RESPONSE DN", "UNFOLDED PROTEIN RESPONSE",
                   "REACTIVE OXYGEN SPECIES PATHWAY", "P53 PATHWAY", "HYPOXIA", "G2M CHECKPOINT", "E2F TARGETS",
                   "DNA REPAIR", "APOPTOSIS")
  Inflammation <- c("TNFA SIGNALING VIA NFKB", "INTERFERON GAMMA RESPONSE", "INTERFERON ALPHA RESPONSE",
                    "IL6 JAK STAT3 SIGNALING", "IL2 STAT5 SIGNALING",
                    "INFLAMMATORY RESPONSE", "COMPLEMENT", "ALLOGRAFT REJECTION")
  Cell_cycle_metabolism <- c("PROTEIN SECRETION", "PEROXISOME", 
                             "OXIDATIVE PHOSPHORYLATION", "MTORC1 SIGNALING",
                             "GLYCOLYSIS", "FATTY ACID METABOLISM", "ADIPOGENESIS")
  Cancer <- c("TGF BETA SIGNALING", "PI3K AKT MTOR SIGNALING", "MYC TARGETS V2", "MYC TARGETS V1", "MITOTIC SPINDLE",
              "KRAS SIGNALING UP", "KRAS SIGNALING DN", 
              "ESTROGEN RESPONSE EARLY", "ESTROGEN RESPONSE LATE", "EPITHELIAL MESENCHYMAL TRANSITION",
              "APICAL SURFACE", "APICAL JUNCTION", "ANDROGEN RESPONSE")
  
  # Subset the Hallmark Pathways to just be those represented in the complete cases
  Miscellaneous <- Miscellaneous[Miscellaneous %in% unique(c(data_master_TNBC2$Pathway, data_master_ER2$Pathway))]
  Cell_stress <- Cell_stress[Cell_stress %in% unique(c(data_master_TNBC2$Pathway, data_master_ER2$Pathway))]
  Cell_cycle_metabolism <- Cell_cycle_metabolism[Cell_cycle_metabolism %in% unique(c(data_master_TNBC2$Pathway, data_master_ER2$Pathway))]
  Inflammation <- Inflammation[Inflammation %in% unique(c(data_master_TNBC2$Pathway, data_master_ER2$Pathway))]
  Cancer <- Cancer[Cancer %in% unique(c(data_master_TNBC2$Pathway, data_master_ER2$Pathway))]
  
  # Make a character string of the list names for looping
  hallmark_list <- list(Miscellaneous, Cell_stress, Inflammation, Cell_cycle_metabolism, Cancer)
  names(hallmark_list) <- c("Miscellaneous", "Cell_stress", "Inflammation", "Cell_cycle_metabolism", "Cancer")

  # Subset the data for pathways within each grouping of pathways
  for(ls in 1:length(hallmark_list)){
    pathway_df <- data[data$Pathway %in% hallmark_list[[ls]],]
    
    # Re-format some cell type labels for plotting
    pathway_df$Celltype <- gsub("CAFs MSC ", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub("CAFs ", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub(" SC", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub("Endothelial ", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub(" cells", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub("-cells", "", pathway_df$Celltype)
    pathway_df$Celltype <- gsub("_", " ", pathway_df$Celltype)
    
    # Re-order the newly formatted cell types to group broader cell types together
    pathway_df$Celltype <- factor(pathway_df$Celltype, levels=c("B Memory", "B Naive", "Plasmablasts",
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
    
   # Plot everything!
    p <- ggplot(pathway_df, mapping = aes(x = Celltype,
                                          y= Pathway, size=abs)) +
      geom_point(alpha = 1, aes(fill = Enrich_Score), pch = 21
      ) +  scale_fill_gradientn(colours = c("blue","white","red"), limits = c(max*-1, max))+
    labs(
        size = "Magnitude GSEA Cor.",
        fill = "GSEA NES",
        x = "",
        y = "") +

      guides(colour = guide_legend(order = 1),
             size = guide_legend(order = 2)) +
      scale_shape_identity() +
      theme(axis.text.x = element_text(size = 6, angle = 90, hjust=1, vjust = 1)) +
      theme(axis.text.y = element_text(size = 8))+
      theme(legend.position = "right")+
      theme(legend.key.size = unit(0.3, 'cm'))+
      scale_size(range = c(0, 2))
    
    assign(paste0(names(hallmark_list[ls])), p)}
  

    
    q <- ggarrange(`Cancer` + 
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank() ),
                   `Inflammation` + 
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank() ),
                   `Cell_stress` + 
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank() ),
                   `Cell_cycle_metabolism` + 
                     theme(axis.text.x = element_blank(),
                           axis.ticks.x = element_blank(),
                           axis.title.x = element_blank() ), 
                   `Miscellaneous`, ncol = 1, common.legend = T, align = "v", heights = 
                     c(length(hallmark_list[["Cancer"]]), length(hallmark_list[["Inflammation"]]),
                       length(hallmark_list[["Cell_stress"]]), length(hallmark_list[["Cell_cycle_metabolism"]]),
                       length(hallmark_list[["Miscellaneous"]]) + 8
                     ), font.label = list(size = 8) 
                   ) 
    ggsave(paste0(subtype_str, "_final.pdf"), q, device = "pdf",
           height = 7, width = 5, units = "in")
    
    
}

setwd("<YOUR RESULTS PATH>")

# Run the bubble plotting functions
bubble_plotting(data_master_TNBC, "TNBC", max = maxboth)
bubble_plotting(data_master_ER, "ER+", max = maxboth)
