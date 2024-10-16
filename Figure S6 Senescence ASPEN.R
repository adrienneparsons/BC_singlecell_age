# Figure S6: ASPEN on senescence-associated pathways
###########################################################
# Adrienne Parsons, 2024-09-10

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

# Set up ------------------------------------------------------------------------------------------
options(Seurat.object.assay.version = "v4")

setwd("<YOUR DATA PATH>")
#setwd("~/DropboxMGB/Projects/Single-cell_breast_cancer/AnalysisAdrienne") # for Peter
rm(list=ls())

# Function to split character string
cutf <- function(x, f=1, d="/") sapply(strsplit(x, d), function(i) paste(i[f], collapse=d))

# Load data ---------------------------------------------------------------------------------------

# Load supplemental table, found at this link: https://www.nature.com/articles/s41588-021-00911-1#Sec39
sup.tib <- read_excel("./wu,swarbrick2021 - Supplement.xlsx", skip = 3)

# Load expression and metadata, downloaded from GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078
folders <- list.files("../Data_Swarbrick", full.names = T)
folders <- folders[!grepl("xlsx", folders)]
data.tib <- tibble(Sample = cutf(folders, d = "/", f = 2),
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

##############################################################################
# Generate a matrix of the donors' ages
age_matrix <- matrix(pull(sup.tib, 3))
rownames(age_matrix) <- gsub("-", "", paste0("CID", sup.tib$`Case ID`))
rownames(age_matrix) <- gsub("CID4290", "CID4290A", rownames(age_matrix))
rownames(age_matrix) <- gsub("CID4530", "CID4530N", rownames(age_matrix))
colnames(age_matrix) <- "age"

# GSEA analysis

# GSEA analysis
# Set the working directory and import the GSEA file for analysis
# gmts can be downloaded here: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
setwd("<YOUR DATA PATH>")
GO_fileC2 <- "c2.all.v2023.2.Hs.symbols.gmt"
GO_fileC5 <- "c5.all.v2023.2.Hs.symbols.gmt"
senescence_genesets <- c("REACTOME_SENESCENCE_ASSOCIATED_SECRETORY_PHENOTYPE_SASP",
                         "GOBP_CELLULAR_SENESCENCE",
                         "CHICAS_RB1_TARGETS_SENESCENT",
                         "FRIDMAN_SENESCENCE_UP",
                         "FRIDMAN_SENESCENCE_DN",
                         "KAMMINGA_SENESCENCE",
                         "REACTOME_CELLULAR_SENESCENCE",
                         "REACTOME_DNA_DAMAGE_TELOMERE_STRESS_INDUCED_SENESCENCE",
                         "REACTOME_ONCOGENE_INDUCED_SENESCENCE",
                         "TANG_SENESCENCE_TP53_TARGETS_DN",
                         "TANG_SENESCENCE_TP53_TARGEETS_UP",
                         "WP_GLYCOLYSIS_IN_SENESCENCE",
                         "WP_NAD_METABOLISM_IN_ONCOGENE_INDUCED_SENESCENCE_AND_MITOCHONDRIAL_DYSFUNCTIONASSOCIATED_SENESCENCE",
                         "WP_PROSTAGLANDIN_AND_LEUKOTRIENE_METABOLISM_IN_SENESCENCE",
                         "WP_TCA_CYCLE_IN_SENESCENCE",
                         "GOBP_NEGATIVE_REGULATION_OF_CELLULAR_SENESCENCE",
                         "GOBP_STRESS_INDUCED_PREMATURE_SENESCENCE"
                         )

# GSEA function (adapted from https://bioinformaticsbreakdown.com/how-to-gsea/)
GSEA <- function(gene_list, pval) {
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
  
  # Get the gene set info for C2 and C5 then subset down to the specific senescence-related ones
  myGOa = fgsea::gmtPathways(GO_fileC2)
  myGOb = fgsea::gmtPathways(GO_fileC5)
  myGO <- c(myGOa, myGOb)
  myGO <- myGO[names(myGO) %in% senescence_genesets]
  
  print("running fgsea")
  
  # Run fgsea
  fgRes <- fgsea::fgsea(pathways = myGO, 
                        stats = gene_list,
                        minSize=15,
                        maxSize=600) %>% 
    as.data.frame() %>% 
    dplyr::filter(padj < !!pval)
  if(dim(fgRes)[1]!=0){
    
    ## Filter FGSEA by using gage results. Must be significant and in same direction to keep 
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
      
      # Optional plotting
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
    }else{
      
      # If there are no significant pathways, generate a dataframe of 0s with the right dimensions of
      # pathways x results
      myGOa = fgsea::gmtPathways(GO_fileC2)
      myGOb = fgsea::gmtPathways(GO_fileC5)
      myGO <- c(myGOa, myGOb)
      myGO <- myGO[names(myGO) %in% senescence_genesets]
      resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 6))
      rownames(resultsnone) <- names(myGO)
      resultsnone[,1] <- names(myGO)
      return(list("Results"= resultsnone))
    }
  }else{
    # If there are no significant pathways, generate a dataframe of 0s with the right dimensions of
    # pathways x results
    myGOa = fgsea::gmtPathways(GO_fileC2)
    myGOb = fgsea::gmtPathways(GO_fileC5)
    myGO <- c(myGOa, myGOb)
    myGO <- myGO[names(myGO) %in% senescence_genesets]
    resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 6))
    rownames(resultsnone) <- names(myGO)
    resultsnone[,1] <- names(myGO)
    return(list("Results"= resultsnone))
  }
}

########################################################################################
# Get the Hallmark gene sets from the msigdbr package
gene_setsC2 = msigdbr(species = "human", category = "C2")
gene_setsC5 = msigdbr(species = "human", category = "C5")
gene_sets <- rbind(gene_setsC2, gene_setsC5)

# make a smaller dataframe of the data we need that will be subsetted out 
gene_set_names <- unique(gene_sets$gs_name)
gene_set_df <- data.frame(gene_sets$gs_name, gene_sets$human_gene_symbol)
gene_set_names <- gene_set_names[gene_set_names %in% senescence_genesets]
gene_set_df <- gene_set_df[gene_set_df$gene_sets.gs_name %in% senescence_genesets,]

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
    gene_list <- list(subset(gene_set_df$gene_sets.human_gene_symbol, gene_set_df$gene_sets.gs_name == i))
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
  
  correlation_df <- as.data.frame(matrix(data = NA, nrow = length(gene_set_names), ncol = length(celltype.ls)))
  rownames(correlation_df) <-gene_set_names
  colnames(correlation_df) <- celltype.ls
  
  
  # Create another data frame that is the p-values for each correlation 
  # each pathway and cell type)
  p_val_df <- as.data.frame(matrix(data = NA, nrow = length(gene_set_names), ncol = length(celltype.ls)))
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
  
  # Get the p-values for the signature scoring correlation
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
  
  # Change NAs to 0 and then return NA values to cell types that have no enrichment score (not enough cells per donor)
  data_master[is.na(data_master)] <- 0
  data_master$Celltype <- as.character(data_master$Celltype)
  
  for(celltype in unique(data_master$Celltype)){
    if(is.numeric(get(paste0(celltype, subtype_str))) == F) {
      data_master$Enrich_Score[data_master$Celltype == celltype] <- NA
    }
  }
  
  # Format the Cell type labels
  data_master$Celltype <- gsub("_", " ", data_master$Celltype)
  
  
# Return the final data frame
  return(data_master)
}


# Function for running fulll ASPEN analysis
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
      # Only use genes that are expressed in >6 samples
      
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
  
  # Run GSEA
  results_df=data.frame(Cell="trash", Pathway="trash", NES="trash", padj = "trash")
  for(j in 1:length(cell_type)){
    temp2 <- list("Results" = NA, "Plot" = NA)
    if(is.numeric(get(paste0(cell_type[j], subtype_str)))){
      print(cell_type[j])
      
      tryCatch({
        temp2 <- GSEA(get(paste0(cell_type[j], subtype_str)), pval=0.05)
      }, error = function(e){
        print("no significant pathways")
        temp2 <- list("Results" = NA, "Plot" = NA)
      })
      
      if(class(temp2$Results)=="data.frame"){
        x=data.frame(cbind(cell_type[j], temp2$Results[,1], temp2$Results[,6], temp2$Results[,3]))
        colnames(x)=c("Cell", "Pathway", "NES", "padj")
        results_df=rbind(results_df, x)} else{
          myGOa = fgsea::gmtPathways(GO_fileC2)
          myGOb = fgsea::gmtPathways(GO_fileC5)
          myGO <- c(myGOa, myGOb)
          myGO <- myGO[names(myGO) %in% senescence_genesets]
          resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 4))
          resultsnone[,2] <- names(myGO)
          resultsnone[,1] <- rep(cell_type[j], length(myGO))
          resultsnone[,3] <- rep(0, length(myGO))
          resultsnone[,4] <- rep("NS", length(myGO))
          colnames(resultsnone) <- c("Cell", "Pathway", "NES", "padj")
          results_df=rbind(results_df, resultsnone)
        }
      
      
    } else{
      print("no vector found")
      myGOa = fgsea::gmtPathways(GO_fileC2)
      myGOb = fgsea::gmtPathways(GO_fileC5)
      myGO <- c(myGOa, myGOb)
      myGO <- myGO[names(myGO) %in% senescence_genesets]
      resultsnone <- data.frame(matrix(0, nrow = length(myGO), ncol = 4))
      resultsnone[,2] <- names(myGO)
      resultsnone[,1] <- rep(cell_type[j], length(myGO))
      resultsnone[,3] <- rep(0, length(myGO))
      resultsnone[,4] <- rep("NS", length(myGO))
      colnames(resultsnone) <- c("Cell", "Pathway", "NES", "padj")
      results_df=rbind(results_df, resultsnone)}
    
  }
  
  
  
  # Prep for bubble plot
  results_df2=results_df
  results_df2=results_df2[-1,]
  results_df3 <- results_df2
  results_df3$padj <- NULL
  data_wide <- spread(results_df3, key = "Pathway", value = "NES", fill = NA)
  
  # Get p-values for data table
  p_wide <- spread(results_df2, key = "Pathway", value = "padj", fill = NA)
  ps <- pivot_longer(p_wide, cols = 3:length(p_wide))
  
  assign(paste0(subtype_str, "_pvals"), na.omit(ps), envir = .GlobalEnv)
  
  # Find the pathways that have all NA values for every cell type
  index <- rowSums(is.na(data_wide)) != ncol(data_wide)-1
  index2 <- which(index == "TRUE")
  
  for(i in index2){
    if(rowSums(is.na(data_wide[1,])) != ncol(data_wide)-1){
      nas <- which(is.na(data_wide[i,]) == TRUE)
      data_wide[i, nas] <- 0
      
    }
  }
  
  # Finish prepping the data for bubble plot
  row.names(data_wide)=data_wide[,1]
  data_wide=data_wide[,-1]
  data_wide2=as.matrix(data_wide)
  data_wide3<- matrix(as.numeric(data_wide2),    # Convert to numeric matrix
                      ncol = ncol(data_wide2))
  row.names(data_wide3)=row.names(data_wide2)
  colnames(data_wide3)=colnames(data_wide2)
  data_wide3=rbind(data_wide3, colSums(abs(data_wide3)))
  data_wide4=data_wide3[,order(-data_wide3[nrow(data_wide3),])]
  
  # Generate the final ASPEN results
  data_master <- signature_scoring(gene_set_names, obj, data_wide4, subtype_str)
  
  # Remove variables from the environment that aren't needed
  for(type in cell_type){
    rm(list = paste0(type, subtype_str))
  }
  return(data_master)
}

# Run the functions for TNBC and ER+
data_master_TNBC <- subset_masterdata("TNBC")
data_master_ER <- subset_masterdata("ER+")

# Remove variables from the environment that aren't needed
for(type in unique(seu.all$celltype_minor)){
  rm(list = paste0(type, "ER+"))
  rm(list = paste0(type, "TNBC"))
}

# Get the maximum NES for both TNBC and ER+
maxboth <- max(c(abs(data_master_TNBC$Enrich_Score), 
                 abs(data_master_ER$Enrich_Score)), na.rm = T)

# Make an empty data frame with the column names from data.master
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

# Bubble plotting function
bubble_plotting <- function(data, subtype_str){
  
  data$Celltype <- gsub("CAFs MSC ", "", data$Celltype)
  data$Celltype <- gsub("CAFs ", "", data$Celltype)
  data$Celltype <- gsub(" SC", "", data$Celltype)
  data$Celltype <- gsub("Endothelial ", "", data$Celltype)
  data$Celltype <- gsub(" cells", "", data$Celltype)
  data$Celltype <- gsub("-cells", "", data$Celltype)
  data$Celltype <- gsub("_", " ", data$Celltype)
  
  data$Celltype <- factor(data$Celltype, levels=c("B Memory", "B Naive", "Plasmablasts",
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
  
  pathway_df <- data
  
  aspect_ratio <- (((length(unique(pathway_df$Pathway))) / 2.5))/ (ncol(pathway_df)+1.5)
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
  
  ggsave(paste0(subtype_str, "_senescence.pdf"), p, device = "pdf",
         height = 3, width = 12, units = "in")
}

setwd("<YOUR RESULTS PATH>")

# Run the bubble plotting functions
bubble_plotting(data_master_TNBC, "TNBC")
bubble_plotting(data_master_ER, "ER+")
