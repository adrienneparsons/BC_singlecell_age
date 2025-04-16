# Figures 6 and 7, Extended Data Figures 7, 8, 9: CellChat RankNet plots and Bubble Plots
###########################################################
# Esther Sauras Colon and Adrienne Parsons, 2024-07-08

# Prerequsites
library(xlsx)
library(tidyr)
library(CellChat)
library(gtools)

## MINOR CELL TYPES
# First, generate the CellChat objects using the Create CellChat Objects.R file (cellgroup = "celltype_minor")
# Load the cellchat objects previously generated

data_path <- "<YOUR DATA PATH>"
options(Seurat.object.assay.version = "v5")
#--------------------------------------------------------------------------------------------------------------------------------
# TNBC, results from "Create CellChat Objects.R"
cellchatTNBC.young <- readRDS(paste0(data_path,"/cellchatTN.young_computeCommunProb_psize_TRUE.Rdata"))
cellchatTNBC.old <- readRDS(paste0(data_path,"/cellchatTN.old_computeCommunProb_psize_TRUE.Rdata"))
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# Merge CellChat object of each dataset together
cellchatTNBC <- mergeCellChat(object.list, add.names = names(object.list))

# ER, results from "Create CellChat Objects.R"
cellchatER.young <- readRDS(paste0(data_path,"/cellchatER.young_computeCommunProb_psize_TRUE.Rdata"))
cellchatER.old <- readRDS(paste0(data_path,"/cellchatER.old_computeCommunProb_psize_TRUE.Rdata"))
object.list <- list(ER.young = cellchatER.young, ER.old = cellchatER.old)

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

# liftCellchat() if necessary 
lapply(object.list, function(x) length(levels(x@idents)))
# ER young has no B cells Naive, so we need to lift up
object.list[[1]] <- liftCellChat(object.list[[1]], levels(object.list[[2]]@idents))

# Merge CellChat object of each dataset together
cellchatER <- mergeCellChat(object.list, add.names = names(object.list))
#--------------------------------------------------------------------------------------------------------------------------------

# Write a function that will be useful for formatting the data later
fun <- function(data){
  
  # Get the subtype for the data
  dataset <- gsub("\\..*","", data$dataset[1])
  
  # Start indexing and set up the final data frame with the correct column names
  i <- 1
  final_data <- data[1,]
  final_data <- final_data[-1,]
  
  # Start with the first line of the input data
  while(i < nrow(data)){
    print(i)
    
    # identify the group name, interaction name, and if the line is representative of old or young
    group.name <- data$group.names[i]
    int.name <- data$interaction_name[i]
    oy <- data$dataset[i]
    
    # Indexing for the next row down
    n <- i+1
    
    # If lines i and n (adjacent lines) have the same group and interactions and the dataset column in line n is not
    # the same value as the dataset value in line i, add both lines to the final data frame and move to the line after
    # n
    if(data$group.names[n] == group.name & data$interaction_name[n] == int.name & data$dataset[n] != oy){
      final_data <- rbind(final_data, data[c(i, n),])
      i <- n+1} else({print(group.name) # if the two adjacent lines are not old and young of the same group & int
                                        # print the group and interaction names and do the following:
        print(int.name)
        print(oy)
        
        # If line i is an "old" line
        if(oy == paste0(dataset, ".old")){
          # Add a "young" line
          print("adding young")
          # Duplicate over the group names, interaction names, etc. to the final data
          # first, add line i
          final_data <- rbind(final_data, data[i,])
          # then duplicate line i but change "old" to young in that newest line
          final_data <- rbind(final_data, data[i,])
          final_data$dataset[nrow(final_data)] <- paste0(dataset, ".young")
          # And change the probability value for that final line to 0
          final_data$prob[nrow(final_data)] <- 0
          # Then move to the next line in the input data and repeat
          i <- i + 1
        }
        
        # If the line i is a "young" line
        if(oy == paste0(dataset, ".young")){
          # Add an "old" line
          print("adding old")
          # Duplicate over the group names, interaction names, etc. to the final data
          # first, add line i
          final_data <- rbind(final_data, data[i,])
          #then duplicate line i but change "young" to old in that newest line
          final_data <- rbind(final_data, data[i,])
          final_data$dataset[nrow(final_data)] <- paste0(dataset, ".old")
          # And change the probability value for that final line to 0
          final_data$prob[nrow(final_data)] <- 0
          # Then move to the next line in the input data and repeat
          i <- i + 1
        }})
    
  }
  # Return the data with included 0 probability lines
  return(final_data)
}



## RankNet TNBC
## iCAF, myCAF, basal, CD4, CD8, macrophages, monocytes

TNBCdir <- "<YOUR DATA PATH FOR TNBC RESULTS>"

# Set up a list of unique pathway names
pathwaysTNBC <- list()

setwd(TNBCdir)
# For each of the seven cell types make a rank net with each of the seven cell types and save
for(type in c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC")){
  for(type2 in c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC")){
    print(paste0(type, "->", type2))
    rn <- rankNet(cellchatTNBC, measure = "weight", mode = "comparison", stacked = F, sources.use = type, 
                  targets.use = type2, do.stat = FALSE, show.raw = FALSE, 
                  title = paste0("Signaling pathways for ", type, "--", type2))+ 
                  theme(axis.text.y = element_text(size=8), title = element_text(size=10))
    
    ggsave(paste0(type, "--", type2, ".pdf"), rn, device = "pdf")
    
    # Get the pathways that are represented in the iterated celltype-celltype interaction
    data <- rankNet(cellchatTNBC, measure = "weight", mode = "comparison", stacked = F, sources.use = type, targets.use = type2, do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
    data <- data$signaling.contribution
    data <- as.character(data$name)
    uniques <- unique(data)
    
    # Add those pathways to the final list
    pathwaysTNBC <- append(pathwaysTNBC, uniques)
  }
}

# Get the number of unique pathways represented in these 49 interactions
length(unique(unlist(pathwaysTNBC)))
pathwaysTNBC <- unlist(pathwaysTNBC)

# Make a table of the frequency of these pathways
table_pathwaysTNBC <- table(pathwaysTNBC)[order(table(pathwaysTNBC), decreasing=T)]
names_pathwaysTNBC <- names(table_pathwaysTNBC[table_pathwaysTNBC != 0])

# TNBC Bubble plots for all pathways
netVisual_bubble(cellchatTNBC, signaling = unique(pathwaysTNBC), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), comparison = c(1, 2), angle.x = 90)
bubble_data_all <- netVisual_bubble(cellchatTNBC, signaling = unique(pathwaysTNBC), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), comparison = c(1, 2), angle.x = 90, return.data = TRUE)
LRpair <- netAnalysis_contribution(object.list[[1]], signaling = unique(pathwaysTNBC), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), return.data = TRUE)

# Keep the ones that have p<0.01 (pval == 3)
bubble_data_all <- bubble_data_all$communication[bubble_data_all$communication$pval == 3,]

# Create the row when the probability is 0
bubble_data_all <- bubble_data_all %>% arrange(interaction_name, desc(group.names))
bubble_data_all_2 <- fun(bubble_data_all)
bubble_dataTNBC <- bubble_data_all_2

# Logistic regression
setwd(TNBCdir)
bubble_data_all_3 <- bubble_data_all_2 %>% arrange(interaction_name, desc(group.names), dataset)
write.xlsx(bubble_data_all_3, "Bubbledata_TNBC.xlsx", sheetName = "Sheet1")

# Generate an empty data frame to populate with logistic regression results
df <- data.frame(pathway_name = character(0),
                 cell_interaction = character(0),
                 TNBC.young = numeric(0),
                 TNBC.old = numeric(0))

# For each row in the ligand-receptor data
for (i in 1:nrow(bubble_data_all_3)){
  # Get the pathway name
  sp <- bubble_data_all_3$pathway_name[i]
  
  # get the interaction name
  int <- bubble_data_all_3$group.names[i]
  
  # Get the mean probability for all LR pairs for that given pathway and interaction for young and old
  prob_TNBC_young <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "TNBC.young" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  prob_TNBC_old <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "TNBC.old" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  
  # Add the probabilities to the blank data frame
  df <- rbind(df, data.frame(pathway_name = sp, cell_interaction = int, TNBC.young = prob_TNBC_young, TNBC.old = prob_TNBC_old))
}

# Reduce redundancies in the final data
df <- unique(df)

# Break up the data into young and old then re-join with annotations
dfyoung <- cbind(df[,1:3],rep(c("Young"),646))
colnames(dfyoung) <- c("pathway_name","cell_interaction","avg_prob","dataset")
dfold <- cbind(df[,c(1:2,4)],rep(c("Old"),646))
colnames(dfold) <- c("pathway_name","cell_interaction","avg_prob","dataset")
df2 <- rbind(dfyoung,dfold)

# Format the data for logistic regression, making na values 0
df3 <- pivot_wider(df2, names_from = pathway_name, values_from = avg_prob)
df3[is.na(df3)] <- 0
write.xlsx(df3, "Bubbledata_TNBC_glm.xlsx", sheetName = "Sheet1")

# Change annotations to 0 and 1 for logistic regression
df3$dataset <- factor(df3$dataset, levels=c("Young","Old"), labels = c(0,1))

# Format the column names for ease
colnames(df3) <- make.names(colnames(df3))

# Univariate logistic regression model
res <- lapply(3:73, function(i) summary(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))
res_confint <- lapply(3:73, function(i) confint(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

# Supplementary table
suppltable <- data.frame(pathway_name = character(0),
                         estimate = numeric(0),
                         pvalue = numeric(0))

ci95 <- data.frame(pathway_name = character(0),
                   "2.5 %" = numeric(0),
                   "97.5 %" = numeric(0)
                   
)
for (x in res) {
  suppltable <- as.data.frame(rbind(suppltable,c(rownames(x$coefficients)[2],x$coefficients[2, "Estimate"], x$coefficients[2, "Pr(>|z|)"])))
}
colnames(suppltable) <- c("pathway_name","estimate","pvalue")

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

write.xlsx(suppltable_values_def,"Bubbledata_TNBC_glm_res.xlsx",sheetName = "Sheet1")

# Bubble plot of TNBC logistic regression pathways
sp_from_lr_ranknet <- c("GALECTIN","CypA","PLAU","FN1","THBS","APP","MHC-I","MPZ","ADGRE")
netVisual_bubble(cellchatTNBC, signaling = sp_from_lr_ranknet, sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), comparison = c(1, 2), angle.x = 45, thresh = 0.01)
databubble <- netVisual_bubble(cellchatTNBC, signaling = sp_from_lr_ranknet, sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Monocyte", "Macrophage", "Cancer Basal SC"), comparison = c(1, 2), angle.x = 45, return.data = TRUE)



## RankNet ER
## iCAF, myCAF, lumA, CD4, CD8, macrophages, PVL Differentiated, ACKR1

ERdir <- "<YOUR DATA PATH FOR ER RESULTS>"

# Set up a list of unique pathway names
pathwaysER <- list()

setwd(ERdir)
# For each of the seven cell types make a rank net with each of the seven cell types and save
for(type in c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC")){
  for(type2 in c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC")){
    print(paste0(type, "->", type2))
    tryCatch({rn <- rankNet(cellchatER, measure = "weight", mode = "comparison", stacked = F, sources.use = type, 
                            targets.use = type2, do.stat = FALSE, show.raw = FALSE, 
                            title = paste0("Signaling pathways for ", type, "--", type2))+ 
      theme(axis.text.y = element_text(size=8), title = element_text(size=10))
    
    ggsave(paste0(type, "--", type2, ".pdf"), rn, device = "pdf")
    
    # Get the pathways that are represented in the iterated celltype-celltype interaction
    data <- rankNet(cellchatER, measure = "weight", mode = "comparison", stacked = F, sources.use = type, targets.use = type2, do.stat = TRUE, show.raw = FALSE, return.data = TRUE)
    data <- data$signaling.contribution
    data <- as.character(data$name)
    uniques <- unique(data)
    
    # Add those pathways to the final list
    pathwaysER <- append(pathwaysER, uniques)}, error = function(e){print("no interaction found!")})
  }
}

# Get the number of unique pathways represented in these 49 interactions
length(unique(unlist(pathwaysER)))
pathwaysER <- unlist(pathwaysER)

# Make a table of the frequency of these pathways
table_pathwaysER <- table(pathwaysER)[order(table(pathwaysER), decreasing=T)]
names_pathwaysER <- names(table_pathwaysER[table_pathwaysER != 0])

# ER Bubble plots for all pathways
netVisual_bubble(cellchatER, signaling = unique(pathwaysER), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), comparison = c(1, 2), angle.x = 90)
bubble_data_all <- netVisual_bubble(cellchatER, signaling = unique(pathwaysER), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), comparison = c(1, 2), angle.x = 90, return.data = TRUE)
LRpair <- netAnalysis_contribution(object.list[[1]], signaling = unique(pathwaysER), sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), return.data = TRUE)

# Keep the ones that have p<0.01 (pval == 3)
bubble_data_all <- bubble_data_all$communication[bubble_data_all$communication$pval == 3,]

# Create the row when the probability is 0
bubble_data_all <- bubble_data_all %>% arrange(interaction_name, desc(group.names))
bubble_data_all_2 <- fun(bubble_data_all)
bubble_dataER <- bubble_data_all_2

# Logistic regression
setwd(ERdir)
bubble_data_all_3 <- bubble_data_all_2 %>% arrange(interaction_name, desc(group.names), dataset)
write.xlsx(bubble_data_all_3, "Bubbledata_ER.xlsx", sheetName = "Sheet1")

# Generate an empty data frame to populate with logistic regression results
df <- data.frame(pathway_name = character(0),
                 cell_interaction = character(0),
                 ER.young = numeric(0),
                 ER.old = numeric(0))

# For each row in the ligand-receptor data
for (i in 1:nrow(bubble_data_all_3)){
  # Get the pathway name
  sp <- bubble_data_all_3$pathway_name[i]
  
  # get the interaction name
  int <- bubble_data_all_3$group.names[i]
  
  # Get the mean probability for all LR pairs for that given pathway and interaction for young and old
  prob_ER_young <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "ER.young" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  prob_ER_old <- mean(bubble_data_all_3$prob[bubble_data_all_3$dataset == "ER.old" & bubble_data_all_3$pathway_name == sp & bubble_data_all_3$group.names == int])
  
  # Add the probabilities to the blank data frame
  df <- rbind(df, data.frame(pathway_name = sp, cell_interaction = int, ER.young = prob_ER_young, ER.old = prob_ER_old))
}

# Reduce redundancies in the final data
df <- unique(df)

# Break up the data into young and old then re-join with annotations
dfyoung <- cbind(df[,1:3],rep(c("Young"),745))
colnames(dfyoung) <- c("pathway_name","cell_interaction","avg_prob","dataset")
dfold <- cbind(df[,c(1:2,4)],rep(c("Old"),745))
colnames(dfold) <- c("pathway_name","cell_interaction","avg_prob","dataset")
df2 <- rbind(dfyoung,dfold)

# Format the data for logistic regression, making na values 0
df3 <- pivot_wider(df2, names_from = pathway_name, values_from = avg_prob)
df3[is.na(df3)] <- 0
write.xlsx(df3, "Bubbledata_ER_glm.xlsx", sheetName = "Sheet1")

# Change annotations to 0 and 1 for logistic regression
df3$dataset <- factor(df3$dataset, levels=c("Young","Old"), labels = c(0,1))

# Format the column names for ease
colnames(df3) <- make.names(colnames(df3))

# Univariate logistic regression model
res <- lapply(3:86, function(i) summary(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))
res_confint <- lapply(3:86, function(i) confint(glm(as.formula(paste0('dataset ~ ', colnames(df3)[i])), data = df3, family = binomial)))

# Supplementary table
suppltable <- data.frame(pathway_name = character(0),
                         estimate = numeric(0),
                         pvalue = numeric(0))

ci95 <- data.frame(pathway_name = character(0),
                   "2.5 %" = numeric(0),
                   "97.5 %" = numeric(0)
                   
)
for (x in res) {
  suppltable <- as.data.frame(rbind(suppltable,c(rownames(x$coefficients)[2],x$coefficients[2, "Estimate"], x$coefficients[2, "Pr(>|z|)"])))
}
colnames(suppltable) <- c("pathway_name","estimate","pvalue")

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

write.xlsx(suppltable_values_def,"Bubbledata_ER_glm_res.xlsx",sheetName = "Sheet1")

# Bubble plot of ER logistic regression pathways
sp_from_lr_ranknet <- c("FN1","THBS","MIF","LAMININ")
netVisual_bubble(cellchatER, signaling = sp_from_lr_ranknet, sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), comparison = c(1, 2), angle.x = 45, thresh = 0.01)
databubble <- netVisual_bubble(cellchatER, signaling = sp_from_lr_ranknet, sources.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), targets.use = c("CAFs MSC iCAF-like", "CAFs myCAF-like", "T cells CD8+", "T cells CD4+", "Endothelial ACKR1", "Macrophage", "PVL Differentiated","Cancer LumA SC"), comparison = c(1, 2), angle.x = 45, return.data = TRUE)

################################################################################
# Get the frequencies of TNBC and ER+ results

# Start wiith ER+
data <- data.frame("Sender" = NA, "Receiver" = NA, "Pathway" = NA, "Probyoung" = NA,
                   "Probold" = NA)
# For the senders and receivers of interest
for(sender in unique(cellchatER@meta$celltype_minor)[c(1,4,5,6,11,12,16,21)]){
  for(receiver in unique(cellchatER@meta$celltype_minor)[c(1,4,5,6,11,12,16,21)]){
    print(paste0(sender, "->", receiver))
    
    # Try to run the RankNet analysis
    tryCatch({
      rn <- rankNet(cellchatER, measure = "weight", mode = "comparison", stacked = F, sources.use = sender, targets.use = receiver, do.stat = FALSE, show.raw = FALSE, return.data = T)
      
      # For each pathway add the sender, receiver, pathway, and probabilities to the final data
      for(pathway in unique(rn$signaling.contribution$name)){
        ls <- c(sender, receiver, pathway, rn$signaling.contribution$contribution.scaled[rn$signaling.contribution$group == "ER.young" & rn$signaling.contribution$name == pathway],
                rn$signaling.contribution$contribution.scaled[rn$signaling.contribution$group == "ER.old" & rn$signaling.contribution$name == pathway])
        data <- rbind(data, ls)
      }
    }, error = function(e){
      print("no communication found")
    })
    
  }
}

# Format the data and get the numbers of times the probability is greater or exclusive to one group
data <- data[-1,]
nrow(data[data$Probyoung > data$Probold,])
nrow(data[data$Probyoung < data$Probold,])
nrow(data[data$Probyoung !=0 & data$Probold == 0,])
nrow(data[data$Probold !=0 & data$Probyoung == 0,])

# Repeat for TNBC
data <- data.frame("Sender" = NA, "Receiver" = NA, "Pathway" = NA, "Probyoung" = NA,
                   "Probold" = NA)

# For each sender and receiver
for(sender in unique(cellchatTNBC@meta$celltype_minor)[c(4,5,8,9,13,14,26)]){
  for(receiver in unique(cellchatTNBC@meta$celltype_minor)[c(4,5,8,9,13,14,26)]){
    print(paste0(sender, "->", receiver))
    
    # Run the ranknet and then add the sender, receiver, pathway, and probabilities to the final data
    tryCatch({
      rn <- rankNet(cellchatTNBC, measure = "weight", mode = "comparison", stacked = F, sources.use = sender, targets.use = receiver, do.stat = FALSE, show.raw = FALSE, return.data = T)
      
      for(pathway in unique(rn$signaling.contribution$name)){
        ls <- c(sender, receiver, pathway, rn$signaling.contribution$contribution.scaled[rn$signaling.contribution$group == "TNBC.young" & rn$signaling.contribution$name == pathway],
                rn$signaling.contribution$contribution.scaled[rn$signaling.contribution$group == "TNBC.old" & rn$signaling.contribution$name == pathway])
        data <- rbind(data, ls)
      }
    }, error = function(e){
      print("no communication found")
    })
    
  }
}

# Get the final data and get the number of times the probability is greater or exclusive to one group
data <- data[-1,]
nrow(data[data$Probyoung > data$Probold,])
nrow(data[data$Probyoung < data$Probold,])
nrow(data[data$Probyoung !=0 & data$Probold == 0,])
nrow(data[data$Probold !=0 & data$Probyoung == 0,])

#################################################################################

# Identify which LR pairs need boxes in Figure 5/6 because the difference is >1.2 fold between old and young
Boxes_ER <- bubble_dataER[bubble_dataER$pathway_name %in% c("MIF", "LAMININ", "THBS", "FN1"),]

# Start index and generate a data frame with the correct column names
n <- 1
Boxes_ER2 <- Boxes_ER[1,]
Boxes_ER2$FC <- NA
Boxes_ER2 <- Boxes_ER2[-1,]

# While the index is less than the number of rows in the data
while(n < nrow(Boxes_ER)){
  print(n)
  subset <- Boxes_ER[c(n,n+1),]
  print(paste0(unique(subset$source), ", ", unique(subset$target), ", ", unique(subset$interaction_name)))
  
  # Get the old and young data for a given source, target, and interaction name
  old <- which(subset$dataset == "ER.old")
  young <- which(subset$dataset == "ER.young")
  
  # If one of the groups has a 0 probability, add the other group to the final data (it gets a box)
  if(0 %in% subset$prob){
    print("0 present")
    to_add <- c(subset[which(subset$prob !=0),], "unique")
    names(to_add) <- c(names(subset[which(subset$prob != 0),]), "FC")
    Boxes_ER2 <- rbind(Boxes_ER2, to_add)
  } else({
    # If they are both nonzero, get the fold change difference
    ovsy <- foldchange(subset$prob[old], subset$prob[young])

    # If the fold change difference is >1.2 or <-1.2 add that data to the final data (one row gets a box)
    if(abs(ovsy) > 1.2){
      print(ovsy)
      to_add <- c(subset[which(subset$prob == max(subset$prob)),], ovsy)
      names(to_add) <- c(names(subset[which(subset$prob == max(subset$prob)),]), "FC")
      Boxes_ER2 <- rbind(Boxes_ER2, to_add)
    } else(Boxes_ER2 <- Boxes_ER2)
  })
  
  # Then move to the next pair
  n <- n+2
}

write.csv(Boxes_ER2, file = "foldchangesERbubble.csv")

# Figure out the balance of young to old boxes in ER+ bubble data
nrow(Boxes_ER2[Boxes_ER2$dataset == "ER.old",])
nrow(Boxes_ER2[Boxes_ER2$dataset == "ER.young",])

# Repeat for TNBC
# Identify which LR pairs need boxes in Figure 5/6 because the difference is >1.2 fold between old and young
Boxes_TNBC <- bubble_dataTNBC[bubble_dataTNBC$pathway_name %in% c("GALECTIN","CypA","PLAU","FN1","THBS","APP","MHC-I","MPZ","ADGRE"),]

# Start index and generate a data frame with the correct column names
n <- 1
Boxes_TNBC2 <- Boxes_TNBC[1,]
Boxes_TNBC2$FC <- NA
Boxes_TNBC2 <- Boxes_TNBC2[-1,]

# While the index is less than the number of rows in the data
while(n < nrow(Boxes_TNBC)){
  print(n)
  subset <- Boxes_TNBC[c(n,n+1),]
  print(paste0(unique(subset$source), ", ", unique(subset$target), ", ", unique(subset$interaction_name)))
  
  # Get the old and young data for a given source, target, and interaction name
  old <- which(subset$dataset == "TNBC.old")
  young <- which(subset$dataset == "TNBC.young")
  
  # If one of the groups has a 0 probability, add the other group to the final data (it gets a box)
  if(0 %in% subset$prob){
    print("0 present")
    to_add <- c(subset[which(subset$prob !=0),], "unique")
    names(to_add) <- c(names(subset[which(subset$prob != 0),]), "FC")
    Boxes_TNBC2 <- rbind(Boxes_TNBC2, to_add)
  } else({
    # If they are both nonzero, get the fold change difference
    ovsy <- foldchange(subset$prob[old], subset$prob[young])

    # If the fold change difference is >1.2 or <-1.2 add that data to the final data (one row gets a box)
    if(abs(ovsy) > 1.2){
      print(ovsy)
      to_add <- c(subset[which(subset$prob == max(subset$prob)),], ovsy)
      names(to_add) <- c(names(subset[which(subset$prob == max(subset$prob)),]), "FC")
      Boxes_TNBC2 <- rbind(Boxes_TNBC2, to_add)
    } else(Boxes_TNBC2 <- Boxes_TNBC2)
  })
  
  # Then move to the next pair
  n <- n+2
}

nrow(Boxes_TNBC2[Boxes_TNBC2$dataset == "TNBC.old",])
nrow(Boxes_TNBC2[Boxes_TNBC2$dataset == "TNBC.young",])

write.csv(Boxes_TNBC2, file = "foldchangesTNBCbubble.csv")

