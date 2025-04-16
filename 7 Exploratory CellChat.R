# Figure 5 and Extended Data Figure 5: Comparison analysis of multiple datasets using CellChat
###########################################################
# Esther Sauras Colon & Adrienne Parsons, 2025-02-28

# Load the required libraries
library(CellChat)
library(patchwork)
library(RColorBrewer)
library(ComplexHeatmap)

## MAJOR CELL TYPES
# First, generate the CellChat objects using the Create CellChat Objects.R file (cellgroup = "celltype_major")
# Load the cellchat objects previously generated
# Create an object list
# DO A SEPARATE ANALYSIS FOR TNBC AND ER

data_path <- "<YOUR DATA PATH>"
options(Seurat.object.assay.version = "v5")
#--------------------------------------------------------------------------------------------------------------------------------
# TNBC, results from "Create CellChat Objects.R"
cellchatTNBC.young <- readRDS(paste0(data_path,"/cellchatTN.young_computeCommunProb_psize_TRUE_celltype_major.Rdata"))
cellchatTNBC.old <- readRDS(paste0(data_path,"/cellchatTN.old_computeCommunProb_psize_TRUE_celltype_major.Rdata"))
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# ER, results from "Create CellChat Objects.R"
cellchatER.young <- readRDS(paste0(data_path,"/cellchatER.young_computeCommunProb_psize_TRUE_celltype_major.Rdata"))
cellchatER.old <- readRDS(paste0(data_path,"/cellchatER.old_computeCommunProb_psize_TRUE_celltype_major.Rdata"))
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
# First, generate the CellChat objects using the Create CellChat Objects.R file (cellgroup = "celltype_minor")
# Load the cellchat objects previously generated

# Create an object list
# DO A SEPARATE ANALYSIS FOR TNBC AND ER

#--------------------------------------------------------------------------------------------------------------------------------
# TNBC, results from "Create CellChat Objects.R"
cellchatTNBC.young <- readRDS(paste0(data_path,"/cellchatTN.young_computeCommunProb_psize_TRUE.Rdata"))
cellchatTNBC.old <- readRDS(paste0(data_path,"/cellchatTN.old_computeCommunProb_psize_TRUE.Rdata"))
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# ER, results from "Create CellChat Objects.R"
cellchatER.young <- readRDS(paste0(data_path,"/cellchatER.young_computeCommunProb_psize_TRUE.Rdata"))
cellchatER.old <- readRDS(paste0(data_path,"/cellchatER.old_computeCommunProb_psize_TRUE.Rdata"))
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
# Remove the following code comment only for ER (ER young has no B cells Naive, so we need to lift up)
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

# --------------------------------------------------------------------------------------------
# Single-donor CellChat analysis

# Heat map of full cohort celtype major (aggregated from multiple patients):
#--------------------------------------------------------------------------------------------------------------------------------
# TNBC, results from "Create CellChat Objects.R
cellchatTNBC.young <- readRDS("<YOUR DATA PATH>/cellchatTN.young_computeCommunProb_psize_TRUE_celltype_major.Rdata")
cellchatTNBC.old <- readRDS("<YOUR DATA PATH>/cellchatTN.old_computeCommunProb_psize_TRUE_celltype_major.Rdata")
object.list <- list(TNBC.young = cellchatTNBC.young, TNBC.old = cellchatTNBC.old)

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

lapply(object.list, function(x) length(levels(x@idents)))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap with the differential interaction strength among different cell populations
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg2

# Extract the communication results from the heatmap data and format for allternative plotting
TNBC_res <- as.data.frame(gg2@matrix)
TNBC_res$Sender <- rownames(TNBC_res)
TNBC_res <- pivot_longer(TNBC_res, cols = 1:ncol(TNBC_res)-1)


# ER, results from "Create CellChat Objects.R"
cellchatER.young <- readRDS("<YOUR DATA PATH>/cellchatER.young_computeCommunProb_psize_TRUE_celltype_major.Rdata")
cellchatER.old <- readRDS("<YOUR DATA PATH>/cellchatER.old_computeCommunProb_psize_TRUE_celltype_major.Rdata")
object.list <- list(ER.young = cellchatER.young, ER.old = cellchatER.old)
#--------------------------------------------------------------------------------------------------------------------------------

# Check min.cells
lapply(object.list, function(x) summary(x@idents))
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
object.list <- lapply(object.list, function(x) filterCommunication(x, min.cells = 3))
object.list <- lapply(object.list, function(x) computeCommunProbPathway(x))
object.list <- lapply(object.list, function(x) aggregateNet(x))

lapply(object.list, function(x) length(levels(x@idents)))

# Merge CellChat object of each dataset together
cellchat <- mergeCellChat(object.list, add.names = names(object.list))

# Heatmap with the differential interaction strength among different cell populations
gg2 <- netVisual_heatmap(cellchat, measure = "weight")
gg2

# Extract the communication results from the heatmap data and format for allternative plotting
ER_res <- as.data.frame(gg2@matrix)
ER_res$Sender <- rownames(ER_res)
ER_res <- pivot_longer(ER_res, cols = 1:ncol(ER_res)-1)

# Subset to just cell types of interest
# The cell types of interest as determined through regression analysis
TNBC_ofinterest <- c("CAFs", "Cancer Epithelial", "Myeloid", "T-cells")
ER_ofinterest <- c("CAFs", "Cancer Epithelial", "Myeloid", "T-cells", "Endothelial", "PVL")

# Subset the heatmap data to just be major cell types of interest
TNBC_res <- TNBC_res[TNBC_res$Sender %in% TNBC_ofinterest & TNBC_res$name %in% TNBC_ofinterest,]
ER_res <- ER_res[ER_res$Sender %in% ER_ofinterest & ER_res$name %in% ER_ofinterest,]

# Make a new column showing sender:receiver
TNBC_res$Commun <- paste0(TNBC_res$Sender, ":", TNBC_res$name)
ER_res$Commun <- paste0(ER_res$Sender, ":", ER_res$name)

# Get the maximum communication probability for plotting
max_TNBC <- max(abs(TNBC_res$value))
max_ER <- max(abs(ER_res$value))

# Order the data so the plot is in descending communication probability
factors_TNBC <- factor(TNBC_res$Commun, levels = levels(fct_reorder(TNBC_res$Commun, TNBC_res$value)))
factors_TNBC2 <- levels(factors_TNBC)

factors_ER <- factor(ER_res$Commun, levels = levels(fct_reorder(ER_res$Commun, ER_res$value)))
factors_ER2 <- levels(factors_ER)

# Plot TNBC
ppp <- ggplot(TNBC_res, aes(x = "Communicaiton", y = fct_reorder(Commun, value), fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick",
                       limits = c(max_TNBC*-1, max_TNBC))+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+theme_minimal()

# Plot ER+
qqq <- ggplot(ER_res, aes(x = "Communicaiton", y = fct_reorder(Commun, value), fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "dodgerblue3", mid = "white", high = "firebrick",
                       limits = c(max_ER*-1, max_ER))+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+theme_minimal()

# View
ppp
qqq

# Save
setwd("<YOUR RESULTS PATH>")
ggsave("TNBC_celltype_major_allcohort_heatmap.pdf", ppp, device = "pdf", width = 6)
ggsave("ER_celltype_major_allcohort_heatmap.pdf", qqq, device = "pdf", width = 6)

#--------------------------------------------------------------------------------------------------------------------------------
# Data path for single-patient CellChat objects
data_path <- "<YOUR DATA PATH>"


# list the files in the data path
files <- list.files(data_path)

# Pull the CellChat objects that include major cell types
majors <- files[grep("major", files)]

# For each objet:
for(obj in majors){
  obj1 <- readRDS(paste0(data_path, "/" ,obj))
  
  # CellChat Pre-processing
  obj1 <- filterCommunication(obj1, min.cells = 3)
  obj1 <- computeCommunProbPathway(obj1)
  obj1 <- aggregateNet(obj1)
  
  # Plot the commmunication network plots for each donor individually
  weight.max <- getMaxWeight(list(obj1), attribute = c("idents","weight"))
  par(mfrow = c(1,2), xpd=TRUE)
  gg <- netVisual_circle(obj1@net$weight, weight.scale = T, label.edge= F, edge.weight.max = weight.max[2], edge.width.max = 12, title.name = paste0("Strength of interactions - ", gsub("_computeCommunProb_psize_T_celltype_major.Rdata", "", obj)))
  
  # Save the communication probability data as a new object with the patient ID and age
  assign(paste0(gsub("_computeCommunProb_psize_T_celltype_major.Rdata", "", obj), "_", unique(obj1@meta$age), "_", unique(obj1@meta$subtype),"_communprobs"), obj1@net$weight)
}

# Initialize a data frame to populate with all of the necessary information
alldata <- data.frame(senders = character(0), 
                      name = character(0),
                      value = numeric(0),
                      commun = character(0),
                      age = numeric(0),
                      subtype = character(0),
                      ncells_total = numeric(0))

# For each object:
for(obj in majors){
  # Read in the object to get a character string to pull in the commmunication data
  print(obj)
  obj1 <- readRDS(paste0(data_path, "/" ,obj))
  data <- paste0(gsub("_computeCommunProb_psize_T_celltype_major.Rdata", "", obj), "_", unique(obj1@meta$age), "_", unique(obj1@meta$subtype),"_communprobs")
  data1 <- data.frame(get(data))
  
  # Format the data for plotting:
  data1$senders <- rownames(data1)
  data2 <- pivot_longer(data1, cols = 1:length(data1)-1)
  
  # Change the values in the "name" column to be consistent with the values in the 
  # "senders" column
  data2$name[data2$name == "B.cells"] <- "B-cells"
  data2$name[data2$name == "T.cells"] <- "T-cells"
  data2$name[data2$name == "Cancer.Epithelial"] <- "Cancer Epithelial"
  data2$name[data2$name == "Normal.Epithelial"] <- "Normal Epithelial"
  
  # Make a new column that lists sender:receiver
  data2$commun <- paste0(data2$senders, ":", data2$name)
  
  # Pull the patient's age from the name of the dataset being formatted
  data2$age = as.numeric(unlist(strsplit(data, split = "_"))[2])
  
  # Pull the patient's molecular subtype from the name of the dataset being formatted
  data2$subtype = unlist(strsplit(data, split = "_"))[3]
  
  # Add a new column to list the total number of cells being analyzed
  data2$ncells_total <- NA
  
  # For each sender/receiver combination
  for(type1 in unique(obj1@meta$celltype_major)){
    for(type2 in unique(obj1@meta$celltype_major)){
      
      # Total the number of cells being analyzed by adding the row numbers in the metadata
      # That represent cells of the sender or receiver type
      ncells <- nrow(obj1@meta[obj1@meta$celltype_major == type1,])+
        nrow(obj1@meta[obj1@meta$celltype_major == type2,])
      
      # If a cell type is missing or the sum ends up being NA, make it 0 instead
      if(is.na(ncells)){
        ncells <- 0
      }
      
      # Add the total number of cells to the appropriate part of the data 
      data2$ncells_total[data2$senders == type1 & data2$name == type2] <- ncells
      
    }
  }
  
  # Aggregate the data into one big data frame
  alldata <- rbind(alldata, data2)
}

# Two TNBC donors are 52; for plotting purposes make one donor 52.1
alldata$age[50:74] <- 52.1

# Subset the data to TNBC patients and TNBC cell types of interest
TNBC <- alldata[alldata$subtype == "TNBC",]
TNBC <- TNBC[TNBC$senders %in% TNBC_ofinterest & TNBC$name %in% TNBC_ofinterest,]
TNBC$commun <- factor(TNBC$commun, levels = factors_TNBC2)

# Subset the data to ER+ patients and ER+ cell types of interest
ER <- alldata[alldata$subtype == "ER",]
ER <- ER[ER$senders %in% ER_ofinterest & ER$name %in% ER_ofinterest,]
ER$commun <- factor(ER$commun, levels = factors_ER2)

# Plot a TNBC heat map
pp <- ggplot(TNBC, aes(x = as.factor(age), y = commun, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "white", mid = "orchid", high = "maroon4",
                       breaks = seq(0, .45, .03), limits = c(0, .45),
                       midpoint = .28)+xlab("Donor Age")+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+
  scale_y_discrete(factors_TNBC2)+
  theme_minimal()

# View
pp

# Plot an ER+ heat map
qq <- ggplot(ER, aes(x = as.factor(age), y = commun, fill = value))+
  geom_tile()+
  scale_fill_gradient2(low = "white", mid = "orchid", high = "maroon4",
                       breaks = seq(0, .18, .02), limits = c(0, .18),
                       midpoint = .09)+xlab("Donor Age")+
  scale_y_discrete(factors_ER2)+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+theme_minimal()

# View
qq

# Save
setwd("<YOUR RESULTS DIRECTORY>")

ggsave("TNBC_communication_strengths_byage_heatmap.pdf", pp, device = "pdf", height = 8, width = 8)
ggsave("ER_communication_strengths_byage_heatmaap.pdf", qq, device = "pdf", width = 8, height = 10)

#-----------------------------------------------------------------------------
# Heat maps averaging probabilities predicted from single donors by age group

# Make a new TNBC data frame for plotting the averages
avg_TNBC <- data_frame(commun = rep(unique(TNBC$commun), 2), age_group = NA, avg_value = NA)
avg_TNBC$age_group[1:16] <- "Young"
avg_TNBC$age_group[17:32] <- "Older"

# For each sender/receiver combination:
for(comm in unique(avg_TNBC$commun)){
  
  # Subset to that sender/receiver combination
  sub <- TNBC[TNBC$commun == comm,]
  
  # Find the average probability by age group
  avg_young <- mean(sub$value[sub$age <= 55])
  avg_old <- mean(sub$value[sub$age > 55])
  
  # Add the averages to the data frame
  avg_TNBC$avg_value[avg_TNBC$commun == comm & avg_TNBC$age_group == "Young"] <-avg_young
  avg_TNBC$avg_value[avg_TNBC$commun == comm & avg_TNBC$age_group == "Older"] <-avg_old
  
}

# Relevel the factors so they are plotted as younger first then older
avg_TNBC$age_group <- factor(avg_TNBC$age_group, levels = c("Young", "Older"))

# Plot
p2 <- ggplot(avg_TNBC, aes(x = age_group, y = commun, fill = avg_value))+
  geom_tile()+
  scale_y_discrete(factors_TNBC2)+
  scale_fill_gradient2(low = "white", mid = "orchid", high = "maroon4",
                       breaks = seq(0, .16, .02), limits = c(0, .16),
                       midpoint = .11)+xlab("Age Group")+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+theme_minimal()

# Repeat for ER+
# Make a new data frame to populate with the average probabilities by age group
avg_ER <- data_frame(commun = rep(unique(ER$commun), 2), age_group = NA, avg_value = NA)
avg_ER$age_group[1:36] <- "Young"
avg_ER$age_group[37:72] <- "Older"

# For each sender/receiver combination:
for(comm in unique(avg_ER$commun)){
  
  # Subset to that sender/receiver combination
  sub <- ER[ER$commun == comm,]
  
  # Find the average probability by age group
  avg_young <- mean(sub$value[sub$age <= 55])
  avg_old <- mean(sub$value[sub$age > 55])
  
  # Add the average probability to the data frame
  avg_ER$avg_value[avg_ER$commun == comm & avg_ER$age_group == "Young"] <-avg_young
  avg_ER$avg_value[avg_ER$commun == comm & avg_ER$age_group == "Older"] <-avg_old
  
}

# Reorder the age groups so it is plotted as younger first then older
avg_ER$age_group <- factor(avg_ER$age_group, levels = c("Young", "Older"))

# Plot
q2 <- ggplot(avg_ER, aes(x = age_group, y = commun, fill = avg_value))+
  geom_tile()+
  scale_fill_gradient2(low = "white", mid = "orchid", high = "maroon4",
                       breaks = seq(0, .09, .02), limits = c(0, .09),
                       midpoint = .055)+xlab("Age Group")+
  scale_y_discrete(factors_ER2)+
  ylab("Source:Target Communication")+labs(fill = "Signaling Probability")+theme_minimal()

# Save both heat maps
ggsave("TNBC_individual_communication_strengths_byagegroup.pdf", p2, device = "pdf", height = 8, width = 6)
ggsave("ER_communication_communication_strengths_byagegroup.pdf", q2, device = "pdf", height = 10,, width = 6)




