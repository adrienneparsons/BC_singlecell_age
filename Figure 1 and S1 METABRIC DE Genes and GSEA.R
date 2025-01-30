# Figure 1: DE gene analysis and GSEA, METABRIC and TCGA
###########################################################
# Adrienne Parsons and Esther Sauras Colon, 2025-01-14
# Get the DE genes between old and young ER+ and TNBC donors in METABRIC and TCGA
# for disease subtype, then run GSEA

# Prerequisite packages
library(ggplot2)
library(forcats)
library(edgeR)
library(limma)
library(Biobase)
library(fgsea)
library(ggpubr)
library(dplyr)
library(ggrepel)

# First, load in the METABRIC clinical data and format, downloaded from cBioPortal: https://www.cbioportal.org/study/summary?id=brca_metabric
clin_metabric <- read.table("<YOUR DATA PATH>/data_clinical_patient.txt",
                            sep = "\t")
colnames(clin_metabric) <- clin_metabric[1,]
clin_metabric <- clin_metabric[-1,]

# Load in the clinical samples data to get stage and format, also downloaded from cBioPortal: https://www.cbioportal.org/study/summary?id=brca_metabric
clin_sample <- read.table("<YOUR DATA PATH>/data_clinical_sample.txt",
                          sep = "\t")
colnames(clin_sample) <- clin_sample[1,]
clin_sample <- clin_sample[-1,]

# Add columns to the clinical data to include tumor stage and make numeric values numeric vectors
clin_metabric$AGE_AT_DIAGNOSIS <- as.numeric(clin_metabric$AGE_AT_DIAGNOSIS)
clin_metabric$RFS_MONTHS <- as.numeric(clin_metabric$RFS_MONTHS)
clin_metabric$TUMOR_STAGE <- clin_sample$TUMOR_STAGE
clin_metabric$TUMOR_STAGE <- as.numeric(clin_metabric$TUMOR_STAGE)

# TNBC donors for each age group
TNBC_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                              clin_metabric$THREEGENE == "ER-/HER2-" &
                              clin_metabric$TUMOR_STAGE <= 3 &
                              clin_metabric$TUMOR_STAGE != 0,]

TNBC_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                            clin_metabric$THREEGENE == "ER-/HER2-"&
                            clin_metabric$TUMOR_STAGE <= 3&
                            clin_metabric$TUMOR_STAGE != 0,]

# ER+ low prolif donors for each age group
ER_lp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                               clin_metabric$THREEGENE == "ER+/HER2- Low Prolif"&
                               clin_metabric$TUMOR_STAGE <= 3&
                               clin_metabric$TUMOR_STAGE != 0,]

ER_lp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- Low Prolif"&
                             clin_metabric$TUMOR_STAGE <= 3&
                             clin_metabric$TUMOR_STAGE != 0,]

# ER+ high prolif donors for each age group
ER_hp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                               clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                               clin_metabric$TUMOR_STAGE <= 3&
                               clin_metabric$TUMOR_STAGE != 0,]

ER_hp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                             clin_metabric$TUMOR_STAGE <= 3&
                             clin_metabric$TUMOR_STAGE != 0,]

# Remove NA values from the subsetted clinical data so as not to include 
# donors with missing values
for(c in c("ER_hp_old", "ER_hp_young", "ER_lp_old", "ER_lp_young",
           "TNBC_old", "TNBC_young")){
  df <- get(c)
  df <- na.omit(df)
  assign(c, df)
}

# We decided to group all ER+ donors together irrespective of proliferation status
ER_young <- rbind(ER_hp_young, ER_lp_young)
ER_old <- rbind(ER_hp_old, ER_lp_old)

# Remove variables to avoid confusion
rm("ER_hp_old", "ER_lp_old", "ER_hp_young", "ER_lp_young")

# Read in and format the gene expression data for TNBC and ER+
# These data can be found on FigShare; see README for links
#these are the resultant files from "METABRIC data wrangling.R"
TNBC <- read.csv("<YOUR DATA PATH>/2024_METABRIC_TNBC.csv")
rownames(TNBC) <- TNBC[,1]
TNBC <- TNBC[,-1]


ER <- read.csv("<YOUR DATA PATH>/2024_METABRIC_ER.csv")
rownames(ER) <- ER[,1]
ER <- ER[,-1]

# Get an annotation of gene names in the correct order
TNBC_genes <- rownames(TNBC)
ER_genes <- rownames(ER)

# Write the donor information tables
TNBC_old2 <- TNBC_old[, c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "THREEGENE", "TUMOR_STAGE")]
TNBC_young2 <- TNBC_young[, c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "THREEGENE", "TUMOR_STAGE")]
ER_old2 <- ER_old[, c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "THREEGENE", "TUMOR_STAGE")]
ER_young2 <- ER_young[, c("PATIENT_ID", "AGE_AT_DIAGNOSIS", "THREEGENE", "TUMOR_STAGE")]


TNBC_Table1 <- rbind(TNBC_old2, TNBC_young2)
ER_Table1 <- rbind(ER_old2, ER_young2)

setwd("<YOUR RESULTS PATH>")
write.csv(file = "METABRIC_TNBC_donors.csv", TNBC_Table1)
write.csv(file = "METABRIC_ER_donors.csv", ER_Table1)

# Start METABRIC DE gene analysis

# Make a data frame of the formatted gene expression, but only ER+ donors
# within the older and younger groups
ER_young2 <- ER[,make.names(colnames(ER)) %in% make.names(ER_young$PATIENT_ID)]
ER_old2 <- ER[,make.names(colnames(ER)) %in% make.names(ER_old$PATIENT_ID)]
df2 <- cbind(ER_young2, ER_old2)

# Make an expressionset of these data
aseset <- ExpressionSet(assayData = as.matrix(df2))

# Find the genes with the highest standard deviation
var <- apply(exprs(aseset), 1, sd)
var <- var[order(var)]

# Plot
plot(1:length(var), var)

# Only run limma on genes with highest stdev
df3 <- df2[rownames(df2) %in% names(var[(length(var)-674):length(var)]),]

# exponentiate the data to raw count values (METABRIC is log2 transformed on cBioPortal)
df3 <- 2^df3
# Create a DGEList object
d0 <- DGEList(df3)

# Calculate the normalization factors
d0 <- calcNormFactors(d0)

# Make the annotations for DE gene expression analysis
annots <- c(rep("young", length(colnames(ER_young2) %in% colnames(df3))),
            rep("old", length(colnames(ER_old2) %in% colnames(df3))))


# Specify the model to be fitted
mm <- model.matrix(~0 + annots)

# Normalize using Voom 
y <- voom(d0, mm, plot = F)

# fit a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)

# specify which groups to compare
contr <- makeContrasts(annotsold - annotsyoung, levels = colnames(coef(fit)))

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# empirical Bayes
tmp <- eBayes(tmp)

# get DE genes
top.table_ER <- topTable(tmp, sort.by = "logFC", n = Inf)

# Remove variables to prevent confusion
rm("df", "df2", "df3", "fit", "mm", "tmp", "y", "d0")

# Repeat for TNBC. First, filter the gene expression data to just be the oldest and youngest TNBC
# donors
TNBC_young2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_young$PATIENT_ID)]
TNBC_old2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_old$PATIENT_ID)]
df2 <- cbind(TNBC_young2, TNBC_old2)

# Make that data an ExpressionSet
aseset <- ExpressionSet(assayData = as.matrix(df2))

# Find the genes with the highest spread of standard deviation expression
var <- apply(exprs(aseset), 1, sd)
var <- var[order(var)]

# plot
plot(1:length(var), var)

# Subset the data to just be the genes with the highest spread of expression
df3 <- df2[rownames(df2) %in% names(var[(length(var)-674):length(var)]),]

# Exponentiate the data to resemble raw counts (cBioPortal log2-normalizes data)
df3 <- 2^df3

# Create a DGEList object
d0 <- DGEList(df3)

# Calculate the normalization factors
d0 <- calcNormFactors(d0)

# Make the annotations for DE gene expression analysis
annots <- c(rep("young", length(colnames(TNBC_young2) %in% colnames(df3))),
            rep("old", length(colnames(TNBC_old2) %in% colnames(df3))))

# Specify the model to be fitted
mm <- model.matrix(~0 + annots)

# Normalize using Voom 
y <- voom(d0, mm, plot = F)

# fit a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)

# specify which groups to compare
contr <- makeContrasts(annotsold - annotsyoung, levels = colnames(coef(fit)))

#estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# empirical Bayes
tmp <- eBayes(tmp)

# get DE genes
top.table_TNBC <- topTable(tmp, sort.by = "logFC", n = Inf)

# Remove variables to prevent confusion
rm("df2", "df3", "fit", "mm", "tmp", "y", "d0")

# Plotting
setwd("<YOUR RESULTS PATH>")

# For the two results tables:

for(df in c("top.table_ER", "top.table_TNBC")){
  print(df)
  dd3 <- get(df)
  
  # Copy the data for volcano plot
  dd3$gene <- rownames(dd3)
  dd3$abs <- abs(dd3$logFC)
  
  # Generate a Volcano plot based on the limma analysis differentially
  # abundant genes
  dd3$minuslog10p <- -log10(dd3$adj.P.Val)
  dd3$direction <- NA
  
  # Only add gene names to the most differentially expressed genes
  # Add the gene names to any gene with adjusted p < 0.05 and magnitude fold change > 0.5
  dd3$label <- NA
  dd3$label[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 0.5] <- dd3$gene[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 0.5]
  
  # Add in enrichment annotations for coloring the dots
  dd3$direction <- "enriched in older"
  dd3$direction[dd3$logFC < 0] <- "enriched in younger"
  dd3$direction[dd3$minuslog10p < -log10(0.05)] <- "0"
  
  # Write a table of the genes
  dd4 <- arrange(dd3, desc(minuslog10p), desc(logFC))
  dd5 <- dd4[,c("gene", "logFC", "minuslog10p", "P.Value", "adj.P.Val")]
  assign(paste0(gsub("top.table_", "", df), "_METABRIC"), dd5)
  
  #write.csv(dd5, file = paste0(gsub("top.table_", "", df), "_METABRIC_volcanocoords.csv"), quote = F,row.names = F)
}

# Clear the environment and prepare for GSEA
#rm(list = ls())
#---------------------------------------------------------------------------------------
wdir <- "<YOUR DATA PATH>"
# First, load in the TCGA clinical data and format, downloaded from cBioPortal: https://www.cbioportal.org/study/summary?id=brca_tcga_pan_can_atlas_2018
clin_tcga <- read.table(paste0(wdir,"/brca_tcga_pan_can_atlas_2018_clinical_data.tsv"),
                        sep = "\t")
colnames(clin_tcga) <- clin_tcga[1,]
clin_tcga <- clin_tcga[-1,]

# Make numeric values numeric vectors and change "Sample ID"
clin_tcga$`Diagnosis Age` <- as.numeric(clin_tcga$`Diagnosis Age`)
clin_tcga$`Progress Free Survival (Months)` <- as.numeric(clin_tcga$`Progress Free Survival (Months)`)
clin_tcga$TUMOR_STAGE <- as.numeric(sub("T([0-9]+).*", "\\1", clin_tcga$`American Joint Committee on Cancer Tumor Stage Code`)) #TX appears as NA
clin_tcga$`Sample ID` <- gsub("-", ".", clin_tcga$`Sample ID`)

# TNBC donors for each age group (BRCA_Basal)
TNBC_young <- clin_tcga[clin_tcga$`Diagnosis Age` < 45 & 
                          !is.na(clin_tcga$`Diagnosis Age`) &
                          clin_tcga$Subtype == "BRCA_Basal" &
                          !is.na(clin_tcga$Subtype) &
                          clin_tcga$TUMOR_STAGE <= 3 &
                          !is.na(clin_tcga$TUMOR_STAGE),] #nrow 30

TNBC_old <- clin_tcga[clin_tcga$`Diagnosis Age` > 65 &
                        !is.na(clin_tcga$`Diagnosis Age`) &
                        clin_tcga$Subtype == "BRCA_Basal"&
                        !is.na(clin_tcga$Subtype) &
                        clin_tcga$TUMOR_STAGE <= 3 &
                        !is.na(clin_tcga$TUMOR_STAGE),] #nrow 37

# ER+ donors for each age group (BRCA_LumA)
ER_young <- clin_tcga[clin_tcga$`Diagnosis Age` < 45 & 
                        !is.na(clin_tcga$`Diagnosis Age`) &
                        clin_tcga$Subtype == "BRCA_LumA" &
                        !is.na(clin_tcga$Subtype) &
                        clin_tcga$TUMOR_STAGE <= 3 &
                        !is.na(clin_tcga$TUMOR_STAGE),] #nrow 68

ER_old <- clin_tcga[clin_tcga$`Diagnosis Age` > 65 &
                      !is.na(clin_tcga$`Diagnosis Age`) &
                      clin_tcga$Subtype == "BRCA_LumA" &
                      !is.na(clin_tcga$Subtype) &
                      clin_tcga$TUMOR_STAGE <= 3 &
                      !is.na(clin_tcga$TUMOR_STAGE),] #nrow 152

# Read in and format the gene expression data for TNBC and ER+
# Download brca_tcga_pan_can_atlas_2018.tar.gz from cBioPortal
# Uncompress data_mrna_seq_v2_rsem.txt
X <- read.table("<YOUR DATA PATH>/data_mrna_seq_v2_rsem.txt", header=TRUE, sep="\t", fill=TRUE, stringsAsFactors=FALSE)

TNBC_col <- colnames(X)[colnames(X) %in% TNBC_old$`Sample ID` | colnames(X) %in% TNBC_young$`Sample ID`]
TNBC_col <- unique(c(colnames(X)[1], TNBC_col))
TNBC_pre <- X[,TNBC_col, drop=FALSE]
TNBC_pre <- TNBC_pre[TNBC_pre$Hugo_Symbol != "",] #Filter the rows with Hugo_Symbol
TNBC_pre[duplicated(TNBC_pre$Hugo_Symbol)| duplicated(TNBC_pre$Hugo_Symbol, fromLast = TRUE),] #There are 7 genes duplicated
TNBC <- TNBC_pre %>%
  group_by(Hugo_Symbol) %>%  
  slice_max(across(where(is.numeric)), with_ties = FALSE) %>% #Select the row with the highest gene expression values
  ungroup() %>%
  as.data.frame()
rownames(TNBC) <- TNBC[,1]
TNBC <- TNBC[,-1]


ER_col <- colnames(X)[colnames(X) %in% ER_old$`Sample ID` | colnames(X) %in% ER_young$`Sample ID`]
ER_col <- unique(c(colnames(X)[1], ER_col))
ER_pre <- X[,ER_col, drop=FALSE]
ER_pre <- ER_pre[ER_pre$Hugo_Symbol != "",] #Filter the rows with Hugo_Symbol
#options(max.print=9999)
ER_pre[duplicated(ER_pre$Hugo_Symbol)| duplicated(ER_pre$Hugo_Symbol, fromLast = TRUE),] #There are 7 genes duplicated
ER <- ER_pre %>%
  group_by(Hugo_Symbol) %>%  
  slice_max(across(where(is.numeric)), with_ties = FALSE) %>% #Select the row with the highest gene expression values
  ungroup() %>%
  as.data.frame()
rownames(ER) <- ER[,1]
ER <- ER[,-1]

setwd("<YOUR RESULTS PATH>")
# Get the necessary donor information for table 1
TNBC_old2 <- TNBC_old[, c("Sample ID", "Diagnosis Age", "Subtype", "TUMOR_STAGE")]
TNBC_young2 <- TNBC_young[, c("Sample ID", "Diagnosis Age", "Subtype", "TUMOR_STAGE")]
ER_old2 <- ER_old[, c("Sample ID", "Diagnosis Age", "Subtype", "TUMOR_STAGE")]
ER_young2 <- ER_young[, c("Sample ID", "Diagnosis Age", "Subtype", "TUMOR_STAGE")]


TNBC_Table1 <- rbind(TNBC_old2, TNBC_young2)
ER_Table1 <- rbind(ER_old2, ER_young2)

write.csv(file = "TCGA_TNBC_donors.csv", TNBC_Table1)
write.csv(file = "TCGA_ER_donors.csv", ER_Table1)


# Get an annotation of gene names in the correct order
TNBC_genes <- rownames(TNBC)
ER_genes <- rownames(ER)


# Start DE gene analysis
# RSEM (batch normalized from Illumina HiSeq_RNASeqV2)

# Make a data frame of the formatted gene expression, but only ER+ donors
# within the older and younger groups
ER_young2 <- ER[,make.names(colnames(ER)) %in% make.names(ER_young$`Sample ID`)]
ER_old2 <- ER[,make.names(colnames(ER)) %in% make.names(ER_old$`Sample ID`)]
df2 <- cbind(ER_young2, ER_old2)

group_ER <- as.data.frame(rbind(cbind(ER_young$`Sample ID`,c("Young")),cbind(ER_old$`Sample ID`,c("Old"))))

# Filter low expressed genes
# Density plot of log-CPM values for raw data and filtered data (based on cpm)
cpmvalue <- 1
repthreshold <- 10
keep <- rowSums(cpm(df2) > cpmvalue) >= repthreshold
df3 <- df2[keep,]
dim(df3) # 16151 x 220 keep

par(mfrow=c(1,2)) 
#plot(density(cpm(df2, log=T)[,1]),lwd=2,ylim=c(0,0.20),las=2,main="",xlab="") 
#title(main="A. Raw data",xlab="Log-cpm") 
#abline(v=log2(cpmvalue),lty=3) 
#for (i in 2:ncol(df2)){  
#  den<-density(cpm(df2, log=T)[,i])  
#  lines(den$x,den$y,lwd=2)
#} 

#plot(density(cpm(df3, log=T)[,1]),lwd=2,ylim=c(0,0.20),las=2,main="",xlab="") 
#title(main="B. Filtered data",xlab="Log-cpm") 
#abline(v=log2(cpmvalue),lty=3) 
#for (i in 2:ncol(df3)){  
#  den<-density(cpm(df3, log=T)[,i])  
#  lines(den$x,den$y,lwd=2) 
#} 

# Unormalized data vs normalized data
# Create a DGEList object
d1 <- DGEList(df3)

par(mfrow=c(1,2)) 
#boxplot(cpm(d1,log=TRUE),las=2,main="") 
#title(main="A. Unnormalized data",ylab="Log-cpm") 

# Calculate the normalization factors
d0 <- calcNormFactors(d1)
#boxplot(cpm(d0,log=TRUE),las=2,main="") 
#title(main="B. Normalized data",ylab="Log-cpm")

# Make the annotations for DE gene expression analysis
annots <- c(rep("young", length(colnames(ER_young2) %in% colnames(df3))),
            rep("old", length(colnames(ER_old2) %in% colnames(df3))))

dev.off()
#plotMDS(d0, col=as.numeric(as.factor(annots)))


# Specify the model to be fitted
mm <- model.matrix(~0 + annots)

# Normalize using Voom 
y <- voom(d0, mm, plot = T)

# Fit a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)

# Specify which groups to compare
contr <- makeContrasts(annotsold - annotsyoung, levels = colnames(coef(fit)))

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# Empirical Bayes
tmp <- eBayes(tmp)

# Get DE genes
top.table_ER <- topTable(tmp, sort.by = "logFC", n = Inf)
nrow(top.table_ER[top.table_ER$adj.P.Val < 0.05,]) #2588 genes with adj.P.Val < 0.05

# Remove variables to prevent confusion
rm("df2", "df3", "fit", "mm", "tmp", "y", "d0", "d1")

# Repeat for TNBC. First, filter the gene expression data to just be the oldest and youngest TNBC
# donors
TNBC_young2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_young$`Sample ID`)]
TNBC_old2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_old$`Sample ID`)]
df2 <- cbind(TNBC_young2, TNBC_old2)

# Filter low expressed genes
# Density plot of log-CPM values for raw data and filtered data (based on cpm)
cpmvalue <- 1
repthreshold <- 10
keep <- rowSums(cpm(df2) > cpmvalue) >= repthreshold
df3 <- df2[keep,]
dim(df3) # 15439 x 67

par(mfrow=c(1,2)) 
#plot(density(cpm(df2, log=T)[,1]),lwd=2,ylim=c(0,0.20),las=2,main="",xlab="") 
#title(main="A. Raw data",xlab="Log-cpm") 
#abline(v=log2(cpmvalue),lty=3) 
#for (i in 2:ncol(df2)){  
#  den<-density(cpm(df2, log=T)[,i])  
#  lines(den$x,den$y,lwd=2)
#} 

#plot(density(cpm(df3, log=T)[,1]),lwd=2,ylim=c(0,0.20),las=2,main="",xlab="") 
#title(main="B. Filtered data",xlab="Log-cpm") 
#abline(v=log2(cpmvalue),lty=3) 
#for (i in 2:ncol(df3)){  
#  den<-density(cpm(df3, log=T)[,i])  
#  lines(den$x,den$y,lwd=2) 
#} 

# Unormalized data vs normalized data
# Create a DGEList object
d1 <- DGEList(df3)

par(mfrow=c(1,2)) 
#boxplot(cpm(d1,log=TRUE),las=2,main="") 
#title(main="A. Unnormalized data",ylab="Log-cpm") 

# Calculate the normalization factors
d0 <- calcNormFactors(d1)
#boxplot(cpm(d0,log=TRUE),las=2,main="")
#title(main="B. Normalized data",ylab="Log-cpm")

# Make the annotations for DE gene expression analysis
annots <- c(rep("young", length(colnames(TNBC_young2) %in% colnames(df3))),
            rep("old", length(colnames(TNBC_old2) %in% colnames(df3))))

dev.off()
#plotMDS(d0, col=as.numeric(as.factor(annots)))

# Specify the model to be fitted
mm <- model.matrix(~0 + annots)

# Normalize using Voom 
y <- voom(d0, mm, plot = T)

# Fit a linear model using weighted least squares for each gene
fit <- lmFit(y, mm)

# Specify which groups to compare
contr <- makeContrasts(annotsold - annotsyoung, levels = colnames(coef(fit)))

# Estimate contrast for each gene
tmp <- contrasts.fit(fit, contr)

# Empirical Bayes
tmp <- eBayes(tmp)

# Get DE genes
top.table_TNBC <- topTable(tmp, sort.by = "logFC", n = Inf)
nrow(top.table_TNBC[top.table_TNBC$adj.P.Val < 0.05,]) #0 genes with adj.P.Val < 0.05
nrow(top.table_TNBC[top.table_TNBC$P.Value < 0.05,]) #1666 genes with P.Val < 0.05

# Remove variables to prevent confusion
rm("df2", "df3", "fit", "mm", "tmp", "y", "d0", "d1")

# Plotting
setwd("YOUR RESULTS PATH>")

# For the ER results table:

for(df in c("top.table_ER")){
  print(df)
  dd3 <- get(df)
  
  # Copy the data for volcano plot
  dd3$gene <- rownames(dd3)
  dd3$abs <- abs(dd3$logFC)
  
  # Generate a Volcano plot based on the limma analysis differentially
  # abundant genes
  dd3$minuslog10p <- -log10(dd3$adj.P.Val)
  dd3$direction <- NA
  
  # Only add gene names to the most differentially expressed genes
  # Add the gene names to any gene with adjusted p < 0.05 and magnitude fold change > 1.4
  dd3$label <- NA
  dd3$label[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 1.4] <- dd3$gene[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 1.4]
  
  # Add in enrichment annotations for coloring the dots
  dd3$direction <- "enriched in older"
  dd3$direction[dd3$logFC < 0] <- "enriched in younger"
  dd3$direction[dd3$minuslog10p < -log10(0.05)] <- "0"
  
  # Write a table of the genes
  dd4 <- arrange(dd3, desc(minuslog10p), desc(logFC))
  dd5 <- dd4[,c("gene", "logFC", "minuslog10p", "P.Value", "adj.P.Val")]
  
  assign(paste0(gsub("top.table_", "", df), "_TCGA"), dd5)
  
  #write.csv(dd5, file = paste0(gsub("top.table_", "", df), "_TCGA_volcanocoords.csv"), quote = F, row.names = F)
}

# For the TNBC results table:

for(df in c("top.table_TNBC")){
  print(df)
  dd3 <- get(df)
  
  # Copy the data for volcano plot
  dd3$gene <- rownames(dd3)
  dd3$abs <- abs(dd3$logFC)
  
  # Generate a Volcano plot based on the limma analysis differentially
  # abundant genes
  dd3$minuslog10p <- -log10(dd3$P.Value)
  dd3$direction <- NA
  
  # Only add gene names to the most differentially expressed genes
  # Add the gene names to any gene with adjusted p < 0.05 and magnitude fold change > 1.4
  dd3$label <- NA
  dd3$label[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 1.4] <- dd3$gene[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 1.4]
  
  # Add in enrichment annotations for coloring the dots
  dd3$direction <- "enriched in older"
  dd3$direction[dd3$logFC < 0] <- "enriched in younger"
  dd3$direction[dd3$minuslog10p < -log10(0.05)] <- "0"
  
  # Write a table of the genes
  dd4 <- arrange(dd3, desc(minuslog10p), desc(logFC))
  dd5 <- dd4[,c("gene", "logFC", "minuslog10p", "P.Value", "adj.P.Val")]
  assign(paste0(gsub("top.table_", "", df), "_TCGA"), dd5)
  #write.csv(dd5, file = paste0(gsub("top.table_", "", df), "_TCGA_volcanocoords.csv"), quote = F, row.names = F)
}

# Clear the environment and prepare for GSEA
#rm(list = ls())

# Volcano plotting
setwd("<YOUR RESULTS PATH>")

# Load the data
TCGA_TNBC <- read.csv("TNBC_TCGA_volcanocoords.csv")
METABRIC_TNBC <- read.csv("TNBC_METABRIC_volcanocoords.csv")
TCGA_ER <- read.csv("ER_TCGA_volcanocoords.csv")
METABRIC_ER <- read.csv("ER_METABRIC_volcanocoords.csv")

# Check the genes that are significant and overlap between datasets (in the same direction)
ER_overlaps_up <- as.data.frame(Reduce(intersect, list(METABRIC_ER$gene[METABRIC_ER$adj.P.Val < 0.05 & METABRIC_ER$logFC > 0], TCGA_ER$gene[TCGA_ER$adj.P.Val < 0.05 & TCGA_ER$logFC > 0])))
TNBC_overlaps_up <- as.data.frame(Reduce(intersect, list(METABRIC_TNBC$gene[METABRIC_TNBC$adj.P.Val < 0.05 & METABRIC_TNBC$logFC > 0], TCGA_TNBC$gene[TCGA_TNBC$P.Value < 0.05 & TCGA_TNBC$logFC > 0])))
ER_overlaps_down <- as.data.frame(Reduce(intersect, list(METABRIC_ER$gene[METABRIC_ER$adj.P.Val < 0.05 & METABRIC_ER$logFC < 0], TCGA_ER$gene[TCGA_ER$adj.P.Val < 0.05 & TCGA_ER$logFC < 0])))
TNBC_overlaps_down <- as.data.frame(Reduce(intersect, list(METABRIC_TNBC$gene[METABRIC_TNBC$adj.P.Val < 0.05 & METABRIC_TNBC$logFC < 0], TCGA_TNBC$gene[TCGA_TNBC$P.Value < 0.05 & TCGA_TNBC$logFC < 0])))

# Format and add a direction annotation
colnames(ER_overlaps_up) <- "gene"
colnames(ER_overlaps_down) <- "gene"
colnames(TNBC_overlaps_up) <- "gene"
colnames(TNBC_overlaps_down) <- "gene"
ER_overlaps_up$direction <- "Up"
ER_overlaps_down$direction <- "Down"
TNBC_overlaps_up$direction <- "Up"
TNBC_overlaps_down$direction <- "Down"

# Create the final data frames and save
ER_overlaps <- rbind(ER_overlaps_up, ER_overlaps_down)
TNBC_overlaps <- rbind(TNBC_overlaps_up, TNBC_overlaps_down)
write.csv(ER_overlaps, "Figure1_commongenes_ER.csv", row.names = F)
write.csv(TNBC_overlaps, "Figure1_commongenes_TNBC.csv", row.names = F)

# Add labels to significant and high fold difference data points
# For METABRIC, take all genes with log fold change >0.5 or <-.5 and padj < 0.05
METABRIC_TNBC$label <- NA
METABRIC_TNBC$label[abs(METABRIC_TNBC$logFC) > 0.5 & METABRIC_TNBC$minuslog10p > 1.75] <- METABRIC_TNBC$gene[abs(METABRIC_TNBC$logFC) > 0.5 & METABRIC_TNBC$minuslog10p > 1.75]
TNBC_labels <- na.omit(METABRIC_TNBC$label)

METABRIC_ER$label <- NA
METABRIC_ER$label[abs(METABRIC_ER$logFC) > 0.5 & METABRIC_ER$minuslog10p > 7.5] <- METABRIC_ER$gene[abs(METABRIC_ER$logFC) > 0.5 & METABRIC_ER$minuslog10p > 7.5]
ER_labels <- na.omit(METABRIC_ER$label)

# Repeat for TCGA
TCGA_TNBC$label <- NA
TCGA_TNBC$label[abs(TCGA_TNBC$logFC) > 0.5 & TCGA_TNBC$minuslog10p > 3.25] <- TCGA_TNBC$gene[abs(TCGA_TNBC$logFC) > 0.5 & TCGA_TNBC$minuslog10p > 3.25]
TNBC_labels_TCGA <- na.omit(TCGA_TNBC$label)

TCGA_ER$label <- NA
TCGA_ER$label[abs(TCGA_ER$logFC) > 0.5 & TCGA_ER$minuslog10p > 7.5] <- TCGA_ER$gene[abs(TCGA_ER$logFC) > 0.5 & TCGA_ER$minuslog10p > 7.5]
ER_labels_TCGA <- na.omit(TCGA_ER$label)

# Check overlaps of labeled genes
Reduce(intersect, list(ER_labels, ER_labels_TCGA))
Reduce(intersect, list(TNBC_labels, TNBC_labels_TCGA))


# Plot
for(df in c("METABRIC_ER", "METABRIC_TNBC", "TCGA_ER", "TCGA_TNBC")){
  print(df)
  df1 <- get(df)
  df1$direction <- "enriched in older"
  df1$direction[df1$logFC < 0] <- "enriched in younger"
  df1$direction[df1$minuslog10p < -log10(0.05)] <- "0"
  df1$abs <- abs(df1$logFC)
  # Make the volcano plot and save
  volcano <- ggplot(data=df1, aes(x=logFC, y=minuslog10p, label = label, fill = direction)) + 
    geom_point(size = .8, shape = 21, stroke = 0)+
    xlim(c(max(df1$abs)*-1-1, max(df1$abs)+1))+theme_bw()+
    scale_fill_manual(values = c("lightgrey", "red", "blue"))+
    ylim(c(0, max(df1$minuslog10p) + 0.75))+
    xlab("Log2 Fold Change")+
    ylab("-log10(FDR)")+
    geom_text_repel(max.overlaps = 30, size = 1.75, box.padding = 0.1)+
    theme(aspect.ratio = 1)
  
  ggsave(paste0(df, "_volcano.pdf"), volcano, device = "pdf", height = 2, width =4)
}



# ---------------------------------------------------------------------------------------
# Running GSEA on results; gmt files obtained from this link: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
setwd("Users/addie/Desktop/GSEA/gmts")
gmts <- list.files("<YOUR DATA PATH>/gmts")

# For the ER+ and TNBC results:
for(res in c("TCGA_ER", "TCGA_TNBC", "METABRIC_ER", "METABRIC_TNBC")){
  
  # Get the data and extract the logFC for ranking
  results <- get(res)
  results$logFC <- as.numeric(results$logFC)
  results$P.Value <- as.numeric(results$P.Value)
  results$adj.P.Val <- as.numeric(results$adj.P.Val)
  results <- results %>% arrange(-logFC)
  genes <- results$logFC
  names(genes) <- results$gene
  
  # Set up the GSEA results data frame and remove the row of NAs
  gsea_data <- data.frame("pathway" = NA, "pval" = NA, "padj" = NA, "log2err" = NA,
                          "ES" = NA, "NES" = NA, "size" = NA, "leadingEdge" = NA, "gmt" = NA)
  gsea_data <- gsea_data[-1,]
  
  
  # For each gmt file (different gene set)
  # Run gsea using fgsea and the ranked gene list
  # if there are any pathways with adjusted p-values < 0.05 add it to the final data
  for(gmt in gmts){
    print(gmt)
    set.seed(42)
    myGO <- fgsea::gmtPathways(gmt)
    gsea_out <- fgsea(pathways = myGO,
                      stats = genes)
    gsea_out2 <- gsea_out[gsea_out$padj < 0.05]
    gsea_out2$gmt <- gmt
    gsea_data <- rbind(gsea_data, gsea_out2)
  }
  
  # Compile all results and write to the environment
  assign(paste0("sig_gsea_", res), gsea_data)
}

setwd("<YOUR RESULTS PATH>")

# Order the significant TNBC GSEA results by decreasing NES and write to computer
sig_gsea_METABRIC_TNBC <- sig_gsea_METABRIC_TNBC[order(sig_gsea_METABRIC_TNBC$NES, decreasing = T),]
sig_gsea_METABRIC_TNBC <- sig_gsea_METABRIC_TNBC[,-c(7,8)]
sig_gsea_METABRIC_TNBC2 <- apply(sig_gsea_METABRIC_TNBC,2,as.character)
write.csv(sig_gsea_METABRIC_TNBC2, file = "METABRIC_TNBC_GSEA.csv", quote = F, row.names = F)

sig_gsea_TCGA_TNBC <- sig_gsea_TCGA_TNBC[order(sig_gsea_TCGA_TNBC$NES, decreasing = T),]
sig_gsea_TCGA_TNBC <- sig_gsea_TCGA_TNBC[,-c(7,8)]
sig_gsea_TCGA_TNBC2 <- apply(sig_gsea_TCGA_TNBC,2,as.character)
write.csv(sig_gsea_TCGA_TNBC2, file = "TCGA_TNBC_GSEA.csv", quote = F, row.names = F)

# Repeat for ER+
sig_gsea_METABRIC_ER <- sig_gsea_METABRIC_ER[order(sig_gsea_METABRIC_ER$NES, decreasing = T),]
sig_gsea_METABRIC_ER <- sig_gsea_METABRIC_ER[,-c(7,8)]
sig_gsea_METABRIC_ER2 <- apply(sig_gsea_METABRIC_ER,2,as.character)
write.csv(sig_gsea_METABRIC_ER2, file = "METABRIC_ER_GSEA.csv", quote = F, row.names = F)

sig_gsea_TCGA_ER <- sig_gsea_TCGA_ER[order(sig_gsea_TCGA_ER$NES, decreasing = T),]
sig_gsea_TCGA_ER <- sig_gsea_TCGA_ER[,-c(7,8)]
sig_gsea_TCGA_ER2 <- apply(sig_gsea_TCGA_ER,2,as.character)
write.csv(sig_gsea_TCGA_ER2, file = "TCGA_ER_GSEA.csv", quote = F, row.names = F)
# Clear the environment
rm(list = ls())

#---------------------------------------------------------------------------------------
# Load the GSEA results from METABRIC and TCGA
sig_GSEA_ER_METABRIC <- read.csv("<YOUR RESULTS PATH>/METABRIC_ER_GSEA.csv")
sig_GSEA_TNBC_METABRIC <- read.csv("<YOUR RESULTS PATH>/METABRIC_TNBC_GSEA.csv")
sig_GSEA_ER_TCGA <- read.csv("<YOUR RESULTS PATH>/TCGA_ER_GSEA.csv")
sig_GSEA_TNBC_TCGA <- read.csv("<YOUR RESULTS PATH>/TCGA_TNBC_GSEA.csv")

ER_pathways <- Reduce(intersect, list(sig_GSEA_ER_METABRIC$pathway, sig_GSEA_ER_TCGA$pathway))
TNBC_pathways <- Reduce(intersect, list(sig_GSEA_TNBC_METABRIC$pathway, sig_GSEA_TNBC_TCGA$pathway))

ER_pathways
TNBC_pathways

# Manually curate lists of significant pathways that are functionally or biologically similar
ER_BC <- c(
  "HOLLERN_SQUAMOUS_BREAST_TUMOR",
  "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN",
  "MCBRYAN_PUBERTAL_BREAST_4_5WK_UP",
  "LIM_MAMMARY_STEM_CELL_UP",
  "MCBRYAN_PUBERTAL_BREAST_6_7WK_DN",
  "SMID_BREAST_CANCER_BASAL_UP",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN",
  "SMID_BREAST_CANCER_NORMAL_LIKE_UP",
  "SMID_BREAST_CANCER_LUMINAL_B_DN"
)

ER_OtherCancer <- c(
  "LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN",
  "LINDGREN_BLADDER_CANCER_CLUSTER_2B",
  "CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_MODERATELY_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_DN",
  "LEI_MYB_TARGETS"
)

ER_RegulatoryAndSignaling <- c(
  "CHICAS_RB1_TARGETS_GROWING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "MARTINEZ_RB1_AND_TP53_TARGETS_DN"
)

ER_Misc <- c(
  "HP_PAIN_IN_HEAD_AND_NECK_REGION",
  "WP_HAIR_FOLLICLE_DEVELOPMENT_CYTODIFFERENTIATION_PART_3_OF_3",
  "GOMF_STRUCTURAL_MOLECULE_ACTIVITY",
  "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE",
  "REACTOME_KERATINIZATION",
  "MA_RAT_AGING_DN",
  "GAUSSMANN_MLL_AF4_FUSION_TARGETS_E_UP",
  "PAL_PRMT5_TARGETS_UP"
)

# Check to make sure all of the pathways are included

setdiff(c(ER_BC, ER_Misc, ER_OtherCancer, ER_RegulatoryAndSignaling), sig_GSEA_ER_METABRIC$pathway)
setdiff(c(ER_BC, ER_Misc, ER_OtherCancer, ER_RegulatoryAndSignaling), sig_GSEA_ER_TCGA$pathway)

# Repeat for TNBC

TNBC_BC <- c(
  "SMID_BREAST_CANCER_BASAL_DN",
  "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP",
  "DOANE_BREAST_CANCER_ESR1_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "SMID_BREAST_CANCER_ERBB2_UP",
  "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN",
  "VANTVEER_BREAST_CANCER_ESR1_DN",
  "TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN",
  "HOLLERN_SQUAMOUS_BREAST_TUMOR",
  "YANG_BREAST_CANCER_ESR1_DN",
  "SMID_BREAST_CANCER_LUMINAL_B_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BRAIN_UP",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN",
  "SMID_BREAST_CANCER_BASAL_UP"
)

TNBC_Immune <- c(
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_POSITIVE_REGULATION_OF_IMMUNE_RESPONSE",
  "GOBP_POSITIVE_REGULATION_OF_IMMUNE_SYSTEM_PROCESS",
  "HALLMARK_ALLOGRAFT_REJECTION",
  "REACTOME_MHC_CLASS_II_ANTIGEN_PRESENTATION",
  "BASSO_CD40_SIGNALING_UP",
  "WUNDER_INFLAMMATORY_RESPONSE_AND_CHOLESTEROL_UP",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "KEGG_MEDICUS_REFERENCE_ANTIGEN_PROCESSING_AND_PRESENTATION_BY_MHC_CLASS_II_MOLECULES",
  "KEGG_ALLOGRAFT_REJECTION",
  "MORI_IMMATURE_B_LYMPHOCYTE_DN",
  "WP_ALLOGRAFT_REJECTION",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
  
)

TNBC_OtherCancer <- c(
  "VECCHI_GASTRIC_CANCER_EARLY_UP",
  "NAKAYAMA_SOFT_TISSUE_TUMORS_PCA1_UP",
  "SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_DN",
  "LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN",
  "DEURIG_T_CELL_PROLYMPHOCYTIC_LEUKEMIA_DN",
  "SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6",
  "ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER",
  "WHITEFORD_PEDIATRIC_CANCER_MARKERS",
  "SABATES_COLORECTAL_ADENOMA_UP"
)

TNBC_DiseaseSignaling <- c(
  "THUM_SYSTOLIC_HEART_FAILURE_UP",
  "KEGG_ASTHMA",
  "KEGG_AUTOIMMUNE_THYROID_DISEASE",
  "KEGG_GRAFT_VERSUS_HOST_DISEASE",
  "KEGG_TYPE_I_DIABETES_MELLITUS",
  "KEGG_MEDICUS_PATHOGEN_HTLV_1_TAX_TO_NFY_MEDIATED_TRANSCRIPTION",
  "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",
  "KEGG_VIRAL_MYOCARDITIS"
  
)

TNBC_RegulatoryNetworks <- c(
  "MARSON_BOUND_BY_E2F4_UNSTIMULATED",
  "SANSOM_APC_TARGETS_DN",
  "FISCHER_DREAM_TARGETS",
  "HALLMARK_E2F_TARGETS"
)

TNBC_Misc <- c(
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",
  "GAURNIER_PSMD4_TARGETS",
  "MCLACHLAN_DENTAL_CARIES_UP",
  "WIELAND_UP_BY_HBV_INFECTION",
  "BENPORATH_ES_CORE_NINE_CORRELATED",
  "LEE_BMP2_TARGETS_DN",
  "GSE33292_WT_VS_TCF1_KO_DN3_THYMOCYTE_DN",
  "GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP",
  "MEBARKI_HCC_PROGENITOR_FZD8CRD_UP",
  "GSE3400_UNTREATED_VS_IFNB_TREATED_MEF_DN",
  "PUJANA_CHEK2_PCC_NETWORK",
  "KRAS.LUNG_UP.V1_DN"
)

setdiff(c(TNBC_BC, TNBC_DiseaseSignaling, TNBC_Immune, TNBC_Misc, TNBC_OtherCancer, TNBC_RegulatoryNetworks), sig_GSEA_TNBC_METABRIC$pathway)
setdiff(c(TNBC_BC, TNBC_DiseaseSignaling, TNBC_Immune, TNBC_Misc, TNBC_OtherCancer, TNBC_RegulatoryNetworks), sig_GSEA_TNBC_TCGA$pathway)

# Visualize for TNBC
setwd("<YOUR RESULTS PATH>")
for(pathways in c("TNBC_BC", "TNBC_DiseaseSignaling", "TNBC_Immune", "TNBC_Misc", "TNBC_OtherCancer", "TNBC_RegulatoryNetworks")){
  # Get the list of pathway names to plot
  ls <- get(pathways)
  
  # get the highest NES in the total dataset
  max <- max(abs(sig_GSEA_TNBC_METABRIC$NES))
  
  # Subset the significant GSEA results pathways to just be the ones in the group you are iterating through
  data <- sig_GSEA_TNBC_METABRIC[sig_GSEA_TNBC_METABRIC$pathway %in% ls,]
  
  # Get the -log10 FDR and format the pathway names for plotting
  data$minuslog10 <- -log10(data$padj)
  data$pathway <- gsub("_", " ", data$pathway)
  
  # Plot, ordering by NES
  p <- ggplot(data, aes(x = as.numeric(NES), y = fct_reorder(pathway, .x = as.numeric(NES), .desc = F), fill = NES, size = minuslog10))+
    geom_point(shape = 21)+
    # Set the x-axis limits to be large enough so the bubbles don't get clipped
    xlim(c(max*-1.5, max+1.5))+
    xlab("Normalized Enrichment Score")+
    ylab("Pathway")+
    ggtitle(as.character(pathways)) +  
    scale_fill_gradientn(colours = c("blue","white","red"), limits = c(max*-1, max))+
    labs(
      size = "-log10FDR",
      fill = "GSEA NES",
      x = "",
      y = "")+
    theme(legend.position = "bottom")
  
  # Add the plot for each functional grouping to the environment
  assign(paste0(pathways, "_plot"), p)
  
}

# Arrange the smaller plots into one figure
q <- ggarrange(`TNBC_BC_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_Immune_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_DiseaseSignaling_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_OtherCancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_RegulatoryNetworks_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_Misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(TNBC_BC), length(TNBC_Immune),
                   length(TNBC_DiseaseSignaling)+1,
                   length(TNBC_OtherCancer)+1,
                   length(TNBC_RegulatoryNetworks)+1,
                   length(TNBC_Misc)
                 ), widths = 12, font.label = list(size = 8) 
) 

ggsave("TNBC_GSEA_METABRIC.pdf", q, height = 20, width = 9)


# Repeat for ER+ data
setwd("<YOUR RESULTS PATH>")

for(pathways in c("ER_BC", "ER_OtherCancer", "ER_RegulatoryAndSignaling", "ER_Misc")){
  # Get the list of pathways to plot
  ls <- get(pathways)
  
  # Find the largest magnitude NES value in the data
  max <- max(abs(sig_GSEA_ER_METABRIC$NES))
  
  # Subset the GSEA data to just be the pathways in the group being iterated over
  data <- sig_GSEA_ER_METABRIC[sig_GSEA_ER_METABRIC$pathway %in% ls,]
  
  # Get the -log10FDR and format pathway names for plotting
  data$minuslog10 <- -log10(data$padj)
  data$pathway <- gsub("_", " ", data$pathway)
  
  # Plot, ordering by NES value
  p <- ggplot(data, aes(x = as.numeric(NES), y = fct_reorder(pathway, .x = as.numeric(NES), .desc = F), fill = NES, size = minuslog10))+
    geom_point(shape = 21)+
    xlim(c(max*-1.5, max*1.5))+
    xlab("Normalized Enrichment Score")+
    ylab("Pathway")+
    ggtitle(as.character(pathways)) +  
    scale_fill_gradientn(colours = c("blue","white","red"), limits = c(max*-1, max))+
    labs(
      size = "-log10FDR",
      fill = "GSEA NES",
      x = "",
      y = "")+
    theme(legend.position = "bottom")
  
  # Add the plot to the environment
  assign(paste0(pathways, "_plot"), p)
}

# Add the plots together into one larger figure
q <- ggarrange(`ER_BC_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_OtherCancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_RegulatoryAndSignaling_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_Misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(ER_BC),
                   length(ER_OtherCancer),
                   length(ER_RegulatoryAndSignaling),
                   length(ER_Misc)
                 ), widths = 4, font.label = list(size = 8) 
) 

# Save the plots
ggsave("ER_GSEA_METABRIC.pdf", q, height = 11)

# Visualize for TNBC
for(pathways in c("TNBC_BC", "TNBC_DiseaseSignaling", "TNBC_Immune", "TNBC_Misc", "TNBC_OtherCancer", "TNBC_RegulatoryNetworks")){
  # Get the list of pathway names to plot
  print(pathways)
  ls <- get(pathways)
  
  # get the highest NES in the total dataset
  max <- max(abs(sig_GSEA_TNBC_TCGA$NES))
  
  # Subset the significant GSEA results pathways to just be the ones in the group you are iterating through
  data <- sig_GSEA_TNBC_TCGA[sig_GSEA_TNBC_TCGA$pathway %in% ls,]
  
  # Get the -log10 FDR and format the pathway names for plotting
  data$minuslog10 <- -log10(data$padj)
  
  # These plots are ordered by METABRIC NES values
  pathway_order <- sig_GSEA_TNBC_METABRIC$pathway[sig_GSEA_TNBC_METABRIC$pathway %in% ls]
  pathway_order <- rev(pathway_order)
  
  data$pathway <- gsub("_", " ", data$pathway)
  pathway_order <- gsub("_", " ", pathway_order)
  pathway_order2 <- match(data$pathway, pathway_order)
  
  data$pathway <- factor(data$pathway, levels = pathway_order)
  
  # Plot, ordering by METABRIC NES
  p <- ggplot(data, aes(x = as.numeric(NES), y = pathway, fill = NES, size = minuslog10))+
    geom_point(shape = 21)+
    # Set the x-axis limits to be large enough so the bubbles don't get clipped
    xlim(c(max*-1.5, max*1.5))+
    xlab("Normalized Enrichment Score")+
    ylab("Pathway")+
    ggtitle(as.character(pathways)) +  
    scale_fill_gradientn(colours = c("blue","white","red"), limits = c(max*-1, max))+
    labs(
      size = "-log10FDR",
      fill = "GSEA NES",
      x = "",
      y = "")+
    theme(legend.position = "bottom")
  
  # Add the plot for each functional grouping to the environment
  assign(paste0(pathways, "_plot"), p, envir = .GlobalEnv)
}


# Arrange the smaller plots into one figure
q <- ggarrange(`TNBC_BC_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_Immune_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_DiseaseSignaling_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_OtherCancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_RegulatoryNetworks_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_Misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(TNBC_BC), length(TNBC_Immune),
                   length(TNBC_DiseaseSignaling)+1,
                   length(TNBC_OtherCancer)+1,
                   length(TNBC_RegulatoryNetworks)+1,
                   length(TNBC_Misc)
                 ), widths = 3, font.label = list(size = 8) 
) 

ggsave("TNBC_GSEA_TCGA.pdf", q, height = 20, width = 9)


# Repeat for ER+ data

for(pathways in c("ER_BC", "ER_OtherCancer", "ER_RegulatoryAndSignaling", "ER_Misc")){
  # Get the list of pathways to plot
  ls <- get(pathways)
  
  # Find the largest magnitude NES value in the data
  max <- max(abs(sig_GSEA_ER_TCGA$NES))
  
  # Subset the GSEA data to just be the pathways in the group being iterated over
  data <- sig_GSEA_ER_TCGA[sig_GSEA_ER_TCGA$pathway %in% ls,]
  
  # Get the -log10FDR and format pathway names for plotting
  data$minuslog10 <- -log10(data$padj)
  
  # The TCGA plots are ordered by the METABRIC NES values so pathwas are in the same order
  pathway_order <- sig_GSEA_ER_METABRIC$pathway[sig_GSEA_ER_METABRIC$pathway %in% ls]
  pathway_order <- rev(pathway_order)
  
  data$pathway <- gsub("_", " ", data$pathway)
  pathway_order <- gsub("_", " ", pathway_order)
  pathway_order2 <- match(data$pathway, pathway_order)
  
  data$pathway <- factor(data$pathway, levels = pathway_order)
  
  
  # Plot, ordering by METABRIC ES value
  p <- ggplot(data, aes(x = as.numeric(NES), y = pathway, fill = NES, size = minuslog10))+
    geom_point(shape = 21)+
    xlim(c(max*-1.5, max*1.5))+
    xlab("Normalized Enrichment Score")+
    ylab("Pathway")+
    ggtitle(as.character(pathways)) +  
    scale_fill_gradientn(colours = c("blue","white","red"), limits = c(max*-1, max))+
    labs(
      size = "-log10FDR",
      fill = "GSEA NES",
      x = "",
      y = "")+
    theme(legend.position = "bottom")
  
  # Add the plot to the environment
  assign(paste0(pathways, "_plot"), p)
}

# Add the plots together into one larger figure
q <- ggarrange(`ER_BC_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_OtherCancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_RegulatoryAndSignaling_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_Misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(ER_BC),
                   length(ER_OtherCancer),
                   length(ER_RegulatoryAndSignaling),
                   length(ER_Misc)
                 ), widths = 3, font.label = list(size = 8) 
) 

# Save the plots
ggsave("ER_GSEA_TCGA.pdf", q, height = 11)