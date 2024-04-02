# Get the DE genes between early and late recurrence donors in METABRIC
# for age groups and disease subtype

# First, load in the METABRIC clinical data and format
clin_metabric <- read.table("/Users/addie/Dropbox (Partners HealthCare)/Single-cell_breast_cancer/AnalysisAdrienne/Cibersort input data/brca_metabric/data_clinical_patient.txt",
                                 sep = "\t")
colnames(clin_metabric) <- clin_metabric[1,]
clin_metabric <- clin_metabric[-1,]

# Load in the clinical samples data to get stage and format
clin_sample <- read.table("/Users/addie/Dropbox (Partners HealthCare)/Single-cell_breast_cancer/AnalysisAdrienne/Cibersort input data/brca_metabric/data_clinical_sample.txt",
                          sep = "\t")
colnames(clin_sample) <- clin_sample[1,]
clin_sample <- clin_sample[-1,]

# Add columns to the clinical data to include tumor stage and make numerics numeric
clin_metabric$AGE_AT_DIAGNOSIS <- as.numeric(clin_metabric$AGE_AT_DIAGNOSIS)
clin_metabric$RFS_MONTHS <- as.numeric(clin_metabric$RFS_MONTHS)
clin_metabric$TUMOR_STAGE <- clin_sample$TUMOR_STAGE
clin_metabric$TUMOR_STAGE <- as.numeric(clin_metabric$TUMOR_STAGE)

# TNBC donors for each age group
TNBC_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                              clin_metabric$THREEGENE == "ER-/HER2-" &
                              clin_metabric$TUMOR_STAGE <= 2,]


TNBC_mid <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS >= 45 & 
                            clin_metabric$AGE_AT_DIAGNOSIS <= 65 &
                            clin_metabric$THREEGENE == "ER-/HER2-"&
                            clin_metabric$TUMOR_STAGE <= 2,]

TNBC_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                            clin_metabric$THREEGENE == "ER-/HER2-"&
                            clin_metabric$TUMOR_STAGE <= 2,]

TNBC_all <- clin_metabric[clin_metabric$THREEGENE == "ER-/HER2-"&
                            clin_metabric$TUMOR_STAGE <= 2,]

# ER+ low prolif donors for each age group
ER_lp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                              clin_metabric$THREEGENE == "ER+/HER2- Low Prolif",]
ER_lp_mid <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS >= 45 & 
                            clin_metabric$AGE_AT_DIAGNOSIS <= 65 &
                            clin_metabric$THREEGENE == "ER+/HER2- Low Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]
ER_lp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                            clin_metabric$THREEGENE == "ER+/HER2- Low Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]

# ER+ high prolif donors for each age group
ER_hp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                               clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                               clin_metabric$TUMOR_STAGE <= 2,]
ER_hp_mid <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS >= 45 & 
                             clin_metabric$AGE_AT_DIAGNOSIS <= 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]
ER_hp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]

ER_all <- clin_metabric[clin_metabric$THREEGENE == "ER+/HER2- High Prolif" |
                          clin_metabric$THREEGENE == "ER+/HER2- Low Prolif",]

ER_all <- ER_all[ER_all$TUMOR_STAGE <= 2,]

# Remove NA values from the subsetted clinical data so as not to include 
# donors with missing values
for(c in c("ER_hp_mid", "ER_hp_old", "ER_hp_young", "ER_lp_mid", "ER_lp_old", "ER_lp_young",
           "TNBC_mid", "TNBC_old", "TNBC_young", "TNBC_all", "ER_all")){
  df <- get(c)
  df <- na.omit(df)
  assign(c, df)
}

# We decided to group all ER+ donors together irrespective of proliferation status
ER_young <- rbind(ER_hp_young, ER_lp_young)
ER_mid <- rbind(ER_hp_mid, ER_lp_mid)
ER_old <- rbind(ER_hp_old, ER_lp_old)

# Read in and format the gene expression data for TNBC and ER+
# THESE DATA ARE TOO LARGE FOR GITHUB
TNBC <- read.csv("/Users/addie/Dropbox (Partners HealthCare)/Single-cell_breast_cancer/Data_METABRIC/2024_METABRIC_TNBC.csv")
rownames(TNBC) <- TNBC[,1]
TNBC <- TNBC[,-1]


ER <- read.csv("/Users/addie/Dropbox (Partners HealthCare)/Single-cell_breast_cancer/Data_METABRIC/2024_METABRIC_ER.csv")
rownames(ER) <- ER[,1]
ER <- ER[,-1]

# Get an annotation of gene names in the correct order for later
TNBC_genes <- rownames(TNBC)
ER_genes <- rownames(ER)

# Prep the matrix for DE gene analysis
TNBC3 <- exp(TNBC)
TNBC4 <- sapply(data.frame(TNBC3), as.integer)
rownames(TNBC4) <- TNBC_genes

# Prep the ER+ matrix for DE gene analysi
ER3 <- exp(ER)
ER4 <- sapply(data.frame(ER3), as.integer)
rownames(ER4) <- ER_genes

library(DESeq2)
# for the three TNBC age groups:
for(df in c("TNBC_old", "TNBC_mid", "TNBC_young", "TNBC_all")){
  print(df)
  dfi <- get(df)
  # Format the data frame for east of analysis (remove hyphens)
  dfi$PATIENT_ID <- make.names(dfi$PATIENT_ID)
  donors_early <- make.names(dfi$PATIENT_ID[dfi$RFS_MONTHS < 24 & dfi$RFS_STATUS == "1:Recurred"])
  donors_late <- make.names(dfi$PATIENT_ID[dfi$RFS_MONTHS > 120])
  
  print("number early recurrence")
  print(length(donors_early))
  print("mean age early")
  print(mean(as.numeric(dfi$AGE_AT_DIAGNOSIS[dfi$PATIENT_ID %in% donors_early])))
  print("number late recurrence")
  print(length(donors_late))
  print("mean age late")
  print(mean(as.numeric(dfi$AGE_AT_DIAGNOSIS[dfi$PATIENT_ID %in% donors_late])))
  print("mean recurrence time early")
  print(mean(as.numeric(dfi$RFS_MONTHS[dfi$PATIENT_ID %in% donors_early])))
  print("mean recurrence time late")
  print(mean(as.numeric(dfi$RFS_MONTHS[dfi$PATIENT_ID %in% donors_late])))
  
  
  # Make a data frame of the formatted gene expression, but only donors
  # within the early or late recurrence groups
  df2 <- data.frame(cbind(TNBC4[,colnames(TNBC4) %in% donors_early],
              TNBC4[,colnames(TNBC4) %in% donors_late]))
  
  # Make a vector of annotations for the data
  annots <- data.frame(c(rep("early", length(data.frame(TNBC4)[,colnames(TNBC4) %in% donors_early])),
                       rep("late", length(data.frame(TNBC4)[,colnames(TNBC4) %in% donors_late]))))
  colnames(annots) <- "annots"
  
  # Run DESeq2
  y <- DESeqDataSetFromMatrix(countData = df2, colData = annots, design = ~annots
  )
  y <- DESeq(y)
  
  # Get the results and filter to just significant geness
  res <- data.frame(results(y))
  #rownames(res) <- rownames(top.table)
  res2 <- subset(res, subset = res$padj < 0.05)
  assign(paste0(df, "_res"), res2)
}

# Repeat for ER+
for(df in c("ER_old", "ER_mid", "ER_young", "ER_all")){
  print(df)
  dfi <- get(df)
  dfi$PATIENT_ID <- make.names(dfi$PATIENT_ID)
  donors_early <- make.names(dfi$PATIENT_ID[dfi$RFS_MONTHS < 24 & dfi$RFS_STATUS == "1:Recurred"])
  donors_late <- make.names(dfi$PATIENT_ID[dfi$RFS_MONTHS > 120])
  
  print("number early recurrence")
  print(length(donors_early))
  print("mean age early")
  print(mean(as.numeric(dfi$AGE_AT_DIAGNOSIS[dfi$PATIENT_ID %in% donors_early])))
  print("number late recurrence")
  print(length(donors_late))
  print("mean age late")
  print(mean(as.numeric(dfi$AGE_AT_DIAGNOSIS[dfi$PATIENT_ID %in% donors_late])))
  print("mean recurrence time early")
  print(mean(as.numeric(dfi$RFS_MONTHS[dfi$PATIENT_ID %in% donors_early])))
  print("mean recurrence time late")
  print(mean(as.numeric(dfi$RFS_MONTHS[dfi$PATIENT_ID %in% donors_late])))
  
  #assign(df, donors_late)
  
  df2 <- data.frame(cbind(ER4[,colnames(ER4) %in% donors_early],
                          ER4[,colnames(ER4) %in% donors_late]))
  
  annots <- data.frame(c(rep("early", length(data.frame(ER4)[,colnames(ER4) %in% donors_early])),
                         rep("late", length(data.frame(ER4)[,colnames(ER4) %in% donors_late]))))
  colnames(annots) <- "annots"
  
  df3 <- na.omit(df2)
  
  y <- DESeqDataSetFromMatrix(countData = df3, colData = annots, design = ~annots
  )
  y <- DESeq(y)
  
  res <- data.frame(results(y))
  #rownames(res) <- rownames(top.table)
  res2 <- subset(res, subset = res$padj < 0.05)
  assign(paste0(df, "_res"), res2)
}


#################################################################################
# Plotting
setwd("/Users/addie/desktop/Figure1")

library(ggplot2)
library(forcats)
for(df in c("ER_old_res", "ER_young_res", "ER_mid_res", "TNBC_old_res", 
            "TNBC_mid_res", "TNBC_young_res", "TNBC_all_res", "ER_all_res")){
  print(df)
  dfi <- get(df)
  
  dfi$abs <- abs(dfi$log2FoldChange)
  dfi$gene <- rownames(dfi)
  
  # Only select the top 20 most differentially expressed genes by abs(logfc)
  dd <- dfi[order(dfi$abs, decreasing = T),]
  if(nrow(dd) > 20){
    dd <- dd[1:20,]
  }
  dd$logFC <- as.numeric(dd$log2FoldChange)
  
  dd$direction <- NA
  dd$direction[dd$logFC < 0] <- "enriched in early"
  dd$direction[dd$logFC > 0] <- "enriched in late"
  
  
  plot <- ggplot(dd, aes(x = logFC, y = fct_reorder(gene, logFC), fill = direction))+
    geom_bar(stat = "identity")+
    xlim(c(max(abs(as.numeric(dd$logFC)))*-1, max(abs(as.numeric(dd$logFC)))))+
    scale_fill_manual(values = c("turquoise4", "red2"))+
    ylab("Gene")+theme(aspect.ratio = 1.75)+xlab("Log2 Fold Change")
  
  assign(paste0(df, "_plot"), plot)
  ggsave(paste0(df, "_plot.pdf"), plot, device = "pdf", height = nrow(dd)/5)
}

# Get the lists of significant genes for the UpSet plot
eryoung <- rownames(ER_young_res)
ermid <- rownames(ER_mid_res)
erold <- rownames(ER_old_res)
erall <- rownames(ER_all_res)
tnbcyoung <- rownames(TNBC_young_res)
tnbcmid <- rownames(TNBC_mid_res)
tnbcold <- rownames(TNBC_old_res)
tnbcall <- rownames(TNBC_all_res)


library(UpSetR)
# Dataset
listinput_TNBC <- list(
                  TNBC_Young = tnbcyoung, TNBC_Mid = tnbcmid, TNBC_Old = tnbcold,
                  TNBC_All = tnbcall)

listinput_ER <- list(ER_Young = eryoung, ER_Mid = ermid, ER_Old = erold,
                       ER_All = erall)

# Plot
p <- upset(fromList(listinput_TNBC), 
      nintersects = 40,
      nsets = 4, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.6, 0.4),
      number.angles = 0, 
      text.scale = 1.1, 
      point.size = 2.8, 
      line.size = 1
)

q <- upset(fromList(listinput_ER), 
           nintersects = 40,
           nsets = 4, 
           order.by = "freq", 
           decreasing = T, 
           mb.ratio = c(0.6, 0.4),
           number.angles = 0, 
           text.scale = 1.1, 
           point.size = 2.8, 
           line.size = 1
)

dev.off()

pdf("UpSet.pdf", p)

dev.off()

#Kaplan-Meier Plots for TNBC and ER+
library(ggsurvfit)
library(survival)

TNBC_young$AGE_GROUP <- "Young"
TNBC_mid$AGE_GROUP <- "Middle-aged"
TNBC_old$AGE_GROUP <- "Old"

# Merge the TNBC young, mid, and old data
TNBC2 <- rbind(TNBC_young, TNBC_mid)
TNBC2 <- rbind(TNBC2, TNBC_old)

# Subset to donors with early or late recurrence and format the annotations
TNBC2$RECURRENCE <- NA
TNBC2$RECURRENCE[TNBC2$RFS_MONTHS < 24 & TNBC2$RFS_STATUS == "1:Recurred"] <- "Early"
TNBC2$RECURRENCE[TNBC2$RFS_MONTHS > 120] <- "Late"
TNBC2$RFS_STATUS[TNBC2$RFS_STATUS == "1:Recurred"] <- 1
TNBC2$RFS_STATUS[TNBC2$RFS_STATUS == "0:Not Recurred"] <- 0
TNBC2$RFS_STATUS <- as.numeric(TNBC2$RFS_STATUS)

# Add a new column that generates each combination of age group and recurrence
# (6 grounps)
TNBC2$OVERALL <- paste(TNBC2$AGE_GROUP, TNBC2$RECURRENCE)

# Remove any NAs if necessary
TNBC2 <- na.omit(TNBC2)

# Plot the survival curve
TNBC_surv <- survfit2(Surv(RFS_MONTHS, RFS_STATUS) ~ c(OVERALL), data = TNBC2) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "RFS Probability"
  )+theme(aspect.ratio = 1)

#Merge the ER+ data
ER_young$AGE_GROUP <- "Young"
ER_mid$AGE_GROUP <- "Middle-aged"
ER_old$AGE_GROUP <- "Old"

ER2 <- rbind(ER_young, ER_mid)
ER2 <- rbind(ER2, ER_old)

# Subset to early and late recurrence ER+ donors
ER2$RECURRENCE <- NA
ER2$RECURRENCE[ER2$RFS_MONTHS < 24 & ER2$RFS_STATUS == "1:Recurred"] <- "Early"
ER2$RECURRENCE[ER2$RFS_MONTHS > 120] <- "Late"
ER2$RFS_STATUS[ER2$RFS_STATUS == "1:Recurred"] <- 1
ER2$RFS_STATUS[ER2$RFS_STATUS == "0:Not Recurred"] <- 0
ER2$RFS_STATUS <- as.numeric(ER2$RFS_STATUS)

# Add a new column that generates each combination of age group and recurrence
# (6 grounps)
ER2$OVERALL <- paste(ER2$AGE_GROUP, ER2$RECURRENCE)

# Remove NAs if necessary
ER2 <- na.omit(ER2)

# Plot

ER_surv <- survfit2(Surv(RFS_MONTHS, RFS_STATUS) ~ c(OVERALL), data = ER2) %>% 
  ggsurvfit() +
  labs(
    x = "Months",
    y = "RFS Probability"
  )+theme(aspect.ratio = 1)+
  xlim(0, 280)

ggsave("TNBC_surv.pdf", TNBC_surv, device = "pdf")
ggsave("ER_surv.pdf", ER_surv, device = "pdf")
