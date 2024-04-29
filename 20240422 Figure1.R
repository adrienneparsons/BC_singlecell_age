# Get the DE genes between early and late recurrence donors and old and young in METABRIC
# for disease subtype

# Prerequisite packages
library(ggplot2)
library(forcats)
library(edgeR)
library(limma)
library(Biobase)
library(fgsea)
library(ggpubr)

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

TNBC_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                            clin_metabric$THREEGENE == "ER-/HER2-"&
                            clin_metabric$TUMOR_STAGE <= 2,]

# ER+ low prolif donors for each age group
ER_lp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                               clin_metabric$THREEGENE == "ER+/HER2- Low Prolif",]

ER_lp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- Low Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]

# ER+ high prolif donors for each age group
ER_hp_young <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS < 45 & 
                               clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                               clin_metabric$TUMOR_STAGE <= 2,]

ER_hp_old <- clin_metabric[clin_metabric$AGE_AT_DIAGNOSIS > 65 &
                             clin_metabric$THREEGENE == "ER+/HER2- High Prolif"&
                             clin_metabric$TUMOR_STAGE <= 2,]

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


# Start DE gene analysis

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
df3 <- df2[rownames(df2) %in% names(var[23500:length(var)]),]

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

# Repeat for TNBC. First, filter the gene expression data to just be the oldest and youngest TNBC
# donors
TNBC_young2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_young$PATIENT_ID)]
TNBC_old2 <- TNBC[,make.names(colnames(TNBC)) %in% make.names(TNBC_old$PATIENT_ID)]
df2 <- cbind(TNBC_young2, TNBC_old2)

# Make that data an ExpressionSet
aseset <- ExpressionSet(assayData = as.matrix(df2))

# Find the genes with the highest spread of standard deviation exprression
var <- apply(exprs(aseset), 1, sd)
var <- var[order(var)]

# plot
plot(1:length(var), var)

# Subset the data to just be the 675 genes with the highest spread of expression
df3 <- df2[rownames(df2) %in% names(var[23500:length(var)]),]

# Exponeentiate the data to resemble raw counts (cBioPortal log2-normalizes data)
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

# Plotting
setwd("/Users/addie/desktop")

# For the two results tables:

for(df in c("top.table_ER", "top.table_TNBC")){
  print(df)
  dfi <- get(df)
  
  # Copy the data for volcano plot
  dd3 <- dfi
  dd3$gene <- rownames(dd3)
  dd3$abs <- abs(dd3$logFC)
  
  # Subset the data to just the genes with low adjusted p-values and print number of genes in each direction
  df_sub <- dfi[dfi$adj.P.Val < 0.15,]
  print("logfc < 0")
  print(nrow(df_sub[df_sub$logFC < 0,]))
  print("logfc > 0")
  print(nrow(df_sub[df_sub$logFC > 0,]))
  
  df_sub$abs <- abs(df_sub$logFC)
  df_sub$gene <- rownames(df_sub)
  
  # Only select the top 20 most differentially expressed genes by abs(logfc), selecting up to 20 enriched in
  # either direction
  dd <- df_sub[order(df_sub$abs, decreasing = T),]
  
  dd_up <- dd[dd$logFC > 0,]
  dd_down <- dd[dd$logFC < 0,]
  
  #Get up to the top 20 most differentially expressed in young or older
  if(nrow(dd_up) > 20){
    dd_up <- dd_up[1:20,]
  }
  
  if(nrow(dd_down) > 20){
    dd_down <- dd_down[1:20,]
  }
  
  dd2 <- rbind(dd_up, dd_down)
  
  # Set up for bar chart plotting
  dd2$logFC <- as.numeric(dd2$logFC)
  
  # Label enrichment direction by whether log fold change > or < 0
  dd2$direction <- NA
  dd2$direction[dd2$logFC > 0] <- "enriched in older"
  dd2$direction[dd2$logFC < 0] <- "enriched in younger"
  
  # Plot the barr charts of most DE genes in both directions
  plot <- ggplot(dd2, aes(x = logFC, y = fct_reorder(gene, logFC), fill = direction))+
    geom_bar(stat = "identity")+
    xlim(c(max(abs(as.numeric(dd2$logFC)))*-1, max(abs(as.numeric(dd2$logFC)))))+
    scale_fill_manual(values = c("turquoise4", "red2"))+
    ylab("Gene")+theme(aspect.ratio = 3)+xlab("Log2 Fold Change")+
    theme(text = element_text(size = 8))
  
  assign(paste0(df, "_plot"), plot)
  ggsave(paste0(df, "_plot.pdf"), plot, device = "pdf", height = nrow(dd2)/12)
  
  
  # Generate a Volcano plot based on the limma analysis differentially
  # abundant genes
  dd3$minuslog10p <- -log10(dd3$adj.P.Val)
  dd3$direction <- NA
  
  # Only add gene names to the most differentially expressed genes
  dd3$label <- NA
  
  # If analyzing TNBC, add the gene names to any gene with adjusted p < 0.15
  # If analyzing ER+, add the gene names to any gene with adjusted p < 0.05 and magnitude fold change > 0.5
  if(df == "top.table_TNBC"){
    dd3$label[dd3$minuslog10p > -log10(0.15)] <- dd3$gene[dd3$minuslog10p > -log10(0.15)]
  } else(dd3$label[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 0.5] <- dd3$gene[dd3$minuslog10p > -log10(0.05) & abs(dd3$logFC) > 0.5])
  
  
  # Add in enrichment annotations for coloring the dots
  dd3$direction <- "enriched in older"
  dd3$direction[dd3$logFC < 0] <- "enriched in younger"
  dd3$direction[is.na(dd3$label)] <- "0"
  
  # Write a table of the genes
  dd4 <- arrange(dd3, desc(minuslog10p), desc(logFC))
  dd5 <- dd4[,c("logFC", "minuslog10p", "gene")]
  
  write.table(dd5, file = paste0(gsub("top.table_", "", df), "_volcanocoords.tsv"), quote = F, sep = "\t", row.names = F)
  
  # Make the volcano plot and save
  volcano <- ggplot(data=dd3, aes(x=logFC, y=minuslog10p, label = label, fill = direction)) + 
    geom_point(size = 3, shape = 21)+
    geom_text(vjust = -0.5, hjust = 1, size = 1)+
    xlim(c(max(dd3$abs)*-1-.25, max(dd3$abs)+.25))+theme_bw()+
    scale_fill_manual(values = c("black", "red", "blue"))+
    ylim(c(0, max(dd3$minuslog10p) + 0.5))+
    xlab("Log2 Fold Change Difference")+
    ylab("-log10(FDR)")
  
  ggsave(paste0(df, "_volcano.pdf"), volcano, device = "pdf", height = 3, width = 8)
}

# Running GSEA on results
setwd("/Users/addie/desktop/GSEA/gmts")
gmts <- list.files("/Users/addie/desktop/GSEA/gmts")

# For the ER+ and TNBC results:
for(res in c("top.table_ER", "top.table_TNBC")){
  
  # Get the data and extract the logFC for ranking
  results <- get(res)
  genes <- results$logFC
  names(genes) <- rownames(results)
  
  # Set up the GSEA results data frame and remove the row of NAs
  gsea_data <- data.frame("pathway" = NA, "pval" = NA, "padj" = NA, "log2err" = NA,
                          "ES" = NA, "NES" = NA, "size" = NA, "leadingEdge" = NA, "gmt" = NA)
  gsea_data <- gsea_data[-1,]
  

  # For each gmt file (different gene set)
  # Run gsea using fgsea and the ranked gene list
  # if there are any pathways with adjusted p-values < 0.05 add it to the final data
  for(gmt in gmts){
    print(gmt)
    myGO <- fgsea::gmtPathways(gmt)
    gsea_out <- fgsea(pathways = myGO,
                      stats = genes)
    gsea_out2 <- gsea_out[gsea_out$padj < 0.05]
    gsea_out2$gmt <- gmt
    gsea_data <- rbind(gsea_data, gsea_out2)
  }
  
  # Compile all results and write to the environment
  group <- gsub("top.table_", "", res)
  
  assign(paste0("sig_gsea_", group), gsea_data)
}

# Order the significant TNBC GSEA results by decreasing NES and write to computer
sig_gsea_TNBC <- sig_gsea_TNBC[order(sig_gsea_TNBC$NES, decreasing = T),]
sig_gsea_TNBC2 <- apply(sig_gsea_TNBC,2,as.character)
write.csv(sig_gsea_TNBC2, file = "20240424_TNBC_GSEA_FINAL.csv", quote = F, row.names = F)

# Repeat for ER+
sig_gsea_ER <- sig_gsea_ER[order(sig_gsea_ER$NES, decreasing = T),]
sig_gsea_ER2 <- apply(sig_gsea_ER,2,as.character)
write.csv(sig_gsea_ER2, file = "20240424_ER_GSEA_FINAL.csv", quote = F, row.names = F)

# Manually curate lists of significant pathways that are functionally or biologically similar
TNBC_BC <- c(
  "SMID_BREAST_CANCER_BASAL_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_UP",
  "SMID_BREAST_CANCER_RELAPSE_IN_LUNG_DN",
  "SMID_BREAST_CANCER_ERBB2_UP",
  "SMID_BREAST_CANCER_LUMINAL_B_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "SMID_BREAST_CANCER_BASAL_UP",
  "LIEN_BREAST_CARCINOMA_METAPLASTIC_VS_DUCTAL_DN",
  "VANTVEER_BREAST_CANCER_ESR1_DN",
  "DOANE_BREAST_CANCER_ESR1_DN",
  "DOANE_BREAST_CANCER_ESR1_UP",
  "SOTIRIOU_BREAST_CANCER_GRADE_1_VS_3_UP",
  "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN",
  "LIM_MAMMARY_LUMINAL_MATURE_UP"
)

TNBC_immune <- c("BASSO_CD40_SIGNALING_UP",
  "RUTELLA_RESPONSE_TO_HGF_VS_CSF2RB_AND_IL4_UP",
  "REACTOME_INTERFERON_GAMMA_SIGNALING",
  "WUNDER_INFLAMMATORY_RESPONSE_AND_CHOLESTEROL_UP",
  "CHEMNITZ_RESPONSE_TO_PROSTAGLANDIN_E2_UP",
  "BIOCARTA_MHC_PATHWAY",
  "NAKAJIMA_EOSINOPHIL",
  "GAVIN_FOXP3_TARGETS_CLUSTER_P6",
  "MORI_IMMATURE_B_LYMPHOCYTE_DN",
  "GSE15750_DAY6_VS_DAY10_TRAF6KO_EFF_CD8_TCELL_UP",
  "GSE15750_DAY6_VS_DAY10_EFF_CD8_TCELL_UP",
  "GSE21927_SPLEEN_C57BL6_VS_4T1_TUMOR_BALBC_MONOCYTES_DN",
  "GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_YOUNG_DN",
  "GSE39110_DAY3_VS_DAY6_POST_IMMUNIZATION_CD8_TCELL_DN",
  "GSE36476_CTRL_VS_TSST_ACT_72H_MEMORY_CD4_TCELL_OLD_DN",
  "GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_YOUNG_DN"
)


TNBC_cellcycle <- c(
  "HALLMARK_G2M_CHECKPOINT",
  "HALLMARK_E2F_TARGETS",
  "FISCHER_G2_M_CELL_CYCLE",
  "CHEN_METABOLIC_SYNDROM_NETWORK",
  "FISCHER_DREAM_TARGETS"
)


TNBC_othercancer <- c(
  "ROSTY_CERVICAL_CANCER_PROLIFERATION_CLUSTER",
  "NAKAYAMA_SOFT_TISSUE_TUMORS_PCA1_UP",
  "HORIUCHI_WTAP_TARGETS_DN",
  "SHEDDEN_LUNG_CANCER_POOR_SURVIVAL_A6",
  "CAIRO_HEPATOBLASTOMA_CLASSES_UP",
  "WHITEFORD_PEDIATRIC_CANCER_MARKERS",
  "LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN",
  "SENGUPTA_NASOPHARYNGEAL_CARCINOMA_WITH_LMP1_DN",
  "VECCHI_GASTRIC_CANCER_EARLY_UP",
  "RODRIGUES_THYROID_CARCINOMA_POORLY_DIFFERENTIATED_UP"
)

TNBC_misc <- c(
  "THUM_SYSTOLIC_HEART_FAILURE_UP",
  "GAURNIER_PSMD4_TARGETS",
  "KEGG_SYSTEMIC_LUPUS_ERYTHEMATOSUS",
  "MCLACHLAN_DENTAL_CARIES_UP",
  "JOHNSTONE_PARVB_TARGETS_3_DN",
  "LEE_BMP2_TARGETS_DN",
  "PUJANA_CHEK2_PCC_NETWORK",
  "BASAKI_YBX1_TARGETS_UP",
  "GSE33292_WT_VS_TCF1_KO_DN3_THYMOCYTE_DN",
  "GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP",
  "MEBARKI_HCC_PROGENITOR_FZD8CRD_UP",
  "MCBRYAN_PUBERTAL_BREAST_6_7WK_DN",
  "GAVIN_FOXP3_TARGETS_CLUSTER_P6",
  "MARSON_BOUND_BY_E2F4_UNSTIMULATED",
  "BENPORATH_ES_CORE_NINE_CORRELATED",
  "BLANCO_MELO_BRONCHIAL_EPITHELIAL_CELLS_INFLUENZA_A_DEL_NS1_INFECTION_DN"
)

ER_BC <- c(
           "SMID_BREAST_CANCER_BASAL_UP",
           "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
           "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN",
           "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN",
           "SMID_BREAST_CANCER_NORMAL_LIKE_UP",
           "SMID_BREAST_CANCER_LUMINAL_B_DN",
           "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN",
           "HOLLERN_SQUAMOUS_BREAST_TUMOR"
  
)

ER_othercancer <- c(
  "LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_MODERATELY_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_DN",
  "CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_DN",
  "WEST_ADRENOCORTICAL_TUMOR_DN",
  "LINDGREN_BLADDER_CANCER_CLUSTER_2B"
)

ER_developmentdiff <- c(
  "LEI_MYB_TARGETS",
  "LIM_MAMMARY_STEM_CELL_UP",
  "BOQUEST_STEM_CELL_CULTURED_VS_FRESH_DN"
)

ER_misc <- c(
  "MA_RAT_AGING_DN",
  "NAGASHIMA_NRG1_SIGNALING_UP",
  "PID_REG_GR_PATHWAY",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "REACTOME_CELLULAR_RESPONSES_TO_STIMULI",
  "GOMF_STRUCTURAL_MOLECULE_ACTIVITY",
  "REACTOME_FORMATION_OF_THE_CORNIFIED_ENVELOPE",
  "REACTOME_KERATINIZATION"
)

# Visualize for TNBC
setwd("/Users/addie/desktop/GSEAbubblecharts/TNBC")
for(pathways in c("TNBC_BC", "TNBC_cellcycle", "TNBC_othercancer", "TNBC_immune", "TNBC_misc")){
  # Get the list of pathway names to plot
  ls <- get(pathways)
  
  # get the highest NES in the total dataset
  max <- max(abs(sig_gsea_TNBC$NES))
  
  # Subset the significant GSEA results pathways to just be the ones in the group you are iterating through
  data <- sig_gsea_TNBC[sig_gsea_TNBC$pathway %in% ls,]
  
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
               `TNBC_immune_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_cellcycle_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_othercancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `TNBC_misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(TNBC_BC), length(TNBC_immune),
                   length(TNBC_cellcycle)+3,
                   length(TNBC_othercancer),
                   length(TNBC_misc)
                 ), widths = 12, font.label = list(size = 8) 
) 

ggsave("TNBC_GSEA.pdf", q, height = 15)


# Repeat for ER+ data
setwd("/Users/addie/desktop/GSEAbarcharts/ER")

for(pathways in c("ER_BC", "ER_developmentdiff", "ER_othercancer", "ER_misc")){
  # Get the list of pathways to plot
  ls <- get(pathways)
  
  # Find the largest magnitude NES value in the data
  max <- max(abs(sig_gsea_ER$NES))
  
  # Subset the GSEA data to just be the pathways in the group being iterated over
  data <- sig_gsea_ER[sig_gsea_ER$pathway %in% ls,]
  
  # Get the -log10FDR and format pathway names for plotting
  data$minuslog10 <- -log10(data$padj)
  data$pathway <- gsub("_", " ", data$pathway)
  
  # Plot, ordering by NES value
  p <- ggplot(data, aes(x = as.numeric(NES), y = fct_reorder(pathway, .x = as.numeric(NES), .desc = F), fill = NES, size = minuslog10))+
    geom_point(shape = 21)+
    xlim(c(max*-1, max))+
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
               `ER_developmentdiff_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_othercancer_plot` + 
                 theme(axis.text.x = element_blank(),
                       axis.ticks.x = element_blank(),
                       axis.title.x = element_blank() ),
               `ER_misc_plot`, ncol = 1, common.legend = T, align = "v", heights = 
                 c(length(ER_BC), length(ER_developmentdiff),
                   length(ER_othercancer),
                   length(ER_misc)
                 ), widths = 4, font.label = list(size = 8) 
) 

# Save the plots
ggsave("ER_GSEA.pdf", q, height = 12)

