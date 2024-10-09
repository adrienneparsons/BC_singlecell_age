# Figure 1: DE gene analysis and GSEA
###########################################################
# Adrienne Parsons, 2024-05-07

# Get the DE genes between early and late recurrence donors and old and young in METABRIC
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
# THESE DATA ARE TOO LARGE FOR GITHUB; these are the resultant files from "METABRIC data wrangling.R"
TNBC <- read.csv("<YOUR DATA PATH>/2024_METABRIC_TNBC.csv")
rownames(TNBC) <- TNBC[,1]
TNBC <- TNBC[,-1]


ER <- read.csv("<YOUR DATA PATH>/2024_METABRIC_ER.csv")
rownames(ER) <- ER[,1]
ER <- ER[,-1]

# Get an annotation of gene names in the correct order
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
df3 <- df2[rownames(df2) %in% names(var[23500:length(var)]),]

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
setwd("<YOUR DATA PATH>")

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
  
  write.table(dd5, file = paste0(gsub("top.table_", "", df), "_volcanocoords.tsv"), quote = F, sep = "\t", row.names = F)
  
  # Make the volcano plot and save
  volcano <- ggplot(data=dd3, aes(x=logFC, y=minuslog10p, label = label, fill = direction)) + 
    geom_point(size = 3, shape = 21, stroke = 0)+
    xlim(c(max(dd3$abs)*-1-.25, max(dd3$abs)+.25))+theme_bw()+
    scale_fill_manual(values = c("lightgrey", "red", "blue"))+
    ylim(c(0, max(dd3$minuslog10p) + 0.5))+
    xlab("Log2 Fold Change")+
    ylab("-log10(FDR)")+
    geom_text_repel(max.overlaps = 200, size = 3)+
    theme(aspect.ratio = 1)
  
  ggsave(paste0(df, "_volcano2.pdf"), volcano, device = "pdf", height = 6, width = 12)
}

# Running GSEA on results; gmt files obtained from this link: https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp
setwd("<YOUR DATA PATH>/gmts")
gmts <- list.files("<YOUR DATA PATH>/gmts")
gmts <- gmts[-grep(".csv", gmts)]

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
    set.seed(42)
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


setwd("<YOUR DATA PATH>")
# Order the significant TNBC GSEA results by decreasing NES and write to computer
sig_gsea_TNBC <- sig_gsea_TNBC[order(sig_gsea_TNBC$NES, decreasing = T),]
sig_gsea_TNBC <- sig_gsea_TNBC[,-c("size", "leadingEdge")]
sig_gsea_TNBC2 <- apply(sig_gsea_TNBC,2,as.character)
write.csv(sig_gsea_TNBC2, file = "TNBC_GSEA_inclStageIII.csv", quote = F, row.names = F)

# Repeat for ER+
sig_gsea_ER <- sig_gsea_ER[order(sig_gsea_ER$NES, decreasing = T),]
sig_gsea_ER <- sig_gsea_ER[,-c("size", "leadingEdge")]
sig_gsea_ER2 <- apply(sig_gsea_ER,2,as.character)
write.csv(sig_gsea_ER2, file = "ER_GSEA_inclStageIII.csv", quote = F, row.names = F)

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
  "SMID_BREAST_CANCER_LUMINAL_B_DN",
  "SMID_BREAST_CANCER_BASAL_UP"
)

ER_OtherCancer <- c(
  "LIU_OVARIAN_CANCER_TUMORS_AND_XENOGRAFTS_XDGS_DN",
  "LINDGREN_BLADDER_CANCER_CLUSTER_2B",
  "WEST_ADRENOCORTICAL_TUMOR_DN",
  "CHIANG_LIVER_CANCER_SUBCLASS_CTNNB1_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_MODERATELY_DN",
  "RICKMAN_TUMOR_DIFFERENTIATED_WELL_VS_POORLY_DN",
  "LEI_MYB_TARGETS"
)

ER_RegulatoryAndSignaling <- c(
  "PID_REG_GR_PATHWAY",
  "CHICAS_RB1_TARGETS_GROWING",
  "HALLMARK_TNFA_SIGNALING_VIA_NFKB",
  "NAGASHIMA_NRG1_SIGNALING_UP",
  "MARTINEZ_RB1_AND_TP53_TARGETS_DN",
  "REACTOME_CELLULAR_RESPONSES_TO_STIMULI"
)

ER_Misc <- c(
  "BOQUEST_STEM_CELL_CULTURED_VS_FRESH_DN",
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
setdiff(sig_gsea_ER$pathway, c(ER_BC, ER_Misc, ER_OtherCancer, ER_RegulatoryAndSignaling))
setdiff(c(ER_BC, ER_Misc, ER_OtherCancer, ER_RegulatoryAndSignaling), sig_gsea_ER$pathway)

# Repeat for TNBC

TNBC_BC <- c(
  "SMID_BREAST_CANCER_BASAL_DN",
  "SCHUETZ_BREAST_CANCER_DUCTAL_INVASIVE_UP",
  "DOANE_BREAST_CANCER_ESR1_UP",
  "LIEN_BREAST_CARCINOMA_METAPLASTIC_VS_DUCTAL_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_UP",
  "SMID_BREAST_CANCER_ERBB2_UP",
  "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_UP",
  "SMID_BREAST_CANCER_RELAPSE_IN_LUNG_DN",
  "SMID_BREAST_CANCER_LUMINAL_A_UP",
  "VANTVEER_BREAST_CANCER_ESR1_UP",
  "DOANE_BREAST_CANCER_CLASSES_UP",
  "TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN",
  "HOLLERN_SQUAMOUS_BREAST_TUMOR",
  "YANG_BREAST_CANCER_ESR1_DN",
  "SMID_BREAST_CANCER_LUMINAL_B_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BRAIN_UP",
  "LIM_MAMMARY_LUMINAL_MATURE_UP",
  "LIM_MAMMARY_LUMINAL_MATURE_DN",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_DUCTAL_NORMAL_DN",
  "CHARAFE_BREAST_CANCER_LUMINAL_VS_BASAL_DN",
  "TURASHVILI_BREAST_DUCTAL_CARCINOMA_VS_LOBULAR_NORMAL_DN",
  "TURASHVILI_BREAST_LOBULAR_CARCINOMA_VS_DUCTAL_NORMAL_DN",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "SMID_BREAST_CANCER_BASAL_UP",
  "SMID_BREAST_CANCER_RELAPSE_IN_BONE_DN",
  "DOANE_BREAST_CANCER_ESR1_DN",
  "VANTVEER_BREAST_CANCER_ESR1_DN"
)

TNBC_Immune <- c(
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_PEPTIDE_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION",
  "GOBP_ANTIGEN_PROCESSING_AND_PRESENTATION_OF_EXOGENOUS_PEPTIDE_ANTIGEN_VIA_MHC_CLASS_II",
  "GOCC_MHC_CLASS_II_PROTEIN_COMPLEX",
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
  "GSE25123_WT_VS_PPARG_KO_MACROPHAGE_DN",
  "GOMF_MHC_CLASS_II_PROTEIN_COMPLEX_BINDING",
  "GSE36476_CTRL_VS_TSST_ACT_40H_MEMORY_CD4_TCELL_YOUNG_DN",
  "KEGG_ANTIGEN_PROCESSING_AND_PRESENTATION"
  
)

TNBC_OtherCancer <- c(
  "VECCHI_GASTRIC_CANCER_EARLY_UP",
  "DURCHDEWALD_SKIN_CARCINOGENESIS_DN",
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
  "FISCHER_G2_M_CELL_CYCLE",
  "MARSON_BOUND_BY_E2F4_UNSTIMULATED",
  "SANSOM_APC_TARGETS_DN",
  "KOINUMA_TARGETS_OF_SMAD2_OR_SMAD3",
  "FISCHER_DREAM_TARGETS",
  "MTOR_UP.V1_UP",
  "HALLMARK_E2F_TARGETS"
)

TNBC_Misc <- c(
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",
  "JAATINEN_HEMATOPOIETIC_STEM_CELL_DN",
  "GAURNIER_PSMD4_TARGETS",
  "MCLACHLAN_DENTAL_CARIES_UP",
  "WIELAND_UP_BY_HBV_INFECTION",
  "HP_EXCESSIVE_DAYTIME_SOMNOLENCE",
  "HP_PROMINENT_NASAL_BRIDGE",
  "BENPORATH_ES_CORE_NINE_CORRELATED",
  "LEE_BMP2_TARGETS_DN",
  "HORIUCHI_WTAP_TARGETS_DN",
  "GSE33292_WT_VS_TCF1_KO_DN3_THYMOCYTE_DN",
  "GOBERT_OLIGODENDROCYTE_DIFFERENTIATION_UP",
  "MEBARKI_HCC_PROGENITOR_FZD8CRD_UP",
  "GSE3400_UNTREATED_VS_IFNB_TREATED_MEF_DN",
  "PUJANA_CHEK2_PCC_NETWORK",
  "KRAS.LUNG_UP.V1_DN"
)

setdiff(sig_gsea_TNBC$pathway, c(TNBC_BC, TNBC_DiseaseSignaling, TNBC_Immune, TNBC_Misc, TNBC_OtherCancer, TNBC_RegulatoryNetworks))
setdiff(c(TNBC_BC, TNBC_DiseaseSignaling, TNBC_Immune, TNBC_Misc, TNBC_OtherCancer, TNBC_RegulatoryNetworks), sig_gsea_TNBC$pathway)

# Visualize for TNBC
setwd("<YOUR DATA PATH>/TNBC")
for(pathways in c("TNBC_BC", "TNBC_DiseaseSignaling", "TNBC_Immune", "TNBC_Misc", "TNBC_OtherCancer", "TNBC_RegulatoryNetworks")){
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

ggsave("TNBC_GSEA.pdf", q, height = 20, width = 9)


# Repeat for ER+ data
setwd("<YOUR DATA PATH>/ER")

for(pathways in c("ER_BC", "ER_OtherCancer", "ER_RegulatoryAndSignaling", "ER_Misc")){
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
ggsave("ER_GSEA.pdf", q, height = 11)

# Save the donor information for a supplemental table
suppl <- clin_metabric[make.names(clin_metabric$PATIENT_ID) %in% c(colnames(ER_old2), colnames(ER_young2), 
                                                                   colnames(TNBC_old2), colnames(TNBC_young2)),]
suppl2 <- data.frame(PATIENT_ID = suppl$PATIENT_ID, AGE_AT_DIAGNOSIS = suppl$AGE_AT_DIAGNOSIS, 
                     THREEGENE = suppl$THREEGENE, TUMOR_STAGE = suppl$TUMOR_STAGE)

setwd("../")
write.csv(file = "METABRIC_donortable.csv", suppl2)
