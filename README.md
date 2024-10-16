# BC_singlecell_age

This repository contains several analysis scripts and files used to assess the role of age in bulk microarray (METABRIC) and a single-cell RNA sequencing atlas (Wu… Swarbrick Nature Genetics, 2021). It is in support of a publication by Parsons et al. 2024. The included scripts have multiple aims: 1) infer biological processes that differ with patient age and molecular disease subtype (ER+ and triple-negative, TNBC) through differential gene expression followed by GSEA of the METABRIC data. 2)  it examines changes of cell type composition with age using the single-cell atlas, 3) using an analysis pipeline designed for the accompanying manuscript we have termed ASPEN (Age-Specific Program ENrichment), transcriptomic trends that correlate with age are distilled into genetic program enrichments with age (using pathways found on MSigDB). 4) through analysis using the CellChat package (Jin et al. Nature Communications 2021), the scripts infer intercellular communication between given cell types in the single cell atlas, comparing older and younger donors within disease subtype groups. The analyses from CellChat include inference of interactome networks that depict the communication probabilities between cell types, an application of logistic regression for determining signaling pathways of interest, and an assessment of the ligand-receptor pairs through which the significant pathways are engaged.

### Set up for the scripts
The 2021 single-cell RNA seq data can be downloaded from GEO at accession number GSE176078. Matrix, barcode, feature, and metadata files for each donor should be in unique folders. METABRIC data and metadata can be found on cBioPortal. For ease, we have included the clinical and sample-specific metadata from cBioPortal in the Data folder of this repository. The raw METABRIC data, subsetted by ER+ and TNBC donors (output from “METABRIC data wrangling.R”) can be found on FigShare at:

TNBC: https://figshare.com/articles/dataset/2024_METABRIC_TNBC_csv/27242253?file=4983452

ER+:
https://figshare.com/articles/dataset/2024_METABRIC_ER_csv/27242256?file=49834524

1.  Additional files needed for analysis include the .gmt files for the c2, c5, c6, c7, and Hallmark pathways from MSigDB, which are also included in the Data folder of this repository. For ease of replicating CellChat analyses, we have also included .Rdata files of the CellChat objects output from “Create CellChat Objects.R”.

All analyses are completed in R v4.4.1 and can be completed on a local machine. Both Seurat v4 and v5 are needed for different parts of the scripts. The version required for each R file is noted in the file’s comments.

### Running the scripts
#### METABRIC data wrangling and Figures 1 and S1
Using the “data_clinical_patient.txt” and “data_clinical_sample.txt” files on cBioPortal or copied into this repository, as well as the raw data from cBioPortal, the “METABRIC data wrangling.R” file subsets the large gene expression files by molecular subtype (ER+ or TNBC) and removes duplicate genes. This script takes >30 minutes to complete; for ease, we have made the outputs of this script available on FigShare:

 TNBC: https://figshare.com/articles/dataset/2024_METABRIC_TNBC_csv/27242253?file=4983452

ER+:
https://figshare.com/articles/dataset/2024_METABRIC_ER_csv/27242256?file=49834524

Once the gene expression files are formatted or obtained from FigShare, “Figure 1 and S1 METABRIC DE Genes and GSEA.R” performs differential gene expression analysis between the ER+ and TNBC donors <45 and >65 years of age with stage I-III tumors. Once tables of the necessary donors’ gene expression are made, the genes with the highest variance across the entire old and young donor pool per subtype are identified, and differential gene expression analysis is performed using limma. The results of this analysis are visualized as a volcano plot:

![Volcano Plots.](/www/VolcanoPlots.png)

The highly variable genes are then ranked by fold change difference (highest positive fold change to lowest negative fold change) and GSEA is performed using the .gmt files found on MSigDB or in the Data folder of this repository. The results, which use the fgsea package to find Normalized Enrichment Scores, are visualized as bubble plots:

![GSEA Bubble Plots.](/www/GSEA_bubbleplots.png)

Table outputs include the donor metadata used for the analysis, as well as the DE gene expression analysis results for TNBC and ER+.

### Figures S2 and S3, Proportion age analysis
Using the single cell atlas’s cell type annotations of multiple granularities, the “Figure S2 and S3 Proportion Age Analysis.R” file assesses cell type composition with age. It generates stacked bar charts of each minor cell type (celltype_minor) as a proportion of each major cell type (celltype_major), for each donor, and orders the donors on the x-axis by age. For additional ease of interpretation, scatterplots of donor age vs. proportion for each minor cell type are generated and Pearson Correlation testing is done to test for significant associations between proportion values and age.

<img src="/www/Proportions.png" alt="proportions example" width="250"/>

### Figures 3, S4, S5, and S6, ASPEN
We developed ASPEN (Age-Specific Program ENrichment) as a means of identifying pathways or gene programs that are enriched with older or younger age in the single-cell atlas data. The details of ASPEN’s functions can be found in the comments of “Figure 3 Hallmark Aspen.R” and the associated manuscript. Briefly, Seurat objects for each donor are subsetted by celltype_minor. For each cell type, the average gene expression of every gene is calculated, and these averages are correlated to donor age (n = 10 for TNBC, n = 11 for ER+), and genes are then ranked by correlation coefficient. Genes with a 0 coefficient are removed to avoid erroneous ranking, and GSEA is performed. Second, Seurat signature scoring commands are used to assign a signature score for each pathway of interest to Seurat objects containing all TNBC or ER+ donors; these scored objects are then further subsetted to find the average signature score per cell type per donor. For each cell type, these scores are correlated to donor age. The results from both analyses are then visualized in a bubble plot:

![Hallmark Pathway ASPEN.](/www/Hallmark_ASPEN.png)

This script also outputs the numerical results from ASPEN analysis. The Normalized Enrichment Scores calculated from ASPEN showed cell type-pathway combinations were predominantly enriched in older patients in TNBC and in younger patients in ER+ breast cancer. We show the distribution of Normalized Enrichment Scores for TNBC vs. ER+ in a box plot, generated by “Figure S4 ASPEN results boxplot.R”.

<img src="/www/boxplot.png" alt="boxplot" width="250"/>

ASPEN can be used to assess age associations of any pathway found on MSigDB and for any cell type annotations. Figure 3 utilizes the celltype_minor annotations and the Hallmark Pathways, but we also demonstrate ASPEN’s capabilities for more granular cell populations (celltype_subset, “Figure S5 Hallmark ASPEN celltype_subset.R”) and for selected pathways related to cellular senescence (“Figure S6 Senescnece ASPEN.R”)

![Celltype_subset ASPEN.](/www/Subset_ASPEN.png)

![Senescence ASPEN.](/www/Senescence.png)

### Cellchat object creation and Figure 4, CellChat interactome analysis
To better understand the activity of specific cell types within young and aged breast tumor microenvironments, these scripts utilized the package CellChat (Jin et al. Nature Communications, 2021). First, CellChat objects are generated from the single-cell atlas in the file “Create CellChat Objects.R”. The objects created include all ER+ donors, all TNBC donors, then all ER+ donors ≤55, all ER+ >55, all TNBC donors ≤55, and all TNBC donors >55. Because these samples are unsorted, we set the argument computeCommunProb() to be TRUE. This script takes >60 minutes to run, so for ease of replication, we have included the resulting CellChat objects used in downstream analyses in the Data folder of this repository.

Once the CellChat objects are made or obtained from the Data folder, the script uses the standard CellChat pipeline to examine cell type to cell type interactions broadly. First, it infers the communication networks between major cell populations (celltype_major) and generate circle plots for the TNBC old and young, then ER+ old and young populations, and then within subtype a differential signaling network to show differences in interactome with age is generated. In these plots, red edges indicate enrichment in the older age group, while blue edges indicate enrichment in the younger group.

![Network Plots.](/www/Networks.png)

The script also looks at communication probability differences between cell types with age using additional standard CellChat pipelines, which create scatter plots of the minor cell types (celltype_minor), where each data point is a cell type, and the axes represent incoming and outgoing communication probabilities. 

![Scatter Plots.](/www/Scatters.png)

Finally, to visualize the age-associated differences in celltype_minor to celltype_minor communication, the script uses standard CellChat commands to generate differential heat maps, where each square represents the communication probability of a given source/target cell type, red color indicates enrichment in the older cohort, and blue color indicates enrichment in the younger cohort.

![Heatmaps.](/www/heatmaps.png)

### Figure S7, Assessment of ESR1 and MHCII genes by cell type in ER+ and TNBC
In TNBC, the increased macrophage-T cell interactions with age, along with signals of increased MHC-II presentation from the METABRIC analysis suggested an age-associated difference in MHC-II activity. We therefore used the single-cell atlas to examine cell type-specific gene expression of MHC-II genes in TNBC, comparing donors ≤55 and >55. We also examined the cell type-specific expression of ESR1 in ER+ donors following the METABRIC analyses from Figure 1. This analysis is completed using the file “Figure S7 EST1 and MHC II gene expression.R,” which generates sina plots for each gene of interest for each minor cell type (celltype_minor) in donors ≤55 and >55 years.

![Sina Plots.](/www/sinas.png)

### Figures 5, 6, S8, S9, and S10, Further CellChat-guided investigation.
The script “Figure 5 6 S8 S9 S10 RankNet and bubble plots.R” completes the analyses done for this manuscript by further investigating the specific pathways through which the cell types are communicating using CellChat commands. As described in the accompanying manuscript, we examined the communication probability differences with age (the numerical values associated with the colors in the CellChat heatmaps) and determined that there were certain cell types of interest within the TNBC and ER+ cohorts. We focused our analyses on these 7 (TNBC) or 8 (ER+) cell types and their interactions with one another. From there, the aim was to identify which signaling pathways these cell types were interacting through, and specifically which ligand-receptor interactions showed enrichment in one of the age cohorts within each subtype.

The script first leads the created CellChat objects and introduces a function used for formatting data tables in downstream analysis. It then uses standard CellChat commands, specifying the senders/receivers of interest, to generate bar charts that show differential enrichment of given pathways within each sender-receiver combination, using the rankNet() function. The pathways identified in each sender/receiver plot are saved, as is the number of times a given pathway from the CellChat database is represented in a plot. 

![Example RankNet Bar Charts.](/www/RankNets.png)

The resultant list of pathways included too many to investigate individually, so from the list of pathways present in every rankNet analysis, the script implements a univariate logistic regression to identify pathways of interest. Each sender-receiver communication probability for a given cell type and age group is an individual data point, and probabilities of association to the >55 group are calculated. Pathways that have significant associations in either the positive or negative direction that also appear >15 times throughout all plots were deemed “of interest” for the final analysis. This final analysis takes the pathways of interest, and using the standard netVisual_bubble() command, generates a bubble plot of ligand-receptor interaction probabilities for each source/target cell combination, each ligand-receptor interaction within the identified pathways, for both the young and aged cohort.

![CellChat Bubble Plots.](/www/cellchat_bubble.png)

 For ease of interpretation in the final figure of the manuscript, we also include code to identify the fold-difference in communication probability between the young and aged cohorts for each ligand-receptor interaction/sender/receiver combination. Any combination with where the probability is significant at p < 0.01 exclusively in either the young or aged cohort, or where the difference between probabilities in both cohorts is >1.2 are identified and saved to a table, which was used to highlight specific differences in the manuscript’s figure.
