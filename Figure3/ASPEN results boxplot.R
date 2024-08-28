#Prerequisites
library(readxl)
library(ggplot2)

setwd("/Users/addie/Partners HealthCare Dropbox/Adrienne Parsons/Single-cell_breast_cancer/BC_Aging_Paper/Figures/Draft 6/")

# Read in the ASPEN results and make a box plot of enrichment scores based on molecular subtype
TNBC_ASPEN_results <- read_excel("/Users/addie/Partners HealthCare Dropbox/Adrienne Parsons/Single-cell_breast_cancer/BC_Aging_Paper/Figures/Draft 6/Suppl tables/Suppl Table 4 ASPENResults.xlsx", sheet = 1)
ER_ASPEN_results <- read_excel("/Users/addie/Partners HealthCare Dropbox/Adrienne Parsons/Single-cell_breast_cancer/BC_Aging_Paper/Figures/Draft 6/Suppl tables/Suppl Table 4 ASPENResults.xlsx", sheet = 2)

TNBC_ASPEN_results$Subtype <- "TNBC"
ER_ASPEN_results$Subtype <- "ER+"

ASPEN_results <- rbind(TNBC_ASPEN_results, ER_ASPEN_results)

ASPEN_results2 <- ASPEN_results[ASPEN_results$Enrich_Score != "NA" & ASPEN_results$Enrich_Score != 0,]

boxplot <- ggplot(ASPEN_results2, aes(x = Subtype, y = as.numeric(Enrich_Score)))+geom_boxplot(outlier.shape = NA)+
  ylab("Normalized Enrichment Score")+
  geom_jitter(alpha = 0.5, colour = "purple4")+
  theme(aspect.ratio = 1)+
  theme_bw()
  
  ggsave("20240724_ASPEN_results_boxplot.pdf", boxplot, device = "pdf")

# Check to ensure the enrichment scores are significantly different
t.test(as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "ER+"]), as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "TNBC"]))
