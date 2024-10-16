# Figure S4: ASPEN results boxplot
###########################################################
# Adrienne Parsons, 2024-09-10

#Prerequisites
library(readxl)
library(ggplot2)

setwd("<YOUR DATA DIRECTORY>")

# Read in the ASPEN results and make a box plot of enrichment scores based on molecular subtype
# These tables are outputs from "Figure 3 Hallmark Aspen.R"
TNBC_ASPEN_results <- read.csv("TNBC_ASPENResults.csv")
ER_ASPEN_results <- read.csv("ER_ASPENResults.csv")

# Merge the two results tables
TNBC_ASPEN_results$Subtype <- "TNBC"
ER_ASPEN_results$Subtype <- "ER+"
ASPEN_results <- rbind(TNBC_ASPEN_results, ER_ASPEN_results)

# Subset to just the results with a numerical NES
ASPEN_results2 <- ASPEN_results[ASPEN_results$Enrich_Score != "NA" & ASPEN_results$Enrich_Score != 0,]
ASPEN_results2 <- na.omit(ASPEN_results2)

# Plot the boxplot and save
boxplot <- ggplot(ASPEN_results2, aes(x = Subtype, y = as.numeric(Enrich_Score)))+geom_boxplot(outlier.shape = NA)+
  ylab("Normalized Enrichment Score")+
  geom_jitter(alpha = 0.5, colour = "purple4")+
  theme(aspect.ratio = 1)+
  theme_bw()
  
  ggsave("ASPEN_results_boxplot.pdf", boxplot, device = "pdf")

# Check to ensure the enrichment scores are significantly different
t.test(as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "ER+"]), as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "TNBC"]))
