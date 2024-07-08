#Prerequisites
library(readxl)
library(ggplot2)

# Read in the ASPEN results and make a box plot of enrichment scores based on molecular subtype
ASPEN_results <- read_excel("/Users/addie/Desktop/Suppl Table 4 ASPENResults.xlsx", sheet = 1)
ASPEN_results2 <- ASPEN_results[ASPEN_results$Enrich_Score != "NA" & ASPEN_results$Enrich_Score != 0,]

boxplot <- ggplot(ASPEN_results2, aes(x = Subtype, y = as.numeric(Enrich_Score)))+geom_boxplot()+
  ylab("Normalized Enrichment Score")

# Check to ensure the enrichment scores are significantly different
t.test(as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "ER"]), as.numeric(ASPEN_results2$Enrich_Score[ASPEN_results2$Subtype == "TNBC"]))
