# Data wrangling: remove duplicated genes and break down by subtype
###########################################################
# Adrienne Parsons, 2024-04-01

# load packages
library(tidyverse)

# read in the raw data from METABRIC, acquired from cBioPortal (https://www.cbioportal.org/study/summary?id=brca_metabric)
clinical_data <- read.table("<YOUR DATA PATH>/data_clinical_patient.txt",
                            sep = "\t")
transcript_data <- read_excel("<YOUR DATA PATH>/data_mrna_agilent_microarray.xlsx")

# Format the data a little so that the first row of the table is the
# column names
colnames(clinical_data) <- clinical_data[1,]
clinical_data <- clinical_data[-1,]

# Filter down the clinical data to TNBC donors who have ages on file
clinical_data_TNBC <- clinical_data %>% filter(., THREEGENE == "ER-/HER2-")
clinical_data_TNBC <- clinical_data_TNBC[is.na(clinical_data_TNBC$AGE_AT_DIAGNOSIS) == F,]

# Filter down the clinical data to ER+ donors who have ages on file
clinical_data_ER <- clinical_data[clinical_data$THREEGENE == "ER+/HER2- High Prolif" |
                                    clinical_data$THREEGENE == "ER+/HER2- Low Prolif",]
clinical_data_ER <- clinical_data_ER[is.na(clinical_data_ER$AGE_AT_DIAGNOSIS) == F,]

# Get a list of the donor IDs in this edited list
donor_IDs_TNBC <- clinical_data_TNBC$PATIENT_ID
donor_IDs_ER <- clinical_data_ER$PATIENT_ID

# Some gene are duplicated in the METABRIC data; average the values
# First, make the data frame values numeric and add the gene symbols back in
td2 <- data.frame(sapply(data.frame(transcript_data[,3:length(transcript_data)]), as.numeric))
numtd <- cbind(transcript_data[,1], td2)

# Make NA values 0
numtd[is.na(numtd)] <- 0

# Check the number of unique genes; this indicates there are some duplicates
length(unique(numtd$Hugo_Symbol))

# Get a list of every gene entry and the unique genes
allgenes <- numtd$Hugo_Symbol
unique_genes <- unique(allgenes)

# Initialize an empty list to populate with the duplicate gene names
dupgenes <- list()

# Iterate through all the genes and add the duplicate gene names to the list
# This takes a long time!
for(gene in unique_genes){
# To track progress, print every a number every 50 genes
  if(which(unique_genes == gene) %% 50 == 00){
    print(which(unique_genes == gene))
  }
  # If duplicate, add to list
  if(nrow(numtd[numtd$Hugo_Symbol == gene,]) > 1){
    dupgenes <- append(dupgenes, gene)
  }
}

# Make the list of duplicate genes a character string
dupgenes <- unlist(dupgenes)

# Initialize a data frame that just has the gene names as row names
avgdup_transcriptdata <- numtd[1,]
avgdup_transcriptdata <- avgdup_transcriptdata[-1,]

# for each gene in the duplicate gene list
for(gene2 in dupgenes){
  # Print the progress
  print(which(dupgenes == gene2))
  # Subset the numerical gene data to just be the rows with the duplicate gene
  transcript_data_sub <- numtd[which(numtd$Hugo_Symbol == gene2), 2:length(numtd)]
  # get the mean expression for duplicates
  avgs <- colMeans(transcript_data_sub)
  # Add the averages and the gene names to the empty data frame
  avgdup_transcriptdata <- rbind(avgdup_transcriptdata, c(gene2, avgs))
}

# Format the averaged data frame so the donor names are included
colnames(avgdup_transcriptdata) <- colnames(numtd)

# append the single gene entries to the averaged data for duplicate genes
final_data <- rbind(avgdup_transcriptdata, numtd[numtd$Hugo_Symbol %in% dupgenes == F,])

# Check to ensure that there are no duplicates in the final data
length(unique(final_data$Hugo_Symbol)) == nrow(final_data)

# Filter down the transcript data to just have the IDs in the 
# edited clinical data, and add the gene symbols back to that data frame
TD_TNBC <- final_data[,colnames(final_data) %in% make.names(donor_IDs_TNBC)]
TD_TNBC$Gene_Symbol <- final_data$Hugo_Symbol
TD_TNBC2 <- cbind(TD_TNBC[,length(TD_TNBC)], TD_TNBC[,1:length(TD_TNBC)-1])
colnames(TD_TNBC2) <- c("Hugo_Symbol", colnames(TD_TNBC2[,2:length(TD_TNBC2)]))

TD_ER <- final_data[,colnames(final_data) %in% make.names(donor_IDs_ER)]
TD_ER$Gene_Symbol <- final_data$Hugo_Symbol
TD_ER2 <- cbind(TD_ER[,length(TD_ER)], TD_ER[,1:length(TD_ER)-1])
colnames(TD_ER2) <- c("Hugo_Symbol", colnames(TD_ER2[,2:length(TD_ER2)]))

# Save the data tables
write.csv(TD_TNBC2, file = "2024_METABRIC_TNBC.csv", quote = F, row.names = F)
write.csv(TD_ER2, file = "2024_METABRIC_ER.csv", quote = F, row.names = F)


