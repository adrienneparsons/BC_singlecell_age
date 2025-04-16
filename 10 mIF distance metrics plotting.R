# Figure 4: Age associations in mIF images
###########################################################
# Hanyun Z and Adrienne Parsons, 2025-03-05

# Prerequisite libraries
library(dplyr)
library(ggplot2)

# Load the results
clin <- read.csv("<YOUR DATA DIRECTORY>/clinical_df_merge_age.csv")
dist <- read.csv("<YOUR DATA DIRECTORY>/distance_trait.csv")
ratio <- read.csv("<YOUR DATA DIRECTORY>/tumour_stromal_ratio.csv")

# Check that all patients in the clinical data are represented in the distance and ratio results
clin$Patient.ID %in% dist$Patient.ID %>% table()
clin$Patient.ID %in% ratio$X %>% table()

# Add patient age and molecular subtype to the data
dist$Age <- NA
dist$subtype <- NA
ratio$Age <- NA
ratio$subtype <- NA

for(id in clin$Patient.ID){
  dist$Age[dist$Patient.ID == id] <- clin$Age[clin$Patient.ID == id]
  dist$subtype[dist$Patient.ID == id] <- clin$subtype[clin$Patient.ID == id]
  ratio$Age[dist$Patient.ID == id] <- clin$Age[clin$Patient.ID == id]
  ratio$subtype[dist$Patient.ID == id] <- clin$subtype[clin$Patient.ID == id]
}

# Check for any NA values
is.na(dist$Age) %>% table()
is.na(dist$subtype) %>% table()
is.na(ratio$Age) %>% table()
is.na(ratio$subtype) %>% table()

# Check the subtype names (should be "TNBC" and "ER+)
unique(dist$subtype)

# Subset the results to be TNBC or ER+ alone
TNBC_dist <- dist[dist$subtype == "TNBC",]
ER_dist <- dist[dist$subtype == "ER+",]
TNBC_ratio <- ratio[ratio$subtype == "TNBC",]
ER_ratio <- ratio[ratio$subtype == "ER+",]

# Different age groups for comparison:
# Make lists of donor IDs that fall into different young and old categories
group1 <- clin$Patient.ID[clin$Age <= 55]
group2 <- clin$Patient.ID[clin$Age > 55]
group3 <- clin$Patient.ID[clin$Age >55 & clin$Age < 75]
group4 <- clin$Patient.ID[clin$Age <= 45]
group5 <- clin$Patient.ID[clin$Age >=65]
group6 <- clin$Patient.ID[clin$Age > 65 & clin$Age < 75]

# Identify which groups are young and which are old; in this case we just want to compare
# Young to Old with a threshold of 55, so we are comparing Group 1 to Group 2
young_groupings <- list("Group1" = group1)
old_groupings <- list("Group2" = group2)

# Metrics in the distance data frame we want to compare
dist_metrics <- c("endothelial_c_general_icaf_mycaf_Reference_perc_r30",
                  "general_icaf_mycaf_endothelial_c_Reference_perc_r30",
                  "endothelial_c_general_CTL_Reference_perc_r30",       
                  "general_CTL_endothelial_c_Reference_perc_r30",       
                  "general_CTL_general_icaf_mycaf_Reference_perc_r30",  
                  "general_CTL_icaf_Reference_perc_r30",                
                  "general_CTL_myCAF_Reference_perc_r30",               
                  "general_icaf_mycaf_general_CTL_Reference_perc_r30",  
                  "icaf_general_CTL_Reference_perc_r30",                
                  "mycaf_general_CTL_Reference_perc_r30",               
                  "general_CTL_epithelial_Reference_perc_r30",          
                  "endothelial_c_epithelial_Reference_perc_r30") 

# Metrics in the ratio data we want to compare
ratio_metrics <- c("cd8_t_general_tumour_perc", "cd8_t_general_stroma_perc", "endothelial_c_tumour_perc",
                   "endothelial_c_stroma_perc")

setwd("<YOUR RESULTS DIRECTORY>/TNBC")

# For the old groupings and the young groupings:
for(old in names(old_groupings)){
  for(young in names(young_groupings)){
    # Initialize a data frame with the comparison made, means and medians and the p-value
    df <- data.frame("comp" = character(0), "mean_old" = numeric(0), "median_old" = numeric(0),
                     "mean_young" = numeric(0), "median_young" = numeric(0), "pval" = numeric(0)
    )
    # Get the patients for the comparison being made
    olds <- na.omit(unlist(old_groupings[old]))
    youngs <- na.omit(unlist(young_groupings[young]))
    
    # Subset the data to just be the patients to be compared
    TNBC <- TNBC_dist[TNBC_dist$Patient.ID %in% c(olds, youngs),]
    
    # Add an annotation if the patient is Young or Old
    TNBC$Group <- NA
    TNBC$Group[TNBC$Patient.ID %in% olds] <- "Older"
    TNBC$Group[TNBC$Patient.ID %in% youngs] <- "Young"
    
    # For each metric to analyze:
    for(metric in dist_metrics){
      # Plot the metric valuues per patient by age group in a bix plot
      p <- ggplot(TNBC, aes(x = Group, y = .data[[metric]], fill = Group))+geom_boxplot()+
        geom_jitter()+
        ylab(metric)+
        theme(aspect.ratio = 1)+
        theme_minimal()
      
      # Save
      ggsave(filename = paste0(old, " vs. ", young, " ", metric, ".pdf"), p, device = "pdf")
      
      # Get a list of the values for the metric of interest for each age group
      o <- na.omit(TNBC[[metric]][TNBC$Patient.ID %in% olds])
      y <- na.omit(TNBC[[metric]][TNBC$Patient.ID %in% youngs])
      
      # Wilcox test for significance
      t <- wilcox.test(o, y)
      
      # Get the results and add to the blank data frame
      ls <- list("comp" = metric, "mean_old" = mean(o), "median_old" = median(o), "mean_young" = 
                   mean(y), "median_young" = median(y), "pval" = t$p.value)
      
      df <- rbind(df, ls)
      
    }
    
    # Save the complete statistics
    assign("TNBC_distmetrics", df)
  }
}

# Endothelial cells weren't of interest in TNBC; remove those comparisons
TNBC_distmetrics[grep("endothelial", TNBC_distmetrics$comp),] <- NA
TNBC_distmetrics <- na.omit(TNBC_distmetrics)

# Remove additional metrics not of interest
TNBC_distmetrics[2:5,] <- NA
TNBC_distmetrics <- na.omit(TNBC_distmetrics)

# Adjust p-values (FDR)
TNBC_distmetrics$padj <- p.adjust(TNBC_distmetrics$pval, method = "BH")

# Repeat for ratio metrics:
# For the old and young groupings:
for(old in names(old_groupings)){
  for(young in names(young_groupings)){
    
    # Initialize an empty data frame to populate with statistics
    df <- data.frame("comp" = character(0), "mean_old" = numeric(0), "median_old" = numeric(0),
                     "mean_young" = numeric(0), "median_young" = numeric(0), "pval" = numeric(0)
    )
    
    # Get the patient IDs for the young and old groups
    olds <- na.omit(unlist(old_groupings[old]))
    youngs <- na.omit(unlist(young_groupings[young]))
    
    # Subset to just the patients to analyze
    TNBC <- TNBC_ratio[TNBC_ratio$X %in% c(olds, youngs),]
    
    # Add an age group annotation
    TNBC$Group <- NA
    TNBC$Group[TNBC$X %in% olds] <- "Older"
    TNBC$Group[TNBC$X %in% youngs] <- "Young"
    
    # For each metric to compare in the ratio metrics:
    for(metric in ratio_metrics){
      
      # Plot a box plot per patient of that metric by age group
      p <- ggplot(TNBC, aes(x = Group, y = .data[[metric]], fill = Group))+geom_boxplot()+
        geom_jitter()+
        ylab(metric)+
        theme(aspect.ratio = 1)+
        theme_minimal()
      
      # Save
      ggsave(filename = paste0(old, " vs. ", young, " ", metric, ".pdf"), p, device = "pdf")
      
      # Get the data for the old and young patients for that metric
      o <- na.omit(TNBC[[metric]][TNBC$X %in% olds])
      y <- na.omit(TNBC[[metric]][TNBC$X %in% youngs])
      
      # Wilcox test for significance
      t <- wilcox.test(o, y)
      
      # Collect relevant data
      ls <- list("comp" = metric, "mean_old" = mean(o), "median_old" = median(o), "mean_young" = 
                   mean(y), "median_young" = median(y), "pval" = t$p.value)
      
      # Add to data frame
      df <- rbind(df, ls)
      
    }
    
    # Save complete data
    assign("TNBC_ratiometrics", df)
  }
}

# Endothelial cells weren't of interest in TNBC; remove
TNBC_ratiometrics[grep("endothelial", TNBC_ratiometrics$comp),] <- NA
TNBC_ratiometrics <- na.omit(TNBC_ratiometrics)

# Adjust p-values (FDR)
TNBC_ratiometrics$padj <- p.adjust(TNBC_ratiometrics$pval, method = "BH")


# Repeat for ER+:
setwd("<YOUR RESULTS DIRECTORY>/ER")

# For the old and young groupings to compare:
for(old in names(old_groupings)){
  for(young in names(young_groupings)){
    
    # Initialize and empty data frame to populatte with statistics
    df <- data.frame("comp" = character(0), "mean_old" = numeric(0), "median_old" = numeric(0),
                     "mean_young" = numeric(0), "median_young" = numeric(0), "pval" = numeric(0)
    )
    
    # Get the old and young patient IDs for the comparison
    olds <- na.omit(unlist(old_groupings[old]))
    youngs <- na.omit(unlist(young_groupings[young]))
    
    # Subset to the data of interest
    ER <- ER_dist[ER_dist$Patient.ID %in% c(olds, youngs),]
    
    # Add an age annotation to the patients
    ER$Group <- NA
    ER$Group[ER$Patient.ID %in% olds] <- "Older"
    ER$Group[ER$Patient.ID %in% youngs] <- "Young"
    
    # For each metric of interest:
    for(metric in dist_metrics){
      
      # Plot a boxplot comparing each patient's value for that metric, grouped by age group
      p <- ggplot(ER, aes(x = Group, y = .data[[metric]], fill = Group))+geom_boxplot()+
        geom_jitter()+
        ylab(metric)+
        theme(aspect.ratio = 1)+
        theme_minimal()
      
      # Save
      ggsave(filename = paste0(old, " vs. ", young, " ", metric, ".pdf"), p, device = "pdf")
      
      # Get the relevant metric data for the comparison being made
      o <- na.omit(ER[[metric]][ER$Patient.ID %in% olds])
      y <- na.omit(ER[[metric]][ER$Patient.ID %in% youngs])
      
      # Wilcox test for significance
      t <- wilcox.test(o, y)
      
      # Collect relevant statistics
      ls <- list("comp" = metric, "mean_old" = mean(o), "median_old" = median(o), "mean_young" = 
                   mean(y), "median_young" = median(y), "pval" = t$p.value)
      
      # Add to data frame
      df <- rbind(df, ls)
      
    }
    
    # Save the complete data
    assign("ER_distmetrics", df)
  }
}

# Remove metrics not of interest
ER_distmetrics[c(3:10, 12), ] <- NA
ER_distmetrics <- na.omit(ER_distmetrics)

# Adjust p-values (FDR)
ER_distmetrics$padj <- p.adjust(ER_distmetrics$pval, method = "BH")

# Repeat for ratio metrics:

# For the groups you'd like to compare:
for(old in names(old_groupings)){
  for(young in names(young_groupings)){
    
    # Initialize a data frame to populate with statistics
    df <- data.frame("comp" = character(0), "mean_old" = numeric(0), "median_old" = numeric(0),
                     "mean_young" = numeric(0), "median_young" = numeric(0), "pval" = numeric(0)
    )
    
    # Get the patients IDs for the groups to compare
    olds <- na.omit(unlist(old_groupings[old]))
    youngs <- na.omit(unlist(young_groupings[young]))
    
    # Subset the data to just be patients of interest
    ER <- ER_ratio[ER_ratio$X %in% c(olds, youngs),]
    
    # Add an age group annotation
    ER$Group <- NA
    ER$Group[ER$X %in% olds] <- "Older"
    ER$Group[ER$X %in% youngs] <- "Young"
    
    # For each of the ratio metrics of interest:
    for(metric in ratio_metrics){
      
      # plot a box plot comparing the metrics values for each patient, grouped by age group
      p <- ggplot(ER, aes(x = Group, y = .data[[metric]], fill = Group))+geom_boxplot()+
        geom_jitter()+
        ylab(metric)+
        theme(aspect.ratio = 1)+
        theme_minimal()
      
      # Save
      ggsave(filename = paste0(old, " vs. ", young, " ", metric, ".pdf"), p, device = "pdf")
      
      # Get the relevant metric data for the patients being compared
      o <- na.omit(ER[[metric]][ER$X %in% olds])
      y <- na.omit(ER[[metric]][ER$X %in% youngs])
      
      # Wilcox test for significance
      t <- wilcox.test(o, y)
      
      # Get the relevant statistics
      ls <- list("comp" = metric, "mean_old" = mean(o), "median_old" = median(o), "mean_young" = 
                   mean(y), "median_young" = median(y), "pval" = t$p.value)
      
      # Add to the data frame
      df <- rbind(df, ls)
      
    }
    
    # Save the final data to environment
    assign("ER_ratiometrics", df)
  }
}

# Adjust p-values (FDR)
ER_ratiometrics$padj <- p.adjust(ER_ratiometrics$pval, method = "BH")