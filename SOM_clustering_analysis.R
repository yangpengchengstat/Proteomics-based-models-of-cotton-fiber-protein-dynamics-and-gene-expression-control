# ---------------------------------------------------------------------------------------
# SOM Clustering Analysis for Cotton Proteomics Data (Merged Version)
# 
# This script reads protein abundance data from three subcellular fractions (APO, P200, S200),
# applies filtering (≥ 3 DPAs, ≥ 2 consecutive DPAs, etc.), interpolates missing values,
# performs Self-Organizing Map (SOM) clustering, handles outlier/unreliable removal, then
# merges the “informative” clusters across fractions for cross-fraction analysis.
#
# It produces CSV files and PDFs documenting each step (raw clustering, outlier removal,
# unreliable profile filtering, final informative clusters, etc.).
#
# Author: Pengcheng Yang (2023). Prepared for GitHub usage.
# ---------------------------------------------------------------------------------------

rm(list = ls())  # Clear the R workspace

# Required libraries
library(kohonen)    # Self-Organizing Maps
library(tempR)      # Additional time-series / correlation functions
library(dplyr)      # Data manipulation
library(ggplot2)    # Plotting
library(openxlsx)   # Reading/writing .xlsx files
library(tibble)     # Tidy data frames
library(tidyverse)  # Comprehensive data manipulation & I/O
require(Rcpp)       # Integrate R and C++ for certain distance measures

# Expose the wcc0 C++ function (WCC-based distance) from the kohonen package
# 'mustWork=TRUE' ensures an error if it's not found
sourceCpp(system.file("Distances", "wcc0.cpp", package = "kohonen", mustWork = TRUE))

# Check current working directory (for debugging)
getwd()

# User-defined working directory for input data
workingDir = "/path/to/your/input_data_directory"
setwd(workingDir)

# Avoid string→factor conversion
options(stringsAsFactors = FALSE)

##########################################################
# 1) Data importation & basic setup
##########################################################

# The three subcellular fractions
fraction_names <- c("APO", "P200", "S200")

# Initialize lists to store data
prot_exp_data         <- list()  # Protein expression table (per fraction)
prot_group_score_data <- list()  # Score/annotation (e.g., which cluster is "informative")
prot_notation_data    <- list()  # Additional ID annotation (gene name, ID, etc.)

# Loop over fractions to read .xlsx files
for (frac in fraction_names){
  # 1) Protein abundance data
  prot_exp_table <- read.xlsx(paste0(frac,".xlsx"))
  prot_exp_table <- prot_exp_table %>%
    remove_rownames %>%
    column_to_rownames(var=names(prot_exp_table)[1]) %>%
    as.data.frame()
  prot_exp_data[[frac]] <- prot_exp_table
  rm(prot_exp_table)
  
  # 2) Group/score table (which proteins are considered "I" = informative, etc.)
  score_table <- read.xlsx(paste0("score_", frac, ".xlsx"))
  prot_group_score_data[[frac]] <- score_table
  rm(score_table)
  
  # 3) Notation info (gene ID, name). The second column is set as rownames
  prots_noation_table <- read.xlsx(paste0(frac,"_name_id.xlsx"))
  prots_noation_table <- prots_noation_table %>%
    remove_rownames %>%
    column_to_rownames(var=names(prots_noation_table)[2]) %>%
    as.data.frame()
  prot_notation_data[[frac]] <- prots_noation_table
  rm(prots_noation_table)
}

# Additional data: protein-mRNA ID matching, mRNA data
prot_mRNA_ID_map_data <- read.xlsx("proteinIDmap.xlsx")
mRNA_data <- read.xlsx("mRNA seq exp data.xlsx")
mRNA_data <- mRNA_data %>%
  remove_rownames %>%
  column_to_rownames(var=names(mRNA_data)[1]) %>%
  as.data.frame()

##########################################################
# 2) SOM clustering analysis on each fraction
##########################################################

# We'll store results in these structures
SOM_clustering_result          <- list()
SOM_codes_result               <- list()
SOM_prot_cluster               <- list()
mean_of_normalized_informative <- list()

# Hist y-limit for some plots
histplot_ylim <- c(0, 25)

# For each fraction: read data, filter, run SOM, remove outliers, etc.
for (frac in fraction_names){
  
  # Adjust these paths or create them as needed
  workingDir = paste0("/path/to/SOM_clustering_results/", frac)
  setwd(workingDir)
  
  # Set SOM grid and abundance cutoffs per fraction
  if (frac == "APO"){
    x_limit <- 5
    y_limit <- 6
    abundance_cutoff <- 20
  } else if (frac == "S200"){
    x_limit <- 5
    y_limit <- 4
    abundance_cutoff <- 20
  } else if (frac == "P200"){
    x_limit <- 5
    y_limit <- 8
    abundance_cutoff <- 19
  } else {
    print("Fraction name not recognized.")
  }
  
  # Possibly adjust histogram scale if fraction is S200
  if(frac=="S200"){
    histplot_ylim <- c(0, 25)
  } else {
    histplot_ylim <- c(0, 70)
  }
  
  # Retrieve protein expression data (fraction-specific) and transpose so rows = samples/time
  protfrac_exp_data <- prot_exp_data[[frac]]
  protfrac_exp_data_T <- as.data.frame(t(protfrac_exp_data))
  
  # Identify proteins with average abundance > certain cutoff
  prots_high_abundance <- c()
  for (p in names(protfrac_exp_data_T)) {
    p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
    if (p_ab >= abundance_cutoff) {
      prots_high_abundance <- c(prots_high_abundance, p)
    }
  }
  
  # Number of DPAs / timepoints
  num_dpa <- nrow(protfrac_exp_data_T)
  
  # Create a small text log
  f_name <- paste(paste("LOG_OF_COMPUTATION_on", frac, sep="_"), "txt", sep=".")
  cat("LOG OF COMPUTATION\n\n", file=f_name)
  cat("1: There are", ncol(protfrac_exp_data_T), "proteins in the data set of", frac, ".\n",
      file=f_name, append=TRUE)
  
  # Filter: proteins present in ≥1 DPA
  prots_data_least1 <- data.frame(a=1:num_dpa)
  count_gt01 <- 0
  for (i in names(protfrac_exp_data_T)){
    df_i <- protfrac_exp_data_T[i]
    df_i_temp <- df_i[,1]
    if (sum(df_i_temp>0)>0){
      count_gt01 <- count_gt01 + 1
      prots_data_least1 <- cbind(prots_data_least1, df_i)
    }
  }
  prots_data_least1 <- prots_data_least1[,-1]
  mean_of_least1dpa_abandance <- c()
  for (p in names(prots_data_least1)){
    profile_p <- prots_data_least1[,p]
    profile_p_non0 <- profile_p[profile_p>0]
    mean_of_least1dpa_abandance <- c(mean_of_least1dpa_abandance, mean(profile_p_non0))
  }
  
  # Filter: proteins present in ≥3 DPAs & at least 2 consecutive DPAs
  prots_data_least3 <- data.frame(a=1:num_dpa)
  count_gt03 <- 0
  for (i in names(protfrac_exp_data_T)){
    df_i <- protfrac_exp_data_T[i]
    df_i_temp <- df_i[,1]
    if (sum(df_i_temp>0)>2){
      count_gt03 <- count_gt03 + 1
      check_num_gt0 <- (df_i>0)
      index_gt0 <- which(check_num_gt0==TRUE)
      # If there's at least one place with consecutive (difference==1)
      if (any((index_gt0[2:length(index_gt0)] - index_gt0[1:(length(index_gt0)-1)])==1)){
        prots_data_least3 <- cbind(prots_data_least3, df_i)
      }
    }
  }
  prots_data_least3 <- prots_data_least3[,-1]
  cat("2: There are", count_gt03, "proteins present in ≥3 DPAs.\n", file=f_name, append=TRUE)
  cat("3: There are", ncol(prots_data_least3),
      "proteins present in ≥3 DPAs and ≥2 consecutive DPAs.\n", 
      file=f_name, append=TRUE)
  
  mean_of_least3dpa_abandance <- c()
  for (p in names(prots_data_least3)){
    profile_p <- prots_data_least3[,p]
    profile_p_non0 <- profile_p[profile_p>0]
    mean_of_least3dpa_abandance <- c(mean_of_least3dpa_abandance, mean(profile_p_non0))
  }
  
  # Filter: proteins present in ≥10 DPAs
  prots_data_least10 <- data.frame(a=1:num_dpa)
  for (i in names(protfrac_exp_data_T)){
    df_i <- protfrac_exp_data_T[i]
    df_i_temp <- df_i[,1]
    if (sum(df_i_temp>0)>9){
      prots_data_least10 <- cbind(prots_data_least10, df_i)
    }
  }
  prots_data_least10 <- prots_data_least10[,-1]
  cat("4: There are", ncol(prots_data_least10), "proteins present in ≥10 DPAs.\n", file=f_name, append=TRUE)
  
  mean_of_least10dpa_abandance <- c()
  for (p in names(prots_data_least10)){
    profile_p <- prots_data_least10[,p]
    profile_p_non0 <- profile_p[profile_p>0]
    mean_of_least10dpa_abandance <- c(mean_of_least10dpa_abandance, mean(profile_p_non0))
  }
  
  # Interpolate data
  data_introp <- data.frame(a=1:num_dpa)
  for (i in names(prots_data_least3)){
    prot_d <- prots_data_least3[, i]
    df_i <- prots_data_least3[i]
    check_0 <- (prot_d==0)
    
    # Single zero between two non-zeros
    for (j in 2:(length(check_0)-1)){
      if (check_0[j]==TRUE && check_0[j-1]==FALSE && check_0[j+1]==FALSE){
        df_i[j,i] <- (df_i[j-1,i] + df_i[j+1,i]) / 2
      }
    }
    # 2 zeros in 4-step
    for (j in 4:length(check_0)){
      if (all(c(check_0[j-3]==FALSE, check_0[j-2]==TRUE, check_0[j-1]==TRUE, check_0[j]==FALSE))){
        slope_j <- (df_i[j,i] - df_i[j-3,i]) / 3
        df_i[j-2,i] <- slope_j + df_i[j-3,i]
        df_i[j-1,i] <- (slope_j*2) + df_i[j-3,i]
      }
    }
    # Edge correction for zeros at start/end
    if (check_0[1]==TRUE && check_0[2]==FALSE){
      df_i[1,i] <- df_i[2,i]
    }
    if (check_0[length(check_0)]==TRUE && check_0[length(check_0)-1]==FALSE){
      df_i[length(check_0),i] <- df_i[length(check_0)-1,i]
    }
    if (all(c(check_0[1]==TRUE, check_0[2]==TRUE, check_0[3]==FALSE, check_0[4]==FALSE))){
      df_i[1,i] <- df_i[3,i]
      df_i[2,i] <- df_i[4,i]
    } else if (all(c(check_0[1]==TRUE, check_0[2]==TRUE, check_0[3]==FALSE, check_0[4]==TRUE))){
      df_i[1,i] <- df_i[3,i]
      df_i[2,i] <- df_i[3,i]
    }
    if (all(c(check_0[length(check_0)], check_0[length(check_0)-1], check_0[length(check_0)-2]==FALSE, check_0[length(check_0)-3]==FALSE))){
      df_i[length(check_0),i] <- df_i[length(check_0)-2,i]
      df_i[length(check_0)-1,i] <- df_i[length(check_0)-3,i]
    }
    data_introp <- cbind(data_introp, df_i)
  }
  data_introp <- data_introp[,-1]
  
  # Build SOM
  Data <- as.matrix(t(data_introp))
  
  # Normalize row-wise (L2 norm)
  normalized_Data <- Data
  for (p in rownames(normalized_Data)){
    normalized_Data[p,] <- normalized_Data[p,] / sqrt( normalized_Data[p,] %*% normalized_Data[p,] )[1,1]
  }
  
  # Correlation matrix for cluster-based outlier removal
  all_corr_matrix <- normalized_Data %*% t(normalized_Data)
  all_corr <- all_corr_matrix[upper.tri(all_corr_matrix, diag = FALSE)]
  
  set.seed(123)
  system.time(
    SOM_wcc <- supersom(
      Data,
      grid = somgrid(xdim=x_limit, ydim=y_limit, topo="rectangular",
                     neighbourhood.fct="gaussian", toroidal=FALSE),
      dist.fcts = c("WCCd0"),
      rlen=500,
      normalizeDataLayers=FALSE
    )
  )
  
  # Summarize cluster assignments
  prot_id <- rownames(Data)
  label_of_sample <- SOM_wcc$unit.classif
  cluster_ID <- levels(factor(label_of_sample))
  
  sample_in_each_unit <- list()
  for (i in seq_len(x_limit * y_limit)){
    postion_index <- (label_of_sample==i)
    name <- paste0('V', as.character(i))
    if (any(postion_index)){
      sample_in_each_unit[[name]] <- prot_id[postion_index]
    } else {
      sample_in_each_unit[[name]] <- NULL
    }
  }
  SOM_prot_cluster[[frac]] <- sample_in_each_unit
  
  # Retrieve codebook vectors
  codes_som_all <- SOM_wcc$codes[[1]]
  normalized_codes_som_all <- codes_som_all
  for (p in rownames(normalized_codes_som_all)){
    normalized_codes_som_all[p,] <- normalized_codes_som_all[p,] / sqrt( normalized_codes_som_all[p,] %*% normalized_codes_som_all[p,] )[1,1]
  }
  
  # Preliminary histogram plotting
  h1_t <- hist(mean_of_least1dpa_abandance,
               xlim=c(as.integer(range(mean_of_least1dpa_abandance)[1]),
                      as.integer(range(mean_of_least1dpa_abandance)[2])+1),
               breaks=seq(as.integer(range(mean_of_least1dpa_abandance)[1]),
                          as.integer(range(mean_of_least1dpa_abandance)[2])+1, 1))
  h2_t <- hist(mean_of_least3dpa_abandance,
               xlim=c(as.integer(range(mean_of_least3dpa_abandance)[1]),
                      as.integer(range(mean_of_least3dpa_abandance)[2])+1),
               breaks=seq(as.integer(range(mean_of_least3dpa_abandance)[1]),
                          as.integer(range(mean_of_least3dpa_abandance)[2])+1, 1))
  h3_t <- hist(mean_of_least10dpa_abandance,
               xlim=c(as.integer(range(mean_of_least10dpa_abandance)[1]),
                      as.integer(range(mean_of_least10dpa_abandance)[2])+1),
               breaks=seq(as.integer(range(mean_of_least10dpa_abandance)[1]),
                          as.integer(range(mean_of_least10dpa_abandance)[2])+1, 1))
  
  # Prepare output PDF
  tables_clusters_SOM               <- NULL
  tables_clusters_after_rmoutlier   <- NULL
  tables_clusters_after_rmunreliable<- NULL
  tables_clusters_after_rmlowabundance <- NULL
  
  pdf(paste(
    paste0(frac,"_proteins_hist_of_avg_abundance_and_interpolated_profile_SOM_",Sys.Date()), 
    "pdf", sep='.'), width=12,height=15)
  
  # Quick overlay hist of # proteins in >=1, >=3, >=10 DPAs
  cluster_data <- list()
  h1 <- hist(mean_of_least1dpa_abandance,
             main="Histogram: avg abundances of non-zero values",
             xlab="avg abundance",
             xlim=c(as.integer(range(mean_of_least1dpa_abandance)[1]),
                    as.integer(range(mean_of_least1dpa_abandance)[2])+1),
             ylim=c(0, max(h1_t$counts, h2_t$counts, h3_t$counts)+50),
             col="blue", border="black",
             breaks=seq(as.integer(range(mean_of_least1dpa_abandance)[1]),
                        as.integer(range(mean_of_least1dpa_abandance)[2])+1,1)
  )
  h2 <- hist(mean_of_least3dpa_abandance,
             xlim=c(as.integer(range(mean_of_least3dpa_abandance)[1]),
                    as.integer(range(mean_of_least3dpa_abandance)[2])+1),
             ylim=c(0, max(h1_t$counts,h2_t$counts,h3_t$counts)+50),
             col="green", border="black",
             breaks=seq(as.integer(range(mean_of_least3dpa_abandance)[1]),
                        as.integer(range(mean_of_least3dpa_abandance)[2])+1,1), add=TRUE
  )
  h3 <- hist(mean_of_least10dpa_abandance,
             xlim=c(as.integer(range(mean_of_least10dpa_abandance)[1]),
                    as.integer(range(mean_of_least10dpa_abandance)[2])+1),
             ylim=c(0,max(h1_t$counts,h2_t$counts,h3_t$counts)+50),
             col="orange", border="black",
             breaks=seq(as.integer(range(mean_of_least10dpa_abandance)[1]),
                        as.integer(range(mean_of_least10dpa_abandance)[2])+1,1), add=TRUE
  )
  
  l <- min(as.integer(range(mean_of_least10dpa_abandance)[1]),
           as.integer(range(mean_of_least1dpa_abandance)[1]),
           as.integer(range(mean_of_least3dpa_abandance)[1]))
  r <- max(as.integer(range(mean_of_least10dpa_abandance)[2])+1,
           as.integer(range(mean_of_least1dpa_abandance)[2])+1,
           as.integer(range(mean_of_least3dpa_abandance)[2])+1)
  
  axis(side=1, at=l:r, labels=l:r)
  text(h1$mids,h1$counts, labels=h1$counts, adj=c(0.5, -0.5), cex=0.8)
  text(h2$mids,h2$counts, labels=h2$counts, adj=c(0.5, 0.95), cex=0.8)
  text(h3$mids,h3$counts, labels=h3$counts, adj=c(0.5, 0.95), cex=0.8)
  legend("topright",
         c("≥1 DPA", "≥3 DPAs & ≥2 consecutive", "≥10 DPAs"),
         fill=c("blue","green","orange"))
  
  # Plot code vectors and mapping
  plot(SOM_wcc, type="codes", main="Codes")
  plot(SOM_wcc, type="mapping", pch=1, main="Count Plot", keepMargins=TRUE)
  
  # Retrieve ID info from fraction table
  prots_ID <- prot_notation_data[[frac]]
  
  # For each SOM cluster node
  for (c in names(sample_in_each_unit)){
    proteins <- sample_in_each_unit[[c]]
    data_c <- NULL
    # Collect the raw and interpolated profiles for each protein
    for (p in proteins){
      id_info <- prots_ID[p, ]
      profile_info_intero <- data_introp[, p]
      profile_info_raw <- protfrac_exp_data_T[, p]
      df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
      colnames(df_p) <- c("cluster_id", colnames(prots_ID),
                          paste(rownames(protfrac_exp_data_T),"initial",sep="_"),
                          paste(rownames(data_introp),"interpolated",sep="_"))
      rownames(df_p) <- p
      data_c <- rbind(data_c, df_p)
    }
    tables_clusters_SOM <- rbind(tables_clusters_SOM, data_c)
    
    if(length(proteins)<2){
      next
    }
    
    # Outlier detection by correlation threshold
    corr_c <- all_corr_matrix[proteins, proteins]
    corr_c <- corr_c[upper.tri(corr_c, diag=FALSE)]
    cut_off_c <- min(0.8, quantile(corr_c, 0.5))
    normcode_v <- normalized_codes_som_all[c,]
    
    # Plot raw data in this cluster
    par(mfrow=c(4,2))
    xVals <- 5:24
    if (length(proteins)>0){
      par(mar=c(5,4,4,8), xpd=TRUE)
      plot(xVals, Data[proteins[1],], type="b",xaxt="n",
           xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35), main=c)
      axis(1, at=xVals, labels=xVals)
      c1=1
      if (length(proteins)>1){
        for (i in 2:length(proteins)){
          c1=c1+1
          lines(xVals, Data[proteins[i],], type="b", col=c1)
        }
      }
    }
    
    # Show a histogram of average abundance for these proteins
    proteins_aboundance <- c()
    for (p in proteins) {
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]]>0])
      proteins_aboundance <- c(proteins_aboundance, p_ab)
    }
    if (length(proteins_aboundance)>0){
      h_ylim <- hist(proteins_aboundance,
                     breaks=seq(as.integer(range(proteins_aboundance)[1]),
                                as.integer(range(proteins_aboundance)[2])+1,1), 
                     plot=FALSE)
      par(mar=c(5,4,4,8), xpd=TRUE)
      hist(proteins_aboundance,
           main="Histogram: avg abundance (non-zero)",
           xlab="avg abundance",
           xlim=c(15,35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(as.integer(range(proteins_aboundance)[1]),
                      as.integer(range(proteins_aboundance)[2])+1,1))
    }
    
    # Corr-based removal
    prots_nonoutlier <- c()
    for (p in proteins){
      corrpc <- normalized_Data[p,] %*% normcode_v
      if (corrpc>cut_off_c){
        prots_nonoutlier <- c(prots_nonoutlier, p)
      }
    }
    
    if (length(prots_nonoutlier)>0){
      par(mar=c(5,4,4,8), xpd=TRUE)
      plot(xVals, Data[prots_nonoutlier[1],], type="b",xaxt="n",
           xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35),
           main=paste0(c,", after remove outlier, cutoff=", cut_off_c))
      axis(1, at=xVals, labels=xVals)
      c2=1
      if (length(prots_nonoutlier)>1){
        for (i in 2:length(prots_nonoutlier)){
          c2=c2+1
          lines(xVals, Data[prots_nonoutlier[i],], type="b", col=c2)
        }
      }
      data_c <- NULL
      for (p in prots_nonoutlier){
        id_info <- prots_ID[p,]
        profile_info_intero <- data_introp[,p]
        profile_info_raw <- protfrac_exp_data_T[,p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID),
                            paste(rownames(protfrac_exp_data_T),"initial",sep="_"),
                            paste(rownames(data_introp),"interoplated",sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmoutlier <- rbind(tables_clusters_after_rmoutlier, data_c)
    }
    
    # Additional histogram for the non-outlier set
    prots_nonoutlier_aboundance <- c()
    for (p in prots_nonoutlier){
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]]>0])
      prots_nonoutlier_aboundance <- c(prots_nonoutlier_aboundance, p_ab)
    }
    if (length(prots_nonoutlier_aboundance)>0){
      h <- hist(prots_nonoutlier_aboundance,
                breaks=seq(as.integer(range(prots_nonoutlier_aboundance)[1]),
                           as.integer(range(prots_nonoutlier_aboundance)[2])+1,1),
                plot=FALSE)
      par(mar=c(5,4,4,8), xpd=TRUE)
      hist(prots_nonoutlier_aboundance,
           main="Histogram: avg abundance (non-outlier proteins)",
           xlab="avg abundance",
           xlim=c(15,35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(as.integer(range(prots_nonoutlier_aboundance)[1]),
                      as.integer(range(prots_nonoutlier_aboundance)[2])+1,1))
    }
    
    # Unreliable pattern check
    prots_unreliable <- c()
    for (p in prots_nonoutlier){
      raw_p <- protfrac_exp_data_T[[p]]
      # A naive pattern check, e.g. F,T,F in first or last 3 points
      if ( all((raw_p[1:3]>0) == c(FALSE,TRUE,FALSE)) || all((raw_p[18:20]>0) == c(FALSE,TRUE,FALSE)) ){
        prots_unreliable <- c(prots_unreliable,p)
      }
      # Checking for 3 consecutive zeros in the interpolated version
      intro_p <- data_introp[[p]]
      stack_p <- c()
      for (idx in seq_along(intro_p)){
        if (intro_p[idx]==0){
          if(length(stack_p)==0){
            stack_p <- c(stack_p, idx)
          } else {
            if ((idx - stack_p[length(stack_p)])==1){
              stack_p <- c(stack_p, idx)
            } else {
              if (length(stack_p)==3){
                prots_unreliable <- c(prots_unreliable,p)
                break
              }
              stack_p <- c(idx)
            }
          }
        }
      }
      if (length(stack_p)==3){
        prots_unreliable <- c(prots_unreliable,p)
      }
    }
    
    # Another pass: 0-non0-0 pattern check
    prots_unreliable_again <- c()
    for (p in prots_nonoutlier){
      intro_p <- c(0, data_introp[[p]], 0)
      for (i in 1:(length(intro_p)-2)){
        if (intro_p[i]==0 && intro_p[i+1]>0 && intro_p[i+2]==0){
          prots_unreliable_again <- c(prots_unreliable_again,p)
          break
        }
      }
    }
    prots_all_unreliable <- unique(c(prots_unreliable, prots_unreliable_again))
    prots_reliable <- setdiff(prots_nonoutlier, prots_all_unreliable)
    
    # Plot the final "reliable" subset
    if (length(prots_reliable)>0){
      par(mar=c(5,4,4,8), xpd=TRUE)
      plot(xVals, Data[prots_reliable[1],], type="b",
           xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35),
           main=paste0(c,", after remove unreliable"))
      axis(1, at=xVals, labels=xVals)
      c3=1
      if (length(prots_reliable)>1){
        for (i in 2:length(prots_reliable)){
          c3=c3+1
          lines(xVals, Data[prots_reliable[i],], type="b", col=c3)
        }
      }
    }
    
    # Collect them
    if (length(prots_reliable)>0){
      data_c <- NULL
      for (p in prots_reliable){
        id_info <- prots_ID[p,]
        profile_info_intero <- data_introp[,p]
        profile_info_raw <- protfrac_exp_data_T[,p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID),
                            paste(rownames(protfrac_exp_data_T),"initial",sep="_"),
                            paste(rownames(data_introp),"interpolated",sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmunreliable <- rbind(tables_clusters_after_rmunreliable, data_c)
    }
    
    # Quick histogram for reliable subset
    prots_reliable_aboundance <- c()
    for (p in prots_reliable){
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]]>0])
      prots_reliable_aboundance <- c(prots_reliable_aboundance, p_ab)
    }
    if (length(prots_reliable_aboundance)>0){
      h <- hist(prots_reliable_aboundance,
                breaks=seq(as.integer(range(prots_reliable_aboundance)[1]),
                           as.integer(range(prots_reliable_aboundance)[2])+1,1),
                plot=FALSE)
      par(mar=c(5,4,4,8), xpd=TRUE)
      hist(prots_reliable_aboundance,
           main="Histogram: avg abundance (reliable proteins)",
           xlab="avg abundance",
           xlim=c(15,35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(as.integer(range(prots_reliable_aboundance)[1]),
                      as.integer(range(prots_reliable_aboundance)[2])+1,1))
    }
    
    # Next, pick out "high abundance" among the reliable
    prots_reliable_high_abundance <- intersect(prots_reliable, prots_high_abundance)
    if (length(prots_reliable_high_abundance)>0){
      par(mar=c(5,4,4,8), xpd=TRUE)
      plot(xVals, Data[prots_reliable_high_abundance[1],],
           type="b", xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35),
           main=paste0(c,", after remove low abundance"))
      axis(1, at=xVals, labels=xVals)
      c3=1
      if (length(prots_reliable_high_abundance)>1){
        for (i in 2:length(prots_reliable_high_abundance)){
          c3=c3+1
          lines(xVals, Data[prots_reliable_high_abundance[i],], type="b", col=c3)
        }
      }
    }
    
    # Collect them
    if (length(prots_reliable_high_abundance)>0){
      data_c <- NULL
      for (p in prots_reliable_high_abundance){
        id_info <- prots_ID[p,]
        profile_info_intero <- data_introp[,p]
        profile_info_raw <- protfrac_exp_data_T[,p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID),
                            paste(rownames(protfrac_exp_data_T),"initial",sep="_"),
                            paste(rownames(data_introp),"interpolated",sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmlowabundance <- rbind(tables_clusters_after_rmlowabundance, data_c)
    }
    
    # A histogram for the high abundance subset
    # ...
  }
  dev.off()  # End PDF output
  
  # Identify the 'informative' cluster IDs from the fraction's score table
  score_frac <- prot_group_score_data[[frac]]
  informative_id <- filter(score_frac, score=="I")$CID
  
  # Subset final 'informative' table
  tables_clusters_informative <- NULL
  for (id in informative_id){
    dt <- filter(tables_clusters_after_rmunreliable, cluster_id==id)
    dt <- mutate(dt, cluster_id=paste(frac, cluster_id, sep="_"))
    tables_clusters_informative <- rbind(tables_clusters_informative, dt)
  }
  
  # Save stepwise results
  SOM_clustering_result[[frac]] <- list()
  SOM_codes_result[[frac]] <- list()
  
  write.csv(tables_clusters_SOM,
            paste(paste0(frac,"_Table1_step1_SOM_clustering_result_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_clustering_result[[frac]][["step1"]] <- tables_clusters_SOM
  
  write.csv(tables_clusters_after_rmoutlier,
            paste(paste0(frac,"_Table2_step2_SOM_clustering_outliers_removed_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_clustering_result[[frac]][["step2"]] <- tables_clusters_after_rmoutlier
  
  write.csv(tables_clusters_after_rmunreliable,
            paste(paste0(frac,"_Table3_step3_SOM_clustering_unreliable_removed_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_clustering_result[[frac]][["step3"]] <- tables_clusters_after_rmunreliable
  
  write.csv(tables_clusters_after_rmlowabundance,
            paste(paste0(frac,"_Table4_step3_SOM_clustering_lowAbundance_removed_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_clustering_result[[frac]][["step4"]] <- tables_clusters_after_rmlowabundance
  
  write.csv(tables_clusters_informative,
            paste(paste0(frac,"_Table5_step4_SOM_informativeGroups_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_clustering_result[[frac]][["step5"]] <- tables_clusters_informative
  
  # High abundance among informative
  tables_clusters_informative_high_abundance <- tables_clusters_informative[prots_high_abundance,]
  tables_clusters_informative_high_abundance <- filter(tables_clusters_informative_high_abundance, !is.na(cluster_id))
  
  # Some indexing logic for cluster_id substring
  indexOfcluster <- ifelse(frac=="APO", 6, 7)
  
  # Sort by cluster ID numeric
  tables_clusters_informative_high_abundance <- arrange(
    tables_clusters_informative_high_abundance,
    as.numeric(substring(cluster_id, indexOfcluster, indexOfcluster+1))
  )
  numOfCol <- ncol(tables_clusters_informative_high_abundance)
  
  # Normalize the last 20 columns
  tables_clusters_informative_temp <- tables_clusters_informative_high_abundance[, c(1,(numOfCol-19):numOfCol)]
  for (i in seq_len(nrow(tables_clusters_informative_temp))){
    temp_i <- tables_clusters_informative_temp[i,2:21]
    for (j in seq(2,21)){
      tables_clusters_informative_temp[i,j] <- tables_clusters_informative_temp[i,j]/max(temp_i)
    }
  }
  tables_clusters_informative_high_abundance_normalized <- tables_clusters_informative_temp
  
  # Compute means across the 20 columns for each cluster
  tables_mean_each_cluster <- NULL
  for (c in levels(factor(tables_clusters_informative_high_abundance_normalized$cluster_id,
                          levels=unique(tables_clusters_informative_high_abundance_normalized$cluster_id)))){
    temp_c  <- filter(tables_clusters_informative_high_abundance_normalized, cluster_id==c)[,-1]
    temp_df <- data.frame(colMeans(temp_c))
    names(temp_df) <- c
    tempdf <- data.frame(t(temp_df))
    tables_mean_each_cluster <- rbind(tables_mean_each_cluster, tempdf)
  }
  mean_of_normalized_informative[[frac]] <- tables_mean_each_cluster
  
  # Save outputs
  write.csv(tables_clusters_informative_high_abundance,
            paste(paste0(frac,"_Table_informativeGroups_highAbundance_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  write.csv(tables_clusters_informative_high_abundance_normalized,
            paste(paste0(frac,"_Table_informativeGroups_highAbundance_normalized_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  write.csv(tables_mean_each_cluster,
            paste(paste0(frac,"_Table_meanOfNormalizedInformative_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  
  # Retrieve SOM codes in a data frame
  code_som      <- data.frame(codes_som_all)
  informative_code <- code_som[informative_id,]
  rownames(informative_code) <- paste(frac, informative_id, sep="_")
  
  write.csv(code_som,
            paste(paste0(frac,"_Table5_SOM_codes_profile_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_codes_result[[frac]][["codes"]] <- code_som
  
  write.csv(informative_code,
            paste(paste0(frac,"_Table6_SOM_codes_profile_informative_",Sys.Date()), "csv", sep='.'), row.names=TRUE)
  SOM_codes_result[[frac]][["informative_codes"]] <- informative_code
}  # End fraction loop

##########################################################
# Merge the "informative" means across fractions
##########################################################
workingDir = "/path/to/your/merge_informative_folder"
setwd(workingDir)

meannorm_informative <- NULL
for (frac in fraction_names){
  meannorm_informative <- rbind(meannorm_informative, mean_of_normalized_informative[[frac]])
}
Data_meannorm <- as.matrix(meannorm_informative)

# Function merges the codes from multiple fractions, runs another SOM to identify cross-fraction patterns
merge_informative_codes <- function(som_grid_x, som_grid_y, Data, codes_informative){
  
  x_limit <- som_grid_x
  y_limit <- som_grid_y
  grid_size <- paste0("(SOM grid size of ", x_limit*y_limit, ")")
  set.seed(123)
  
  system.time(
    SOM_wcc <- supersom(
      Data,
      grid = somgrid(xdim=x_limit, ydim=y_limit, topo="rectangular",
                     neighbourhood.fct="gaussian", toroidal=FALSE),
      dist.fcts = c("WCCd0"),
      rlen=500,
      normalizeDataLayers=FALSE
    )
  )
  
  prot_id <- rownames(Data)
  label_of_sample <- SOM_wcc$unit.classif
  cluster_ID <- levels(factor(label_of_sample))
  
  sample_in_each_unit <- list()
  for (i in seq_len(x_limit*y_limit)){
    postion_index <- (label_of_sample==i)
    name <- paste0('V', as.character(i))
    if (any(postion_index)){
      sample_in_each_unit[[name]] <- prot_id[postion_index]
    } else {
      sample_in_each_unit[[name]] <- NULL
    }
  }
  
  # Compile cluster-coded table
  clustering_inform_code <- NULL
  for (c in names(sample_in_each_unit)){
    df_c <- codes_informative[sample_in_each_unit[[c]],]
    df_c[["cluster_ID"]] <- rep(c, nrow(df_c))
    clustering_inform_code <- rbind(clustering_inform_code, df_c)
  }
  
  # Export cluster membership
  write.csv(
    clustering_inform_code,
    paste(
      paste0("Table_SOM_clustering_of_informative_meannorm_APO_S200_P200_", grid_size),
      Sys.Date(),
      "csv", sep='_'
    ),
    row.names=TRUE
  )
  
  # Save codebook vectors
  codes_merge <- SOM_wcc$codes[[1]]
  write.csv(
    codes_merge,
    paste(
      paste0("Table_SOM_codes_of_merging_informative_meannorm_APO_S200_P200_",grid_size),
      Sys.Date(),"csv",sep='_'
    ),
    row.names=TRUE
  )
  
  # Build data for line plots
  plot_data <- list()
  for (v in names(sample_in_each_unit)){
    df <- codes_informative[sample_in_each_unit[[v]],]
    value <- c()
    type_ <- c()
    DPA <- c()
    names_ <- c()
    for (r in rownames(df)){
      zx <- t(df[r,])
      cc <- zx[, colnames(zx)]
      value <- c(value, cc)
      type_ <- c(type_, rep(strsplit(r, "_")[[1]][1], 20))
      names_ <- c(names_, rep(r, 20))
      DPA <- c(DPA, 5:24)
    }
    # Also append the SOM code
    value <- c(value, SOM_wcc$codes[[1]][v,])
    type_ <- c(type_, rep(v, 20))
    names_ <- c(names_, rep(v, 20))
    DPA <- c(DPA, 5:24)
    plot_data[[v]] <- data.frame(value=value, DPA=DPA, Fraction=type_, Code=names_)
  }
  
  # Plot to PDF
  pdf(
    paste(
      paste0("Merge_informative_meannorm_plot_APO_P200_S200_", grid_size,"_", Sys.Date()),
      "pdf", sep='.'
    ),
    width=12, height=8
  )
  
  plot(SOM_wcc, type="codes", main="Codes")
  plot(SOM_wcc, type="mapping", pch=1, main="Count Plot", keepMargins=TRUE)
  
  z = c("APO","P200","S200","V")
  color_set <- c("#ed5151","#149ece","#a7c636","black")
  
  for (V in names(plot_data)){
    df_t <- plot_data[[V]]
    colors_t <- levels(factor(df_t$Fraction))
    # Move "V" last if it appears
    if (substr(colors_t[length(colors_t)],1,1)=="V"){
      colors_t <- c(colors_t[-length(colors_t)], "V")
    }
    colors <- color_set[z %in% colors_t]
    lineso <- c(2:(length(levels(factor(df_t$Code)))),1)
    plot_p <- ggplot(df_t, aes(x=DPA, y=value, linetype=Code, color=Fraction,
                               group=interaction(Code, Fraction)))+
      geom_line()+
      scale_linetype_manual(values=lineso)+
      scale_colour_manual(values=colors)+
      labs(
        title=V,
        subtitle=paste0("{",paste(sample_in_each_unit[[V]],collapse=";"),"}")
      )+
      scale_x_continuous(breaks=seq(5,24,1), limits=c(5,24))+
      scale_y_continuous(breaks=seq(0,1,0.1), limits=c(0,1))+
      theme(
        plot.title=element_text(hjust=0.5, vjust=12, margin=margin(t=40,b=-30)),
        title=element_text(size=14),
        axis.text.x=element_text(size=10),
        axis.text.y=element_text(size=10),
        axis.title.x=element_text(size=12),
        axis.title.y=element_text(size=12)
      )+
      theme(
        legend.title=element_text(size=8, face=8),
        legend.text=element_text(size=8),
        legend.position="bottom"
      )
    print(plot_p)
  }
  dev.off()
}

# Example usage
x_lim <- 3
y_lim <- 3
merge_informative_codes(som_grid_x=x_lim, som_grid_y=y_lim,
                        Data=Data_meannorm,
                        codes_informative=meannorm_informative)

# Another set of grids (5 x 4, 5 x 6, 5 x 8)
x_lim <- 5
for (y_lim in c(4,6,8)){
  merge_informative_codes(som_grid_x=x_lim,
                          som_grid_y=y_lim,
                          Data=Data_meannorm,
                          codes_informative=meannorm_informative)
}

# Additional example: combining codes_informative from each fraction
# (If needed for further final merges)
# ...
# End of Script
# ---------------------------------------------------------------------------------------
