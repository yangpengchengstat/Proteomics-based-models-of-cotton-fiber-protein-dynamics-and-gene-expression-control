rm(list = ls())#clean the working space
library(kohonen) # self organizing map
library(tempR)
library(dplyr)
library(ggplot2)
library(openxlsx)
library(tibble)
library(tidyverse) # read CVS
require(Rcpp) # help to integrate R and C++ via R functions
sourceCpp(system.file("Distances", "wcc0.cpp", package = "kohonen",  mustWork = TRUE)) # expose C++ function wcc0 to R level
getwd()
workingDir = "/Users/pengchengyang/cotton pj/input/Input Data"
setwd(workingDir)
options(stringsAsFactors = FALSE)


########################################
# 1 Data importation
########################################
histplot_ylim <- c(0, 25)
# Fraction data importation
fraction_names <- c("APO", "P200", "S200") # Three subcellular fractions
prot_exp_data <- list()
prot_group_score_data <- list()
prot_notation_data <- list()
for (frac in fraction_names){
  prot_exp_table <- read.xlsx(paste0(frac,".xlsx")) # read in the protein expression (Protein Abundances) data
  prot_exp_table <- prot_exp_table %>% remove_rownames %>% column_to_rownames(var=names(prot_exp_table)[1]) %>% as.data.frame()
  prot_exp_data[[frac]] <- prot_exp_table
  rm(prot_exp_table)
  score_table <- read.xlsx(paste0(paste0("score_",frac),".xlsx"))
  prot_group_score_data[[frac]] <- score_table # read in the evaluation of the 
  rm(score_table)
  prots_noation_table <- read.xlsx(paste0(frac,"_name_id.xlsx")) # read in the relevant information: like gene ID, name
  prots_noation_table <- prots_noation_table %>% remove_rownames %>% column_to_rownames(var=names(prots_noation_table)[2]) %>% as.data.frame()
  prot_notation_data[[frac]] <- prots_noation_table
  rm(prots_noation_table)
}
# Prot and mRNA ID matching data
prot_mRNA_ID_map_data <- read.xlsx("proteinIDmap.xlsx")
# mRNA data
mRNA_data <- read.xlsx("mRNA seq exp data.xlsx")
mRNA_data <- mRNA_data %>% remove_rownames %>% column_to_rownames(var=names(mRNA_data)[1]) %>% as.data.frame()
#######################################
# 2 SOM clustering analysis on fractions
#######################################
SOM_clustering_result <- list()
SOM_codes_result <- list()
SOM_prot_cluster <- list()
mean_of_normalized_informative <- list()
histplot_ylim <- c(0, 25)
for (frac in fraction_names){
  workingDir = paste0("/Users/pengchengyang/cotton pj/output/SOM clustering on data of 3 fractions/",frac)
  setwd(workingDir)
  if (frac == "APO"){
    x_limit <- 5
    y_limit <- 6
    abundance_cutoff <- 20
  }else if (frac == "S200"){
    x_limit <- 5
    y_limit <- 4
    abundance_cutoff <- 20
  }else if (frac == "P200"){
    x_limit <- 5
    y_limit <- 8
    abundance_cutoff <- 19
  }else{
    print("Wrong names.")
  }
  if(frac=="S200"){
    histplot_ylim <- c(0, 25)
  }else{
    histplot_ylim <- c(0, 70)
  }
  protfrac_exp_data <- prot_exp_data[[frac]]
  protfrac_exp_data_T <- as.data.frame(t(protfrac_exp_data))
  # prots with abundance > cutoff
  prots_high_abundance <- c()
  for (p in names(protfrac_exp_data_T)) {
    p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
    if (p_ab >= abundance_cutoff) {
      prots_high_abundance <- c(prots_high_abundance, p)
    }
  }
  num_dpa <- dim(protfrac_exp_data_T)[1]
  f_name <- paste(paste("LOG OF COMPUTATION on", frac, sep=" "),
                  "txt", sep=".")
  cat("LOG OF COMPUTATION\n\n", file=f_name)
  cat("1: There are",  dim(protfrac_exp_data_T)[2], " proteins in data set of ", frac, ".\n",
      file=f_name, append=TRUE)
  
  prots_data_least1 <- data.frame(a=1:num_dpa) #initialize space
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
  
  prots_data_least3 <- data.frame(a=1:num_dpa) #initialize space
  count_gt03 <- 0
  for (i in names(protfrac_exp_data_T)){
    df_i <- protfrac_exp_data_T[i]
    df_i_temp <- df_i[,1]
    if (sum(df_i_temp>0)>2){
      count_gt03 <- count_gt03 + 1
      check_num_gt0 <- df_i>0
      index_gt0 <- which(check_num_gt0==T)
      if (any((index_gt0[2:length(index_gt0)]-index_gt0[1:(length(index_gt0)-1)])==1)){
        prots_data_least3 <- cbind(prots_data_least3, df_i)
      }
    }
  }
  prots_data_least3 <- prots_data_least3[,-1]
  cat("2: There are",  count_gt03, " proteins present in ≥ 3 different DPAs.\n",
      file=f_name, append=TRUE)
  cat("3: There are", dim(prots_data_least3)[2],
      " proteins present in ≥ 3 different DPAs and at least 2 consecutive DPAs.\n", 
      file=f_name, append=TRUE)
  
  
  mean_of_least3dpa_abandance <- c()
  for (p in names(prots_data_least3)){
    profile_p <- prots_data_least3[,p]
    profile_p_non0 <- profile_p[profile_p>0]
    mean_of_least3dpa_abandance <- c(mean_of_least3dpa_abandance, mean(profile_p_non0))
  }
  
  prots_data_least10 <- data.frame(a=1:num_dpa)
  for (i in names(protfrac_exp_data_T)){
    df_i <- protfrac_exp_data_T[i]
    df_i_temp <- df_i[,1]
    if (sum(df_i_temp>0)>9){
      prots_data_least10 <- cbind(prots_data_least10, df_i)
    }
  }
  prots_data_least10 <- prots_data_least10[,-1]
  cat("3: There are", dim(prots_data_least10)[2],
      " proteins present in 10 3 different DPAs.\n", 
      file=f_name, append=TRUE)
  
  mean_of_least10dpa_abandance <- c()
  for (p in names(prots_data_least10)){
    profile_p <- prots_data_least10[,p]
    profile_p_non0 <- profile_p[profile_p>0]
    mean_of_least10dpa_abandance <- c(mean_of_least10dpa_abandance, mean(profile_p_non0))
  }
  
  # data interpolation
  data_introp <- data.frame(a=1:num_dpa)
  interpolation_mid <- c()
  interpolation4 <- c()
  interpolation_0_1 <- c()
  interpolation_0_2 <- c()
  for (i in names(prots_data_least3)){
    prot_d <- prots_data_least3[, i]
    df_i <- prots_data_least3[i]
    check_0 <- prot_d==0
    for (j in 2:(length(check_0)-1)){
      if (all(c(check_0[j]==T, check_0[j-1]==F, check_0[j+1]==F))){
        df_i[j,i] <- (df_i[j-1,i] + df_i[j+1,i])/2
        interpolation_mid <- c(interpolation_mid, i)
      }
    }
    for (j in 4:length(check_0)){
      if (all(c(check_0[j-3]==F, check_0[j-1]==T, check_0[j-2]==T, check_0[j]==F))){
        slope_j <- (df_i[j,i]-df_i[j-3,i])/3
        df_i[j-2,i] <- slope_j+df_i[j-3,i]
        df_i[j-1,i] <- slope_j*2+df_i[j-3,i]
        interpolation4 <- c(interpolation4, i)
      }
    }
    check_0 <- df_i[,i]==0
    if (check_0[1]==T && check_0[2]==F){
      df_i[1,i] <- df_i[2,i]
      interpolation_0_1 <- c(interpolation_0_1, i)
    }
    if (check_0[length(check_0)]==T && check_0[length(check_0)-1]==F){
      df_i[length(check_0),i] <- df_i[length(check_0)-1,i]
    }
    if (all(c(check_0[1]==T, check_0[2]==T, check_0[3]==F, check_0[4]==F))){
      df_i[1,i] <- df_i[3,i]
      df_i[2,i] <- df_i[4,i]
      interpolation_0_2 <- c(interpolation_0_2, i)
    }else if (all(c(check_0[1]==T, check_0[2]==T, check_0[3]==F, check_0[4]==T))){
      df_i[1,i] <- df_i[3,i]
      df_i[2,i] <- df_i[3,i]
    }
    if (all(c(check_0[length(check_0)]==T, check_0[length(check_0)-1]==T, check_0[length(check_0)-2]==F, check_0[length(check_0)-3]==F))){
      df_i[length(check_0),i] <- df_i[length(check_0)-2,i]
      df_i[length(check_0)-1,i] <- df_i[length(check_0)-3,i]
      #print(paste0(2,i))
    }else if (all(c(check_0[length(check_0)]==T, check_0[length(check_0)-1]==T, check_0[length(check_0)-2]==F, check_0[length(check_0)-3]==T))){
      df_i[length(check_0),i] <- df_i[length(check_0)-2,i]
      df_i[length(check_0)-1,i] <- df_i[length(check_0)-2,i]
      #print(paste0(1,i))
    }
    data_introp <- cbind(data_introp, df_i)
  }
  data_introp <- data_introp[,-1]
  
  # SOM Clustering
  Data <- as.matrix(t(data_introp))
  normalized_Data <- Data
  for (p in rownames(normalized_Data)){
    normalized_Data[p,] <- normalized_Data[p,]/sqrt(normalized_Data[p,] %*% normalized_Data[p,])[1,1]
  }
  all_corr_matrix <- normalized_Data %*% t(normalized_Data)
  all_corr <- all_corr_matrix[upper.tri(all_corr_matrix, diag = FALSE)]
  
  
  set.seed(123)
  
  system.time(SOM_wcc <- supersom(Data,
                                  grid=somgrid(xdim = x_limit, ydim = y_limit, topo = c("rectangular"),
                                               neighbourhood.fct = "gaussian", toroidal = F),
                                  dist.fcts = c("WCCd0"), rlen=500, normalizeDataLayers = FALSE))
  
  prot_id <- rownames(Data)
  label_of_sample <- SOM_wcc$unit.classif
  cluster_ID <- levels(factor(label_of_sample))
  
  sample_in_each_unit <- list()
  for (i in c(1:(x_limit*y_limit))){
    postion_index <- c(label_of_sample==i)
    name <- paste0('V',as.character(i), sep="")
    if (length(label_of_sample[postion_index])>0){
      sample_in_each_unit[[name]] <- prot_id[postion_index]
    }else{
      sample_in_each_unit[[name]] <- NULL
    }
  }
  SOM_prot_cluster[[frac]] <- sample_in_each_unit
  codes_som_all <- SOM_wcc$codes[[1]]
  normalized_codes_som_all <- codes_som_all
  for (p in rownames(normalized_codes_som_all)){
    normalized_codes_som_all[p,] <- normalized_codes_som_all[p,]/sqrt(normalized_codes_som_all[p,] %*% normalized_codes_som_all[p,])[1,1]
  }
  
  h1_t <- hist(mean_of_least1dpa_abandance,
               xlim=c(as.integer(range(mean_of_least1dpa_abandance)[1]), 
                      as.integer(range(mean_of_least1dpa_abandance)[2])+1),
               breaks=seq(from =as.integer(range(mean_of_least1dpa_abandance)[1]),
                          to =as.integer(range(mean_of_least1dpa_abandance)[2])+1, by = 1))
  h2_t <- hist(mean_of_least3dpa_abandance,
               xlim=c(as.integer(range(mean_of_least3dpa_abandance)[1]), 
                      as.integer(range(mean_of_least3dpa_abandance)[2])+1),
               breaks=seq(from =as.integer(range(mean_of_least3dpa_abandance)[1]),
                          to =as.integer(range(mean_of_least3dpa_abandance)[2])+1, by = 1))
  h3_t <- hist(mean_of_least10dpa_abandance,
               xlim=c(as.integer(range(mean_of_least10dpa_abandance)[1]), 
                      as.integer(range(mean_of_least10dpa_abandance)[2])+1),
               breaks=seq(from=as.integer(range(mean_of_least10dpa_abandance)[1]),
                          to=as.integer(range(mean_of_least10dpa_abandance)[2])+1, by = 1))
  
  tables_clusters_SOM <- NULL
  tables_clusters_after_rmoutlier <- NULL
  tables_clusters_after_rmunreliable <- NULL
  tables_clusters_after_rmlowabundance <- NULL
  pdf(paste(paste(paste0(frac," proteins ", "hist of avg abundance and interplated profile plot for SOM clusters"), Sys.Date(), sep='_'),
            "pdf", sep='.'), width=12,height=15)
  
  cluster_data <- list()
  h1 <- hist(mean_of_least1dpa_abandance,
             main="Histogram of the average abundances of non-zero values ",
             xlab="average abundance",
             xlim=c(as.integer(range(mean_of_least1dpa_abandance)[1]), 
                    as.integer(range(mean_of_least1dpa_abandance)[2])+1),
             ylim=c(0, max(h1_t$counts, h2_t$counts, h3_t$counts)+50),
             col="blue",
             border="black",
             breaks=seq(from =as.integer(range(mean_of_least1dpa_abandance)[1]),
                        to =as.integer(range(mean_of_least1dpa_abandance)[2])+1, by = 1)
  )
  
  h2 <- hist(mean_of_least3dpa_abandance,
             main="Histogram of the average abundances of non-zero values ",
             xlab="average abundance",
             xlim=c(as.integer(range(mean_of_least3dpa_abandance)[1]), 
                    as.integer(range(mean_of_least3dpa_abandance)[2])+1),
             ylim=c(0, max(h1_t$counts, h2_t$counts, h3_t$counts)+50),
             col="green",
             border="black",
             breaks=seq(from =as.integer(range(mean_of_least3dpa_abandance)[1]),
                        to =as.integer(range(mean_of_least3dpa_abandance)[2])+1, by = 1), add=T
  )
  
  h3 <- hist(mean_of_least10dpa_abandance,
             main="Histogram of the average abundances of non-zero values ",
             xlab="average abundance",
             xlim=c(as.integer(range(mean_of_least10dpa_abandance)[1]), 
                    as.integer(range(mean_of_least10dpa_abandance)[2])+1),
             ylim=c(0, max(h1_t$counts, h2_t$counts, h3_t$counts)+50),
             col="orange",
             border="black",
             breaks=seq(from=as.integer(range(mean_of_least10dpa_abandance)[1]),
                        to=as.integer(range(mean_of_least10dpa_abandance)[2])+1, by = 1), add=T
  )
  l <- min(as.integer(range(mean_of_least10dpa_abandance)[1]),
          as.integer(range(mean_of_least1dpa_abandance)[1]),
          as.integer(range(mean_of_least3dpa_abandance)[1]))
  r <- max(as.integer(range(mean_of_least10dpa_abandance)[2])+1,
         as.integer(range(mean_of_least1dpa_abandance)[2])+1,
         as.integer(range(mean_of_least3dpa_abandance)[2])+1)
  axis(side=1, at = l:r, labels = l:r)
  text(h1$mids,h1$counts,labels=h1$counts, adj=c(0.5, -0.5), cex=0.8)
  text(h2$mids,h2$counts,labels=h2$counts, adj=c(0.5, 0.95), cex=0.8)
  text(h3$mids,h3$counts,labels=h3$counts, adj=c(0.5, 0.95), cex=0.8)
  legend("topright", c("at least 1 DPA", ">=3 and >= 2 consecutive DPAs", ">=10 and >= 2 consecutive DPAs"),
         fill=c("blue", "green", "orange"))
  
  plot(SOM_wcc, type = "codes", main="Codes")
  plot(SOM_wcc, type = "mapping",
       pch = 1, main = "Count plot", keepMargins = TRUE)
  prots_ID <- prot_notation_data[[frac]]
  for (c in names(sample_in_each_unit)){
    proteins <- sample_in_each_unit[[c]]
    data_c <- NULL
    for (p in proteins){
      id_info <- prots_ID[p, ]
      profile_info_intero <- data_introp[, p]
      profile_info_raw <- protfrac_exp_data_T[, p]
      df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
      colnames(df_p) <- c("cluster_id", colnames(prots_ID), 
                          paste(rownames(protfrac_exp_data_T), "initial", sep="_"), 
                          paste(rownames(data_introp), "interoplated", sep="_"))
      rownames(df_p) <- p
      data_c <- rbind(data_c, df_p)
    }
    tables_clusters_SOM <- rbind(tables_clusters_SOM, data_c)
    
    if(length(proteins)<2){
      next
    }
    
    corr_c <- all_corr_matrix[proteins, proteins]
    corr_c <- corr_c[upper.tri(corr_c, diag = FALSE)]
    #cut_off_c <- quantile(corr_c, 0.5)
    cut_off_c <- min(0.8, quantile(corr_c, 0.5))
    normcode_v <- normalized_codes_som_all[c,]
    #par(cex=0.7, mai=c(0.1,0.1,0.2,0.1))
    par(mfrow=c(4,2))
    x = 5:24
    if (length(proteins)>0){
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      plot(x, Data[proteins[1], ], type="b",xaxt="n",xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35), main=c)
      axis(1, at=x,labels=x)
      c1 = 1
      if (length(proteins)>1){
        for (i in 2:length(proteins)){
          c1 <- c1 + 1
          lines(x, Data[proteins[i], ], type="b", col=c1)
        }
      }
    }
    proteins_aboundance <- c()#prots_high_abundance <- c()
    for (p in proteins) {
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
      proteins_aboundance <- c(proteins_aboundance, p_ab)
    }
    if (length(proteins_aboundance) > 0) {
      h_ylim <- hist(proteins_aboundance, breaks=seq(from=as.integer(range(proteins_aboundance)[1]),
                                                to=as.integer(range(proteins_aboundance)[2])+1, by = 1), plot=FALSE)
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      hist(proteins_aboundance,main="Histogram of the average abundances of non-zero values ",
           xlab="average abundance",
           xlim=c(15, 35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(from=as.integer(range(proteins_aboundance)[1]),
                      to=as.integer(range(proteins_aboundance)[2])+1, by = 1))
    }
    
    
    prots_nonoutlier <- c()
    for (p in proteins){
      corrpc <- normalized_Data[p,] %*% normcode_v
      if (corrpc >cut_off_c){
        prots_nonoutlier <- c(prots_nonoutlier, p)
      }
    }
    if (length(prots_nonoutlier)>0){
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      plot(x, Data[prots_nonoutlier[1], ], type="b",xaxt="n", xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35), 
           main=paste(paste(c,"after remove outlier", sep = ", "), cut_off_c, sep=" and cutoff is "))
      axis(1, at=x,labels=x)
      c2 = 1
      if (length(prots_nonoutlier)>1){
        for (i in 2:length(prots_nonoutlier)){
          c2 <- c2 + 1
          lines(x, Data[prots_nonoutlier[i], ], type="b", col=c2)
        }
      }
      data_c <- NULL
      for (p in prots_nonoutlier){
        id_info <- prots_ID[p, ]
        profile_info_intero <- data_introp[, p]
        profile_info_raw <- protfrac_exp_data_T[, p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID), 
                            paste(rownames(protfrac_exp_data_T), "initial", sep="_"), 
                            paste(rownames(data_introp), "interoplated", sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmoutlier <- rbind(tables_clusters_after_rmoutlier, data_c)
    }
    prots_nonoutlier_aboundance <- c()#prots_high_abundance <- c()
    for (p in prots_nonoutlier) {
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
      prots_nonoutlier_aboundance <- c(prots_nonoutlier_aboundance, p_ab)
    }
    if (length(prots_nonoutlier_aboundance) > 0) {
      h <- hist(prots_nonoutlier_aboundance,breaks=seq(from=as.integer(range(prots_nonoutlier_aboundance)[1]),
                                                       to=as.integer(range(prots_nonoutlier_aboundance)[2])+1, by = 1), plot=FALSE)
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      hist(prots_nonoutlier_aboundance, main="Histogram of the average abundances of non-outlier proteins ",
           xlab="average abundance",
           xlim=c(15, 35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(from=as.integer(range(prots_nonoutlier_aboundance)[1]),
                      to=as.integer(range(prots_nonoutlier_aboundance)[2])+1, by = 1))
    }
    
    prots_unreliable <- c()
    for (p in prots_nonoutlier){
      raw_p <- protfrac_exp_data_T[[p]]
      if (all((raw_p[1:3]>0) == c(F,T,F)) | all((raw_p[18:20]>0) == c(F,T,F))){
        prots_unreliable <- c(prots_unreliable,p)
      }
      intro_p <- data_introp[[p]]
      stack_p <- c()
      for (i in 1:num_dpa){
        if (intro_p[i] == 0){
          if(length(stack_p)==0){
            stack_p <- c(stack_p, i)
          }else{
            if (i-stack_p[length(stack_p)]==1){
              stack_p <- c(stack_p, i)
            }else{
              if (length(stack_p)==3){
                prots_unreliable <- c(prots_unreliable,p)
                break
              }
              stack_p <- c(i)
            }
          }
        }
      }
      if (length(stack_p)==3){
        prots_unreliable <- c(prots_unreliable,p)
        
      }
    }
    prots_unreliable_again <- c()
    for (p in prots_nonoutlier){
      intro_p <- c(0, data_introp[[p]], 0)
      for (i in 1:(length(intro_p)-2)){
        if (all(c(intro_p[i]==0, intro_p[i+1]>0, intro_p[i+2]==0))){
          prots_unreliable_again <- c(prots_unreliable_again, p)
          break
        }
      }
    }
    prots_all_unreliable <- levels(factor(c(prots_unreliable, prots_unreliable_again)))
    prots_reliable <- prots_nonoutlier[!(prots_nonoutlier %in% prots_all_unreliable)]
    if (length(prots_reliable)>0){
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      plot(x, Data[prots_reliable[1], ], type="b",xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35), 
           main=paste(c,"after remove unreliable", sep = ", "))
      axis(1, at=x,labels=x)
      c3 = 1
      if (length(prots_reliable)>1){
        for (i in 2:length(prots_reliable)){
          c3 <- c3 + 1
          lines(x, Data[prots_reliable[i], ], type="b", col=c3)
        }
      }
    }
    if(length(prots_reliable) > 0){
      data_c <- NULL
      for (p in prots_reliable){
        id_info <- prots_ID[p, ]
        profile_info_intero <- data_introp[, p]
        profile_info_raw <- protfrac_exp_data_T[, p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID), 
                            paste(rownames(protfrac_exp_data_T), "initial", sep="_"), 
                            paste(rownames(data_introp), "interoplated", sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmunreliable <- rbind(tables_clusters_after_rmunreliable, data_c)
      
    }
    prots_reliable_aboundance <- c()#prots_high_abundance <- c()
    for (p in prots_reliable) {
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
      prots_reliable_aboundance <- c(prots_reliable_aboundance, p_ab)
    }
    if (length(prots_reliable_aboundance) > 0) {
      h <- hist(prots_reliable_aboundance,breaks=seq(from=as.integer(range(prots_reliable_aboundance)[1]),
                                                     to=as.integer(range(prots_reliable_aboundance)[2])+1, by = 1),plot=FALSE)
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      hist(prots_reliable_aboundance, main="Histogram of the average abundances of reliable proteins ",
           xlab="average abundance",
           xlim=c(15, 35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(from=as.integer(range(prots_reliable_aboundance)[1]),
                      to=as.integer(range(prots_reliable_aboundance)[2])+1, by = 1))
    }
    
    ######
    prots_reliable_high_abundance <- prots_reliable[prots_reliable %in% prots_high_abundance]
    if (length(prots_reliable_high_abundance)>0){
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      plot(x, Data[prots_reliable_high_abundance[1], ], type="b",xlab="DPA", ylab="Log Intensity", col=1, ylim=range(0,35), 
           main=paste(c,"after remove low abundance", sep = ", "))
      axis(1, at=x,labels=x)
      c3 = 1
      if (length(prots_reliable_high_abundance)>1){
        for (i in 2:length(prots_reliable_high_abundance)){
          c3 <- c3 + 1
          lines(x, Data[prots_reliable_high_abundance[i], ], type="b", col=c3)
        }
      }
    }
    if(length(prots_reliable_high_abundance) > 0){
      data_c <- NULL
      for (p in prots_reliable_high_abundance){
        id_info <- prots_ID[p, ]
        profile_info_intero <- data_introp[, p]
        profile_info_raw <- protfrac_exp_data_T[, p]
        df_p <- data.frame(c(c, id_info, profile_info_raw, profile_info_intero))
        colnames(df_p) <- c("cluster_id", colnames(prots_ID), 
                            paste(rownames(protfrac_exp_data_T), "initial", sep="_"), 
                            paste(rownames(data_introp), "interoplated", sep="_"))
        rownames(df_p) <- p
        data_c <- rbind(data_c, df_p)
      }
      tables_clusters_after_rmlowabundance <- rbind(tables_clusters_after_rmlowabundance, data_c)
      
    }
    prots_high_aboundance <- c()#prots_high_abundance <- c()
    for (p in prots_reliable_high_abundance) {
      p_ab <- mean(protfrac_exp_data_T[[p]][protfrac_exp_data_T[[p]] > 0])
      prots_high_aboundance <- c(prots_high_aboundance, p_ab)
    }
    if (length(prots_high_aboundance) > 0) {
      h <- hist(prots_high_aboundance,breaks=seq(from=as.integer(range(prots_high_aboundance)[1]),
                                                 to=as.integer(range(prots_high_aboundance)[2])+1, by = 1),plot=FALSE)
      par(mar=c(5, 4, 4, 8), xpd=TRUE)
      hist(prots_high_aboundance, main="Histogram of the average abundances of high abundance proteins ",
           xlab="average abundance",
           xlim=c(15, 35),
           ylim=c(0, max(h_ylim$counts)+5),
           breaks=seq(from=as.integer(range(prots_high_aboundance)[1]),
                      to=as.integer(range(prots_high_aboundance)[2])+1, by = 1))
    }
  }
  dev.off()
  score_frac <- prot_group_score_data[[frac]]
  informative_id <- filter(score_frac, score=="I")$CID
  tables_clusters_informative <- NULL
  for (id in informative_id){
    dt <- filter(tables_clusters_after_rmunreliable, cluster_id==id)
    dt <- mutate(dt, cluster_id=paste(frac, cluster_id, sep="_"))
    tables_clusters_informative <- rbind(tables_clusters_informative, dt)
  }
  
  SOM_clustering_result[[frac]] <- list()
  SOM_codes_result[[frac]] <- list()
  
  write.csv(tables_clusters_SOM, 
            paste(paste(paste0(frac, " Table 1 step1 SOM clustering result"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_clustering_result[[frac]][["step1"]] <- tables_clusters_SOM
  write.csv(tables_clusters_after_rmoutlier, 
            paste(paste(paste0(frac, " Table 2 step2 SOM clustering result (outlier removed)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_clustering_result[[frac]][["step2"]] <- tables_clusters_after_rmoutlier
  write.csv(tables_clusters_after_rmunreliable, 
            paste(paste(paste0(frac, " Table 3 step3 SOM clustering result (unreliable profile removed)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_clustering_result[[frac]][["step3"]] <- tables_clusters_after_rmunreliable
  write.csv(tables_clusters_after_rmlowabundance, 
            paste(paste(paste0(frac, " Table 4 step3 SOM clustering result (kow abundance profile removed)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_clustering_result[[frac]][["step4"]] <- tables_clusters_after_rmlowabundance
  write.csv(tables_clusters_informative, 
            paste(paste(paste0(frac, " Table 5 step4 SOM clustering result (informative groups)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_clustering_result[[frac]][["step5"]] <- tables_clusters_informative
  tables_clusters_informative_high_abundance <- tables_clusters_informative[prots_high_abundance, ]
  tables_clusters_informative_high_abundance <- filter(tables_clusters_informative_high_abundance, !is.na(cluster_id))
  indexOfcluster <- 7
  if (frac == "APO") {
    indexOfcluster <- 6
  }
  tables_clusters_informative_high_abundance <- arrange(tables_clusters_informative_high_abundance, as.numeric(substring(cluster_id, indexOfcluster,indexOfcluster+1)))
  numOfCol <- dim(tables_clusters_informative_high_abundance)[2]
  tables_clusters_informative_temp <- tables_clusters_informative_high_abundance[,c(1, (numOfCol-19):numOfCol)]
  for (i in 1:dim(tables_clusters_informative_temp)) {
    temp_i <- tables_clusters_informative_temp[i, 2:21]
    for (j in 2:21) {
      tables_clusters_informative_temp[i,j] <- tables_clusters_informative_temp[i,j] / max(temp_i)
    }
  }
  tables_clusters_informative_high_abundance_normalized <- tables_clusters_informative_temp
  #tables_clusters_informative_high_abundance_normalized_temp <- tables_clusters_informative_high_abundance
  tables_mean_each_cluster <- NULL
  for (c in levels(factor(tables_clusters_informative_high_abundance_normalized$cluster_id, levels = unique(tables_clusters_informative_high_abundance_normalized$cluster_id)))) {
    temp_c <- filter(tables_clusters_informative_high_abundance_normalized, cluster_id==c)[,-1]
    temp_c_df <- data.frame(colMeans(temp_c))
    names(temp_c_df) <- c
    tempdf <- data.frame(t(temp_c_df))
    tables_mean_each_cluster <- rbind(tables_mean_each_cluster, tempdf)
  }
  mean_of_normalized_informative[[frac]] <- tables_mean_each_cluster
  write.csv(tables_clusters_informative_high_abundance, 
            paste(paste(paste0(frac, " Table informative groups containing high abundance"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  write.csv(tables_clusters_informative_high_abundance_normalized, 
            paste(paste(paste0(frac, " Table informative groups containing high abundance (normalized)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  write.csv(tables_mean_each_cluster, 
            paste(paste(paste0(frac, " Table mean of normalized informative groups"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  code_som <- data.frame(codes_som_all)
  informative_code <- code_som[informative_id,]
  rownames(informative_code) <- paste(frac, informative_id, sep="_")
  write.csv(code_som, 
            paste(paste(paste0(frac, " Table 5 SOM codes profile"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_codes_result[[frac]][["codes"]] <- code_som
  write.csv(informative_code, 
            paste(paste(paste0(frac, " Table 6 SOM codes profile (informative group)"), Sys.Date(), sep='_'), "csv", sep='.'), row.names = T)
  SOM_codes_result[[frac]][["informative_codes"]] <- informative_code
}
