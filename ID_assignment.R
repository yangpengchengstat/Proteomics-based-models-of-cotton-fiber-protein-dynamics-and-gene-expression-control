# -----------------------------------------------------------------------------------------
# Protein-mRNA ID Assignment
# 
# This script assigns protein IDs from the input table (Prot_ID_map) by applying certain 
# filtering criteria (e.g., removing IDs below 1/4 maximum peptide count), determining 
# the "majority" or best isoform for ortholog/homoeolog pairs, and merging them with 
# relevant fraction-based clustering outputs from SOM_clustering_result. 
#
# Author: Pengcheng Yang (2023), for GitHub usage.
# -----------------------------------------------------------------------------------------

############################
# Protein-mRNA ID assignment
############################

workingDir = "C:/Users/User/Desktop/My Files/cotton pj/Prot mRNA ID assignment"
setwd(workingDir)

# 'Prot_ID_map' typically references the 'proteinIDmap.xlsx' already loaded earlier
# data: protein-mRNA ID matching
Prot_ID_map <- read.xlsx("proteinIDmap.xlsx")

# Define column names for references
prot_id <- "Protein.IDs"
count   <- "Peptide.counts.(unique)"

prot_set <- strsplit(Prot_ID_map[1, prot_id], split=";")[[1]]
Majority_ID <- c()  # Will store the final chosen isoform(s)
Reason      <- c()  # Explanations for each assignment
prots_map   <- list()

# Iterate through each row in 'Prot_ID_map' to assign best isoform(s)
for (i in 1:(dim(Prot_ID_map)[1])){
  protMap <- c()
  
  # Retrieve protein set and peptide counts from the row
  prot_set  <- strsplit(Prot_ID_map[i, prot_id], split=";")[[1]]
  count_set <- as.numeric(strsplit(Prot_ID_map[i, count], split=";")[[1]])
  
  # Basic check that the number of IDs matches number of counts
  if (length(prot_set) != length(count_set)){
    print("problem")
  }
  
  # Remove IDs with less than 1/4 the max count (a threshold approach)
  count_set_frac <- count_set / max(count_set)
  prot_set  <- prot_set [count_set_frac > 0.25]
  count_set <- count_set[count_set_frac > 0.25]
  
  # We'll parse the protein ID to see ortholog/homolog
  ortholog_check <- c()
  homolog_check  <- c()
  for (j in 1:length(prot_set)){
    # Typically, Ghir_A... or Ghir_D... pattern
    # strsplit(...)[[1]][ c(1,2) ] might be "Ghir" + "someGeneID"
    # index [4] might be "A" or "D" etc.
    ortholog_check <- c(ortholog_check,
                        paste(strsplit(prot_set[j], split='.', fixed=TRUE)[[1]][c(1,2)], 
                              collapse='.'))
    homolog_check <- c(homolog_check,
                       strsplit(prot_set[j], split='.', fixed=TRUE)[[1]][4])
  }
  
  # CASE 1: Single ortholog & single homoeolog → pick the isoform w/ max count
  if (length(levels(factor(ortholog_check))) == 1 & length(levels(factor(homolog_check))) == 1){
    prot <- prot_set[which.max(count_set)]
    why  <- "case 1: unique ortholog and unique homoeolog, choose the isoform with max count."
    protMap <- c(protMap, prot)
    Majority_ID <- c(Majority_ID, prot)
    Reason <- c(Reason, why)
  }
  
  # CASE 2: Single ortholog & multiple homoeologs (A vs. D for example)
  # choose isoform w/ max count in each homoeolog set
  if (length(levels(factor(ortholog_check))) == 1 &
      length(levels(factor(homolog_check))) > 1){
    prot_homo_A    <- c()
    prot_homo_B    <- c()
    count_homo_A   <- c()
    count_homo_B   <- c()
    
    for (k in 1:length(prot_set)){
      p <- prot_set[k]
      if (strsplit(prot_set[k], split='.', fixed=TRUE)[[1]][4] == "A"){
        prot_homo_A  <- c(prot_homo_A, p)
        count_homo_A <- c(count_homo_A, count_set[k])
      } else {
        prot_homo_B  <- c(prot_homo_B, p)
        count_homo_B <- c(count_homo_B, count_set[k])
      }
    }
    homo_A_D <- c(prot_homo_A[which.max(count_homo_A)], 
                  prot_homo_B[which.max(count_homo_B)])
    Majority_ID <- c(Majority_ID, paste(homo_A_D, collapse=';'))
    protMap     <- c(protMap, homo_A_D)
    
    why <- "case 2: unique orthologs and multiple homeologs, choose the homeolog A and D with isoform with max count respectively."
    Reason <- c(Reason, why)
  }
  
  # CASE 3: Multiple orthologs & possibly multiple homoeologs
  if (length(levels(factor(ortholog_check))) > 1){
    orth_set      <- levels(factor(ortholog_check))
    orth_sub_prot <- list()
    orth_sub_count<- list()
    
    # Build sub-lists for each orth_set
    for (orth in orth_set){
      orth_sub_prot[[orth]]  <- c()
      orth_sub_count[[orth]] <- c()
      
      for (m in 1:length(prot_set)){
        p_orth <- paste(strsplit(prot_set[m], split='.', fixed=TRUE)[[1]][c(1,2)], 
                        collapse='.')
        if (p_orth == orth){
          orth_sub_prot[[orth]]  <- c(orth_sub_prot[[orth]],  prot_set[m])
          orth_sub_count[[orth]] <- c(orth_sub_count[[orth]], count_set[m])
        }
      }
    }
    orht_count <- c()
    for (orth in orth_set){
      count_o <- count_set[ prot_set %in% orth_sub_prot[[orth]] ]
      orht_count <- c(orht_count, sum(count_o))
    }
    orth_id_set <- orth_set
    
    if (length(orth_id_set) > 1){
      multi_orth_p <- c()
      
      for (orth in orth_id_set){
        orth_homo_check <- c()
        
        for (p in orth_sub_prot[[orth]]){
          orth_homo_check <- c(orth_homo_check,
                               strsplit(p, split='.', fixed=TRUE)[[1]][4])
        }
        if (length(orth_homo_check) == 1){
          # if there's only one homolog, choose isoform with max
          count_h <- c()
          for (p in orth_sub_prot[[orth]]){
            count_h <- c(count_h, count_set[ prot_set==p ])
          }
          multi_orth_p <- c(multi_orth_p, orth_sub_prot[[orth]][ which.max(count_h) ])
        }
        if (length(orth_homo_check) > 1){
          # multiple homoeolog forms
          prot_homo_A   <- c()
          prot_homo_B   <- c()
          count_homo_A  <- c()
          count_homo_B  <- c()
          
          for (p in orth_sub_prot[[orth]]){
            if (strsplit(p, split='.', fixed=TRUE)[[1]][4] == "A"){
              prot_homo_A  <- c(prot_homo_A, p)
              count_homo_A <- c(count_homo_A, count_set[ prot_set==p ])
            } else {
              prot_homo_B  <- c(prot_homo_B, p)
              count_homo_B <- c(count_homo_B, count_set[ prot_set==p ])
            }
          }
          homo_A_D <- c(prot_homo_A[ which.max(count_homo_A) ], 
                        prot_homo_B[ which.max(count_homo_B) ])
          multi_orth_p <- c(multi_orth_p, paste(homo_A_D, collapse=';'))
          protMap_p    <- c(multi_orth_p, homo_A_D)
        }
      }
      Majority_ID <- c(Majority_ID, paste(multi_orth_p, collapse=';'))
      protMap     <- c(protMap, protMap_p)
      why <- "case 3: multiple orthologs and multiple homeologs, choose the homolog A and D with isoform with max count respectively."
      Reason <- c(Reason, why)
    }
  }
  prots_map[[ Prot_ID_map[i, "Leading.protein.IDs"] ]] <- protMap
}

# data_f aggregates relevant columns from 'Prot_ID_map'
data_f <- Prot_ID_map[, c(2,3,4,8)]
data_f[["Prots_assignment"]] <- Majority_ID
data_f[["Explanation"]]      <- Reason

# Define textual patterns for quick referencing
c1 <- "case 1: unique ortholog and unique homoeolog, choose the isoform with max count."
c2 <- "case 2: unique orthologs and multiple homeologs, choose the homeolog A and D with isoform with max count respectively."
c3 <- "case 3: multiple orthologs and multiple homeologs, choose the homolog A and D with isoform with max count respectively."

# Subset the data based on these textual explanations
case1 <- filter(data_f, Explanation==c1)
case2 <- filter(data_f, Explanation==c2)
case3 <- filter(data_f, Explanation==c3)

# Output the assigned ID table
write.csv(
  data_f, 
  paste(paste("Table - Cotton protein ID assignment", Sys.Date(), sep='_'), "csv", sep='.'),
  row.names=FALSE
)

# Additional merging steps with the P200/S200/APO further clustering outputs (step4)
P200_further_clustering <- SOM_clustering_result$P200$step4
P200_further_clustering <- P200_further_clustering[, 30:49]

S200_further_clustering <- SOM_clustering_result$S200$step4
S200_further_clustering <- S200_further_clustering[, 30:49]

APO_further_clustering  <- SOM_clustering_result$APO$step4
APO_further_clustering  <- APO_further_clustering[, 29:48]

# Combine leading IDs with the 'case1' scenario
leading_IDs       <- rownames(P200_further_clustering)
case1_IDs         <- case1$Leading.protein.IDs
leading_IDs_in_case1 <- leading_IDs[ leading_IDs %in% case1_IDs ]

# Retrieve the assigned mRNA IDs from the 'case1' scenario
# (though it's commented out in the original code, we keep it for reference)
#rmna_IDs_in_case1 <- case1$Prots_assignment

rmna_IDs_in_case1 <- (case1$Prots_assignment)[ case1_IDs %in% leading_IDs_in_case1 ]
leading_IDs_in_case1[1:5]
rmna_IDs_in_case1[1:5]

# Convert protein IDs to an example of a new mRNA ID (e.g., set the 3rd index to "1")
mrna_for_leading_IDs <- c()
for (i in leading_IDs_in_case1){
  i_ <- strsplit(i, split=".", fixed=TRUE)[[1]]
  i_[3] <- "1"
  mrna_for_leading_IDs <- c(mrna_for_leading_IDs, paste(i_, collapse="."))
}

# End of the ID assignment script
# ---------------------------------------------------------------------------------------
