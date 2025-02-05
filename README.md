# Proteomics-based-models-of-cotton-fiber-protein-dynamics-and-gene-expression-control

### SOM Clustering Analysis for the Cotton Project

This repository contains three R scripts for cotton proteomics data analysis:

1. **`SOM_clustering_analysis.R`**  
   - Performs Self‐Organizing Map (SOM) clustering of protein expression data.
   - Handles outlier removal and interpolates missing profiles.
   - Generates CSV and PDF outputs summarizing clustering results.
   - Merges the informative data from the SOM clustering steps.
   - Performs additional cross-fraction or cross-condition analyses.
   - Produces combined or summarized CSV files for further interpretation.

2. **`ID_assignment.R`**  
   - Assigns protein IDs and merges them with relevant mRNA transcript identifiers.
   - Helps align homoeolog IDs, resolving multiple isoforms.
   - Outputs consistent ID mappings (CSV) used in the clustering pipeline.

<br />

## Requirements

- **R version 4.2 or higher**  
- The following packages:
  - `kohonen` (for SOM functionality)
    
    1.once installation of ‘kohonen’ completed, find path in your computer:
     path...\R\win-library\R version...\kohonen\Distances
    
    2. copy and paste the file “wcc3.cpp” into “Distances” folder found in step 1.
  - `tempR` (time-series / correlation functions)
  - `dplyr`, `tibble`, `tidyverse` (data manipulation & I/O)
  - `ggplot2` (plotting)
  - `openxlsx` (reading/writing XLSX)
  - `Rcpp` (interfaces to C++ for certain distance measures)
