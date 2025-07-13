# Load libraries
library(tidyverse)
library(here)
library(Seurat)
library(hdf5r)
library(cowplot)
library(glue)
library(ggExtra)
library(DT)
library(kableExtra)
library(glmGamPoi)
library(harmony)
library(presto)
library(openxlsx)
library(SingleR)
library(celldex)
library(scrapper)
library(ggplotify)
library(viridis)
library(gridExtra)
library(ggrepel)
library(dittoSeq)
library(ggpubr)
library(httpgd)
library(Nebulosa)
library(scClassify)
library(pheatmap)
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(tidyr)
library(ggplot2)
library(RColorBrewer)
library(DESeq2)
library(edgeR)
library(edge)
library(DT)
library(patchwork)
library(muscat)
library(scater)
library(EnhancedVolcano)
library(Rvisdiff)

# Additional packages for trajectory analysis
suppressWarnings({
  if (!require("monocle3", quietly = TRUE)) {
    cat("monocle3 not found. Install with: BiocManager::install('monocle3')\n")
  } else {
    library(monocle3)
  }
  
  if (!require("slingshot", quietly = TRUE)) {
    cat("slingshot not found. Install with: BiocManager::install('slingshot')\n")
  } else {
    library(slingshot)
  }
  
  if (!require("tradeSeq", quietly = TRUE)) {
    cat("tradeSeq not found. Install with: BiocManager::install('tradeSeq')\n")
  } else {
    library(tradeSeq)
  }
  
  if (!require("SingleCellExperiment", quietly = TRUE)) {
    cat("SingleCellExperiment not found. Install with: BiocManager::install('SingleCellExperiment')\n")
  } else {
    library(SingleCellExperiment)
  }
  
  if (!require("ggridges", quietly = TRUE)) {
    cat("ggridges not found. Install with: install.packages('ggridges')\n")
  } else {
    library(ggridges)
  }
  
  if (!require("future", quietly = TRUE)) {
    cat("future not found. Install with: install.packages('future')\n")
  } else {
    library(future)
  }
  
  if (!require("future.apply", quietly = TRUE)) {
    cat("future.apply not found. Install with: install.packages('future.apply')\n")
  } else {
    library(future.apply)
  }
  
  if (!require("progressr", quietly = TRUE)) {
    cat("progressr not found. Install with: install.packages('progressr')\n")
  } else {
    library(progressr)
  }
})

# Set up parallel processing
suppressWarnings({
  if (require("future", quietly = TRUE)) {
    plan(multisession, workers = min(4, parallel::detectCores() - 1))
  }
})

cat("Package initialization complete!\n")
cat("Note: If any trajectory analysis packages are missing, run the installation script first.\n")
