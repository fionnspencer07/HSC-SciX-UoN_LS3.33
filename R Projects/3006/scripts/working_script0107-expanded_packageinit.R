# Initialise Packages
## Install Packages install.packages("tidyverse", dependencies = TRUE)
# CRAN packages
install.packages("here", dependencies = TRUE)
install.packages("hdf5r", dependencies = TRUE)
install.packages("cowplot", dependencies = TRUE)
install.packages("glue", dependencies = TRUE)
install.packages("ggExtra", dependencies = TRUE)
install.packages("DT", dependencies = TRUE)
install.packages("kableExtra", dependencies = TRUE)
install.packages("harmony", dependencies = TRUE)
install.packages("openxlsx", dependencies = TRUE)
install.packages("ggplotify", dependencies = TRUE)
install.packages("viridis", dependencies = TRUE)
install.packages("gridExtra", dependencies = TRUE)
install.packages("ggrepel", dependencies = TRUE)
install.packages("ggpubr", dependencies = TRUE)

# Bioconductor setup
install.packages("BiocManager", dependencies = TRUE)
BiocManager::install("glmGamPoi", dependencies = TRUE)
BiocManager::install("SingleR", dependencies = TRUE)
BiocManager::install("celldex", dependencies = TRUE)
BiocManager::install("scrapper", dependencies = TRUE)
BiocManager::install("dittoSeq", dependencies = TRUE)
BiocManager::install("NeuCA", dependencies = TRUE)
BiocManager::install("ComplexHeatmap", dependencies = TRUE)
BiocManager::install("edge", dependencies = TRUE)
BiocManager::install("edgeR", dependencies = TRUE)
BiocManager::install("DESeq2", dependencies = TRUE)

# GitHub package
install.packages("devtools", dependencies = TRUE)
devtools::install_github("immunogenomics/presto", dependencies = TRUE)
remotes::install_github("nx10/httpgd", dependencies = TRUE)
remotes::install_github("immunogenomics/presto", dependencies = TRUE)


# Seurat (CRAN)
install.packages("Seurat", dependencies = TRUE)


renv::install.packages("presto")

## Init Packages 

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



# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs",
                     "plasma_cells", "progenitor_cells", "NKT_cells")


seu <- readRDS("./Data/seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))
seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "Naive", "Sham"))
