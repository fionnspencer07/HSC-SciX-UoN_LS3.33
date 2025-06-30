install.packages("tidyverse", dependencies = TRUE)
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

# GitHub package
install.packages("devtools", dependencies = TRUE)
devtools::install_github("immunogenomics/presto", dependencies = TRUE)
remotes::install_github("nx10/httpgd", dependencies = TRUE)
remotes::install_github("immunogenomics/presto", dependencies = TRUE)


# Seurat (CRAN)
install.packages("Seurat", dependencies = TRUE)


renv::install.packages("presto")
