# Initialise Packages
source("C:/Users/fionn/Newcastle Grammar School/NGS SciX-2025 UoN DMG Cancer Signalling - Bioinformatics and Proteomics/Fionn Git/HSC-SciX-UoN_LS3.33/Custom R Functions and Scripts/init_packages.R")

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs",
                     "plasma_cells", "progenitor_cells", "NKT_cells")


seu <- readRDS("./Data/seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))
seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "Naive", "Sham"))
