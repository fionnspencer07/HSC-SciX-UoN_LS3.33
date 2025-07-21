library(here)
#init packages
source(here("Custom R Functions and Scripts", "init_packages.R"))  # will only work if the rproj opened is HSC SciX UoN LS3.33

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs",
                     "plasma_cells", "progenitor_cells", "NKT_cells")
#importing RDS & Init as a Seurat Object
seu <- readRDS(here("Data",
                    "seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds"))

#Subsetting Seurat Object

## NK, T Cells
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))

## NK, T Cells - No Treatment [Excluding these Conditions (ONC201), (Dexmethasone), (ONC201 & Dexmethasone)]
seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "NAIVE", "SHAM"))

saveRDS(seu_NKT_focused, file = here("Data", "seu_NKT_focused.RDS")) 

# Define Genes of Interest
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")


# Create Output Folder ====
dir.create(here("R Projects", "1807DGE"))
dir.create(here("R Projects", "1807DGE", "output"))
dir.create(here("R Projects", "1807DGE", "output", "modifiedCS-outs"))

?presto

pairwise_comparisons <- list(
  
  "Naive_vs_UT_WB" = c("Naive_WB", "UT_WB"),
  "Naive_vs_UT_BM" = c("Naive_BM", "UT_BM")
  
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- wilcoxauc(seu_NKT_focused, group_by = "sample_name", groups_use = pairwise_comparisons[[i]])
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(group ==  pairwise_comparisons[[i]][[1]]) %>%
    arrange(group, -logFC)
  
  writeData(wb, sheet = i, x = deg_T)
  
  
}

saveWorkbook(wb, here("R Projects", "1807DGE", "output", "modifiedCS-outs", "wilcoxon.xlsx"))
, overwrite = TRUE)