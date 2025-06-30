source("install_packages.R")
source("initialise_packages.R")

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")

seu@meta.data$cell_types <- as.character(seu@meta.data$cell_types)
seu@meta.data$cell_types[seu@meta.data$sub31 == "31_1"] <- "NKT_cells"

DimPlot(seu, group.by = "cell_types", reduction = "umap_harmony_SCT", label = TRUE, repel = TRUE) +
  scale_colour_manual(values = man_cols)

sce <- as.SingleCellExperiment(seu, assay = "RNA", slot = "data")

p <- dittoBarPlot(sce, var = "cell_types_new", group.by = "sample_name")+
  theme_classic()+
  theme(axis.title.x = element_blank())+
  scale_x_discrete(limits = levels(seu@meta.data[["sample_name"]]))+
  rotate_x_text(angle = 90)+
  scale_fill_manual(values = man_cols)


