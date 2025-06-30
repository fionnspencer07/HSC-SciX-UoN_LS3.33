source("initialise_packages.R")

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")

DimPlot(seu, group.by = "cell_types", reduction = "pca", label = FALSE, repel = FALSE)

DimPlot(seu, group.by = "sample_name", reduction = "umap_SCT", label = FALSE, repel = FALSE)

Nebulosa::plot_density(seu, "Cd3e") + 
  facet_wrap(.~seu$sample_name, ncol = 6) + theme_void()

FeaturePlot(seu, features = "Cd3e", reduction = "umap_harmony_SCT", split.by = "sample_name") +
  patchwork::plot_layout(ncol = 6, nrow = 2)

Nebulosa::plot_density(seu, "Cd3g") + 
  facet_wrap(.~seu$sample_name, ncol = 6) + theme_void()

FeaturePlot(seu, features = "Cd3g", reduction = "umap_harmony_SCT", split.by = "sample_name") +
  patchwork::plot_layout(ncol = 6, nrow = 2)

