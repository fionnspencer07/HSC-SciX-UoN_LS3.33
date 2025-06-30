source("initialise_packages.R")

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")

seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))

conditions <- c("Naive", "Sham", "UT", "ONC", "DEX", "DEX_ONC")
sample_types <- c("WB", "BM")

# Build all present sample names
all_samples <- as.vector(outer(conditions, sample_types, paste, sep = "_"))

# Confirm these match exactly with groups in your data
groups_in_data <- unique(seu_NKT$sample_name)

if (!all(all_samples %in% groups_in_data)) {
  warning("Some expected groups are missing from the data:")
  print(setdiff(all_samples, groups_in_data))
}

# Generate all unique pairwise combos of present samples
combos <- t(combn(all_samples, 2))

deg_results <- list()

for (i in seq_len(nrow(combos))) {
  pair <- combos[i, ]
  comp_name <- paste(pair, collapse = "_vs_")
  
  message("Running comparison: ", comp_name)
  
  res <- wilcoxauc(seu_NKT, group_by = "sample_name", groups_use = pair) %>%
    filter(group == pair[1]) %>%
    arrange(group, -logFC)
  
  deg_results[[comp_name]] <- res
}  

