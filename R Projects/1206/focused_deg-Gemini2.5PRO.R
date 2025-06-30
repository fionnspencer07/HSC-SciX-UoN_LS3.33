source("initialise_packages.R")

# Custom colors (can keep this or remove if not used for plotting later)
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")

# --- MODIFICATION 1: Define your focused conditions ---
# Only include the conditions relevant to your core question
conditions <- c("Naive", "Sham", "UT")
sample_types <- c("WB", "BM")

# Build the sample names for the focused analysis
focused_samples <- as.vector(outer(conditions, sample_types, paste, sep = "_"))

# --- MODIFICATION 2: Subset the Seurat object ---
# This creates a smaller object for more efficient analysis
seu_focused <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells") & sample_name %in% focused_samples)

# Check that the groups in your new subset object are as expected
print("Groups present in the focused data:")
print(unique(seu_focused$sample_name))

# Generate all unique pairwise combos of the focused samples
combos <- t(combn(focused_samples, 2))

deg_results <- list()

for (i in seq_len(nrow(combos))) {
  pair <- combos[i, ]
  comp_name <- paste(pair, collapse = "_vs_")
  
  message("Running comparison: ", comp_name)
  
  # Run the test on the smaller 'seu_focused' object
  res <- wilcoxauc(seu_focused, group_by = "sample_name", groups_use = pair) %>%
    filter(group == pair[1]) %>%
    arrange(group, -logFC)
  
  deg_results[[comp_name]] <- res
}

# You are now ready to analyze the results!
message("Focused DEG analysis complete.")
