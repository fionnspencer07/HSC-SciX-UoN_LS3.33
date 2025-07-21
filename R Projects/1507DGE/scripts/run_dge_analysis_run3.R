# ============================================================================ #
# Comprehensive DGE and Visualization Pipeline for NK/T Cell Analysis (Run 2 v8)
#
# Author: GitHub Copilot, based on specifications from fionnspencer07
# Date: 2025-07-16
#
# Description:
# This script performs a multi-method DGE analysis using Wilcoxon, DESeq2 (on
# pseudobulk data), and MAST for ALL pairwise comparisons between treatment-tissue
# combinations. It generates extensive visualizations using ComplexHeatmap and 
# detailed, multi-sheet Excel reports for each method.
#
# v8 Update: Extended to perform all 15 pairwise comparisons between treatment-tissue
# combinations as requested.
# ============================================================================ #

# --- 1. INITIAL SETUP & DATA LOADING ---

# Set seed for reproducibility across all analyses
set.seed(42)

# Load package libraries from the user-defined init script
library(here)
source(here("Custom R Functions and Scripts", "init_packages.R"))
library(DESeq2) # Explicitly load for DESeq2 analysis
library(SummarizedExperiment) # Required for pseudobulking
library(scuttle) # Explicitly load for modern pseudobulking

# Define output directory for Run 2 using here()
output_dir <- here("R Projects", "1507DGE", "output", "Run 3")
# Create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the focused Seurat object
message("Loading focused Seurat object...")
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))
message("Object loaded successfully.")

# Ensure 'treatment' is a factor with a logical order for comparisons
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, levels = c("NAIVE", "SHAM", "UT"))

# Create a combined treatment_tissue identifier
seu_NKT_focused$treatment_tissue <- paste(seu_NKT_focused$treatment, seu_NKT_focused$tissue, sep = "_")

# Define all treatment-tissue combinations
treatment_tissue_combinations <- c(
  "NAIVE_WB", "NAIVE_BM", 
  "SHAM_WB", "SHAM_BM", 
  "UT_WB", "UT_BM"
)

# Define all 15 pairwise comparisons
pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM" = c("NAIVE_WB", "NAIVE_BM"),
  "NAIVE_WB_vs_SHAM_WB" = c("NAIVE_WB", "SHAM_WB"),
  "NAIVE_WB_vs_SHAM_BM" = c("NAIVE_WB", "SHAM_BM"),
  "NAIVE_WB_vs_UT_WB" = c("NAIVE_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM" = c("NAIVE_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB" = c("NAIVE_BM", "SHAM_WB"),
  "NAIVE_BM_vs_SHAM_BM" = c("NAIVE_BM", "SHAM_BM"),
  "NAIVE_BM_vs_UT_WB" = c("NAIVE_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM" = c("NAIVE_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM" = c("SHAM_WB", "SHAM_BM"),
  "SHAM_WB_vs_UT_WB" = c("SHAM_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM" = c("SHAM_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB" = c("SHAM_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM" = c("SHAM_BM", "UT_BM"),
  "UT_WB_vs_UT_BM" = c("UT_WB", "UT_BM")
)

# Define Genes of Interest
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")
all_genes_of_interest <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

# --- 2. HELPER FUNCTION FOR EXCEL EXPORT ---
create_dge_excel <- function(dge_results_list, method_name, output_dir) {
  wb <- createWorkbook()
  
  # Create a combined dataframe with all contrasts
  all_contrasts_df <- bind_rows(lapply(names(dge_results_list), function(contrast_name) {
    dge_results_list[[contrast_name]] %>% 
      mutate(contrast = contrast_name) %>%
      dplyr::select(contrast, everything())
  }))
  
  # Add individual contrast sheets (first 10 to avoid Excel sheet limit issues)
  contrast_names <- names(dge_results_list)
  for (i in seq_along(contrast_names)) {
    if (i <= 10) {  # Excel has sheet name limitations
      sheet_name <- paste0("Contrast_", i)
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, dge_results_list[[contrast_names[i]]])
    }
  }
  
  # Add gene-specific sheets
  addWorksheet(wb, "S1PR_Genes_All_Contrasts")
  writeData(wb, "S1PR_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% s1pr_genes))
  
  addWorksheet(wb, "Trafficking_Genes_All")
  writeData(wb, "Trafficking_Genes_All", all_contrasts_df %>% filter(gene %in% trafficking_genes))
  
  addWorksheet(wb, "Retention_Genes_All")
  writeData(wb, "Retention_Genes_All", all_contrasts_df %>% filter(gene %in% retention_genes))
  
  addWorksheet(wb, "Egress_Genes_All")
  writeData(wb, "Egress_Genes_All", all_contrasts_df %>% filter(gene %in% egress_genes))
  
  # Summary statistics
  summary_stats <- all_contrasts_df %>% 
    group_by(contrast) %>% 
    summarise(
      total_genes_tested = n(),
      significant_upregulated = sum(p_val_adj < 0.05 & avg_log2FC > 0.25, na.rm = TRUE),
      significant_downregulated = sum(p_val_adj < 0.05 & avg_log2FC < -0.25, na.rm = TRUE),
      .groups = 'drop'
    )
  
  addWorksheet(wb, "Summary_Stats")
  writeData(wb, "Summary_Stats", summary_stats)
  
  # Top genes across all contrasts
  top_up <- all_contrasts_df %>% 
    filter(p_val_adj < 0.05) %>% 
    group_by(contrast) %>% 
    arrange(desc(avg_log2FC), .by_group = TRUE) %>% 
    slice_head(n = 50)  # Reduced to 50 to manage file size
  
  top_down <- all_contrasts_df %>% 
    filter(p_val_adj < 0.05) %>% 
    group_by(contrast) %>% 
    arrange(avg_log2FC, .by_group = TRUE) %>% 
    slice_head(n = 50)
  
  addWorksheet(wb, "Top50_Upregulated")
  writeData(wb, "Top50_Upregulated", top_up)
  
  addWorksheet(wb, "Top50_Downregulated")
  writeData(wb, "Top50_Downregulated", top_down)
  
  # Save the workbook
  file_path <- file.path(output_dir, glue("DGE_{method_name}_All_Pairwise_Comparisons.xlsx"))
  saveWorkbook(wb, file_path, overwrite = TRUE)
  message(glue("SUCCESS: Saved Excel report to {file_path}"))
}

# --- 3. DGE ANALYSIS ---

# METHOD 1: Wilcoxon (via Presto)
message("\n--- Starting DGE Analysis: Wilcoxon (presto) for all pairwise comparisons ---")
Idents(seu_NKT_focused) <- "treatment_tissue"
wilcox_results <- list()

for (comparison_name in names(pairwise_comparisons)) {
  message(glue("  Running Wilcoxon for {comparison_name}..."))
  ident1 <- pairwise_comparisons[[comparison_name]][1]
  ident2 <- pairwise_comparisons[[comparison_name]][2]
  
  result <- presto::wilcoxauc(seu_NKT_focused, 
                              group_by = 'treatment_tissue', 
                              seurat_assay = 'RNA', 
                              ident.1 = ident1, 
                              ident.2 = ident2) %>% 
    as_tibble(rownames = "gene") %>% 
    rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)
  
  wilcox_results[[comparison_name]] <- result
}

create_dge_excel(wilcox_results, "Wilcoxon", output_dir)

# METHOD 2: DESeq2 on Pseudobulk Data using muscat::pbDS
message("\n--- Starting DGE Analysis: DESeq2 on Pseudobulk Data using muscat::pbDS ---")

library(SingleCellExperiment)
library(muscat)

# Convert Seurat object to SingleCellExperiment
sce <- as.SingleCellExperiment(seu_NKT_focused)

# Add group and sample metadata for muscat
colData(sce)$group_id <- sce$treatment_tissue   # Experimental groups (treatment_tissue)
colData(sce)$sample_id <- sce$sample_name       # Sample replicates

# Aggregate single-cell counts into pseudobulk samples (grouped by sample_id and group_id)
pb_sce <- muscat::aggregateData(sce, assay = "counts", by = c("sample_id", "group_id"))

# Now run differential expression on pseudobulk data using DESeq2
pb_res <- muscat::pbDS(pb_sce, method = "DESeq2", design = counts_matrix, verbose = TRUE)

# Extract results and reformat for Excel export
deseq2_results <- list()
for (comparison_name in names(pb_res)) {
  res_df <- as.data.frame(pb_res[[comparison_name]]) %>%
    dplyr::rename(
      gene = gene,
      avg_log2FC = logFC,
      p_val = pval,
      p_val_adj = padj
    ) %>%
    dplyr::mutate(contrast = comparison_name) %>%
    dplyr::select(contrast, everything())
  
  deseq2_results[[comparison_name]] <- res_df
}

create_dge_excel(deseq2_results, "DESeq2_Pseudobulk_muscat", output_dir)

# METHOD 3: MAST
message("\n--- Starting DGE Analysis: MAST for all pairwise comparisons ---")
Idents(seu_NKT_focused) <- "treatment_tissue"
mast_results <- list()

for (comparison_name in names(pairwise_comparisons)) {
  message(glue("  Running MAST for {comparison_name}..."))
  ident1 <- pairwise_comparisons[[comparison_name]][1]
  ident2 <- pairwise_comparisons[[comparison_name]][2]
  
  result <- FindMarkers(seu_NKT_focused, 
                        ident.1 = ident1, 
                        ident.2 = ident2, 
                        test.use = "MAST", 
                        logfc.threshold = 0) %>%
    rownames_to_column("gene") %>%
    as_tibble()
  
  mast_results[[comparison_name]] <- result
}

create_dge_excel(mast_results, "MAST", output_dir)

# --- 4. VISUALIZATION ---
message("\n--- Generating Visualizations for Run 2 (All Pairwise Comparisons) ---")

save_svg <- function(plot, filename, width = 10, height = 8) {
  path <- file.path(output_dir, filename)
  if (inherits(plot, "Heatmap")) { 
    svg(path, width = width, height = height)
    draw(plot)
    dev.off() 
  } else { 
    ggsave(path, plot, width = width, height = height, device = "svg") 
  }
  message(glue("  Saved plot: {filename}"))
}

treatment_colors <- c("NAIVE" = "#619CFF", "SHAM" = "#00BA38", "UT" = "#F8766D")
tissue_colors <- c("BM" = "#7CAE00", "WB" = "#C77CFF")
cell_type_colors <- c("T_cells" = "#f3e65d", "NK_cells" = "#86b0cc", "NKT_cells" = "#eeb84c")

# Volcano plot function
plot_volcano <- function(dge_df, title) {
  data_to_plot <- dge_df %>% 
    mutate(
      highlight = case_when(
        gene %in% s1pr_genes ~ "S1P Pathway",
        gene %in% trafficking_genes ~ "Trafficking",
        p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 ~ "Significant",
        TRUE ~ "Not Significant"
      ),
      label = if_else(gene %in% all_genes_of_interest, gene, "")
    )
  
  p <- ggplot(data_to_plot, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highlight, label = label)) +
    geom_point(alpha = 0.7, size = 2) +
    geom_text_repel(max.overlaps = 15, size = 3, bg.color = "white", bg.r = 0.1) +
    scale_color_manual(name = "Gene Group", 
                       values = c("S1P Pathway" = "#E41A1C", 
                                  "Trafficking" = "#377EB8", 
                                  "Significant" = "black", 
                                  "Not Significant" = "grey80")) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  
  return(p)
}

# Generate volcano plots for select comparisons
key_comparisons <- c("NAIVE_WB_vs_SHAM_WB", "NAIVE_WB_vs_UT_WB", "NAIVE_BM_vs_UT_BM")

for (comp in key_comparisons) {
  if (comp %in% names(wilcox_results)) {
    save_svg(plot_volcano(wilcox_results[[comp]], glue("Volcano: Wilcoxon ({comp})")), 
             glue("Volcano_Wilcoxon_{comp}.svg"))
  }
  if (comp %in% names(deseq2_results)) {
    save_svg(plot_volcano(deseq2_results[[comp]], glue("Volcano: DESeq2 ({comp})")), 
             glue("Volcano_DESeq2_{comp}.svg"))
  }
  if (comp %in% names(mast_results)) {
    save_svg(plot_volcano(mast_results[[comp]], glue("Volcano: MAST ({comp})")), 
             glue("Volcano_MAST_{comp}.svg"))
  }
}

# Average expression heatmap by treatment-tissue combinations
avg_exp_grouped <- AverageExpression(seu_NKT_focused, 
                                     assays = "RNA", 
                                     features = all_genes_of_interest, 
                                     group.by = "treatment_tissue")
avg_exp_grouped <- avg_exp_grouped$RNA

# Create annotation for the heatmap
avg_exp_annot_df <- data.frame(col_names = colnames(avg_exp_grouped)) %>%
  tidyr::separate(col_names, into = c("treatment", "tissue"), sep = "_", remove = FALSE)

ha_avg <- ComplexHeatmap::HeatmapAnnotation(
  Treatment = avg_exp_annot_df$treatment,
  Tissue = avg_exp_annot_df$tissue,
  col
  