# ============================================================================ #
# Comprehensive DGE and Visualization Pipeline for NK/T Cell Analysis
#
# Author: GitHub Copilot, based on specifications from fionnspencer07
# Date: 2025-07-15
#
# Description:
# This script performs a multi-method DGE analysis (Wilcoxon, Muscat+DESeq2,
# Muscat+edgeR) on a Seurat object of NK/T cells. It generates extensive
# visualizations and detailed, multi-sheet Excel reports as per the user's
# comprehensive requirements.
# ============================================================================ #


# --- 1. INITIAL SETUP & DATA LOADING ---

# Set seed for reproducibility across all analyses
set.seed(42)

# Load package libraries from the user-defined init script
# This assumes your R project root is 'HSC SciX UoN LS3.33'
library(here)
source(here("Custom R Functions and Scripts", "init_packages.R"))

# Define output directory using here()
output_dir <- here("R Projects", "1507DGE", "output", "Run 1")
# Create the directory if it doesn't exist, suppressing warnings if it does
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the subsetted and focused Seurat object
# This assumes you have already run your setup script to generate this file
message("Loading focused Seurat object...")
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))
message("Object loaded successfully.")

# Verify metadata columns needed for analysis
required_cols <- c("treatment", "Compartment", "cell_types", "sample_name")
if (!all(required_cols %in% colnames(seu_NKT_focused@meta.data))) {
  stop("Error: One or more required metadata columns are missing. Please check: ", 
       paste(required_cols, collapse=", "))
}
# Ensure 'treatment' is a factor with a logical order for comparisons
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, levels = c("NAIVE", "SHAM", "UT"))


# Define Genes of Interest (using user's lists)
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")

# Create a master list of all unique genes of interest for highlighting
all_genes_of_interest <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))


# --- 2. HELPER FUNCTION FOR EXCEL EXPORT ---

#' Create a comprehensive, multi-sheet Excel workbook for DGE results.
#'
#' @param dge_results_list A list of data frames, one for each contrast.
#' @param method_name A string for the DGE method used (e.g., "Wilcoxon").
#' @param output_dir The path to the output directory.
create_dge_excel <- function(dge_results_list, method_name, output_dir) {
  wb <- createWorkbook()
  
  # Combine all contrast results into a single dataframe for easy filtering
  all_contrasts_df <- bind_rows(
    dge_results_list[["UT_vs_SHAM"]] %>% mutate(contrast = "UT_vs_SHAM"),
    dge_results_list[["UT_vs_NAIVE"]] %>% mutate(contrast = "UT_vs_NAIVE"),
    dge_results_list[["SHAM_vs_NAIVE"]] %>% mutate(contrast = "SHAM_vs_NAIVE")
  ) %>% select(contrast, everything()) # Ensure contrast is the first column
  
  # Sheet 1-3: Full, unfiltered results for each main contrast
  addWorksheet(wb, "All_Genes_UT_vs_SHAM")
  writeData(wb, "All_Genes_UT_vs_SHAM", dge_results_list[["UT_vs_SHAM"]])
  
  addWorksheet(wb, "All_Genes_UT_vs_NAIVE")
  writeData(wb, "All_Genes_UT_vs_NAIVE", dge_results_list[["UT_vs_NAIVE"]])
  
  addWorksheet(wb, "All_Genes_SHAM_vs_NAIVE")
  writeData(wb, "All_Genes_SHAM_vs_NAIVE", dge_results_list[["SHAM_vs_NAIVE"]])
  
  # Sheet 4-7: Genes of Interest subsets
  addWorksheet(wb, "S1PR_Genes_All_Contrasts")
  writeData(wb, "S1PR_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% s1pr_genes))
  
  addWorksheet(wb, "Trafficking_Genes_All_Contrasts")
  writeData(wb, "Trafficking_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% trafficking_genes))
  
  addWorksheet(wb, "Retention_Genes_All_Contrasts")
  writeData(wb, "Retention_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% retention_genes))
  
  addWorksheet(wb, "Egress_Genes_All_Contrasts")
  writeData(wb, "Egress_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% egress_genes))
  
  # Sheet 8: Summary Statistics
  summary_stats <- all_contrasts_df %>%
    group_by(contrast) %>%
    summarise(
      total_genes_tested = n(),
      significant_upregulated = sum(p_val_adj < 0.05 & avg_log2FC > 0.25, na.rm = TRUE),
      significant_downregulated = sum(p_val_adj < 0.05 & avg_log2FC < -0.25, na.rm = TRUE),
      .groups = 'drop'
    )
  addWorksheet(wb, "Summary_Stats_By_Contrast")
  writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
  
  # Sheet 9-10: Top 100 Up- and Down-regulated genes for each contrast
  top_up <- all_contrasts_df %>% filter(p_val_adj < 0.05) %>% group_by(contrast) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 100)
  top_down <- all_contrasts_df %>% filter(p_val_adj < 0.05) %>% group_by(contrast) %>% arrange(avg_log2FC, .by_group = TRUE) %>% slice_head(n = 100)
  
  addWorksheet(wb, "Top100_Upregulated_by_Contrast")
  writeData(wb, "Top100_Upregulated_by_Contrast", top_up)
  
  addWorksheet(wb, "Top100_Downregulated_by_Contrast")
  writeData(wb, "Top100_Downregulated_by_Contrast", top_down)
  
  # Save the workbook
  file_path <- file.path(output_dir, glue("DGE_{method_name}_Complete_Analysis.xlsx"))
  saveWorkbook(wb, file_path, overwrite = TRUE)
  message(glue("SUCCESS: Saved Excel report to {file_path}"))
}

# --- 3. DGE ANALYSIS: METHOD 1 - WILCOXON (via Presto) ---
message("\n--- Starting DGE Analysis: Wilcoxon (presto) ---")
# Set Idents to the 'treatment' column for comparison
Idents(seu_NKT_focused) <- "treatment"
wilcox_results <- list()

# Perform pairwise comparisons
message("  Comparing UT vs SHAM...")
wilcox_results[["UT_vs_SHAM"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'UT', ident.2 = 'SHAM') %>%
  as_tibble(rownames = "gene") %>%
  rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)

message("  Comparing UT vs NAIVE...")
wilcox_results[["UT_vs_NAIVE"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'UT', ident.2 = 'NAIVE') %>%
  as_tibble(rownames = "gene") %>%
  rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)

message("  Comparing SHAM vs NAIVE...")
wilcox_results[["SHAM_vs_NAIVE"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'SHAM', ident.2 = 'NAIVE') %>%
  as_tibble(rownames = "gene") %>%
  rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)

# Create the detailed Excel file for Wilcoxon results
create_dge_excel(wilcox_results, "Wilcoxon", output_dir)


# --- 4. DGE ANALYSIS: METHOD 2 & 3 - MUSCAT (Pseudobulk) ---
message("\n--- Starting DGE Analysis: Muscat (DESeq2 & edgeR) ---")

# Create a SingleCellExperiment object required by muscat
sce <- as.SingleCellExperiment(seu_NKT_focused, assay = "RNA")

# Define study design variables in the SCE object
sce$sample_id <- sce$sample_name
sce$cluster_id <- sce$cell_types # Using cell_types for pseudobulking
sce$group_id <- sce$treatment

# Aggregate counts to pseudobulk samples
message("  Aggregating data to pseudobulk...")
pb <- aggregateData(sce,
                    assay = "counts",
                    by = c("cluster_id", "sample_id"),
                    fun = "sum")

# Define the design formula for the GLM
design <- model.matrix(~ 0 + group_id, data = metadata(pb)$experiment_info)

# Define contrasts for pairwise comparisons
contrast_matrix <- makeContrasts(
  UT_vs_NAIVE = group_idUT - group_idNAIVE,
  UT_vs_SHAM = group_idUT - group_idSHAM,
  SHAM_vs_NAIVE = group_idSHAM - group_idNAIVE,
  levels = design
)

# Helper function to process and average muscat results for the Excel sheet
average_muscat_results <- function(res) {
  # This extracts results for all contrasts and cell types
  tbl <- res$table
  
  # Combine results from all cell types and contrasts into one table
  all_res <- lapply(names(tbl), function(contrast_name) {
    contrast_res <- tbl[[contrast_name]]
    bind_rows(lapply(names(contrast_res), function(cluster_name) {
      as.data.frame(contrast_res[[cluster_name]]) %>%
        rownames_to_column("gene") %>%
        mutate(cluster_id = cluster_name)
    })) %>% mutate(contrast = contrast_name)
  }) %>% bind_rows()
  
  # Create a simplified, averaged result list for the main Excel function
  # NOTE: Averaging p-values is not statistically rigorous but is done here
  # for a high-level overview as requested.
  averaged_list <- list()
  for (contrast_name in names(tbl)) {
    averaged_list[[contrast_name]] <- all_res %>%
      filter(contrast == contrast_name) %>%
      group_by(gene) %>%
      summarise(
        avg_log2FC = mean(logFC, na.rm = TRUE),
        p_val = mean(p_val, na.rm = TRUE),
        p_val_adj = mean(p_adj.loc, na.rm = TRUE),
        .groups = 'drop'
      ) %>%
      arrange(p_val_adj)
  }
  return(averaged_list)
}

# Run DGE with DESeq2
message("  Running Muscat with DESeq2...")
pb_res_deseq <- pbDS(pb, design = design, contrast = contrast_matrix, method = "DESeq2")
deseq_results <- average_muscat_results(pb_res_deseq)
create_dge_excel(deseq_results, "Muscat_DESeq2", output_dir)

# Run DGE with edgeR
message("  Running Muscat with edgeR...")
pb_res_edger <- pbDS(pb, design = design, contrast = contrast_matrix, method = "edgeR")
edger_results <- average_muscat_results(pb_res_edger)
create_dge_excel(edger_results, "Muscat_edgeR", output_dir)


# --- 5. VISUALIZATION ---
message("\n--- Generating Visualizations ---")

# Helper function to save plots as SVG
save_svg <- function(plot, filename, width = 10, height = 8) {
  path <- file.path(output_dir, filename)
  ggsave(path, plot, width = width, height = height, device = "svg")
  message(glue("  Saved plot: {filename}"))
}

# 5.1 Volcano Plots (for each DGE method, main contrast)
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
    scale_color_manual(
      name = "Gene Group",
      values = c("S1P Pathway" = "#E41A1C", "Trafficking" = "#377EB8", "Significant" = "black", "Not Significant" = "grey80")
    ) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") +
    labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") +
    theme_bw(base_size = 14) +
    theme(legend.position = "bottom")
  return(p)
}

save_svg(plot_volcano(wilcox_results$UT_vs_NAIVE, "Volcano: Wilcoxon (UT vs NAIVE)"), "Volcano_Wilcoxon_UT_vs_NAIVE.svg")
save_svg(plot_volcano(deseq_results$UT_vs_NAIVE, "Volcano: Muscat+DESeq2 (UT vs NAIVE)"), "Volcano_DESeq2_UT_vs_NAIVE.svg")
save_svg(plot_volcano(edger_results$UT_vs_NAIVE, "Volcano: edgeR (UT vs NAIVE)"), "Volcano_edgeR_UT_vs_NAIVE.svg")

# 5.2 Heatmaps
# Get average expression per group (treatment x compartment) for plotting
avg_exp_grouped <- AverageExpression(seu_NKT_focused, assays = "RNA", features = all_genes_of_interest, group.by = c("treatment", "Compartment"))$RNA
scaled_avg_exp <- t(scale(t(avg_exp_grouped)))

# Heatmap for S1P genes
s1p_heatmap <- pheatmap(
  scaled_avg_exp[intersect(s1pr_genes, rownames(scaled_avg_exp)), ],
  main = "Expression of S1P Pathway Genes",
  cluster_rows = TRUE, cluster_cols = TRUE,
  scale = "none", # Already scaled
  color = rev(RColorBrewer::brewer.pal(n = 9, name = "RdBu")),
  border_color = "grey60",
  fontsize = 10
)
save_svg(ggplotify::as.ggplot(s1p_heatmap), "Heatmap_S1P_Genes.svg")

# Heatmap for top 50 DEGs from Wilcoxon
top50_wilcox <- wilcox_results$UT_vs_NAIVE %>% arrange(p_val_adj) %>% head(50)
top50_data <- GetAssayData(seu_NKT_focused, slot = "scale.data", assay = "RNA")[top50_wilcox$gene, ]
# Create annotation for heatmap columns
col_annot <- seu_NKT_focused@meta.data[, c("treatment", "Compartment", "cell_types")]

top50_heatmap <- pheatmap(
  top50_data, 
  main = "Top 50 DEGs (Wilcoxon, UT vs NAIVE)",
  show_colnames = FALSE,
  annotation_col = col_annot,
  scale = "none" # Already scaled
)
save_svg(ggplotify::as.ggplot(top50_heatmap), "Heatmap_Top50_Wilcoxon.svg", width = 14, height = 12)

# 5.3 Violin and Box Plots
# Using dittoSeq for powerful and easy plotting
save_svg(dittoPlot(seu_NKT_focused, "S1pr1", group.by = "treatment", plots = c("vlnplot", "boxplot"), 
                   color.panel = c("NAIVE" = "#619CFF", "SHAM" = "#00BA38", "UT" = "#F8766D")) +
           labs(title = "S1pr1 Expression by Treatment") +
           stat_compare_means(comparisons = list(c("UT", "NAIVE"), c("UT", "SHAM")), label = "p.signif"),
         "VlnBox_S1pr1_by_Treatment.svg")

# Split by compartment
save_svg(dittoPlot(seu_NKT_focused, "S1pr1", group.by = "treatment", split.by = "Compartment",
                   plots = c("vlnplot", "jitter"), jitter.size = 0.1) +
           labs(title = "S1pr1 Expression by Treatment and Compartment"),
         "Vln_S1pr1_by_Compartment.svg", width = 12)

# 5.4 Nebulosa Density Plots
save_svg(plot_density(seu_NKT_focused, features = "S1pr1", reduction = "harmony", group.by = "treatment") +
           labs(title = "S1pr1 Expression Density by Treatment (UMAP)"),
         "Nebulosa_S1pr1_by_Treatment.svg")

save_svg(plot_density(seu_NKT_focused, features = "S1pr1", reduction = "harmony", group.by = "Compartment") +
           labs(title = "S1pr1 Expression Density by Compartment (UMAP)"),
         "Nebulosa_S1pr1_by_Compartment.svg")

# 5.5 GSEA Pathway Analysis
# Prepare ranked list from Wilcoxon results
ranks <- wilcox_results[["UT_vs_NAIVE"]] %>%
  select(gene, avg_log2FC) %>%
  na.omit() %>%
  distinct(gene, .keep_all = TRUE) %>%
  deframe()

# Define pathways from user's gene lists
pathways <- list(
  S1P_Pathway = s1pr_genes,
  Trafficking = trafficking_genes,
  BM_Retention = retention_genes,
  T_Cell_Egress = egress_genes
)

fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 3, maxSize = 500)

# GSEA Enrichment Plot
save_svg(plotEnrichment(pathways[["S1P_Pathway"]], ranks) + labs(title = "GSEA: S1P Pathway (UT vs NAIVE)"),
         "GSEA_Enrichment_S1P_Pathway.svg")

# GSEA Dot Plot for all tested pathways
save_svg(ggplot(fgsea_res, aes(x = NES, y = reorder(pathway, NES), size = size, color = pval)) +
           geom_point() +
           scale_color_gradient(low = "red", high = "blue") +
           theme_minimal(base_size = 12) +
           labs(title = "GSEA Pathway Analysis (UT vs NAIVE)", x = "Normalized Enrichment Score", y = "Pathway"),
         "GSEA_DotPlot_All_Pathways.svg")

message("\n--- Analysis pipeline complete! ---")
message(glue("All outputs have been saved to: {output_dir}"))