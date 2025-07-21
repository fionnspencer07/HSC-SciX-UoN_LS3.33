# ============================================================================ #
# Comprehensive DGE and Visualization Pipeline for NK/T Cell Analysis (Run 2 v7)
#
# Author: GitHub Copilot, based on specifications from fionnspencer07
# Date: 2025-07-15
#
# Description:
# This script performs a multi-method DGE analysis using Wilcoxon, DESeq2 (on
# pseudobulk data), and MAST. It generates extensive visualizations using
# ComplexHeatmap and detailed, multi-sheet Excel reports for each method.
#
# v2 Fix: Used `dplyr::select()`.
# v3 Fix: Shortened Excel sheet names.
# v4/v5 Fix: Updated pseudobulking to use `scuttle::summarizeAssayByGroup`.
# v6 Fix: Corrected metadata column name to `tissue`.
# v7 Fix: Corrected the `column_to_rownames` error by creating the metadata
#         data frame manually, which robustly handles pre-existing rownames.
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
output_dir <- here("R Projects", "1507DGE", "output", "Run 2")
# Create the directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Load the focused Seurat object
message("Loading focused Seurat object...")
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))
message("Object loaded successfully.")

# Ensure 'treatment' is a factor with a logical order for comparisons
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, levels = c("NAIVE", "SHAM", "UT"))

# Define Genes of Interest
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")
all_genes_of_interest <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))


# --- 2. HELPER FUNCTION FOR EXCEL EXPORT ---
create_dge_excel <- function(dge_results_list, method_name, output_dir) {
  wb <- createWorkbook()
  all_contrasts_df <- bind_rows(
    dge_results_list[["UT_vs_SHAM"]] %>% mutate(contrast = "UT_vs_SHAM"),
    dge_results_list[["UT_vs_NAIVE"]] %>% mutate(contrast = "UT_vs_NAIVE"),
    dge_results_list[["SHAM_vs_NAIVE"]] %>% mutate(contrast = "SHAM_vs_NAIVE")
  ) %>% dplyr::select(contrast, everything())
  
  addWorksheet(wb, "All_Genes_UT_vs_SHAM"); writeData(wb, "All_Genes_UT_vs_SHAM", dge_results_list[["UT_vs_SHAM"]])
  addWorksheet(wb, "All_Genes_UT_vs_NAIVE"); writeData(wb, "All_Genes_UT_vs_NAIVE", dge_results_list[["UT_vs_NAIVE"]])
  addWorksheet(wb, "All_Genes_SHAM_vs_NAIVE"); writeData(wb, "All_Genes_SHAM_vs_NAIVE", dge_results_list[["SHAM_vs_NAIVE"]])
  
  addWorksheet(wb, "S1PR_Genes_All_Contrasts"); writeData(wb, "S1PR_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% s1pr_genes))
  addWorksheet(wb, "Trafficking_Genes_All"); writeData(wb, "Trafficking_Genes_All", all_contrasts_df %>% filter(gene %in% trafficking_genes))
  addWorksheet(wb, "Retention_Genes_All_Contrasts"); writeData(wb, "Retention_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% retention_genes))
  addWorksheet(wb, "Egress_Genes_All_Contrasts"); writeData(wb, "Egress_Genes_All_Contrasts", all_contrasts_df %>% filter(gene %in% egress_genes))
  
  summary_stats <- all_contrasts_df %>% group_by(contrast) %>% summarise(total_genes_tested = n(), significant_upregulated = sum(p_val_adj < 0.05 & avg_log2FC > 0.25, na.rm = TRUE), significant_downregulated = sum(p_val_adj < 0.05 & avg_log2FC < -0.25, na.rm = TRUE), .groups = 'drop')
  addWorksheet(wb, "Summary_Stats_By_Contrast"); writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
  
  top_up <- all_contrasts_df %>% filter(p_val_adj < 0.05) %>% group_by(contrast) %>% arrange(desc(avg_log2FC), .by_group = TRUE) %>% slice_head(n = 100)
  top_down <- all_contrasts_df %>% filter(p_val_adj < 0.05) %>% group_by(contrast) %>% arrange(avg_log2FC, .by_group = TRUE) %>% slice_head(n = 100)
  
  addWorksheet(wb, "Top100_Upregulated"); writeData(wb, "Top100_Upregulated", top_up)
  addWorksheet(wb, "Top100_Downregulated"); writeData(wb, "Top100_Downregulated", top_down)
  
  file_path <- file.path(output_dir, glue("DGE_{method_name}_Complete_Analysis.xlsx"))
  saveWorkbook(wb, file_path, overwrite = TRUE)
  message(glue("SUCCESS: Saved Excel report to {file_path}"))
}


# --- 3. DGE ANALYSIS ---

# METHOD 1: Wilcoxon (via Presto)
message("\n--- Starting DGE Analysis: Wilcoxon (presto) ---")
Idents(seu_NKT_focused) <- "treatment"
wilcox_results <- list()
wilcox_results[["UT_vs_SHAM"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'UT', ident.2 = 'SHAM') %>% as_tibble(rownames = "gene") %>% rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)
wilcox_results[["UT_vs_NAIVE"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'UT', ident.2 = 'NAIVE') %>% as_tibble(rownames = "gene") %>% rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)
wilcox_results[["SHAM_vs_NAIVE"]] <- presto::wilcoxauc(seu_NKT_focused, group_by = 'treatment', seurat_assay = 'RNA', ident.1 = 'SHAM', ident.2 = 'NAIVE') %>% as_tibble(rownames = "gene") %>% rename(avg_log2FC = logFC, p_val_adj = padj, p_val = pval)
create_dge_excel(wilcox_results, "Wilcoxon", output_dir)

# METHOD 2: DESeq2 on Pseudobulk Data
message("\n--- Starting DGE Analysis: DESeq2 on Pseudobulk Data ---")
counts_matrix <- GetAssayData(seu_NKT_focused, assay = "RNA", slot = "counts")
cell_metadata <- seu_NKT_focused@meta.data

message("  Performing pseudobulking by sample...")
pb <- scuttle::summarizeAssayByGroup(x = counts_matrix, ids = cell_metadata$sample_name, statistics = "sum")
pb_counts <- assay(pb, "sum")

# Create sample-level metadata, using `tissue` column
# ** FIX APPLIED HERE: This new method is robust to pre-existing rownames **
sample_info <- cell_metadata %>% 
  dplyr::select(sample_name, treatment, tissue) %>% 
  distinct()

sample_metadata <- as.data.frame(sample_info)
rownames(sample_metadata) <- sample_metadata$sample_name

pb_counts <- pb_counts[, rownames(sample_metadata)]

dds <- DESeqDataSetFromMatrix(countData = pb_counts, colData = sample_metadata, design = ~ treatment)
message("  Running DESeq2...")
dds <- DESeq(dds)
deseq2_results <- list()
deseq2_results[["UT_vs_SHAM"]] <- results(dds, contrast=c("treatment", "UT", "SHAM")) %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble() %>% rename(avg_log2FC = log2FoldChange, p_val = pvalue, p_val_adj = padj)
deseq2_results[["UT_vs_NAIVE"]] <- results(dds, contrast=c("treatment", "UT", "NAIVE")) %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble() %>% rename(avg_log2FC = log2FoldChange, p_val = pvalue, p_val_adj = padj)
deseq2_results[["SHAM_vs_NAIVE"]] <- results(dds, contrast=c("treatment", "SHAM", "NAIVE")) %>% as.data.frame() %>% rownames_to_column("gene") %>% as_tibble() %>% rename(avg_log2FC = log2FoldChange, p_val = pvalue, p_val_adj = padj)
create_dge_excel(deseq2_results, "DESeq2_Pseudobulk", output_dir)


# METHOD 3: MAST
message("\n--- Starting DGE Analysis: MAST ---")
Idents(seu_NKT_focused) <- "treatment"
mast_results <- list()
mast_results[["UT_vs_SHAM"]] <- FindMarkers(seu_NKT_focused, ident.1 = "UT", ident.2 = "SHAM", test.use = "MAST", logfc.threshold = 0) %>% rownames_to_column("gene") %>% as_tibble()
mast_results[["UT_vs_NAIVE"]] <- FindMarkers(seu_NKT_focused, ident.1 = "UT", ident.2 = "NAIVE", test.use = "MAST", logfc.threshold = 0) %>% rownames_to_column("gene") %>% as_tibble()
mast_results[["SHAM_vs_NAIVE"]] <- FindMarkers(seu_NKT_focused, ident.1 = "SHAM", ident.2 = "NAIVE", test.use = "MAST", logfc.threshold = 0) %>% rownames_to_column("gene") %>% as_tibble()
create_dge_excel(mast_results, "MAST", output_dir)


# --- 4. VISUALIZATION ---
message("\n--- Generating Visualizations for Run 2 ---")
save_svg <- function(plot, filename, width = 10, height = 8) {
  path <- file.path(output_dir, filename)
  if (inherits(plot, "Heatmap")) { svg(path, width = width, height = height); draw(plot); dev.off() } 
  else { ggsave(path, plot, width = width, height = height, device = "svg") }
  message(glue("  Saved plot: {filename}"))
}
treatment_colors <- c("NAIVE" = "#619CFF", "SHAM" = "#00BA38", "UT" = "#F8766D")
tissue_colors <- c("BM" = "#7CAE00", "WB" = "#C77CFF")
cell_type_colors <- c("T_cells" = "#f3e65d", "NK_cells" = "#86b0cc", "NKT_cells" = "#eeb84c")

plot_volcano <- function(dge_df, title) {
  data_to_plot <- dge_df %>% mutate(highlight = case_when(gene %in% s1pr_genes ~ "S1P Pathway", gene %in% trafficking_genes ~ "Trafficking", p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 ~ "Significant", TRUE ~ "Not Significant"), label = if_else(gene %in% all_genes_of_interest, gene, ""))
  p <- ggplot(data_to_plot, aes(x = avg_log2FC, y = -log10(p_val_adj), color = highlight, label = label)) + geom_point(alpha = 0.7, size = 2) + geom_text_repel(max.overlaps = 15, size = 3, bg.color = "white", bg.r = 0.1) + scale_color_manual(name = "Gene Group", values = c("S1P Pathway" = "#E41A1C", "Trafficking" = "#377EB8", "Significant" = "black", "Not Significant" = "grey80")) + geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50") + geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50") + labs(title = title, x = "Log2 Fold Change", y = "-log10(Adjusted P-value)") + theme_bw(base_size = 14) + theme(legend.position = "bottom")
  return(p)
}
save_svg(plot_volcano(wilcox_results$UT_vs_NAIVE, "Volcano: Wilcoxon (UT vs NAIVE)"), "Volcano_Wilcoxon_UT_vs_NAIVE.svg")
save_svg(plot_volcano(deseq2_results$UT_vs_NAIVE, "Volcano: DESeq2 Pseudobulk (UT vs NAIVE)"), "Volcano_DESeq2_UT_vs_NAIVE.svg")
save_svg(plot_volcano(mast_results$UT_vs_NAIVE, "Volcano: MAST (UT vs NAIVE)"), "Volcano_MAST_UT_vs_NAIVE.svg")

avg_exp_grouped <- AverageExpression(seu_NKT_focused, assays = "RNA", features = all_genes_of_interest, group.by = c("treatment", "tissue"))
avg_exp_grouped <- avg_exp_grouped$RNA
avg_exp_annot_df <- data.frame(col_names = colnames(avg_exp_grouped)) %>% separate(col_names, into = c("treatment", "tissue"), sep = "_", remove = FALSE)
ha_avg <- HeatmapAnnotation(Treatment = avg_exp_annot_df$treatment, Tissue = avg_exp_annot_df$tissue, col = list(Treatment = treatment_colors, Tissue = tissue_colors), annotation_name_side = "left")
ht_all_goi <- Heatmap(t(scale(t(avg_exp_grouped))), name = "Scaled Avg Exp", column_title = "Average Expression of All Genes of Interest", top_annotation = ha_avg, cluster_columns = TRUE, row_split = case_when(rownames(avg_exp_grouped) %in% s1pr_genes ~ "S1P Pathway", rownames(avg_exp_grouped) %in% egress_genes ~ "Egress Genes", rownames(avg_exp_grouped) %in% retention_genes ~ "Retention Genes", TRUE ~ "Other Trafficking"), row_title_rot = 0, show_row_dend = FALSE)
save_svg(ht_all_goi, "ComplexHeatmap_AvgExp_GenesOfInterest.svg", width = 8, height = 10)

plot_top50_heatmap <- function(dge_res, method_name) {
  top50_genes <- dge_res %>% arrange(p_val_adj, desc(abs(avg_log2FC))) %>% head(50) %>% pull(gene)
  top50_data <- GetAssayData(seu_NKT_focused, slot = "scale.data", assay = "RNA")[top50_genes, ]
  ha_cells <- HeatmapAnnotation(Treatment = seu_NKT_focused$treatment, Tissue = seu_NKT_focused$tissue, `Cell Type` = seu_NKT_focused$cell_types, col = list(Treatment = treatment_colors, Tissue = tissue_colors, `Cell Type` = cell_type_colors))
  ht <- Heatmap(top50_data, name = "Scaled Exp", column_title = glue("Top 50 DEGs ({method_name}, UT vs NAIVE)"), show_column_names = FALSE, top_annotation = ha_cells, column_split = seu_NKT_focused$treatment, use_raster = TRUE, raster_quality = 4)
  return(ht)
}
save_svg(plot_top50_heatmap(wilcox_results$UT_vs_NAIVE, "Wilcoxon"), "ComplexHeatmap_Top50_Wilcoxon.svg", width = 12, height = 10)
save_svg(plot_top50_heatmap(deseq2_results$UT_vs_NAIVE, "DESeq2"), "ComplexHeatmap_Top50_DESeq2.svg", width = 12, height = 10)
save_svg(plot_top50_heatmap(mast_results$UT_vs_NAIVE, "MAST"), "ComplexHeatmap_Top50_MAST.svg", width = 12, height = 10)

save_svg(dittoPlot(seu_NKT_focused, "S1pr1", group.by = "treatment", plots = c("vlnplot", "boxplot"), color.panel = treatment_colors) + labs(title = "S1pr1 Expression by Treatment") + stat_compare_means(comparisons = list(c("UT", "NAIVE"), c("UT", "SHAM")), label = "p.signif"), "VlnBox_S1pr1_by_Treatment.svg")
save_svg(dittoPlot(seu_NKT_focused, "S1pr1", group.by = "treatment", split.by = "tissue", plots = c("vlnplot", "jitter"), jitter.size = 0.1) + labs(title = "S1pr1 Expression by Treatment and Tissue"), "Vln_S1pr1_by_Tissue.svg", width = 12)
save_svg(Nebulosa::plot_density(seu_NKT_focused, "S1pr1") + facet_wrap(.~seu_NKT_focused$sample_name, ncol = 3) + theme_void(), "Nebulosa_S1pr1_by_Treatment.svg")

ranks <- wilcox_results[["UT_vs_NAIVE"]] %>% dplyr::select(gene, avg_log2FC) %>% na.omit() %>% distinct(gene, .keep_all = TRUE) %>% deframe()
pathways <- list(S1P_Pathway = s1pr_genes, Trafficking = trafficking_genes, BM_Retention = retention_genes, T_Cell_Egress = egress_genes)
fgsea_res <- fgsea(pathways = pathways, stats = ranks, minSize = 3, maxSize = 500)
save_svg(plotEnrichment(pathways[["S1P_Pathway"]], ranks) + labs(title = "GSEA: S1P Pathway (UT vs NAIVE)"), "GSEA_Enrichment_S1P_Pathway.svg")
save_svg(ggplot(fgsea_res, aes(x = NES, y = reorder(pathway, NES), size = size, color = pval)) + geom_point() + scale_color_gradient(low = "red", high = "blue") + theme_minimal(base_size = 12) + labs(title = "GSEA Pathway Analysis (UT vs NAIVE)", x = "Normalized Enrichment Score", y = "Pathway"), "GSEA_DotPlot_All_Pathways.svg")

message("\n--- Analysis pipeline for Run 2 complete! ---")
message(glue("All outputs have been saved to: {output_dir}"))
