# ============================================================================ #
# NK/T Cell DGE Comprehensive Pipeline with All Pairwise Contrasts & Reporting
# Author: fionnspencer07 | Date: 2025-07-16
# ============================================================================ #

set.seed(42)
plan("multicore", workers = 4)

library(here)
source(here("Custom R Functions and Scripts", "init_packages.R"))

# ---- 1. Load Focused Seurat Object ----
output_dir <- here("R Projects", "1607DGE", "output", "Run 1")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
message("Loading focused Seurat object...")
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))
message("Object loaded successfully.")

# ---- 2. Metadata and Genes of Interest ----
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, levels = c("NAIVE", "SHAM", "UT")) # UT = DMG
seu_NKT_focused$Compartment <- factor(ifelse(seu_NKT_focused$tissue == "WB", "WB", "BM"))
seu_NKT_focused$group <- paste(seu_NKT_focused$treatment, seu_NKT_focused$Compartment, sep = "_")

contrast_df <- tribble(
  ~A,           ~B,
  "NAIVE_WB",   "NAIVE_BM",
  "NAIVE_WB",   "SHAM_WB",
  "NAIVE_WB",   "SHAM_BM",
  "NAIVE_WB",   "UT_WB",
  "NAIVE_WB",   "UT_BM",
  "NAIVE_BM",   "SHAM_WB",
  "NAIVE_BM",   "SHAM_BM",
  "NAIVE_BM",   "UT_WB",
  "NAIVE_BM",   "UT_BM",
  "SHAM_WB",    "SHAM_BM",
  "SHAM_WB",    "UT_WB",
  "SHAM_WB",    "UT_BM",
  "SHAM_BM",    "UT_WB",
  "SHAM_BM",    "UT_BM",
  "UT_WB",      "UT_BM"
)

s1p_genes <- c("S1pr1", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
bm_retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")
activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
chemokines <- c("Ccl5", "Ccl4", "Cxcl10", "Cxcl9")
cytokines <- c("Ifng", "Il2", "Tnfa", "Il10")
all_goi <- unique(c(s1p_genes, trafficking_genes, bm_retention_genes, activation_genes, chemokines, cytokines))

# ---- 3. Excel Export Helper ----
create_dge_excel <- function(dge_results, method_name, output_dir, all_goi, s1p_genes, trafficking_genes) {
  wb <- createWorkbook()
  for (i in 1:nrow(contrast_df)) {
    comp <- glue("{contrast_df$A[i]}_vs_{contrast_df$B[i]}")
    dat <- dge_results[[comp]]
    addWorksheet(wb, glue("All_Genes_{comp}"))
    writeData(wb, glue("All_Genes_{comp}"), dat)
  }
  all_contrasts <- bind_rows(lapply(names(dge_results), function(nm) dge_results[[nm]] %>% mutate(contrast = nm)), .id = "contrast_idx")
  addWorksheet(wb, "S1P_Genes_All_Contrasts")
  writeData(wb, "S1P_Genes_All_Contrasts", filter(all_contrasts, gene %in% s1p_genes))
  addWorksheet(wb, "Trafficking_Genes_All_Contrasts")
  writeData(wb, "Trafficking_Genes_All_Contrasts", filter(all_contrasts, gene %in% trafficking_genes))
  summary_stats <- all_contrasts %>%
    group_by(contrast) %>%
    summarise(total_genes = n(), sig_up = sum(p_val_adj < 0.05 & avg_log2FC > 0.25, na.rm=TRUE),
              sig_down = sum(p_val_adj < 0.05 & avg_log2FC < -0.25, na.rm=TRUE), .groups = "drop")
  addWorksheet(wb, "Summary_Stats_By_Contrast")
  writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
  addWorksheet(wb, "GOI_by_Compartment")
  writeData(wb, "GOI_by_Compartment", filter(all_contrasts, gene %in% all_goi))
  top_up <- all_contrasts %>% filter(p_val_adj < 0.05) %>%
    group_by(contrast) %>% arrange(desc(avg_log2FC), .by_group=TRUE) %>% slice_head(n=50)
  addWorksheet(wb, "Top_Upregulated_Each_Contrast")
  writeData(wb, "Top_Upregulated_Each_Contrast", top_up)
  top_down <- all_contrasts %>% filter(p_val_adj < 0.05) %>%
    group_by(contrast) %>% arrange(avg_log2FC, .by_group=TRUE) %>% slice_head(n=50)
  addWorksheet(wb, "Top_Downregulated_Each_Contrast")
  writeData(wb, "Top_Downregulated_Each_Contrast", top_down)
  fname <- file.path(output_dir, glue("DGE_{method_name}_Complete_Analysis.xlsx"))
  saveWorkbook(wb, fname, overwrite=TRUE)
  message(glue("Saved Excel report: {fname}"))
}

# ---- 4. Wilcoxon DGE (Presto) ----
Idents(seu_NKT_focused) <- "group"
wilcox_results <- list()
message("Wilcoxon DGE pairwise...")
for (i in 1:nrow(contrast_df)) {
  g1 <- contrast_df$A[i]; g2 <- contrast_df$B[i]
  comp <- glue("{g1}_vs_{g2}")
  tryCatch({
    res <- presto::wilcoxauc(seu_NKT_focused, group_by = "group", seurat_assay = "RNA", ident.1=g1, ident.2=g2) %>%
      as_tibble(rownames="gene") %>% rename(avg_log2FC = logFC, p_val_adj=padj, p_val=pval)
    wilcox_results[[comp]] <- res
    message(glue("  {comp} complete."))
  }, error=function(e) {
    warning(glue("Wilcoxon failed for {comp}: {e$message}"))
    wilcox_results[[comp]] <- tibble()
  })
}
create_dge_excel(wilcox_results, "Wilcoxon", output_dir, all_goi, s1p_genes, trafficking_genes)

# ---- 5. muscat+DESeq2 and muscat+edgeR (Pseudobulk) ----
sce <- as.SingleCellExperiment(seu_NKT_focused)
colData(sce)$sample_id <- seu_NKT_focused$sample_name
colData(sce)$Compartment <- seu_NKT_focused$tissue

message("muscat pseudobulk DGE: DESeq2...")
pb_dge_deseq2 <- pbDS(sce, method="DESeq2")
muscat_deseq2_results <- list()
for (i in 1:nrow(contrast_df)) {
  g1 <- contrast_df$A[i]; g2 <- contrast_df$B[i]; comp <- glue("{g1}_vs_{g2}")
  tryCatch({
    res <- pb_dge_deseq2$table[[glue("{g1}_vs_{g2}")]]
    muscat_deseq2_results[[comp]] <- res %>% as_tibble(rownames="gene") %>%
      rename(avg_log2FC = logFC, p_val=pval, p_val_adj=p_adj.loc)
    message(glue("  muscat-DESeq2 {comp} done."))
  }, error=function(e) {
    warning(glue("muscat-DESeq2 failed: {comp}: {e$message}"))
    muscat_deseq2_results[[comp]] <- tibble()
  })
}
create_dge_excel(muscat_deseq2_results, "Muscat_DESeq2", output_dir, all_goi, s1p_genes, trafficking_genes)

message("muscat pseudobulk DGE: edgeR...")
pb_dge_edger <- pbDS(sce, group_col="group", method="edgeR", min_cells=10, min_samples=2)
muscat_edger_results <- list()
for (i in 1:nrow(contrast_df)) {
  g1 <- contrast_df$A[i]; g2 <- contrast_df$B[i]; comp <- glue("{g1}_vs_{g2}")
  tryCatch({
    res <- pb_dge_edger$table[[glue("{g1}_vs_{g2}")]]
    muscat_edger_results[[comp]] <- res %>% as_tibble(rownames="gene") %>%
      rename(avg_log2FC = logFC, p_val=pval, p_val_adj=p_adj.loc)
    message(glue("  muscat-edgeR {comp} done."))
  }, error=function(e) {
    warning(glue("muscat-edgeR failed: {comp}: {e$message}"))
    muscat_edger_results[[comp]] <- tibble()
  })
}
create_dge_excel(muscat_edger_results, "Muscat_edgeR", output_dir, all_goi, s1p_genes, trafficking_genes)

# ---- 6. Visualizations ----

save_svg <- function(plot, fname, w=10, h=8) {
  svg(file.path(output_dir, fname), width=w, height=h)
  print(plot)
  dev.off()
  message(glue("Saved: {fname}"))
}
treatment_colors <- c("NAIVE" = "#619CFF", "SHAM" = "#00BA38", "UT" = "#F8766D")
compartment_colors <- c("BM" = "#7CAE00", "WB" = "#C77CFF")

# -- 6.1. Heatmaps for top genes (Wilcoxon, DESeq2, edgeR, selected contrasts) --
for (method in c("Wilcoxon", "Muscat_DESeq2", "Muscat_edgeR")) {
  dge <- get(paste0(tolower(method), "_results"))
  for (i in c(1,5,10,15)) {
    comp <- names(dge)[i]
    top_genes <- dge[[comp]] %>% arrange(p_val_adj, desc(abs(avg_log2FC))) %>% filter(!is.na(gene)) %>% head(30) %>% pull(gene)
    mat <- GetAssayData(seu_NKT_focused, assay="RNA", slot="scale.data")[top_genes, ]
    ha <- HeatmapAnnotation(Treatment=seu_NKT_focused$treatment, Compartment=seu_NKT_focused$Compartment, col=list(Treatment=treatment_colors, Compartment=compartment_colors))
    ht <- Heatmap(mat, name="Scaled Exp", row_names_gp=gpar(fontsize=8), show_column_names=FALSE, top_annotation=ha)
    save_svg(ht, glue("Heatmap_{method}_{comp}.svg"), w=12, h=8)
  }
}

# -- 6.2. S1P pathway heatmaps --
s1p_mat <- GetAssayData(seu_NKT_focused, assay="RNA", slot="scale.data")[s1p_genes, ]
ha <- HeatmapAnnotation(Treatment=seu_NKT_focused$treatment, Compartment=seu_NKT_focused$Compartment, col=list(Treatment=treatment_colors, Compartment=compartment_colors))
ht <- Heatmap(s1p_mat, name="S1P Pathway", row_names_gp=gpar(fontsize=12), show_column_names=FALSE, top_annotation=ha)
save_svg(ht, "Heatmap_S1P_Pathway.svg", w=10, h=6)

# -- 6.3. Nebulosa density plots for S1P genes --
for (gene in s1p_genes) {
  p <- plot_density(seu_NKT_focused, features=gene, group.by="treatment", reduction="umap")
  save_svg(p, glue("Nebulosa_{gene}_by_treatment.svg"), w=8, h=6)
}

# -- 6.4. Volcano plots for top contrast for each method --
plot_volcano <- function(dge_df, title) {
  dge_df <- dge_df %>%
    mutate(highlight = case_when(
      gene %in% s1p_genes ~ "S1P",
      gene %in% trafficking_genes ~ "Trafficking",
      p_val_adj < 0.05 & abs(avg_log2FC) > 0.5 ~ "Significant",
      TRUE ~ "NS"),
      label = if_else(gene %in% c(s1p_genes, trafficking_genes), gene, ""))
  ggplot(dge_df, aes(x=avg_log2FC, y=-log10(p_val_adj), color=highlight, label=label)) +
    geom_point(alpha=0.6, size=1.7) +
    geom_text_repel(size=2.9, max.overlaps=12) +
    scale_color_manual(values=c("S1P"="red", "Trafficking"="blue", "Significant"="black", "NS"="grey80")) +
    geom_hline(yintercept=-log10(0.05), lty=2) +
    geom_vline(xintercept=c(-0.25,0.25), lty=2) +
    labs(title=title, x="Log2 Fold Change", y="-log10(adj P)") +
    theme_bw(base_size=13)
}
for (method in c("Wilcoxon", "Muscat_DESeq2", "Muscat_edgeR")) {
  dge <- get(paste0(tolower(method), "_results"))
  comp <- "NAIVE_WB_vs_UT_WB"
  if (comp %in% names(dge)) {
    p <- plot_volcano(dge[[comp]], glue("Volcano: {method} {comp}"))
    save_svg(p, glue("Volcano_{method}_{comp}.svg"), w=9, h=7)
  }
}

# -- 6.5. Violin/box plots for S1P, trafficking by treatment/compartment --
for (gene in c(s1p_genes, trafficking_genes)) {
  p <- dittoPlot(seu_NKT_focused, gene, group.by="treatment", split.by="Compartment", plots=c("vlnplot","jitter"), color.panel=treatment_colors) +
    labs(title=glue("{gene} by treatment/compartment"))
  save_svg(p, glue("Violin_{gene}_by_group.svg"), w=8, h=6)
}

# -- 6.6. GSEA for S1P/trafficking (Wilcoxon top contrast) --
ranks <- wilcox_results[["NAIVE_WB_vs_UT_WB"]] %>% select(gene, avg_log2FC) %>% deframe()
fgsea_res <- fgsea(pathways=list(S1P=s1p_genes, Trafficking=trafficking_genes), stats=ranks)
p <- ggplot(fgsea_res, aes(NES, reorder(pathway, NES), size=size, color=pval)) +
  geom_point() + scale_color_gradient(low="red", high="blue") +
  labs(title="GSEA: S1P/Trafficking", x="NES", y="Pathway") + theme_minimal()
save_svg(p, "GSEA_S1P_Trafficking.svg", w=8, h=6)

# ---- 7. Completion ----
message(glue("All outputs saved to: {output_dir}"))
message("Pipeline complete. Please check the output directory for Excel and SVG files.")