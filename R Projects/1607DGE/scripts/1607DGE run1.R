## Title: Complete S1P Receptor DGE Analysis Pipeline
## Research Question: Is the loss of Sphingosine-1-Phosphate (S1P) receptor externalization 
## on T-cells and impaired lymphocyte trafficking from bone marrow a direct or indirect 
## consequence of DMG presence?

# Record start time for duration calculation
start_time <- Sys.time()

set.seed(42)
library(here)
library(Seurat)
library(dplyr)
library(ggplot2)
library(ComplexHeatmap)
library(EnhancedVolcano)
library(openxlsx)
library(presto)
library(MAST)
library(muscat)
library(DESeq2)
library(edgeR)
library(monocle3)
library(Nebulosa)
library(viridis)
library(RColorBrewer)
library(pheatmap)
library(corrplot)
library(ggridges)
library(cowplot)
library(patchwork)
library(fgsea)
library(msigdbr)
library(clusterProfiler)
library(enrichplot)
library(tidyr)
library(stringr)
library(scales)
library(ggsci)
library(ggpubr)
library(gridExtra)
library(circlize)
library(ggrepel) # Added missing library for volcano plot labels

# Set up progress tracking
cat("Starting S1P Receptor DGE Analysis Pipeline...\n")

# Load the focused Seurat object
cat("Loading Seurat object...\n")
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))

# Create output directories
output_dir <- here("R Projects", "1607DGE", "output", "Run 1")
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "Excel_Files"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Heatmaps"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Volcano_Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Violin_Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Nebulosa_Plots"), showWarnings = FALSE)
dir.create(file.path(output_dir, "Trajectory_Analysis"), showWarnings = FALSE)

# Define comprehensive genes of interest
s1p_genes <- c("S1pr1", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")
activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
chemokine_genes <- c("Ccl5", "Ccl4", "Cxcl10", "Cxcl9")
cytokine_genes <- c("Ifng", "Il2", "Tnfa", "Il10")

all_genes_of_interest <- unique(c(s1p_genes, trafficking_genes, retention_genes, 
                                  activation_genes, chemokine_genes, cytokine_genes))

# Prepare metadata and set factor levels
cat("Preparing metadata and factor levels...\n")
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, 
                                    levels = c("NAIVE", "SHAM", "UT"))
seu_NKT_focused$tissue <- factor(seu_NKT_focused$tissue, 
                                 levels = c("WB", "BM"))

# Create combined treatment-tissue groups
seu_NKT_focused$group <- paste(seu_NKT_focused$treatment, seu_NKT_focused$tissue, sep = "_")
seu_NKT_focused$group <- factor(seu_NKT_focused$group, 
                                levels = c("NAIVE_WB", "NAIVE_BM", "SHAM_WB", 
                                           "SHAM_BM", "UT_WB", "UT_BM"))

# Define all 15 pairwise comparisons
contrasts_list <- list(
  "naive_WB_vs_naive_BM" = c("NAIVE_WB", "NAIVE_BM"),
  "naive_WB_vs_sham_WB" = c("NAIVE_WB", "SHAM_WB"),
  "naive_WB_vs_sham_BM" = c("NAIVE_WB", "SHAM_BM"),
  "naive_WB_vs_UT_WB" = c("NAIVE_WB", "UT_WB"),
  "naive_WB_vs_UT_BM" = c("NAIVE_WB", "UT_BM"),
  "naive_BM_vs_sham_WB" = c("NAIVE_BM", "SHAM_WB"),
  "naive_BM_vs_sham_BM" = c("NAIVE_BM", "SHAM_BM"),
  "naive_BM_vs_UT_WB" = c("NAIVE_BM", "UT_WB"),
  "naive_BM_vs_UT_BM" = c("NAIVE_BM", "UT_BM"),
  "sham_WB_vs_sham_BM" = c("SHAM_WB", "SHAM_BM"),
  "sham_WB_vs_UT_WB" = c("SHAM_WB", "UT_WB"),
  "sham_WB_vs_UT_BM" = c("SHAM_WB", "UT_BM"),
  "sham_BM_vs_UT_WB" = c("SHAM_BM", "UT_WB"),
  "sham_BM_vs_UT_BM" = c("SHAM_BM", "UT_BM"),
  "UT_WB_vs_UT_BM" = c("UT_WB", "UT_BM")
)

# Function to run Wilcoxon test
run_wilcoxon <- function(seu_obj, group1, group2, contrast_name) {
  cat(paste("Running Wilcoxon for", contrast_name, "...\n"))
  
  # Subset cells for comparison
  cells_group1 <- WhichCells(seu_obj, expression = group == group1)
  cells_group2 <- WhichCells(seu_obj, expression = group == group2)
  
  if(length(cells_group1) == 0 || length(cells_group2) == 0) {
    cat(paste("Skipping Wilcoxon for", contrast_name, "- one or both groups have zero cells.\n"))
    return(NULL)
  }
  
  # Set identity for comparison
  seu_obj$temp_ident <- "other"
  seu_obj$temp_ident[colnames(seu_obj) %in% cells_group1] <- group1
  seu_obj$temp_ident[colnames(seu_obj) %in% cells_group2] <- group2
  seu_obj$temp_ident <- factor(seu_obj$temp_ident, levels = c(group1, group2, "other"))
  
  # Run Wilcoxon test
  results <- FindMarkers(seu_obj, 
                         ident.1 = group1, 
                         ident.2 = group2,
                         group.by = "temp_ident",
                         test.use = "wilcox",
                         min.pct = 0.1,
                         logfc.threshold = 0.1)
  
  # Add gene names and additional statistics
  results$gene <- rownames(results)
  results$contrast <- contrast_name
  results$method <- "Wilcoxon"
  results$significant <- results$p_val_adj < 0.05
  results$effect_size <- abs(results$avg_log2FC)
  
  return(results)
}

# Function to run MAST
run_mast <- function(seu_obj, group1, group2, contrast_name) {
  cat(paste("Running MAST for", contrast_name, "...\n"))
  
  # Subset cells for comparison
  cells_group1 <- WhichCells(seu_obj, expression = group == group1)
  cells_group2 <- WhichCells(seu_obj, expression = group == group2)
  
  if(length(cells_group1) == 0 || length(cells_group2) == 0) {
    cat(paste("Skipping MAST for", contrast_name, "- one or both groups have zero cells.\n"))
    return(NULL)
  }
  
  # Set identity for comparison
  seu_obj$temp_ident <- "other"
  seu_obj$temp_ident[colnames(seu_obj) %in% cells_group1] <- group1
  seu_obj$temp_ident[colnames(seu_obj) %in% cells_group2] <- group2
  seu_obj$temp_ident <- factor(seu_obj$temp_ident, levels = c(group1, group2, "other"))
  
  # Run MAST
  results <- FindMarkers(seu_obj, 
                         ident.1 = group1, 
                         ident.2 = group2,
                         group.by = "temp_ident",
                         test.use = "MAST",
                         min.pct = 0.1,
                         logfc.threshold = 0.1)
  
  # Add gene names and additional statistics
  results$gene <- rownames(results)
  results$contrast <- contrast_name
  results$method <- "MAST"
  results$significant <- results$p_val_adj < 0.05
  results$effect_size <- abs(results$avg_log2FC)
  
  return(results)
}

# Function to run pseudobulk DESeq2
run_pseudobulk_deseq2 <- function(seu_obj, group1, group2, contrast_name) {
  cat(paste("Running pseudobulk DESeq2 for", contrast_name, "...\n"))
  
  tryCatch({
    # Subset cells for comparison
    cells_group1 <- WhichCells(seu_obj, expression = group == group1)
    cells_group2 <- WhichCells(seu_obj, expression = group == group2)
    
    if(length(cells_group1) == 0 || length(cells_group2) == 0) {
      cat(paste("Skipping DESeq2 for", contrast_name, "- one or both groups have zero cells.\n"))
      return(NULL)
    }
    
    # Subset Seurat object
    seu_subset <- seu_obj[, c(cells_group1, cells_group2)]
    seu_subset$comparison_group <- ifelse(colnames(seu_subset) %in% cells_group1, group1, group2)
    
    # Create SingleCellExperiment object
    sce <- as.SingleCellExperiment(seu_subset)
    
    # Aggregate by sample
    sce$sample_id <- paste(sce$sample_name, sce$comparison_group, sep = "_")
    
    # Create pseudobulk
    pb <- aggregateData(sce, assay = "counts", fun = "sum", by = "sample_id")
    
    # Prepare for DESeq2
    counts_matrix <- assay(pb)
    coldata <- data.frame(
      sample_id = colnames(counts_matrix),
      group = ifelse(grepl(group1, colnames(counts_matrix)), group1, group2),
      stringsAsFactors = FALSE
    )
    rownames(coldata) <- coldata$sample_id
    coldata$group <- factor(coldata$group, levels = c(group2, group1)) # Set reference level correctly
    
    # Run DESeq2
    dds <- DESeqDataSetFromMatrix(countData = counts_matrix, 
                                  colData = coldata, 
                                  design = ~ group)
    dds <- DESeq(dds, quiet = TRUE)
    
    results <- results(dds, contrast = c("group", group1, group2))
    results_df <- as.data.frame(results)
    results_df$gene <- rownames(results_df)
    results_df$contrast <- contrast_name
    results_df$method <- "pseudobulk_DESeq2"
    results_df$significant <- results_df$padj < 0.05 & !is.na(results_df$padj)
    results_df$effect_size <- abs(results_df$log2FoldChange)
    
    # Rename columns to match other methods
    colnames(results_df)[colnames(results_df) == "log2FoldChange"] <- "avg_log2FC"
    colnames(results_df)[colnames(results_df) == "pvalue"] <- "p_val"
    colnames(results_df)[colnames(results_df) == "padj"] <- "p_val_adj"
    
    return(results_df)
  }, error = function(e) {
    cat(paste("Error in DESeq2 for", contrast_name, ":", e$message, "\n"))
    return(NULL)
  })
}

# Function to run pseudobulk edgeR
run_pseudobulk_edger <- function(seu_obj, group1, group2, contrast_name) {
  cat(paste("Running pseudobulk edgeR for", contrast_name, "...\n"))
  
  tryCatch({
    # Subset cells for comparison
    cells_group1 <- WhichCells(seu_obj, expression = group == group1)
    cells_group2 <- WhichCells(seu_obj, expression = group == group2)
    
    if(length(cells_group1) == 0 || length(cells_group2) == 0) {
      cat(paste("Skipping edgeR for", contrast_name, "- one or both groups have zero cells.\n"))
      return(NULL)
    }
    
    # Subset Seurat object
    seu_subset <- seu_obj[, c(cells_group1, cells_group2)]
    seu_subset$comparison_group <- ifelse(colnames(seu_subset) %in% cells_group1, group1, group2)
    
    # Create SingleCellExperiment object
    sce <- as.SingleCellExperiment(seu_subset)
    
    # Aggregate by sample
    sce$sample_id <- paste(sce$sample_name, sce$comparison_group, sep = "_")
    
    # Create pseudobulk
    pb <- aggregateData(sce, assay = "counts", fun = "sum", by = "sample_id")
    
    # Prepare for edgeR
    counts_matrix <- assay(pb)
    group_info <- ifelse(grepl(group1, colnames(counts_matrix)), group1, group2)
    group_info <- factor(group_info, levels = c(group2, group1)) # Set reference level correctly
    
    # Run edgeR
    dge <- DGEList(counts = counts_matrix, group = group_info)
    dge <- calcNormFactors(dge)
    design <- model.matrix(~group_info)
    dge <- estimateDisp(dge, design)
    
    fit <- glmQLFit(dge, design)
    qlf <- glmQLFTest(fit, coef = 2)
    
    results_df <- topTags(qlf, n = Inf)$table
    results_df$gene <- rownames(results_df)
    results_df$contrast <- contrast_name
    results_df$method <- "pseudobulk_edgeR"
    results_df$significant <- results_df$FDR < 0.05
    results_df$effect_size <- abs(results_df$logFC)
    
    # Rename columns to match other methods
    colnames(results_df)[colnames(results_df) == "logFC"] <- "avg_log2FC"
    colnames(results_df)[colnames(results_df) == "PValue"] <- "p_val"
    colnames(results_df)[colnames(results_df) == "FDR"] <- "p_val_adj"
    
    return(results_df)
  }, error = function(e) {
    cat(paste("Error in edgeR for", contrast_name, ":", e$message, "\n"))
    return(NULL)
  })
}

# Run all DGE methods for all contrasts
cat("Running DGE analysis for all methods and contrasts...\n")
all_results <- list()

methods <- c("wilcoxon", "mast", "pseudobulk_deseq2", "pseudobulk_edger")

for(method in methods) {
  cat(paste("Running", method, "method...\n"))
  method_results <- list()
  
  for(contrast_name in names(contrasts_list)) {
    group1 <- contrasts_list[[contrast_name]][1]
    group2 <- contrasts_list[[contrast_name]][2]
    
    result <- switch(method,
                     "wilcoxon" = run_wilcoxon(seu_NKT_focused, group1, group2, contrast_name),
                     "mast" = run_mast(seu_NKT_focused, group1, group2, contrast_name),
                     "pseudobulk_deseq2" = run_pseudobulk_deseq2(seu_NKT_focused, group1, group2, contrast_name),
                     "pseudobulk_edger" = run_pseudobulk_edger(seu_NKT_focused, group1, group2, contrast_name))
    
    if(!is.null(result)) {
      method_results[[contrast_name]] <- result
    }
  }
  
  all_results[[method]] <- method_results
}

# Function to create comprehensive Excel output
create_excel_output <- function(results_list, method_name) {
  cat(paste("Creating Excel output for", method_name, "...\n"))
  
  wb <- createWorkbook()
  
  # Sheet 1: All contrasts combined
  all_contrasts_df <- do.call(rbind, results_list)
  addWorksheet(wb, "All_Contrasts_Combined")
  writeData(wb, "All_Contrasts_Combined", all_contrasts_df)
  
  # Sheets 2-16: Individual contrasts (all genes)
  for(contrast_name in names(results_list)) {
    sheet_name <- paste("All_Genes", contrast_name, sep = "_")
    sheet_name <- substr(sheet_name, 1, 31)  # Excel sheet name limit
    addWorksheet(wb, sheet_name)
    writeData(wb, sheet_name, results_list[[contrast_name]])
  }
  
  # Sheet: S1P genes across all contrasts
  s1p_results <- all_contrasts_df[all_contrasts_df$gene %in% s1p_genes, ]
  addWorksheet(wb, "S1P_Genes_All_Contrasts")
  writeData(wb, "S1P_Genes_All_Contrasts", s1p_results)
  
  # Sheet: Trafficking genes across all contrasts
  trafficking_results <- all_contrasts_df[all_contrasts_df$gene %in% trafficking_genes, ]
  addWorksheet(wb, "Trafficking_Genes_All_Contrasts")
  writeData(wb, "Trafficking_Genes_All_Contrasts", trafficking_results)
  
  # Sheet: Summary statistics
  summary_stats <- data.frame(
    Contrast = names(results_list),
    Total_Genes = sapply(results_list, nrow),
    Significant_Genes = sapply(results_list, function(x) sum(x$significant, na.rm = TRUE)),
    Upregulated = sapply(results_list, function(x) sum(x$avg_log2FC > 0 & x$significant, na.rm = TRUE)),
    Downregulated = sapply(results_list, function(x) sum(x$avg_log2FC < 0 & x$significant, na.rm = TRUE)),
    S1P_Genes_Detected = sapply(results_list, function(x) sum(x$gene %in% s1p_genes)),
    Trafficking_Genes_Detected = sapply(results_list, function(x) sum(x$gene %in% trafficking_genes))
  )
  addWorksheet(wb, "Summary_Stats_By_Contrast")
  writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
  
  # Sheet: Genes of interest by compartment
  wb_results <- all_contrasts_df[grepl("_WB", all_contrasts_df$contrast) & 
                                   all_contrasts_df$gene %in% all_genes_of_interest, ]
  bm_results <- all_contrasts_df[grepl("_BM", all_contrasts_df$contrast) & 
                                   all_contrasts_df$gene %in% all_genes_of_interest, ]
  
  addWorksheet(wb, "Genes_of_Interest_WB")
  writeData(wb, "Genes_of_Interest_WB", wb_results)
  
  addWorksheet(wb, "Genes_of_Interest_BM")
  writeData(wb, "Genes_of_Interest_BM", bm_results)
  
  # Sheet: Top upregulated genes by contrast
  top_up <- all_contrasts_df %>%
    filter(avg_log2FC > 0 & significant) %>%
    group_by(contrast) %>%
    slice_max(order_by = avg_log2FC, n = 50) %>%
    ungroup()
  addWorksheet(wb, "Top_Upregulated_By_Contrast")
  writeData(wb, "Top_Upregulated_By_Contrast", top_up)
  
  # Sheet: Top downregulated genes by contrast
  top_down <- all_contrasts_df %>%
    filter(avg_log2FC < 0 & significant) %>%
    group_by(contrast) %>%
    slice_min(order_by = avg_log2FC, n = 50) %>%
    ungroup()
  addWorksheet(wb, "Top_Downregulated_By_Contrast")
  writeData(wb, "Top_Downregulated_By_Contrast", top_down)
  
  # Save Excel file
  excel_file <- file.path(output_dir, "Excel_Files", paste0(method_name, "_DGE_Results.xlsx"))
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  
  return(excel_file)
}

# Create Excel outputs for all methods
cat("Creating Excel outputs...\n")
for(method in names(all_results)) {
  if(length(all_results[[method]]) > 0) {
    create_excel_output(all_results[[method]], method)
  }
}

# Visualization functions
create_heatmap <- function(seu_obj, genes, title, filename) {
  tryCatch({
    # Check which genes are present
    genes_present <- genes[genes %in% rownames(seu_obj)]
    
    if(length(genes_present) == 0) {
      cat(paste("No genes found for heatmap:", title, "\n"))
      return(NULL)
    }
    
    # Get expression data
    avg_exp <- AverageExpression(seu_obj, 
                                 features = genes_present, 
                                 group.by = "group")$RNA
    
    # Create heatmap
    p <- pheatmap(avg_exp, 
                  cluster_rows = TRUE, 
                  cluster_cols = TRUE,
                  scale = "row",
                  color = colorRampPalette(c("blue", "white", "red"))(100),
                  main = title,
                  filename = file.path(output_dir, "Heatmaps", filename),
                  width = 10, 
                  height = 8)
    
    return(p)
  }, error = function(e) {
    cat(paste("Error creating heatmap", title, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create heatmaps
cat("Creating heatmaps...\n")
create_heatmap(seu_NKT_focused, s1p_genes, "S1P Pathway Genes", "S1P_pathway_heatmap.png")
create_heatmap(seu_NKT_focused, trafficking_genes, "Trafficking Genes", "trafficking_genes_heatmap.png")
create_heatmap(seu_NKT_focused, all_genes_of_interest, "All Genes of Interest", "all_genes_of_interest_heatmap.png")

# Create volcano plots
create_volcano_plot <- function(results_df, contrast_name, method_name) {
  tryCatch({
    if(is.null(results_df) || nrow(results_df) == 0) {
      return(NULL)
    }
    
    # Prepare data for plotting
    plot_data <- results_df
    plot_data$log10_pval <- -log10(plot_data$p_val_adj + 1e-300)
    
    # Highlight genes of interest
    plot_data$highlight <- ifelse(plot_data$gene %in% all_genes_of_interest, 
                                  as.character(plot_data$gene), "")
    
    # Create volcano plot
    p <- ggplot(plot_data, aes(x = avg_log2FC, y = log10_pval)) +
      geom_point(aes(color = significant), alpha = 0.6) +
      geom_point(data = subset(plot_data, gene %in% all_genes_of_interest),
                 aes(x = avg_log2FC, y = log10_pval), 
                 color = "red", size = 2) +
      geom_text_repel(aes(label = highlight), 
                      size = 3, 
                      max.overlaps = 20,
                      box.padding = 0.5) +
      scale_color_manual(values = c("FALSE" = "grey", "TRUE" = "red"), name = "Significant (p_adj < 0.05)") +
      geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "red") +
      geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "blue") +
      labs(title = paste("Volcano Plot:", method_name, "-", contrast_name),
           x = "Average Log2 Fold Change",
           y = "-Log10 Adjusted P-value") +
      theme_minimal() +
      theme(legend.position = "bottom")
    
    filename <- file.path(output_dir, "Volcano_Plots", 
                          paste0(method_name, "_", contrast_name, "_volcano.png"))
    ggsave(filename, p, width = 12, height = 8, dpi = 300)
    
    return(p)
  }, error = function(e) {
    cat(paste("Error creating volcano plot for", contrast_name, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create volcano plots for all methods and contrasts
cat("Creating volcano plots...\n")
for(method in names(all_results)) {
  for(contrast_name in names(all_results[[method]])) {
    create_volcano_plot(all_results[[method]][[contrast_name]], contrast_name, method)
  }
}

# Create violin plots
create_violin_plots <- function(seu_obj, genes, title_prefix) {
  tryCatch({
    genes_present <- genes[genes %in% rownames(seu_obj)]
    
    if(length(genes_present) == 0) {
      return(NULL)
    }
    
    plots <- list()
    
    for(gene in genes_present) {
      p <- VlnPlot(seu_obj, 
                   features = gene, 
                   group.by = "group",
                   split.by = "tissue",
                   pt.size = 0.1) +
        labs(title = paste(title_prefix, gene)) +
        theme(axis.text.x = element_text(angle = 45, hjust = 1))
      
      plots[[gene]] <- p
    }
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = 2)
    
    filename <- file.path(output_dir, "Violin_Plots", 
                          paste0(gsub(" ", "_", title_prefix), "_violin_plots.png"))
    ggsave(filename, combined_plot, width = 16, height = 4 * ceiling(length(plots)/2), dpi = 300, limitsize = FALSE)
    
    return(combined_plot)
  }, error = function(e) {
    cat(paste("Error creating violin plots for", title_prefix, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create violin plots
cat("Creating violin plots...\n")
create_violin_plots(seu_NKT_focused, s1p_genes, "S1P Pathway")
create_violin_plots(seu_NKT_focused, trafficking_genes, "Trafficking Genes")

# Create nebulosa plots
create_nebulosa_plots <- function(seu_obj, genes, title_prefix) {
  tryCatch({
    genes_present <- genes[genes %in% rownames(seu_obj)]
    
    if(length(genes_present) == 0) {
      return(NULL)
    }
    
    plots <- list()
    
    for(gene in genes_present) {
      p <- plot_density(seu_obj, gene, slot = "data") +
        ggtitle(paste(title_prefix, gene)) +
        theme_minimal()
      
      plots[[gene]] <- p
    }
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = 2)
    
    filename <- file.path(output_dir, "Nebulosa_Plots", 
                          paste0(gsub(" ", "_", title_prefix), "_nebulosa_plots.png"))
    ggsave(filename, combined_plot, width = 16, height = 4 * ceiling(length(plots)/2), dpi = 300, limitsize = FALSE)
    
    return(combined_plot)
  }, error = function(e) {
    cat(paste("Error creating nebulosa plots for", title_prefix, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create nebulosa plots
cat("Creating nebulosa plots...\n")
create_nebulosa_plots(seu_NKT_focused, s1p_genes, "S1P Pathway")
create_nebulosa_plots(seu_NKT_focused, trafficking_genes, "Trafficking Genes")

# Trajectory analysis with Monocle3
cat("Running Monocle3 trajectory analysis...\n")
tryCatch({
  # Convert to cell_data_set
  cds <- as.cell_data_set(seu_NKT_focused)
  
  # Add gene and cell metadata
  cds <- estimate_size_factors(cds)
  cds <- preprocess_cds(cds, num_dim = 50)
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  # Cluster cells and learn trajectory
  cds <- cluster_cells(cds, resolution = 1e-5)
  cds <- learn_graph(cds)
  
  # Order cells in pseudotime
  cds <- order_cells(cds)
  
  # Create UMAP plots colored by pseudotime and treatment
  p1 <- plot_cells(cds, color_cells_by = "pseudotime", 
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 1.5) +
    ggtitle("UMAP colored by Pseudotime") +
    theme_minimal()
  
  p2 <- plot_cells(cds, color_cells_by = "treatment", 
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 1.5) +
    ggtitle("UMAP colored by Treatment") +
    theme_minimal()
  
  p3 <- plot_cells(cds, color_cells_by = "tissue", 
                   label_cell_groups = FALSE,
                   label_leaves = FALSE,
                   label_branch_points = FALSE,
                   graph_label_size = 1.5) +
    ggtitle("UMAP colored by Tissue") +
    theme_minimal()
  
  # Combine UMAP plots
  combined_umap <- (p1 | p2 | p3)
  
  # Save UMAP plots
  ggsave(file.path(output_dir, "Trajectory_Analysis", "UMAP_pseudotime_treatment.png"), 
         combined_umap, width = 18, height = 6, dpi = 300)
  
  # Create pseudotime plots for S1P genes
  s1p_genes_present <- s1p_genes[s1p_genes %in% rownames(cds)]
  
  if(length(s1p_genes_present) > 0) {
    pseudotime_plots <- list()
    
    for(gene in s1p_genes_present) {
      p <- plot_genes_in_pseudotime(cds[gene,], 
                                    color_cells_by = "treatment",
                                    min_expr = 0.1) +
        ggtitle(paste("S1P Gene:", gene, "- Pseudotime Expression")) +
        theme_minimal()
      
      pseudotime_plots[[gene]] <- p
    }
    
    # Combine S1P pseudotime plots
    combined_s1p_pseudotime <- wrap_plots(pseudotime_plots, ncol = 2)
    
    ggsave(file.path(output_dir, "Trajectory_Analysis", "S1P_genes_pseudotime.png"), 
           combined_s1p_pseudotime, width = 16, height = 4 * ceiling(length(pseudotime_plots)/2), dpi = 300, limitsize = FALSE)
  }
  
  # Create pseudotime plots for trafficking genes
  trafficking_genes_present <- trafficking_genes[trafficking_genes %in% rownames(cds)]
  
  if(length(trafficking_genes_present) > 0) {
    trafficking_pseudotime_plots <- list()
    
    for(gene in trafficking_genes_present) {
      p <- plot_genes_in_pseudotime(cds[gene,], 
                                    color_cells_by = "treatment",
                                    min_expr = 0.1) +
        ggtitle(paste("Trafficking Gene:", gene, "- Pseudotime Expression")) +
        theme_minimal()
      
      trafficking_pseudotime_plots[[gene]] <- p
    }
    
    # Combine trafficking pseudotime plots
    combined_trafficking_pseudotime <- wrap_plots(trafficking_pseudotime_plots, ncol = 2)
    
    ggsave(file.path(output_dir, "Trajectory_Analysis", "Trafficking_genes_pseudotime.png"), 
           combined_trafficking_pseudotime, width = 16, height = 4 * ceiling(length(trafficking_pseudotime_plots)/2), dpi = 300, limitsize = FALSE)
  }
  
  cat("Trajectory analysis completed successfully!\n")
  
}, error = function(e) {
  cat(paste("Error in trajectory analysis:", e$message, "\n"))
})

# Create correlation plots between S1P and trafficking genes
cat("Creating correlation plots...\n")
create_correlation_plots <- function(seu_obj, genes1, genes2, title, filename) {
  tryCatch({
    genes1_present <- genes1[genes1 %in% rownames(seu_obj)]
    genes2_present <- genes2[genes2 %in% rownames(seu_obj)]
    
    if(length(genes1_present) == 0 || length(genes2_present) == 0) {
      cat(paste("Insufficient genes for correlation plot:", title, "\n"))
      return(NULL)
    }
    
    # Get expression data
    exp_data <- GetAssayData(seu_obj, slot = "data")
    exp_matrix <- as.matrix(exp_data[c(genes1_present, genes2_present), ])
    
    # Calculate correlation
    cor_matrix <- cor(t(exp_matrix), method = "pearson")
    
    # Create correlation plot
    png(file.path(output_dir, "Plots", filename), width = 800, height = 800)
    corrplot(cor_matrix, 
             method = "color", 
             type = "upper", 
             order = "hclust", 
             tl.cex = 0.8,
             tl.col = "black",
             title = title,
             mar = c(0,0,2,0))
    dev.off()
    
    return(cor_matrix)
  }, error = function(e) {
    cat(paste("Error creating correlation plot", title, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create correlation plots
create_correlation_plots(seu_NKT_focused, s1p_genes, trafficking_genes, 
                         "S1P vs Trafficking Genes Correlation", "S1P_trafficking_correlation.png")

# Create ridge plots for key genes
cat("Creating ridge plots...\n")
create_ridge_plots <- function(seu_obj, genes, title_prefix) {
  tryCatch({
    genes_present <- genes[genes %in% rownames(seu_obj)]
    
    if(length(genes_present) == 0) {
      return(NULL)
    }
    
    plots <- list()
    
    for(gene in genes_present) {
      # Get expression data
      exp_data <- FetchData(seu_obj, vars = c(gene, "group"))
      colnames(exp_data) <- c("expression", "group")
      
      p <- ggplot(exp_data, aes(x = expression, y = group, fill = group)) +
        geom_density_ridges(alpha = 0.7) +
        labs(title = paste(title_prefix, gene),
             x = "Expression Level",
             y = "Group") +
        theme_minimal() +
        theme(legend.position = "none")
      
      plots[[gene]] <- p
    }
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = 2)
    
    filename <- file.path(output_dir, "Plots", 
                          paste0(gsub(" ", "_", title_prefix), "_ridge_plots.png"))
    ggsave(filename, combined_plot, width = 16, height = 4 * ceiling(length(plots)/2), dpi = 300, limitsize = FALSE)
    
    return(combined_plot)
  }, error = function(e) {
    cat(paste("Error creating ridge plots for", title_prefix, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create ridge plots
create_ridge_plots(seu_NKT_focused, s1p_genes, "S1P Pathway")
create_ridge_plots(seu_NKT_focused, trafficking_genes, "Trafficking Genes")

# Create feature plots
cat("Creating feature plots...\n")
create_feature_plots <- function(seu_obj, genes, title_prefix) {
  tryCatch({
    genes_present <- genes[genes %in% rownames(seu_obj)]
    
    if(length(genes_present) == 0) {
      return(NULL)
    }
    
    plots <- list()
    
    for(gene in genes_present) {
      p <- FeaturePlot(seu_obj, features = gene, 
                       split.by = "tissue",
                       cols = c("lightgrey", "red")) +
        ggtitle(paste(title_prefix, gene))
      
      plots[[gene]] <- p
    }
    
    # Combine plots
    combined_plot <- wrap_plots(plots, ncol = 1) # ncol=1 is better for split.by plots
    
    filename <- file.path(output_dir, "Plots", 
                          paste0(gsub(" ", "_", title_prefix), "_feature_plots.png"))
    ggsave(filename, combined_plot, width = 12, height = 6 * length(plots), dpi = 300, limitsize = FALSE)
    
    return(combined_plot)
  }, error = function(e) {
    cat(paste("Error creating feature plots for", title_prefix, ":", e$message, "\n"))
    return(NULL)
  })
}

# Create feature plots
create_feature_plots(seu_NKT_focused, s1p_genes, "S1P Pathway")
create_feature_plots(seu_NKT_focused, trafficking_genes, "Trafficking Genes")

# Create proportion plots
cat("Creating proportion plots...\n")
create_proportion_plots <- function(seu_obj) {
  tryCatch({
    # Check if 'cell_types' column exists in metadata
    if (!"cell_types" %in% colnames(seu_obj@meta.data)) {
      cat("Error: 'cell_types' column not found in Seurat object metadata. Skipping proportion plots.\n")
      return(NULL)
    }
    
    # Calculate proportions
    prop_data <- seu_obj@meta.data %>%
      group_by(treatment, tissue, cell_types) %>%
      summarise(count = n(), .groups = "drop") %>%
      group_by(treatment, tissue) %>%
      mutate(proportion = count / sum(count))
    
    # Create proportion plot
    p <- ggplot(prop_data, aes(x = treatment, y = proportion, fill = cell_types)) +
      geom_bar(stat = "identity", position = "stack") +
      facet_wrap(~tissue) +
      labs(title = "Cell Type Proportions by Treatment and Tissue",
           x = "Treatment",
           y = "Proportion",
           fill = "Cell Type") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave(file.path(output_dir, "Plots", "cell_type_proportions.png"), 
           p, width = 12, height = 8, dpi = 300)
    
    return(p)
  }, error = function(e) {
    cat(paste("Error creating proportion plots:", e$message, "\n"))
    return(NULL)
  })
}

# Create proportion plots
create_proportion_plots(seu_NKT_focused)

# GSEA analysis
cat("Running GSEA analysis...\n")
run_gsea_analysis <- function(results_df, contrast_name, method_name) {
  tryCatch({
    if(is.null(results_df) || nrow(results_df) == 0) {
      return(NULL)
    }
    
    # Prepare gene list for GSEA
    gene_list <- results_df$avg_log2FC
    names(gene_list) <- results_df$gene
    gene_list <- gene_list[!is.na(gene_list) & !is.na(names(gene_list))]
    gene_list <- gene_list[!duplicated(names(gene_list))]
    gene_list <- sort(gene_list, decreasing = TRUE)
    
    # Get mouse gene sets
    m_df <- msigdbr(species = "Mus musculus", category = "H") # Hallmark gene sets
    m_list <- split(x = m_df$gene_symbol, f = m_df$gs_name)
    
    # Run GSEA
    gsea_results <- fgsea(pathways = m_list, 
                          stats = gene_list,
                          minSize = 15,
                          maxSize = 500)
    
    # Filter significant results
    gsea_results_sig <- gsea_results[gsea_results$padj < 0.05, ]
    
    if(nrow(gsea_results_sig) > 0) {
      # Create enrichment plot for top pathway
      top_pathway <- gsea_results_sig[order(gsea_results_sig$padj)[1], ]$pathway
      
      p <- plotEnrichment(m_list[[top_pathway]], gene_list) +
        ggtitle(paste("GSEA:", method_name, "-", contrast_name, "\n", top_pathway)) +
        theme_minimal()
      
      filename <- file.path(output_dir, "Plots", 
                            paste0("GSEA_", method_name, "_", contrast_name, ".png"))
      ggsave(filename, p, width = 10, height = 6, dpi = 300)
    }
    
    return(gsea_results_sig)
  }, error = function(e) {
    cat(paste("Error in GSEA analysis for", contrast_name, ":", e$message, "\n"))
    return(NULL)
  })
}

# Run GSEA for selected contrasts
gsea_results_list <- list()
selected_contrasts <- c("naive_WB_vs_UT_WB", "naive_BM_vs_UT_BM", "UT_WB_vs_UT_BM")

for(method in names(all_results)) {
  for(contrast in selected_contrasts) {
    if(contrast %in% names(all_results[[method]])) {
      gsea_result <- run_gsea_analysis(all_results[[method]][[contrast]], contrast, method)
      if(!is.null(gsea_result)) {
        gsea_results_list[[paste(method, contrast, sep = "_")]] <- gsea_result
      }
    }
  }
}

# Statistical power analysis
cat("Running statistical power analysis...\n")
perform_power_analysis <- function(seu_obj) {
  tryCatch({
    # Calculate power for detecting fold changes
    power_results <- data.frame(
      Contrast = character(),
      Group1_Cells = numeric(),
      Group2_Cells = numeric(),
      Estimated_Power = numeric(),
      stringsAsFactors = FALSE
    )
    
    for(contrast_name in names(contrasts_list)) {
      group1 <- contrasts_list[[contrast_name]][1]
      group2 <- contrasts_list[[contrast_name]][2]
      
      cells_group1 <- length(WhichCells(seu_obj, expression = group == group1))
      cells_group2 <- length(WhichCells(seu_obj, expression = group == group2))
      
      # Simple power estimation based on cell numbers
      min_cells <- min(cells_group1, cells_group2)
      estimated_power <- pmin(1, min_cells / 100)  # Rough estimate
      
      power_results <- rbind(power_results, data.frame(
        Contrast = contrast_name,
        Group1_Cells = cells_group1,
        Group2_Cells = cells_group2,
        Estimated_Power = estimated_power
      ))
    }
    
    # Save power analysis results
    write.csv(power_results, 
              file.path(output_dir, "power_analysis_results.csv"), 
              row.names = FALSE)
    
    return(power_results)
  }, error = function(e) {
    cat(paste("Error in power analysis:", e$message, "\n"))
    return(NULL)
  })
}

# Run power analysis
power_results <- perform_power_analysis(seu_NKT_focused)

# Create summary report
cat("Creating summary report...\n")
create_summary_report <- function() {
  tryCatch({
    # Check if 'cell_types' column exists
    cell_types_info <- if ("cell_types" %in% colnames(seu_NKT_focused@meta.data)) {
      paste(unique(seu_NKT_focused$cell_types), collapse = ", ")
    } else {
      "Not available in metadata"
    }
    
    report_lines <- c(
      "# S1P Receptor DGE Analysis Summary Report",
      paste("Analysis completed on:", Sys.time()),
      "",
      "## Dataset Information",
      paste("Total cells analyzed:", ncol(seu_NKT_focused)),
      paste("Total genes:", nrow(seu_NKT_focused)),
      paste("Cell types:", cell_types_info),
      paste("Treatments:", paste(levels(seu_NKT_focused$treatment), collapse = ", ")),
      paste("Tissues:", paste(levels(seu_NKT_focused$tissue), collapse = ", ")),
      "",
      "## Analysis Methods",
      "- Wilcoxon rank-sum test",
      "- MAST (Model-based Analysis of Single-cell Transcriptomics)",
      "- Pseudobulk DESeq2",
      "- Pseudobulk edgeR",
      "",
      "## Key Findings",
      paste("Total contrasts analyzed:", length(contrasts_list)),
      paste("S1P pathway genes analyzed:", paste(s1p_genes, collapse = ", ")),
      paste("Trafficking genes analyzed:", paste(trafficking_genes, collapse = ", ")),
      "",
      "## Output Files Generated",
      "- Excel files with comprehensive results for each method",
      "- Volcano plots for all contrasts",
      "- Heatmaps for gene expression patterns",
      "- Violin plots showing expression distributions",
      "- Nebulosa density plots",
      "- Trajectory analysis with Monocle3",
      "- Correlation plots between S1P and trafficking genes",
      "- GSEA enrichment analysis",
      "- Statistical power analysis",
      "",
      "## Files Located In:",
      paste("Output directory:", output_dir),
      "- Excel_Files/: Comprehensive DGE results",
      "- Plots/: Various visualization outputs",
      "- Heatmaps/: Expression heatmaps",
      "- Volcano_Plots/: Volcano plots for all contrasts",
      "- Violin_Plots/: Expression distribution plots",
      "- Nebulosa_Plots/: Density plots",
      "- Trajectory_Analysis/: Monocle3 pseudotime analysis"
    )
    
    writeLines(report_lines, file.path(output_dir, "Analysis_Summary_Report.txt"))
    
    return(report_lines)
  }, error = function(e) {
    cat(paste("Error creating summary report:", e$message, "\n"))
    return(NULL)
  })
}

# Create summary report
summary_report <- create_summary_report()

# Final cleanup and session info
cat("Saving session information...\n")
sessionInfo_output <- capture.output(sessionInfo())
writeLines(sessionInfo_output, file.path(output_dir, "sessionInfo.txt"))

cat("Analysis pipeline completed successfully!\n")
cat(paste("All results saved to:", output_dir, "\n"))
cat("Check the Analysis_Summary_Report.txt for a complete overview of outputs.\n")

# Print completion message
cat(paste0("\n", paste(rep("=", 50), collapse = ""), "\n"))
cat("S1P RECEPTOR DGE ANALYSIS PIPELINE COMPLETED\n")
cat(paste0(paste(rep("=", 50), collapse = ""), "\n"))
cat("Research Question: Is the loss of S1P receptor externalization on T-cells\n")
cat("and impaired lymphocyte trafficking from bone marrow a direct or indirect\n")
cat("consequence of DMG presence?\n")
cat("\n")
cat("Key outputs generated:\n")
cat("1. Comprehensive Excel files with DGE results for all 4 methods\n")
cat("2. 15 pairwise comparisons across all treatment-tissue combinations\n")
cat("3. Extensive visualizations including volcano plots, heatmaps, and trajectory analysis\n")
cat("4. Statistical power analysis and GSEA enrichment results\n")
cat("5. Focused analysis on S1P pathway and trafficking genes\n")
cat("\n")
cat(paste("Total analysis time:", round(difftime(Sys.time(), start_time, units = "mins"), 2), "minutes\n"))
cat(paste("Results directory:", output_dir, "\n"))