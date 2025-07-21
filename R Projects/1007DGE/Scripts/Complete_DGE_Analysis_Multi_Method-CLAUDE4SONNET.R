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

print("=== COMPREHENSIVE DGE ANALYSIS ===")
print("Data structure confirmed:")
print(table(seu_NKT_focused$sample_id, seu_NKT_focused$treatment))

# Set default assay
DefaultAssay(seu_NKT_focused) <- "RNA"

# Create results directory
if(!dir.exists(here("R Projects", "1007DGE", "claude4results"))) {
  dir.create(here("R Projects", "1007DGE", "claude4results"))
}

# ============================================================================
# 1. SEURAT WILCOXON TEST (Standard approach)
# ============================================================================

print("=== 1. SEURAT WILCOXON ANALYSIS ===")

wilcox_results <- list()

for(cell_type in unique(seu_NKT_focused$cell_types)) {
  print(paste("Processing", cell_type))
  
  # Subset to specific cell type
  seu_subset <- subset(seu_NKT_focused, subset = cell_types == cell_type)
  
  # Check cell counts
  treatment_counts <- table(seu_subset$treatment)
  print(treatment_counts)
  
  # Set identity to treatment
  Idents(seu_subset) <- "treatment"
  
  # All possible pairwise comparisons
  treatments <- names(treatment_counts)[treatment_counts >= 20]  # Minimum 20 cells
  
  if(length(treatments) >= 2) {
    comparisons <- list(
      c("NAIVE", "UT"),
      c("SHAM", "UT"),
      c("SHAM", "NAIVE")
    )
    
    for(comp in comparisons) {
      if(all(comp %in% treatments)) {
        comparison_name <- paste(comp[1], "vs", comp[2])
        print(paste("  Comparing", comparison_name))
        
        tryCatch({
          markers <- FindMarkers(seu_subset,
                               ident.1 = comp[1],
                               ident.2 = comp[2],
                               min.pct = 0.1,
                               logfc.threshold = 0.1,
                               test.use = "wilcox",
                               verbose = FALSE)
          
          if(nrow(markers) > 0) {
            markers$gene <- rownames(markers)
            markers$comparison <- comparison_name
            markers$cell_type <- cell_type
            markers$method <- "Wilcoxon"
            
            result_name <- paste(cell_type, comparison_name, "wilcox", sep = "_")
            wilcox_results[[result_name]] <- markers
            
            print(paste("    Found", nrow(markers), "differential genes"))
          }
        }, error = function(e) {
          print(paste("    Error:", e$message))
        })
      }
    }
  }
}

# ============================================================================
# 2. PSEUDOBULK ANALYSIS WITH MUSCAT
# ============================================================================

print("=== 2. MUSCAT PSEUDOBULK ANALYSIS ===")

# Convert to SingleCellExperiment
sce <- as.SingleCellExperiment(seu_NKT_focused)

# Add required metadata for muscat
colData(sce)$cluster_id <- sce$cell_types
colData(sce)$sample_id <- sce$sample_id
colData(sce)$group_id <- sce$treatment

# Prepare SCE for muscat
sce <- prepSCE(sce, 
               kid = "cluster_id",
               sid = "sample_id", 
               gid = "group_id",
               drop = TRUE)

print("Sample summary for pseudobulk:")
print(metadata(sce)$experiment_info)

# Create pseudobulks
pb <- aggregateData(sce, assay = "counts", fun = "sum")

# Run differential expression
muscat_results <- pbDS(pb, 
                      method = c("edgeR", "DESeq2"),
                      filter = "both",
                      min_cells = 20)

# Extract results
muscat_edgeR <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cpm = FALSE)
muscat_DESeq2 <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cpm = FALSE, method = "DESeq2")

# Add method column
if(nrow(muscat_edgeR) > 0) {
  muscat_edgeR$method <- "muscat_edgeR"
}
if(nrow(muscat_DESeq2) > 0) {
  muscat_DESeq2$method <- "muscat_DESeq2"
}

# ============================================================================
# 3. MANUAL PSEUDOBULK WITH DESEQ2
# ============================================================================

print("=== 3. MANUAL PSEUDOBULK DESeq2 ===")

# Function to create pseudobulk counts
create_pseudobulk <- function(seu_obj, group_by = c("sample_id", "cell_types")) {
  counts <- GetAssayData(seu_obj, slot = "counts")
  meta <- seu_obj@meta.data
  group_factor <- interaction(meta[, group_by], drop = TRUE)
  
  pb_counts <- sapply(levels(group_factor), function(x) {
    cells <- names(group_factor)[group_factor == x]
    if(length(cells) > 1) {
      Matrix::rowSums(counts[, cells])
    } else if(length(cells) == 1) {
      counts[, cells]
    } else {
      rep(0, nrow(counts))
    }
  })
  
  return(pb_counts)
}

deseq2_results <- list()

for(cell_type in unique(seu_NKT_focused$cell_types)) {
  print(paste("Manual DESeq2 for", cell_type))
  
  # Subset to cell type
  seu_subset <- subset(seu_NKT_focused, subset = cell_types == cell_type)
  
  # Create pseudobulk
  pb_counts <- create_pseudobulk(seu_subset)
  
  # Sample metadata
  sample_info <- seu_subset@meta.data %>%
    group_by(sample_id, treatment) %>%
    summarise(n_cells = n(), .groups = "drop") %>%
    filter(n_cells >= 20) %>%
    as.data.frame()
  
  if(nrow(sample_info) >= 3 && length(unique(sample_info$treatment)) >= 2) {
    # Filter counts
    pb_counts_filtered <- pb_counts[, sample_info$sample_id]
    keep_genes <- rowSums(pb_counts_filtered > 1) >= 2
    pb_counts_filtered <- pb_counts_filtered[keep_genes, ]
    
    # Create DESeq2 object
    rownames(sample_info) <- sample_info$sample_id
    sample_info$treatment <- factor(sample_info$treatment, levels = c("UT", "NAIVE", "SHAM"))
    
    dds <- DESeqDataSetFromMatrix(countData = pb_counts_filtered,
                                 colData = sample_info,
                                 design = ~ treatment)
    
    # Run DESeq2
    dds <- DESeq(dds, quiet = TRUE)
    
    # Get results for each comparison
    comparisons <- list(
      c("treatment", "NAIVE", "UT"),
      c("treatment", "SHAM", "UT"),
      c("treatment", "SHAM", "NAIVE")
    )
    
    for(comp in comparisons) {
      comp_name <- paste(comp[2], "vs", comp[3])
      
      tryCatch({
        res <- results(dds, contrast = comp)
        res_df <- as.data.frame(res)
        res_df$gene <- rownames(res_df)
        res_df$cell_type <- cell_type
        res_df$comparison <- comp_name
        res_df$method <- "manual_DESeq2"
        
        result_name <- paste(cell_type, comp_name, "DESeq2", sep = "_")
        deseq2_results[[result_name]] <- res_df
        
        print(paste("  ", comp_name, ":", sum(res_df$padj < 0.05, na.rm = TRUE), "significant genes"))
      }, error = function(e) {
        print(paste("  Error in", comp_name, ":", e$message))
      })
    }
  }
}

# ============================================================================
# 4. ANALYZE GENES OF INTEREST ACROSS ALL METHODS
# ============================================================================

print("=== 4. GENES OF INTEREST ANALYSIS ===")

# Your genes of interest
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")

all_goi <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

# Function to extract GOI from results
extract_goi <- function(results_list, gene_col = "gene") {
  goi_results <- list()
  
  for(result_name in names(results_list)) {
    result_df <- results_list[[result_name]]
    
    if(gene_col %in% colnames(result_df)) {
      goi_subset <- result_df[result_df[[gene_col]] %in% all_goi, ]
    } else {
      goi_subset <- result_df[rownames(result_df) %in% all_goi, ]
      goi_subset[[gene_col]] <- rownames(goi_subset)
    }
    
    if(nrow(goi_subset) > 0) {
      goi_results[[result_name]] <- goi_subset
    }
  }
  
  return(goi_results)
}

# Extract GOI from all methods
goi_wilcox <- extract_goi(wilcox_results)
goi_deseq2 <- extract_goi(deseq2_results)

# For muscat results
goi_muscat_edgeR <- list()
goi_muscat_DESeq2 <- list()

if(nrow(muscat_edgeR) > 0) {
  muscat_edgeR_goi <- muscat_edgeR[muscat_edgeR$gene %in% all_goi, ]
  if(nrow(muscat_edgeR_goi) > 0) {
    goi_muscat_edgeR[["muscat_edgeR"]] <- muscat_edgeR_goi
  }
}

if(nrow(muscat_DESeq2) > 0) {
  muscat_DESeq2_goi <- muscat_DESeq2[muscat_DESeq2$gene %in% all_goi, ]
  if(nrow(muscat_DESeq2_goi) > 0) {
    goi_muscat_DESeq2[["muscat_DESeq2"]] <- muscat_DESeq2_goi
  }
}

# Combine all GOI results
all_goi_results <- c(goi_wilcox, goi_deseq2, goi_muscat_edgeR, goi_muscat_DESeq2)

# Create comprehensive GOI summary
goi_summary <- data.frame()
for(result_name in names(all_goi_results)) {
  df <- all_goi_results[[result_name]]
  goi_summary <- rbind(goi_summary, df)
}

print("Genes of interest found across all methods:")
if(nrow(goi_summary) > 0) {
  print(table(goi_summary$gene, goi_summary$method))
}

# ============================================================================
# 5. COMPREHENSIVE VISUALIZATION
# ============================================================================

print("=== 5. CREATING COMPREHENSIVE VISUALIZATIONS ===")

# Method comparison for genes of interest
if(nrow(goi_summary) > 0) {
  # Standardize column names for comparison
  goi_summary_std <- goi_summary %>%
    mutate(
      log2FC = case_when(
        !is.na(avg_log2FC) ~ avg_log2FC,
        !is.na(logFC) ~ logFC,
        !is.na(log2FoldChange) ~ log2FoldChange,
        TRUE ~ NA_real_
      ),
      pvalue = case_when(
        !is.na(p_val_adj) ~ p_val_adj,
        !is.na(FDR) ~ FDR,
        !is.na(padj) ~ padj,
        TRUE ~ NA_real_
      )
    ) %>%
    filter(!is.na(log2FC), !is.na(pvalue)) %>%
    mutate(
      significant = pvalue < 0.05,
      comparison_method = paste(comparison, method, sep = "_")
    )
  
  # Create comprehensive heatmap
  if(nrow(goi_summary_std) > 0) {
    p_comprehensive <- ggplot(goi_summary_std, 
                             aes(x = comparison_method, y = gene)) +
      geom_tile(aes(fill = log2FC), color = "white") +
      geom_text(aes(label = ifelse(significant, "*", "")), size = 4) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                          midpoint = 0, name = "Log2FC") +
      facet_wrap(~cell_type, scales = "free_x") +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            plot.title = element_text(hjust = 0.5)) +
      labs(title = "Genes of Interest - All Methods Comparison\n(* indicates p < 0.05)",
           x = "Comparison_Method",
           y = "Gene")
    
    print(p_comprehensive)
    ggsave(here("Results", "comprehensive_GOI_heatmap.pdf"), 
           p_comprehensive, width = 20, height = 12)
  }
}

# Individual gene expression plots
for(gene in all_goi[all_goi %in% rownames(seu_NKT_focused)]) {
  print(paste("Creating plots for", gene))
  
  # Violin plot by treatment and cell type
  p_violin <- VlnPlot(seu_NKT_focused, 
                     features = gene, 
                     group.by = "treatment", 
                     split.by = "cell_types",
                     pt.size = 0) +
    ggtitle(paste("Expression of", gene))
  
  # Box plot by sample
  plot_data <- FetchData(seu_NKT_focused, 
                        vars = c(gene, "treatment", "cell_types", "sample_id"))
  
  p_box <- ggplot(plot_data, aes_string(x = "treatment", y = gene, fill = "treatment")) +
    geom_boxplot(outlier.size = 0.5) +
    facet_wrap(~cell_types, scales = "free_y") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    ggtitle(paste("Expression of", gene, "by Sample"))
  
  # Combine plots
  combined_plot <- p_violin / p_box
  
  ggsave(here("Results", paste0(gene, "_expression_plots.pdf")), 
         combined_plot, width = 14, height = 10)
}

# ============================================================================
# 6. SAVE ALL RESULTS
# ============================================================================

print("=== 6. SAVING RESULTS ===")

# Save individual method results
saveRDS(wilcox_results, here("Results", "wilcoxon_DGE_results.rds"))
saveRDS(deseq2_results, here("Results", "manual_DESeq2_results.rds"))
saveRDS(muscat_results, here("Results", "muscat_DGE_results.rds"))

# Save combined results
write.csv(muscat_edgeR, here("Results", "muscat_edgeR_results.csv"), row.names = FALSE)
write.csv(muscat_DESeq2, here("Results", "muscat_DESeq2_results.csv"), row.names = FALSE)

# Save GOI results
saveRDS(all_goi_results, here("Results", "all_GOI_results.rds"))
if(nrow(goi_summary) > 0) {
  write.csv(goi_summary, here("Results", "comprehensive_GOI_summary.csv"), row.names = FALSE)
}

# Create method comparison summary
method_summary <- data.frame()
for(method_name in c("Wilcoxon", "manual_DESeq2", "muscat_edgeR", "muscat_DESeq2")) {
  if(method_name == "Wilcoxon") {
    n_results <- length(wilcox_results)
    n_genes <- sum(sapply(wilcox_results, nrow))
  } else if(method_name == "manual_DESeq2") {
    n_results <- length(deseq2_results)
    n_genes <- sum(sapply(deseq2_results, nrow))
  } else if(method_name == "muscat_edgeR") {
    n_results <- ifelse(nrow(muscat_edgeR) > 0, 1, 0)
    n_genes <- nrow(muscat_edgeR)
  } else if(method_name == "muscat_DESeq2") {
    n_results <- ifelse(nrow(muscat_DESeq2) > 0, 1, 0)
    n_genes <- nrow(muscat_DESeq2)
  }
  
  method_summary <- rbind(method_summary, 
                         data.frame(Method = method_name,
                                   N_Comparisons = n_results,
                                   Total_Gene_Results = n_genes))
}

write.csv(method_summary, here("Results", "method_comparison_summary.csv"), row.names = FALSE)

print("=== ANALYSIS COMPLETE ===")
print("Generated files:")
print("- comprehensive_GOI_heatmap.pdf (comparison across all methods)")
print("- Individual gene expression plots for each GOI")
print("- wilcoxon_DGE_results.rds")
print("- manual_DESeq2_results.rds")
print("- muscat_DGE_results.rds")
print("- muscat_edgeR_results.csv")
print("- muscat_DESeq2_results.csv")
print("- comprehensive_GOI_summary.csv")
print("- method_comparison_summary.csv")

print("\nMethod Summary:")
print(method_summary)

if(nrow(goi_summary) > 0) {
  print("\nSignificant GOI results summary:")
  sig_summary <- goi_summary_std %>%
    filter(significant) %>%
    group_by(gene, method) %>%
    summarise(n_significant = n(), .groups = "drop")
  print(sig_summary)
}