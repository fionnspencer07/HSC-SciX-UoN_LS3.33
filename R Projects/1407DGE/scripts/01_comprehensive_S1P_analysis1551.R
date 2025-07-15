# ==============================================================================
# COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS PIPELINE
# Date: 2025-07-14
# User: fionnspencer07
# ==============================================================================
library(here)
# Load all required libraries
source(here("Custom R Functions and Scripts", "init_packages.R"))
# Set seed for reproducibility
set.seed(42)

# Create comprehensive output directory structure
base_output <- here("R Projects", "1407DGE", "output", "Run 1")
dir.create(base_output, recursive = TRUE, showWarnings = FALSE)

# Create subdirectories for organization
subdirs <- c("Excel_Outputs", "SVG_Plots", "Heatmaps", "Trajectory_Analysis", 
             "Volcano_Plots", "Feature_Plots", "Density_Plots", "Violin_Plots",
             "Ridge_Plots", "Correlation_Plots", "Summary_Stats", "Pathway_Analysis")

for(subdir in subdirs) {
  dir.create(file.path(base_output, subdir), recursive = TRUE, showWarnings = FALSE)
}

# Define comprehensive color palette
man_cols <- c("#86b0cc", "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252",
              "#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", 
                     "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

condition_colors <- c("NAIVE" = "#2E86AB", "SHAM" = "#A23B72", "DMG" = "#F18F01")
cell_type_colors <- c("T_cells" = "#eeb84c", "NK_cells" = "#86b0cc", "NKT_cells" = "#e9944b")

# ==============================================================================
# 1. DATA LOADING AND COMPREHENSIVE PREPROCESSING
# ==============================================================================

cat("=== LOADING AND PREPROCESSING DATA ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")

# Load the focused dataset
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.RDS"))

# Comprehensive metadata setup
seu_NKT_focused$condition <- seu_NKT_focused$treatment
seu_NKT_focused$condition[seu_NKT_focused$condition == "UT"] <- "DMG"
seu_NKT_focused$condition <- factor(seu_NKT_focused$condition, levels = c("NAIVE", "SHAM", "DMG"))

# Add compartment information (adjust based on your actual metadata)
if(!"compartment" %in% colnames(seu_NKT_focused@meta.data)) {
  seu_NKT_focused$compartment <- "WB"  # Default to whole blood
}
seu_NKT_focused$compartment <- factor(seu_NKT_focused$compartment)

# Create comprehensive sample metadata
seu_NKT_focused$sample_condition <- paste0(seu_NKT_focused$sample_name, "_", seu_NKT_focused$condition)
seu_NKT_focused$cell_condition <- paste0(seu_NKT_focused$cell_types, "_", seu_NKT_focused$condition)

# Define comprehensive gene sets
gene_sets <- list(
  s1p_pathway = c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1"),
  trafficking = c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1"),
  bone_marrow_retention = c("Cxcr4", "Vcam1", "Vla4", "Sdf1", "Cxcl12"),
  t_cell_activation = c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2", "Cd137", "Cd154"),
  chemokines = c("Ccl5", "Ccl4", "Cxcl10", "Cxcl9", "Ccl3", "Ccl2", "Cxcl1"),
  cytokines = c("Ifng", "Il2", "Tnfa", "Il10", "Il4", "Il17a", "Il6"),
  adhesion_molecules = c("Itgae", "Itgal", "Itgam", "Itgax", "Itgb1", "Itgb2", "Itgb7"),
  transcription_factors = c("Tbx21", "Gata3", "Rorc", "Foxp3", "Bcl6", "Prdm1"),
  apoptosis = c("Bcl2", "Bax", "Fas", "Fasl", "Casp3", "Casp8", "Bid"),
  metabolism = c("Hif1a", "Pfkfb3", "Ldha", "Glut1", "Cpt1a", "Fasn")
)

# Get all genes of interest
all_genes_of_interest <- unique(unlist(gene_sets))

# Check gene availability
available_genes <- intersect(all_genes_of_interest, rownames(seu_NKT_focused))
missing_genes <- setdiff(all_genes_of_interest, rownames(seu_NKT_focused))

cat("Total genes of interest:", length(all_genes_of_interest), "\n")
cat("Available genes:", length(available_genes), "\n")
cat("Missing genes:", length(missing_genes), "\n")

# Save gene availability summary
gene_availability <- data.frame(
  gene_set = rep(names(gene_sets), sapply(gene_sets, length)),
  gene = unlist(gene_sets),
  available = unlist(gene_sets) %in% rownames(seu_NKT_focused)
)

write.xlsx(gene_availability, file.path(base_output, "Excel_Outputs", "Gene_Availability_Summary.xlsx"))

# ==============================================================================
# 2. COMPREHENSIVE HELPER FUNCTIONS
# ==============================================================================

# Enhanced Excel workbook creation with formatting
create_comprehensive_excel <- function(results_list, filename, highlight_genes = NULL) {
  wb <- createWorkbook()
  
  # Define styles
  header_style <- createStyle(textDecoration = "bold", fgFill = "#D3D3D3")
  significant_style <- createStyle(fgFill = "#FFE6E6")
  highlight_style <- createStyle(fgFill = "#E6F3FF", textDecoration = "bold")
  
  for(sheet_name in names(results_list)) {
    cat("  Creating sheet:", sheet_name, "\n")
    
    # Add worksheet
    addWorksheet(wb, sheet_name)
    
    # Write data
    if(is.data.frame(results_list[[sheet_name]]) && nrow(results_list[[sheet_name]]) > 0) {
      writeData(wb, sheet_name, results_list[[sheet_name]], rowNames = TRUE)
      
      # Apply header style
      addStyle(wb, sheet_name, header_style, rows = 1, cols = 1:ncol(results_list[[sheet_name]]))
      
      # Highlight significant genes
      if("p_val_adj" %in% colnames(results_list[[sheet_name]])) {
        p_col <- which(colnames(results_list[[sheet_name]]) == "p_val_adj")
        conditionalFormatting(wb, sheet_name, cols = p_col,
                             rows = 2:(nrow(results_list[[sheet_name]])+1),
                             rule = "<0.05", style = significant_style)
      }
      
      # Highlight genes of interest
      if(!is.null(highlight_genes) && "gene" %in% colnames(results_list[[sheet_name]])) {
        gene_col <- which(colnames(results_list[[sheet_name]]) == "gene")
        highlight_rows <- which(results_list[[sheet_name]]$gene %in% highlight_genes) + 1
        if(length(highlight_rows) > 0) {
          addStyle(wb, sheet_name, highlight_style, rows = highlight_rows, cols = 1:ncol(results_list[[sheet_name]]))
        }
      }
    } else {
      writeData(wb, sheet_name, data.frame(message = "No data available"))
    }
  }
  
  saveWorkbook(wb, filename, overwrite = TRUE)
  cat("Saved Excel file:", filename, "\n")
}

# Enhanced gene extraction function
extract_genes_comprehensive <- function(results, gene_sets, contrast_name = "") {
  extracted_results <- list()
  
  for(set_name in names(gene_sets)) {
    available_genes <- intersect(gene_sets[[set_name]], rownames(results))
    
    if(length(available_genes) > 0) {
      subset_results <- results[available_genes, , drop = FALSE]
      subset_results$gene <- rownames(subset_results)
      subset_results$gene_set <- set_name
      subset_results$contrast <- contrast_name
      
      extracted_results[[set_name]] <- subset_results
    }
  }
  
  return(extracted_results)
}

# Enhanced volcano plot function
create_enhanced_volcano <- function(results, title, highlight_genes = NULL, 
                                  gene_sets = NULL, save_path = NULL) {
  if(is.null(results) || nrow(results) == 0) {
    p <- ggplot() + 
      ggtitle(paste("No data for", title)) +
      theme_minimal()
    return(p)
  }
  
  # Prepare data
  results$gene <- rownames(results)
  results$log10_padj <- -log10(pmax(results$p_val_adj, 1e-300))
  
  # Define significance categories
  results$significance <- "NS"
  results$significance[results$p_val_adj < 0.05 & results$avg_log2FC > 0.5] <- "Upregulated"
  results$significance[results$p_val_adj < 0.05 & results$avg_log2FC < -0.5] <- "Downregulated"
  results$significance[results$p_val_adj < 0.01 & abs(results$avg_log2FC) > 1] <- "Highly Significant"
  
  # Highlight genes of interest
  if(!is.null(highlight_genes)) {
    results$highlight <- ifelse(results$gene %in% highlight_genes, "Genes of Interest", "Other")
  }
  
  # Create base plot
  p <- ggplot(results, aes(x = avg_log2FC, y = log10_padj)) +
    geom_point(aes(color = significance), alpha = 0.6, size = 0.8) +
    scale_color_manual(values = c("NS" = "gray70", "Upregulated" = "red", 
                                 "Downregulated" = "blue", "Highly Significant" = "darkred")) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
    labs(title = title, x = "Average log2(Fold Change)", y = "-log10(Adjusted P-value)") +
    theme_minimal() +
    theme(legend.position = "bottom",
          plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
  
  # Add gene labels for significant genes
  if(!is.null(highlight_genes)) {
    label_data <- results[results$gene %in% highlight_genes & results$p_val_adj < 0.05, ]
    if(nrow(label_data) > 0) {
      p <- p + geom_text_repel(data = label_data, 
                              aes(label = gene), 
                              max.overlaps = 20,
                              box.padding = 0.5,
                              point.padding = 0.5,
                              size = 3)
    }
  }
  
  # Save if path provided
  if(!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, device = "svg")
  }
  
  return(p)
}

# Enhanced heatmap function
create_comprehensive_heatmap <- function(seu_object, genes, group_by = "condition", 
                                       title = "Gene Expression Heatmap", 
                                       save_path = NULL) {
  # Get available genes
  available_genes <- intersect(genes, rownames(seu_object))
  
  if(length(available_genes) == 0) {
    cat("No genes available for heatmap:", title, "\n")
    return(NULL)
  }
  
  # Get expression data
  expr_data <- AverageExpression(seu_object, features = available_genes, 
                                group.by = group_by, assays = "RNA")$RNA
  
  # Scale the data
  expr_scaled <- t(scale(t(expr_data)))
  
  # Create annotation
  col_anno <- data.frame(
    condition = colnames(expr_scaled),
    row.names = colnames(expr_scaled)
  )
  
  # Create heatmap
  p <- pheatmap(expr_scaled,
                cluster_rows = TRUE,
                cluster_cols = TRUE,
                scale = "none",
                color = colorRampPalette(c("blue", "white", "red"))(100),
                main = title,
                fontsize = 10,
                angle_col = 45,
                silent = TRUE)
  
  # Save if path provided
  if(!is.null(save_path)) {
    svg(save_path, width = 10, height = 8)
    print(p)
    dev.off()
  }
  
  return(p)
}

# ==============================================================================
# 3. COMPREHENSIVE DIFFERENTIAL GENE EXPRESSION ANALYSIS
# ==============================================================================

cat("\n=== COMPREHENSIVE DIFFERENTIAL GENE EXPRESSION ANALYSIS ===\n")

# Define all possible contrasts
contrasts <- list(
  "DMG_vs_NAIVE" = c("DMG", "NAIVE"),
  "DMG_vs_SHAM" = c("DMG", "SHAM"),
  "SHAM_vs_NAIVE" = c("SHAM", "NAIVE")
)

# ==============================================================================
# 3.1 WILCOXON RANK SUM TEST (COMPREHENSIVE)
# ==============================================================================

cat("\n--- WILCOXON RANK SUM TEST ---\n")

DefaultAssay(seu_NKT_focused) <- "RNA"

# Initialize results storage
wilcox_results <- list()
wilcox_genes_of_interest <- list()
wilcox_summary_stats <- list()

# Perform comprehensive Wilcoxon tests
for(contrast_name in names(contrasts)) {
  cat("Processing contrast:", contrast_name, "\n")
  
  group1 <- contrasts[[contrast_name]][1]
  group2 <- contrasts[[contrast_name]][2]
  
  # Overall contrast
  tryCatch({
    cells_subset <- seu_NKT_focused$condition %in% c(group1, group2)
    seu_subset <- seu_NKT_focused[, cells_subset]
    
    # Run Wilcoxon test
    dge_result <- FindMarkers(seu_subset, 
                             ident.1 = group1,
                             ident.2 = group2,
                             group.by = "condition",
                             test.use = "wilcox",
                             min.pct = 0.05,
                             logfc.threshold = 0.05,
                             verbose = FALSE)
    
    dge_result$gene <- rownames(dge_result)
    dge_result$contrast <- contrast_name
    
    wilcox_results[[paste0("All_Genes_", contrast_name)]] <- dge_result
    
    # Extract genes of interest
    genes_extracted <- extract_genes_comprehensive(dge_result, gene_sets, contrast_name)
    wilcox_genes_of_interest <- c(wilcox_genes_of_interest, genes_extracted)
    
    # Summary statistics
    summary_stats <- data.frame(
      contrast = contrast_name,
      total_genes = nrow(dge_result),
      significant_genes = sum(dge_result$p_val_adj < 0.05),
      upregulated = sum(dge_result$p_val_adj < 0.05 & dge_result$avg_log2FC > 0),
      downregulated = sum(dge_result$p_val_adj < 0.05 & dge_result$avg_log2FC < 0),
      highly_significant = sum(dge_result$p_val_adj < 0.01 & abs(dge_result$avg_log2FC) > 1)
    )
    
    wilcox_summary_stats[[contrast_name]] <- summary_stats
    
    cat("  Found", nrow(dge_result), "genes (", sum(dge_result$p_val_adj < 0.05), "significant)\n")
    
  }, error = function(e) {
    cat("  Error:", e$message, "\n")
    wilcox_results[[paste0("All_Genes_", contrast_name)]] <- data.frame()
  })
  
  # Cell type specific contrasts
  for(cell_type in unique(seu_NKT_focused$cell_types)) {
    cat("  Processing cell type:", cell_type, "\n")
    
    tryCatch({
      cells_subset <- seu_NKT_focused$condition %in% c(group1, group2) & 
                     seu_NKT_focused$cell_types == cell_type
      seu_subset <- seu_NKT_focused[, cells_subset]
      
      if(ncol(seu_subset) > 10) {  # Minimum cells required
        dge_result <- FindMarkers(seu_subset, 
                                 ident.1 = group1,
                                 ident.2 = group2,
                                 group.by = "condition",
                                 test.use = "wilcox",
                                 min.pct = 0.05,
                                 logfc.threshold = 0.05,
                                 verbose = FALSE)
        
        dge_result$gene <- rownames(dge_result)
        dge_result$contrast <- contrast_name
        dge_result$cell_type <- cell_type
        
        wilcox_results[[paste0(cell_type, "_", contrast_name)]] <- dge_result
        
        cat("    Found", nrow(dge_result), "genes (", sum(dge_result$p_val_adj < 0.05), "significant)\n")
      }
      
    }, error = function(e) {
      cat("    Error:", e$message, "\n")
    })
  }
}

# Combine summary statistics
wilcox_summary_combined <- do.call(rbind, wilcox_summary_stats)

# Get top regulated genes
wilcox_top_genes <- list()
for(result_name in names(wilcox_results)) {
  if(nrow(wilcox_results[[result_name]]) > 0) {
    # Top upregulated
    up_genes <- wilcox_results[[result_name]][wilcox_results[[result_name]]$p_val_adj < 0.05 & 
                                              wilcox_results[[result_name]]$avg_log2FC > 0, ]
    if(nrow(up_genes) > 0) {
      up_genes <- up_genes[order(up_genes$avg_log2FC, decreasing = TRUE), ]
      wilcox_top_genes[[paste0(result_name, "_Top_Upregulated")]] <- head(up_genes, 50)
    }
    
    # Top downregulated
    down_genes <- wilcox_results[[result_name]][wilcox_results[[result_name]]$p_val_adj < 0.05 & 
                                                wilcox_results[[result_name]]$avg_log2FC < 0, ]
    if(nrow(down_genes) > 0) {
      down_genes <- down_genes[order(down_genes$avg_log2FC, decreasing = FALSE), ]
      wilcox_top_genes[[paste0(result_name, "_Top_Downregulated")]] <- head(down_genes, 50)
    }
  }
}

# Create comprehensive Wilcoxon Excel output
wilcox_excel_data <- c(
  wilcox_results,
  wilcox_genes_of_interest,
  wilcox_top_genes,
  list("Summary_Statistics" = wilcox_summary_combined)
)

create_comprehensive_excel(wilcox_excel_data, 
                          file.path(base_output, "Excel_Outputs", "DGE_Wilcoxon_Complete_Analysis.xlsx"),
                          highlight_genes = available_genes)

# ==============================================================================
# 3.2 MUSCAT + DESEQ2 (COMPREHENSIVE PSEUDOBULK)
# ==============================================================================

cat("\n--- MUSCAT + DESEQ2 PSEUDOBULK ANALYSIS ---\n")

deseq2_results <- list()
deseq2_genes_of_interest <- list()
deseq2_summary_stats <- list()

tryCatch({
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu_NKT_focused)
  
  # Set up for muscat
  sce$cluster_id <- sce$cell_types
  sce$sample_id <- sce$sample_name
  sce$group_id <- sce$condition
  
  # Remove samples with too few cells
  sample_counts <- table(sce$sample_id, sce$group_id)
  keep_samples <- rownames(sample_counts)[rowSums(sample_counts >= 10) > 0]
  sce <- sce[, sce$sample_id %in% keep_samples]
  
  # Create pseudobulk
  pb <- aggregateData(sce, assay = "counts", fun = "sum", by = c("cluster_id", "sample_id"))
  
  # Run DESeq2 for each contrast
  for(contrast_name in names(contrasts)) {
    cat("Processing contrast:", contrast_name, "\n")
    
    group1 <- contrasts[[contrast_name]][1]
    group2 <- contrasts[[contrast_name]][2]
    
    # Filter samples for this contrast
    samples_to_keep <- pb$group_id %in% c(group1, group2)
    pb_subset <- pb[, samples_to_keep]
    
    # Ensure minimum samples per group
    group_counts <- table(pb_subset$group_id)
    if(all(group_counts >= 3)) {
      
      tryCatch({
        # Run DESeq2
        res <- pbDS(pb_subset, method = "DESeq2", filter = "both", verbose = FALSE)
        
        # Process results
        if(length(res$res) > 0) {
          for(cell_type in names(res$res)) {
            ct_res <- res$res[[cell_type]]
            
            if(nrow(ct_res) > 0) {
              # Standardize column names
              ct_res$gene <- rownames(ct_res)
              ct_res$contrast <- contrast_name
              ct_res$cell_type <- cell_type
              
              # Rename columns to match Seurat format
              if("logFC" %in% colnames(ct_res)) {
                ct_res$avg_log2FC <- ct_res$logFC
              }
              if("p_adj.loc" %in% colnames(ct_res)) {
                ct_res$p_val_adj <- ct_res$p_adj.loc
              }
              
              deseq2_results[[paste0(cell_type, "_", contrast_name)]] <- ct_res
              
              # Extract genes of interest
              genes_extracted <- extract_genes_comprehensive(ct_res, gene_sets, contrast_name)
              deseq2_genes_of_interest <- c(deseq2_genes_of_interest, genes_extracted)
              
              cat("  ", cell_type, ":", nrow(ct_res), "genes (", 
                  sum(ct_res$p_val_adj < 0.05, na.rm = TRUE), "significant)\n")
            }
          }
        }
        
      }, error = function(e) {
        cat("  Error:", e$message, "\n")
      })
    } else {
      cat("  Insufficient samples for", contrast_name, "\n")
    }
  }
  
  # Create comprehensive DESeq2 Excel output
  deseq2_excel_data <- c(
    deseq2_results,
    deseq2_genes_of_interest
  )
  
  create_comprehensive_excel(deseq2_excel_data, 
                            file.path(base_output, "Excel_Outputs", "DGE_Muscat_DESeq2_Complete_Analysis.xlsx"),
                            highlight_genes = available_genes)
  
}, error = function(e) {
  cat("Error in Muscat + DESeq2:", e$message, "\n")
})

# ==============================================================================
# 3.3 MUSCAT + EDGER (COMPREHENSIVE PSEUDOBULK)
# ==============================================================================

cat("\n--- MUSCAT + EDGER PSEUDOBULK ANALYSIS ---\n")

edger_results <- list()
edger_genes_of_interest <- list()
edger_summary_stats <- list()

tryCatch({
  # Use the same SCE object from DESeq2
  if(exists("pb")) {
    
    # Run edgeR for each contrast
    for(contrast_name in names(contrasts)) {
      cat("Processing contrast:", contrast_name, "\n")
      
      group1 <- contrasts[[contrast_name]][1]
      group2 <- contrasts[[contrast_name]][2]
      
      # Filter samples for this contrast
      samples_to_keep <- pb$group_id %in% c(group1, group2)
      pb_subset <- pb[, samples_to_keep]
      
      # Ensure minimum samples per group
      group_counts <- table(pb_subset$group_id)
      if(all(group_counts >= 3)) {
        
        tryCatch({
          # Run edgeR
          res <- pbDS(pb_subset, method = "edgeR", filter = "both", verbose = FALSE)
          
          # Process results
          if(length(res$res) > 0) {
            for(cell_type in names(res$res)) {
              ct_res <- res$res[[cell_type]]
              
              if(nrow(ct_res) > 0) {
                # Standardize column names
                ct_res$gene <- rownames(ct_res)
                ct_res$contrast <- contrast_name
                ct_res$cell_type <- cell_type
                
                # Rename columns to match Seurat format
                if("logFC" %in% colnames(ct_res)) {
                  ct_res$avg_log2FC <- ct_res$logFC
                }
                if("p_adj.loc" %in% colnames(ct_res)) {
                  ct_res$p_val_adj <- ct_res$p_adj.loc
                }
                
                edger_results[[paste0(cell_type, "_", contrast_name)]] <- ct_res
                
                # Extract genes of interest
                genes_extracted <- extract_genes_comprehensive(ct_res, gene_sets, contrast_name)
                edger_genes_of_interest <- c(edger_genes_of_interest, genes_extracted)
                
                cat("  ", cell_type, ":", nrow(ct_res), "genes (", 
                    sum(ct_res$p_val_adj < 0.05, na.rm = TRUE), "significant)\n")
              }
            }
          }
          
        }, error = function(e) {
          cat("  Error:", e$message, "\n")
        })
      } else {
        cat("  Insufficient samples for", contrast_name, "\n")
      }
    }
    
    # Create comprehensive edgeR Excel output
    edger_excel_data <- c(
      edger_results,
      edger_genes_of_interest
    )
    
    create_comprehensive_excel(edger_excel_data, 
                              file.path(base_output, "Excel_Outputs", "DGE_Muscat_edgeR_Complete_Analysis.xlsx"),
                              highlight_genes = available_genes)
    
  } else {
    cat("Pseudobulk object not available for edgeR analysis\n")
  }
  
}, error = function(e) {
  cat("Error in Muscat + edgeR:", e$message, "\n")
})

# ==============================================================================
# 4. COMPREHENSIVE VISUALIZATIONS
# ==============================================================================

cat("\n=== COMPREHENSIVE VISUALIZATIONS ===\n")

# ==============================================================================
# 4.1 VOLCANO PLOTS FOR ALL METHODS
# ==============================================================================

cat("\n--- CREATING VOLCANO PLOTS ---\n")

# Wilcoxon volcano plots
for(result_name in names(wilcox_results)) {
  if(nrow(wilcox_results[[result_name]]) > 0) {
    volcano_path <- file.path(base_output, "Volcano_Plots", paste0("Volcano_Wilcoxon_", result_name, ".svg"))
    create_enhanced_volcano(wilcox_results[[result_name]], 
                           title = paste("Wilcoxon:", result_name),
                           highlight_genes = available_genes,
                           save_path = volcano_path)
  }
}

# DESeq2 volcano plots
for(result_name in names(deseq2_results)) {
  if(nrow(deseq2_results[[result_name]]) > 0) {
    volcano_path <- file.path(base_output, "Volcano_Plots", paste0("Volcano_DESeq2_", result_name, ".svg"))
    create_enhanced_volcano(deseq2_results[[result_name]], 
                           title = paste("DESeq2:", result_name),
                           highlight_genes = available_genes,
                           save_path = volcano_path)
  }
}

# edgeR volcano plots
for(result_name in names(edger_results)) {
  if(nrow(edger_results[[result_name]]) > 0) {
    volcano_path <- file.path(base_output, "Volcano_Plots", paste0("Volcano_edgeR_", result_name, ".svg"))
    create_enhanced_volcano(edger_results[[result_name]], 
                           title = paste("edgeR:", result_name),
                           highlight_genes = available_genes,
                           save_path = volcano_path)
  }
}

# ==============================================================================
# 4.2 COMPREHENSIVE HEATMAPS
# ==============================================================================

cat("\n--- CREATING HEATMAPS ---\n")

# Heatmaps for each gene set
for(set_name in names(gene_sets)) {
  cat("Creating heatmap for:", set_name, "\n")
  
  # Overall heatmap
  heatmap_path <- file.path(base_output, "Heatmaps", paste0("Heatmap_", set_name, "_by_condition.svg"))
  create_comprehensive_heatmap(seu_NKT_focused, gene_sets[[set_name]], 
                              group_by = "condition", 
                              title = paste(set_name, "Expression by Condition"),
                              save_path = heatmap_path)
  
  # Cell type specific heatmaps
  heatmap_path <- file.path(base_output, "Heatmaps", paste0("Heatmap_", set_name, "_by_cell_type.svg"))
  create_comprehensive_heatmap(seu_NKT_focused, gene_sets[[set_name]], 
                              group_by = "cell_types", 
                              title = paste(set_name, "Expression by Cell Type"),
                              save_path = heatmap_path)
}

# Combined heatmap with all genes of interest
combined_heatmap_path <- file.path(base_output, "Heatmaps", "Heatmap_All_Genes_Combined.svg")
create_comprehensive_heatmap(seu_NKT_focused, available_genes, 
                            group_by = "condition", 
                            title = "All Genes of Interest Expression",
                            save_path = combined_heatmap_path)

# ==============================================================================
# 4.3 FEATURE PLOTS AND DENSITY PLOTS
# ==============================================================================

cat("\n--- CREATING FEATURE AND DENSITY PLOTS ---\n")

# Feature plots for each gene set
for(set_name in names(gene_sets)) {
  available_set_genes <- intersect(gene_sets[[set_name]], rownames(seu_NKT_focused))
  
  if(length(available_set_genes) > 0) {
    cat("Creating feature plots for:", set_name, "\n")
    
    # Limit to top 6 genes to avoid overcrowding
    genes_to_plot <- available_set_genes[1:min(6, length(available_set_genes))]
    
    # Standard feature plots
    feature_plot <- FeaturePlot(seu_NKT_focused, 
                               features = genes_to_plot,
                               split.by = "condition",
                               ncol = 3,
                               pt.size = 0.5) +
      plot_annotation(title = paste(set_name, "Expression by Condition"))
    
    feature_path <- file.path(base_output, "Feature_Plots", paste0("Feature_", set_name, ".svg"))
    ggsave(feature_path, feature_plot, width = 15, height = 10, device = "svg")
    
    # Nebulosa density plots
    tryCatch({
      density_plots <- list()
      for(gene in genes_to_plot) {
        density_plots[[gene]] <- plot_density(seu_NKT_focused, gene, size = 0.5) +
          facet_wrap(~condition, ncol = 3) +
          ggtitle(gene) +
          theme(plot.title = element_text(hjust = 0.5))
      }
      
      density_combined <- wrap_plots(density_plots, ncol = 2)
      density_path <- file.path(base_output, "Density_Plots", paste0("Density_", set_name, ".svg"))
      ggsave(density_path, density_combined, width = 12, height = 8 * length(genes_to_plot)/2, device = "svg")
      
    }, error = function(e) {
      cat("Error creating density plots for", set_name, ":", e$message, "\n")
    })
  }
}

# ==============================================================================
# 4.4 VIOLIN AND RIDGE PLOTS
# ==============================================================================

cat("\n--- CREATING VIOLIN AND RIDGE PLOTS ---\n")

# Violin plots for key genes
for(set_name in names(gene_sets)) {
  available_set_genes <- intersect(gene_sets[[set_name]], rownames(seu_NKT_focused))
  
  if(length(available_set_genes) > 0) {
    genes_to_plot <- available_set_genes[1:min(6, length(available_set_genes))]
    
    # Violin plots
    violin_plot <- VlnPlot(seu_NKT_focused, 
                          features = genes_to_plot,
                          group.by = "condition",
                          ncol = 3,
                          pt.size = 0.1) +
      plot_annotation(title = paste(set_name, "Expression Distribution"))
    
    violin_path <- file.path(base_output, "Violin_Plots", paste0("Violin_", set_name, ".svg"))
    ggsave(violin_path, violin_plot, width = 15, height = 10, device = "svg")
    
    # Ridge plots
    ridge_plots <- list()
    for(gene in genes_to_plot) {
      expr_data <- FetchData(seu_NKT_focused, vars = c(gene, "condition"))
      colnames(expr_data) <- c("expression", "condition")
      
      ridge_plots[[gene]] <- ggplot(expr_data, aes(x = expression, y = condition, fill = condition)) +
        geom_density_ridges(alpha = 0.7) +
        scale_fill_manual(values = condition_colors) +
        ggtitle(gene) +
        theme_minimal() +
        theme(legend.position = "none")
    }
    
    ridge_combined <- wrap_plots(ridge_plots, ncol = 2)
    ridge_path <- file.path(base_output, "Ridge_Plots", paste0("Ridge_", set_name, ".svg"))
    ggsave(ridge_path, ridge_combined, width = 12, height = 8, device = "svg")
  }
}

# ==============================================================================
# 4.5 CORRELATION ANALYSIS
# ==============================================================================

cat("\n--- CREATING CORRELATION PLOTS ---\n")

# Correlation between S1P receptors and trafficking genes
s1p_available <- intersect(gene_sets$s1p_pathway, rownames(seu_NKT_focused))
trafficking_available <- intersect(gene_sets$trafficking, rownames(seu_NKT_focused))

if(length(s1p_available) > 0 && length(trafficking_available) > 0) {
  # Get expression data
  expr_data <- FetchData(seu_NKT_focused, vars = c(s1p_available, trafficking_available))
  
  # Calculate correlation
  cor_matrix <- cor(expr_data, use = "complete.obs")
  
  # Create correlation plot
  svg(file.path(base_output, "Correlation_Plots", "S1P_Trafficking_Correlation.svg"), 
      width = 10, height = 8)
  corrplot(cor_matrix, method = "color", type = "upper", order = "hclust",
           tl.cex = 0.8, tl.col = "black", tl.srt = 45,
           title = "Correlation: S1P Pathway vs Trafficking Genes",
           mar = c(0,0,2,0))
  dev.off()
}

# ==============================================================================
# 5. TRAJECTORY ANALYSIS
# ==============================================================================

cat("\n--- TRAJECTORY ANALYSIS ---\n")

tryCatch({
  # Convert to cell_data_set for monocle3
  cds <- as.cell_data_set(seu_NKT_focused)
  
  # Add metadata
  cds@colData$condition <- seu_NKT_focused$condition
  cds@colData$cell_type <- seu_NKT_focused$cell_types
  
  # Preprocess the data
  cds <- preprocess_cds(cds, num_dim = 30)
  cds <- align_cds(cds, alignment_group = "condition")
  cds <- reduce_dimension(cds, reduction_method = "UMAP")
  
  # Cluster cells and learn trajectory
  cds <- cluster_cells(cds, verbose = FALSE)
  cds <- learn_graph(cds, verbose = FALSE)
  
  # Order cells in pseudotime
  cds <- order_cells(cds)
  
  # Create trajectory plots
  traj_condition <- plot_cells(cds, color_cells_by = "condition", 
                              label_cell_groups = FALSE,
                              label_leaves = FALSE,
                              label_branch_points = FALSE) +
    scale_color_manual(values = condition_colors) +
    ggtitle("Trajectory Analysis: Colored by Condition") +
    theme_minimal()
  
  ggsave(file.path(base_output, "Trajectory_Analysis", "Trajectory_by_Condition.svg"), 
         traj_condition, width = 10, height = 8, device = "svg")
  
  traj_cell_type <- plot_cells(cds, color_cells_by = "cell_type", 
                              label_cell_groups = FALSE,
                              label_leaves = FALSE,
                              label_branch_points = FALSE) +
    scale_color_manual(values = cell_type_colors) +
    ggtitle("Trajectory Analysis: Colored by Cell Type") +
    theme_minimal()
  
  ggsave(file.path(base_output, "Trajectory_Analysis", "Trajectory_by_Cell_Type.svg"), 
         traj_cell_type, width = 10, height = 8, device = "svg")
  
  # Pseudotime plot
  traj_pseudotime <- plot_cells(cds, color_cells_by = "pseudotime",
                               label_cell_groups = FALSE,
                               label_leaves = FALSE,
                               label_branch_points = FALSE) +
    ggtitle("Trajectory Analysis: Pseudotime") +
    theme_minimal()
  
  ggsave(file.path(base_output, "Trajectory_Analysis", "Trajectory_Pseudotime.svg"), 
         traj_pseudotime, width = 10, height = 8, device = "svg")
  
  # Gene expression along trajectory
  if(length(s1p_available) > 0) {
    for(gene in s1p_available[1:min(4, length(s1p_available))]) {
      traj_gene <- plot_cells(cds, genes = gene,
                             label_cell_groups = FALSE,
                             label_leaves = FALSE,
                             label_branch_points = FALSE) +
        ggtitle(paste("Trajectory Expression:", gene)) +
        theme_minimal()
      
      ggsave(file.path(base_output, "Trajectory_Analysis", paste0("Trajectory_", gene, ".svg")), 
             traj_gene, width = 10, height = 8, device = "svg")
    }
  }
  
}, error = function(e) {
  cat("Error in trajectory analysis:", e$message, "\n")
})

# ==============================================================================
# 6. FINAL SUMMARY AND SESSION INFO
# ==============================================================================

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Timestamp:", format(Sys.time(), "%Y-%m-%d %H:%M:%S"), "\n")
cat("Output directory:", base_output, "\n")

# Create summary report
summary_info <- list(
  analysis_date = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  user = "fionnspencer07",
  total_cells = ncol(seu_NKT_focused),
  total_genes = nrow(seu_NKT_focused),
  conditions = levels(seu_NKT_focused$condition),
  cell_types = unique(seu_NKT_focused$cell_types),
  genes_of_interest_available = length(available_genes),
  genes_of_interest_missing = length(missing_genes),
  contrasts_analyzed = names(contrasts),
  methods_used = c("Wilcoxon", "Muscat+DESeq2", "Muscat+edgeR")
)

# Save summary
write.xlsx(summary_info, file.path(base_output, "Summary_Stats", "Analysis_Summary.xlsx"))

# Save session info
writeLines(capture.output(sessionInfo()), file.path(base_output, "session_info.txt"))

# List all created files
cat("\nFiles created:\n")
for(subdir in subdirs) {
  files <- list.files(file.path(base_output, subdir), recursive = TRUE)
  if(length(files) > 0) {
    cat(paste0("  ", subdir, "/:\n"))
    for(file in files) {
      cat(paste0("    ", file, "\n"))
    }
  }
}

cat("\nComprehensive analysis pipeline completed successfully!\n")

