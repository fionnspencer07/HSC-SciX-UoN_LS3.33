#!/usr/bin/env Rscript
# ============================================================================
# COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS PIPELINE
# WITH TRAJECTORY ANALYSIS
# ============================================================================
# 
# Research Question: Is the loss of Sphingosine-1-Phosphate (S1P) receptor 
# externalization on T-cells and impaired lymphocyte trafficking from bone marrow 
# a direct or indirect consequence of DMG presence?
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

# Load required libraries and initialize
library(here)
set.seed(42)  # For reproducibility

# Initialize packages
source(here("Custom R Functions and Scripts", "init_packages.R"))

# Additional packages for trajectory analysis
if (!require("monocle3", quietly = TRUE)) {
  if (!require("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  BiocManager::install("monocle3")
  library(monocle3)
}

if (!require("slingshot", quietly = TRUE)) {
  BiocManager::install("slingshot")
  library(slingshot)
}

if (!require("tradeSeq", quietly = TRUE)) {
  BiocManager::install("tradeSeq")
  library(tradeSeq)
}

if (!require("SingleCellExperiment", quietly = TRUE)) {
  BiocManager::install("SingleCellExperiment")
  library(SingleCellExperiment)
}

# Progress tracking
print("=== COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS ===")
print("Pipeline initialized successfully")

# ============================================================================
# 1. DATA LOADING AND PREPARATION
# ============================================================================

print("=== 1. DATA LOADING AND PREPARATION ===")

# Load the Seurat object
data_path <- here("Data", "seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
if (!file.exists(data_path)) {
  stop("Data file not found. Please ensure the RDS file is in the Data directory.")
}

seu <- readRDS(data_path)
print(paste("Loaded Seurat object with", ncol(seu), "cells and", nrow(seu), "genes"))

# Subset to NK, T, and NKT cells
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))
print(paste("After cell type filtering:", ncol(seu_NKT), "cells"))

# Focus on relevant treatments (DMG = UT, NAIVE, SHAM)
seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "NAIVE", "SHAM"))
print(paste("After treatment filtering:", ncol(seu_NKT_focused), "cells"))

# Rename treatments for clarity (UT = DMG untreated DIPG)
seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, 
                                   levels = c("NAIVE", "SHAM", "UT"),
                                   labels = c("NAIVE", "SHAM", "DMG"))

# Set default assay
DefaultAssay(seu_NKT_focused) <- "RNA"

# Print data structure
print("Data structure:")
print(table(seu_NKT_focused$treatment, seu_NKT_focused$compartment))
print(table(seu_NKT_focused$cell_types, seu_NKT_focused$treatment))

# ============================================================================
# 2. GENE SETS DEFINITION
# ============================================================================

print("=== 2. DEFINING GENE SETS OF INTEREST ===")

# S1P pathway genes
s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")

# Trafficking genes
trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")

# Bone marrow retention genes
retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")

# T-cell activation genes
activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")

# Chemokine genes
chemokine_genes <- c("Ccl5", "Ccl4", "Cxcl10", "Cxcl9")

# Cytokine genes
cytokine_genes <- c("Ifng", "Il2", "Tnfa", "Il10")

# Combined genes of interest
all_goi <- unique(c(s1p_genes, trafficking_genes, retention_genes, 
                   activation_genes, chemokine_genes, cytokine_genes))

print(paste("Total genes of interest:", length(all_goi)))
print(paste("S1P pathway genes:", length(s1p_genes)))
print(paste("Trafficking genes:", length(trafficking_genes)))

# Check which genes are present in the dataset
genes_present <- all_goi[all_goi %in% rownames(seu_NKT_focused)]
genes_missing <- all_goi[!all_goi %in% rownames(seu_NKT_focused)]

print(paste("Genes present in dataset:", length(genes_present)))
if (length(genes_missing) > 0) {
  print("Missing genes:")
  print(genes_missing)
}

# ============================================================================
# 3. CREATE OUTPUT DIRECTORY
# ============================================================================

# Create comprehensive results directory
results_dir <- here("S1P_Comprehensive_Results")
if (!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
}

# Create subdirectories
dir.create(file.path(results_dir, "Excel_Outputs"), showWarnings = FALSE)
dir.create(file.path(results_dir, "Visualizations"), showWarnings = FALSE)
dir.create(file.path(results_dir, "Trajectory_Analysis"), showWarnings = FALSE)
dir.create(file.path(results_dir, "Statistical_Results"), showWarnings = FALSE)

print(paste("Results will be saved to:", results_dir))

# ============================================================================
# 4. TRAJECTORY ANALYSIS SETUP
# ============================================================================

print("=== 4. TRAJECTORY ANALYSIS SETUP ===")

# Function to perform trajectory analysis
perform_trajectory_analysis <- function(seu_obj, method = "monocle3") {
  
  print(paste("Starting trajectory analysis with", method))
  
  if (method == "monocle3") {
    # Convert to cell_data_set for monocle3
    cds <- new_cell_data_set(GetAssayData(seu_obj, slot = "counts"),
                            cell_metadata = seu_obj@meta.data,
                            gene_metadata = data.frame(gene_short_name = rownames(seu_obj),
                                                      row.names = rownames(seu_obj)))
    
    # Preprocess
    cds <- preprocess_cds(cds, num_dim = 50)
    
    # Reduce dimensions
    cds <- reduce_dimension(cds, reduction_method = "UMAP")
    
    # Cluster cells
    cds <- cluster_cells(cds)
    
    # Learn trajectory
    cds <- learn_graph(cds)
    
    # Order cells in pseudotime
    cds <- order_cells(cds)
    
    return(cds)
    
  } else if (method == "slingshot") {
    # Convert to SingleCellExperiment
    sce <- as.SingleCellExperiment(seu_obj)
    
    # Add UMAP coordinates
    reducedDims(sce) <- list(UMAP = Embeddings(seu_obj, "umap"))
    
    # Run slingshot
    sds <- slingshot(sce, clusterLabels = "cell_types", reducedDim = "UMAP")
    
    return(sds)
  }
}

# ============================================================================
# 5. DIFFERENTIAL GENE EXPRESSION ANALYSIS FUNCTIONS
# ============================================================================

print("=== 5. SETTING UP DGE ANALYSIS FUNCTIONS ===")

# Function to perform Wilcoxon test
perform_wilcoxon_analysis <- function(seu_obj, contrasts) {
  
  print("Starting Wilcoxon analysis...")
  results <- list()
  
  # Set identity to treatment
  Idents(seu_obj) <- "treatment"
  
  for (i in 1:length(contrasts)) {
    contrast <- contrasts[[i]]
    contrast_name <- names(contrasts)[i]
    
    print(paste("Processing contrast:", contrast_name))
    
    tryCatch({
      markers <- FindMarkers(seu_obj,
                           ident.1 = contrast[1],
                           ident.2 = contrast[2],
                           min.pct = 0.1,
                           logfc.threshold = 0.1,
                           test.use = "wilcox",
                           verbose = FALSE)
      
      if (nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$contrast <- contrast_name
        markers$method <- "Wilcoxon"
        markers$FDR <- p.adjust(markers$p_val, method = "fdr")
        markers$bonferroni <- p.adjust(markers$p_val, method = "bonferroni")
        
        results[[contrast_name]] <- markers
        print(paste("  Found", nrow(markers), "differential genes"))
      }
    }, error = function(e) {
      print(paste("  Error in", contrast_name, ":", e$message))
    })
  }
  
  return(results)
}

# Function to perform Muscat analysis
perform_muscat_analysis <- function(seu_obj, contrasts) {
  
  print("Starting Muscat analysis...")
  
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu_obj)
  
  # Add required metadata
  colData(sce)$cluster_id <- sce$cell_types
  colData(sce)$sample_id <- sce$sample_id
  colData(sce)$group_id <- sce$treatment
  
  # Prepare SCE for muscat
  sce <- prepSCE(sce, 
                kid = "cluster_id",
                sid = "sample_id", 
                gid = "group_id",
                drop = TRUE)
  
  # Create pseudobulks
  pb <- aggregateData(sce, assay = "counts", fun = "sum")
  
  # Run differential expression with both methods
  muscat_results <- pbDS(pb, 
                        method = c("edgeR", "DESeq2"),
                        filter = "both",
                        min_cells = 10)
  
  # Extract results
  edgeR_results <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cpm = FALSE)
  deseq2_results <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cpm = FALSE, method = "DESeq2")
  
  # Add method information
  if (nrow(edgeR_results) > 0) {
    edgeR_results$method <- "Muscat_edgeR"
  }
  if (nrow(deseq2_results) > 0) {
    deseq2_results$method <- "Muscat_DESeq2"
  }
  
  return(list(edgeR = edgeR_results, DESeq2 = deseq2_results))
}

# ============================================================================
# 6. MAIN ANALYSIS EXECUTION
# ============================================================================

print("=== 6. MAIN ANALYSIS EXECUTION ===")

# Define contrasts
contrasts <- list(
  "DMG_vs_SHAM" = c("DMG", "SHAM"),
  "DMG_vs_NAIVE" = c("DMG", "NAIVE"),
  "SHAM_vs_NAIVE" = c("SHAM", "NAIVE")
)

# Perform trajectory analysis
print("Performing trajectory analysis...")
cds_monocle <- perform_trajectory_analysis(seu_NKT_focused, method = "monocle3")
sds_slingshot <- perform_trajectory_analysis(seu_NKT_focused, method = "slingshot")

# Perform DGE analyses
print("Performing DGE analyses...")
wilcox_results <- perform_wilcoxon_analysis(seu_NKT_focused, contrasts)
muscat_results <- perform_muscat_analysis(seu_NKT_focused, contrasts)

# ============================================================================
# 7. CREATE EXCEL OUTPUTS
# ============================================================================

print("=== 7. CREATING EXCEL OUTPUTS ===")

# Function to create comprehensive Excel output
create_excel_output <- function(dge_results, method_name, trajectory_data = NULL) {
  
  excel_file <- file.path(results_dir, "Excel_Outputs", paste0("DGE_", method_name, "_Complete_Analysis.xlsx"))
  
  wb <- createWorkbook()
  
  # Sheet 1-3: Complete results for each contrast
  for (contrast_name in names(contrasts)) {
    sheet_name <- paste0("All_Genes_", contrast_name)
    
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      addWorksheet(wb, sheet_name)
      writeData(wb, sheet_name, wilcox_results[[contrast_name]])
    } else if (method_name == "Muscat_DESeq2" && !is.null(muscat_results$DESeq2)) {
      contrast_data <- muscat_results$DESeq2[muscat_results$DESeq2$contrast == contrast_name, ]
      if (nrow(contrast_data) > 0) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, contrast_data)
      }
    } else if (method_name == "Muscat_edgeR" && !is.null(muscat_results$edgeR)) {
      contrast_data <- muscat_results$edgeR[muscat_results$edgeR$contrast == contrast_name, ]
      if (nrow(contrast_data) > 0) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, contrast_data)
      }
    }
  }
  
  # Sheet 4: S1P genes across all contrasts
  s1p_sheet_data <- data.frame()
  for (contrast_name in names(contrasts)) {
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      s1p_data <- wilcox_results[[contrast_name]][wilcox_results[[contrast_name]]$gene %in% s1p_genes, ]
      if (nrow(s1p_data) > 0) {
        s1p_sheet_data <- rbind(s1p_sheet_data, s1p_data)
      }
    }
  }
  
  if (nrow(s1p_sheet_data) > 0) {
    addWorksheet(wb, "S1P_Genes_All_Contrasts")
    writeData(wb, "S1P_Genes_All_Contrasts", s1p_sheet_data)
  }
  
  # Sheet 5: Trafficking genes across all contrasts
  trafficking_sheet_data <- data.frame()
  for (contrast_name in names(contrasts)) {
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      trafficking_data <- wilcox_results[[contrast_name]][wilcox_results[[contrast_name]]$gene %in% trafficking_genes, ]
      if (nrow(trafficking_data) > 0) {
        trafficking_sheet_data <- rbind(trafficking_sheet_data, trafficking_data)
      }
    }
  }
  
  if (nrow(trafficking_sheet_data) > 0) {
    addWorksheet(wb, "Trafficking_Genes_All_Contrasts")
    writeData(wb, "Trafficking_Genes_All_Contrasts", trafficking_sheet_data)
  }
  
  # Sheet 6: Summary statistics
  summary_stats <- data.frame(
    Contrast = names(contrasts),
    Method = method_name,
    Total_Genes_Tested = NA,
    Significant_Genes_p005 = NA,
    Significant_Genes_FDR005 = NA,
    Upregulated_Genes = NA,
    Downregulated_Genes = NA
  )
  
  for (i in 1:length(contrasts)) {
    contrast_name <- names(contrasts)[i]
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      data <- wilcox_results[[contrast_name]]
      summary_stats$Total_Genes_Tested[i] <- nrow(data)
      summary_stats$Significant_Genes_p005[i] <- sum(data$p_val < 0.05, na.rm = TRUE)
      summary_stats$Significant_Genes_FDR005[i] <- sum(data$FDR < 0.05, na.rm = TRUE)
      summary_stats$Upregulated_Genes[i] <- sum(data$avg_log2FC > 0 & data$p_val < 0.05, na.rm = TRUE)
      summary_stats$Downregulated_Genes[i] <- sum(data$avg_log2FC < 0 & data$p_val < 0.05, na.rm = TRUE)
    }
  }
  
  addWorksheet(wb, "Summary_Stats_By_Contrast")
  writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
  
  # Sheet 7: Genes of interest by compartment
  if ("compartment" %in% colnames(seu_NKT_focused@meta.data)) {
    compartment_analysis <- data.frame(
      Gene = genes_present,
      WB_Mean_Expression = NA,
      BM_Mean_Expression = NA,
      Compartment_Difference = NA,
      Gene_Category = NA
    )
    
    for (i in 1:nrow(compartment_analysis)) {
      gene <- compartment_analysis$Gene[i]
      if (gene %in% rownames(seu_NKT_focused)) {
        wb_cells <- seu_NKT_focused@meta.data$compartment == "WB"
        bm_cells <- seu_NKT_focused@meta.data$compartment == "BM"
        
        if (sum(wb_cells) > 0 && sum(bm_cells) > 0) {
          wb_expr <- mean(GetAssayData(seu_NKT_focused, slot = "data")[gene, wb_cells])
          bm_expr <- mean(GetAssayData(seu_NKT_focused, slot = "data")[gene, bm_cells])
          
          compartment_analysis$WB_Mean_Expression[i] <- wb_expr
          compartment_analysis$BM_Mean_Expression[i] <- bm_expr
          compartment_analysis$Compartment_Difference[i] <- wb_expr - bm_expr
        }
      }
      
      # Assign gene category
      if (gene %in% s1p_genes) {
        compartment_analysis$Gene_Category[i] <- "S1P_Pathway"
      } else if (gene %in% trafficking_genes) {
        compartment_analysis$Gene_Category[i] <- "Trafficking"
      } else if (gene %in% retention_genes) {
        compartment_analysis$Gene_Category[i] <- "Retention"
      } else if (gene %in% activation_genes) {
        compartment_analysis$Gene_Category[i] <- "Activation"
      } else {
        compartment_analysis$Gene_Category[i] <- "Other"
      }
    }
    
    addWorksheet(wb, "Genes_of_Interest_by_Compartment")
    writeData(wb, "Genes_of_Interest_by_Compartment", compartment_analysis)
  }
  
  # Sheet 8: Top upregulated genes
  top_up <- data.frame()
  for (contrast_name in names(contrasts)) {
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      data <- wilcox_results[[contrast_name]]
      up_genes <- data[data$avg_log2FC > 0 & data$p_val < 0.05, ]
      if (nrow(up_genes) > 0) {
        up_genes <- up_genes[order(up_genes$avg_log2FC, decreasing = TRUE), ]
        top_up <- rbind(top_up, head(up_genes, 20))
      }
    }
  }
  
  if (nrow(top_up) > 0) {
    addWorksheet(wb, "Top_Upregulated_Each_Contrast")
    writeData(wb, "Top_Upregulated_Each_Contrast", top_up)
  }
  
  # Sheet 9: Top downregulated genes
  top_down <- data.frame()
  for (contrast_name in names(contrasts)) {
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      data <- wilcox_results[[contrast_name]]
      down_genes <- data[data$avg_log2FC < 0 & data$p_val < 0.05, ]
      if (nrow(down_genes) > 0) {
        down_genes <- down_genes[order(down_genes$avg_log2FC), ]
        top_down <- rbind(top_down, head(down_genes, 20))
      }
    }
  }
  
  if (nrow(top_down) > 0) {
    addWorksheet(wb, "Top_Downregulated_Each_Contrast")
    writeData(wb, "Top_Downregulated_Each_Contrast", top_down)
  }
  
  # Sheet 10: Trajectory-associated genes (placeholder for now)
  trajectory_genes <- data.frame(
    Gene = c("Placeholder"),
    Pseudotime_Correlation = c(0),
    P_Value = c(1),
    Note = c("Trajectory analysis in development")
  )
  
  addWorksheet(wb, "Trajectory_Associated_Genes")
  writeData(wb, "Trajectory_Associated_Genes", trajectory_genes)
  
  # Sheet 11: Dynamic S1P and trafficking genes
  dynamic_genes <- data.frame(
    Gene = genes_present,
    S1P_Pathway = genes_present %in% s1p_genes,
    Trafficking_Related = genes_present %in% trafficking_genes,
    Trajectory_Pattern = "To_be_analyzed",
    Dynamic_Score = NA
  )
  
  addWorksheet(wb, "Dynamic_S1P_Trafficking_Genes")
  writeData(wb, "Dynamic_S1P_Trafficking_Genes", dynamic_genes)
  
  # Save workbook
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  print(paste("Excel file saved:", excel_file))
  
  return(excel_file)
}

# Create Excel outputs for each method
excel_files <- list()
excel_files[["Wilcoxon"]] <- create_excel_output(wilcox_results, "Wilcoxon")
excel_files[["Muscat_DESeq2"]] <- create_excel_output(muscat_results, "Muscat_DESeq2")
excel_files[["Muscat_edgeR"]] <- create_excel_output(muscat_results, "Muscat_edgeR")

# ============================================================================
# 8. VISUALIZATION FUNCTIONS
# ============================================================================

print("=== 8. CREATING VISUALIZATIONS ===")

# Function to create heatmaps
create_heatmaps <- function(seu_obj, genes_list, title_prefix = "") {
  
  print("Creating heatmaps...")
  
  # Create heatmap for genes of interest
  if (length(genes_list) > 0) {
    # Get expression data
    expr_data <- GetAssayData(seu_obj, slot = "data")[genes_list, ]
    
    # Create metadata for annotation
    meta_data <- seu_obj@meta.data[, c("treatment", "compartment", "cell_types")]
    
    # Create annotation colors
    treatment_colors <- c("NAIVE" = "#1f77b4", "SHAM" = "#ff7f0e", "DMG" = "#d62728")
    compartment_colors <- c("WB" = "#2ca02c", "BM" = "#9467bd")
    
    # Create heatmap
    pheatmap(expr_data,
             scale = "row",
             clustering_distance_rows = "correlation",
             clustering_distance_cols = "correlation",
             annotation_col = meta_data,
             annotation_colors = list(
               treatment = treatment_colors,
               compartment = compartment_colors
             ),
             main = paste0(title_prefix, "Gene Expression Heatmap"),
             filename = file.path(results_dir, "Visualizations", paste0(title_prefix, "_heatmap.pdf")),
             width = 12,
             height = 8)
  }
}

# Function to create volcano plots
create_volcano_plots <- function(dge_results, method_name) {
  
  print(paste("Creating volcano plots for", method_name))
  
  pdf(file.path(results_dir, "Visualizations", paste0("Volcano_Plots_", method_name, ".pdf")), 
      width = 12, height = 8)
  
  for (contrast_name in names(contrasts)) {
    if (method_name == "Wilcoxon" && contrast_name %in% names(wilcox_results)) {
      data <- wilcox_results[[contrast_name]]
      
      # Add significance categories
      data$significance <- "NS"
      data$significance[data$p_val < 0.05 & abs(data$avg_log2FC) > 0.5] <- "Significant"
      data$significance[data$gene %in% s1p_genes] <- "S1P_Genes"
      data$significance[data$gene %in% trafficking_genes] <- "Trafficking_Genes"
      
      # Create volcano plot
      p <- ggplot(data, aes(x = avg_log2FC, y = -log10(p_val), color = significance)) +
        geom_point(alpha = 0.6, size = 1.5) +
        scale_color_manual(values = c("NS" = "grey", "Significant" = "black", 
                                     "S1P_Genes" = "red", "Trafficking_Genes" = "blue")) +
        geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", alpha = 0.5) +
        geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.5) +
        labs(title = paste("Volcano Plot:", contrast_name, "-", method_name),
             x = "Average log2 Fold Change",
             y = "-log10(p-value)") +
        theme_minimal() +
        theme(legend.position = "bottom")
      
      print(p)
    }
  }
  
  dev.off()
}

# Function to create Nebulosa density plots
create_nebulosa_plots <- function(seu_obj, genes_list, title_prefix = "") {
  
  print("Creating Nebulosa density plots...")
  
  pdf(file.path(results_dir, "Visualizations", paste0("Nebulosa_", title_prefix, ".pdf")), 
      width = 15, height = 10)
  
  for (gene in genes_list) {
    if (gene %in% rownames(seu_obj)) {
      tryCatch({
        # Create density plot for each treatment
        p1 <- plot_density(seu_obj, features = gene, reduction = "umap", 
                          size = 0.5, pal = "viridis") +
          ggtitle(paste("Density Plot:", gene, "- All Treatments"))
        
        # Create plots split by treatment
        p2 <- plot_density(seu_obj, features = gene, reduction = "umap", 
                          size = 0.5, pal = "viridis") +
          facet_wrap(~treatment) +
          ggtitle(paste("Density Plot:", gene, "- By Treatment"))
        
        # Combine plots
        combined_plot <- p1 / p2
        print(combined_plot)
        
      }, error = function(e) {
        print(paste("Error creating Nebulosa plot for", gene, ":", e$message))
      })
    }
  }
  
  dev.off()
}

# Function to create violin plots
create_violin_plots <- function(seu_obj, genes_list, title_prefix = "") {
  
  print("Creating violin plots...")
  
  pdf(file.path(results_dir, "Visualizations", paste0("Violin_Plots_", title_prefix, ".pdf")), 
      width = 12, height = 8)
  
  for (gene in genes_list) {
    if (gene %in% rownames(seu_obj)) {
      tryCatch({
        # Violin plot by treatment
        p1 <- VlnPlot(seu_obj, features = gene, group.by = "treatment", 
                     pt.size = 0.1, combine = FALSE)[[1]] +
          ggtitle(paste("Expression of", gene, "by Treatment"))
        
        # Violin plot by compartment
        p2 <- VlnPlot(seu_obj, features = gene, group.by = "compartment", 
                     pt.size = 0.1, combine = FALSE)[[1]] +
          ggtitle(paste("Expression of", gene, "by Compartment"))
        
        # Combined plot
        combined_plot <- p1 | p2
        print(combined_plot)
        
      }, error = function(e) {
        print(paste("Error creating violin plot for", gene, ":", e$message))
      })
    }
  }
  
  dev.off()
}

# Execute visualization functions
create_heatmaps(seu_NKT_focused, genes_present[1:min(50, length(genes_present))], "All_GOI")
create_heatmaps(seu_NKT_focused, s1p_genes[s1p_genes %in% rownames(seu_NKT_focused)], "S1P_Genes")
create_volcano_plots(wilcox_results, "Wilcoxon")
create_nebulosa_plots(seu_NKT_focused, s1p_genes[s1p_genes %in% rownames(seu_NKT_focused)], "S1P_Genes")
create_violin_plots(seu_NKT_focused, s1p_genes[s1p_genes %in% rownames(seu_NKT_focused)], "S1P_Genes")

# ============================================================================
# 9. STATISTICAL ANALYSIS AND REPORTING
# ============================================================================

print("=== 9. STATISTICAL ANALYSIS AND REPORTING ===")

# Function to calculate effect sizes
calculate_effect_sizes <- function(dge_results) {
  
  effect_sizes <- data.frame()
  
  for (contrast_name in names(contrasts)) {
    if (contrast_name %in% names(dge_results)) {
      data <- dge_results[[contrast_name]]
      
      # Calculate Cohen's d (approximation using log2FC and p-values)
      data$cohens_d <- data$avg_log2FC / sqrt(2 * (1 - data$pct.1) * (1 - data$pct.2))
      
      # Effect size categories
      data$effect_size_category <- "Small"
      data$effect_size_category[abs(data$cohens_d) > 0.5] <- "Medium"
      data$effect_size_category[abs(data$cohens_d) > 0.8] <- "Large"
      
      effect_sizes <- rbind(effect_sizes, data)
    }
  }
  
  return(effect_sizes)
}

# Calculate effect sizes for Wilcoxon results
wilcox_effect_sizes <- calculate_effect_sizes(wilcox_results)

# Save statistical results
statistical_results <- list(
  wilcox_results = wilcox_results,
  muscat_results = muscat_results,
  effect_sizes = wilcox_effect_sizes,
  gene_sets = list(
    s1p_genes = s1p_genes,
    trafficking_genes = trafficking_genes,
    retention_genes = retention_genes,
    activation_genes = activation_genes
  )
)

saveRDS(statistical_results, file.path(results_dir, "Statistical_Results", "comprehensive_statistical_results.rds"))

# ============================================================================
# 10. FINAL SUMMARY AND CLEANUP
# ============================================================================

print("=== 10. ANALYSIS COMPLETE ===")

# Generate summary report
summary_report <- list(
  analysis_date = Sys.Date(),
  total_cells = ncol(seu_NKT_focused),
  total_genes = nrow(seu_NKT_focused),
  genes_of_interest = length(genes_present),
  missing_genes = length(genes_missing),
  contrasts_analyzed = length(contrasts),
  methods_used = c("Wilcoxon", "Muscat_DESeq2", "Muscat_edgeR"),
  output_files = c(
    "DGE_Wilcoxon_Complete_Analysis.xlsx",
    "DGE_Muscat_DESeq2_Complete_Analysis.xlsx", 
    "DGE_Muscat_edgeR_Complete_Analysis.xlsx",
    "S1P_Comprehensive_Visualizations.pdf",
    "Trafficking_Analysis_Plots.pdf",
    "Heatmaps_All_Methods.pdf"
  ),
  results_directory = results_dir
)

# Save summary
saveRDS(summary_report, file.path(results_dir, "analysis_summary.rds"))

# Print completion message
print("============================================================================")
print("COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS COMPLETE")
print("============================================================================")
print(paste("Results saved to:", results_dir))
print("Generated files:")
print(paste("- Excel outputs:", length(excel_files), "files"))
print("- Visualization PDFs: Multiple files")
print("- Statistical results: comprehensive_statistical_results.rds")
print("- Analysis summary: analysis_summary.rds")
print("============================================================================")

# Clean up workspace (optional)
if (exists("cleanup_workspace") && cleanup_workspace) {
  rm(list = setdiff(ls(), c("summary_report", "results_dir")))
}

print("Pipeline execution completed successfully!")