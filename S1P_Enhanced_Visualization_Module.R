#!/usr/bin/env Rscript
# ============================================================================
# ENHANCED VISUALIZATION MODULE FOR S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING
# ============================================================================
# 
# This module provides comprehensive visualization functions for the S1P
# receptor and lymphocyte trafficking analysis, including heatmaps, volcano
# plots, Nebulosa density plots, violin plots, and trajectory visualizations.
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

# Load required libraries
suppressMessages({
  library(Seurat)
  library(ggplot2)
  library(pheatmap)
  library(ComplexHeatmap)
  library(circlize)
  library(viridis)
  library(RColorBrewer)
  library(Nebulosa)
  library(EnhancedVolcano)
  library(patchwork)
  library(ggrepel)
  library(ggridges)
  library(dplyr)
  library(tidyr)
  library(here)
  library(grid)
  library(gridExtra)
})

# Set theme for consistent plotting
theme_set(theme_minimal())

# ============================================================================
# HEATMAP FUNCTIONS
# ============================================================================

#' Create comprehensive heatmaps with fixed implementation
#'
#' @param seu_obj Seurat object
#' @param genes_list Vector of genes to include in heatmap
#' @param group_by Column name to group cells by
#' @param split_by Column name to split heatmap by
#' @param title Title for the heatmap
#' @param output_file Path to save the heatmap
#' @param width Width of the output file
#' @param height Height of the output file
#' @return ComplexHeatmap object
create_enhanced_heatmap <- function(seu_obj, genes_list, group_by = "treatment", 
                                   split_by = "compartment", title = "Gene Expression Heatmap",
                                   output_file = NULL, width = 12, height = 8) {
  
  cat("Creating enhanced heatmap...\n")
  
  # Filter genes present in the dataset
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in dataset\n")
    return(NULL)
  }
  
  # Get expression data
  expr_data <- GetAssayData(seu_obj, slot = "data")[genes_present, , drop = FALSE]
  
  # Get metadata
  meta_data <- seu_obj@meta.data
  
  # Create sample averages for each group
  if (group_by %in% colnames(meta_data)) {
    
    # Create grouping factor
    if (!is.null(split_by) && split_by %in% colnames(meta_data)) {
      grouping_factor <- paste(meta_data[[group_by]], meta_data[[split_by]], sep = "_")
    } else {
      grouping_factor <- meta_data[[group_by]]
    }
    
    # Calculate average expression for each group
    avg_expr <- sapply(unique(grouping_factor), function(group) {
      cells_in_group <- names(grouping_factor)[grouping_factor == group]
      if (length(cells_in_group) > 1) {
        rowMeans(expr_data[, cells_in_group, drop = FALSE])
      } else if (length(cells_in_group) == 1) {
        expr_data[, cells_in_group]
      } else {
        rep(0, nrow(expr_data))
      }
    })
    
    # Ensure we have a matrix
    if (is.vector(avg_expr)) {
      avg_expr <- matrix(avg_expr, ncol = 1)
      colnames(avg_expr) <- unique(grouping_factor)
    }
    
    # Scale by row (z-score)
    avg_expr_scaled <- t(scale(t(avg_expr)))
    
    # Create annotation for columns
    col_annotation <- data.frame(
      Group = unique(grouping_factor),
      stringsAsFactors = FALSE
    )
    
    # Extract treatment and compartment information
    if (!is.null(split_by) && split_by %in% colnames(meta_data)) {
      col_annotation$Treatment <- sapply(strsplit(col_annotation$Group, "_"), "[", 1)
      col_annotation$Compartment <- sapply(strsplit(col_annotation$Group, "_"), "[", 2)
    } else {
      col_annotation$Treatment <- col_annotation$Group
    }
    
    rownames(col_annotation) <- col_annotation$Group
    
    # Create color schemes
    treatment_colors <- c("NAIVE" = "#1f77b4", "SHAM" = "#ff7f0e", "DMG" = "#d62728")
    compartment_colors <- c("WB" = "#2ca02c", "BM" = "#9467bd")
    
    # Create annotation object
    if (!is.null(split_by) && split_by %in% colnames(meta_data)) {
      col_ann <- HeatmapAnnotation(
        Treatment = col_annotation$Treatment,
        Compartment = col_annotation$Compartment,
        col = list(
          Treatment = treatment_colors,
          Compartment = compartment_colors
        )
      )
    } else {
      col_ann <- HeatmapAnnotation(
        Treatment = col_annotation$Treatment,
        col = list(Treatment = treatment_colors)
      )
    }
    
    # Create row annotation for gene categories
    row_annotation <- data.frame(
      Gene = rownames(avg_expr_scaled),
      Category = "Other",
      stringsAsFactors = FALSE
    )
    
    # Define gene categories
    s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
    trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
    activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
    
    for (i in 1:nrow(row_annotation)) {
      gene <- row_annotation$Gene[i]
      if (gene %in% s1p_genes) {
        row_annotation$Category[i] <- "S1P_Pathway"
      } else if (gene %in% trafficking_genes) {
        row_annotation$Category[i] <- "Trafficking"
      } else if (gene %in% activation_genes) {
        row_annotation$Category[i] <- "Activation"
      }
    }
    
    rownames(row_annotation) <- row_annotation$Gene
    
    # Create row annotation object
    category_colors <- c("S1P_Pathway" = "#e377c2", "Trafficking" = "#8c564b", 
                        "Activation" = "#17becf", "Other" = "#7f7f7f")
    
    row_ann <- rowAnnotation(
      Category = row_annotation$Category,
      col = list(Category = category_colors)
    )
    
    # Create the heatmap
    ht <- Heatmap(
      avg_expr_scaled,
      name = "Expression\n(Z-score)",
      col = colorRamp2(c(-2, 0, 2), c("blue", "white", "red")),
      top_annotation = col_ann,
      left_annotation = row_ann,
      row_names_gp = gpar(fontsize = 8),
      column_names_gp = gpar(fontsize = 8),
      column_title = title,
      column_title_gp = gpar(fontsize = 12, fontface = "bold"),
      clustering_distance_rows = "pearson",
      clustering_distance_columns = "pearson",
      show_row_dend = TRUE,
      show_column_dend = TRUE,
      row_names_side = "left",
      column_names_side = "bottom",
      heatmap_legend_param = list(
        title_position = "topcenter",
        legend_direction = "horizontal",
        legend_width = unit(4, "cm")
      )
    )
    
    # Save if output file is specified
    if (!is.null(output_file)) {
      pdf(output_file, width = width, height = height)
      draw(ht)
      dev.off()
      cat(paste("Heatmap saved to:", output_file, "\n"))
    }
    
    return(ht)
    
  } else {
    cat(paste("Column", group_by, "not found in metadata\n"))
    return(NULL)
  }
}

#' Create pathway-specific heatmaps
#'
#' @param seu_obj Seurat object
#' @param output_dir Directory to save heatmaps
#' @return None (saves heatmaps to files)
create_pathway_heatmaps <- function(seu_obj, output_dir = ".") {
  
  cat("Creating pathway-specific heatmaps...\n")
  
  # Define gene sets
  s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
  trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
  retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")
  activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
  
  # Create S1P pathway heatmap
  create_enhanced_heatmap(seu_obj, s1p_genes, 
                         title = "S1P Pathway Genes",
                         output_file = file.path(output_dir, "S1P_Pathway_Heatmap.pdf"))
  
  # Create trafficking genes heatmap
  create_enhanced_heatmap(seu_obj, trafficking_genes,
                         title = "Trafficking Genes", 
                         output_file = file.path(output_dir, "Trafficking_Genes_Heatmap.pdf"))
  
  # Create retention genes heatmap
  create_enhanced_heatmap(seu_obj, retention_genes,
                         title = "Bone Marrow Retention Genes",
                         output_file = file.path(output_dir, "Retention_Genes_Heatmap.pdf"))
  
  # Create activation genes heatmap
  create_enhanced_heatmap(seu_obj, activation_genes,
                         title = "T-cell Activation Genes",
                         output_file = file.path(output_dir, "Activation_Genes_Heatmap.pdf"))
  
  cat("Pathway heatmaps created successfully!\n")
}

# ============================================================================
# VOLCANO PLOT FUNCTIONS
# ============================================================================

#' Create enhanced volcano plots
#'
#' @param dge_results DGE results data frame
#' @param title Title for the plot
#' @param output_file Path to save the plot
#' @param highlight_genes Vector of genes to highlight
#' @param width Width of the output file
#' @param height Height of the output file
#' @return ggplot object
create_enhanced_volcano <- function(dge_results, title = "Volcano Plot", 
                                   output_file = NULL, highlight_genes = NULL,
                                   width = 10, height = 8) {
  
  cat("Creating enhanced volcano plot...\n")
  
  if (nrow(dge_results) == 0) {
    cat("No data provided for volcano plot\n")
    return(NULL)
  }
  
  # Ensure required columns exist
  if (!all(c("avg_log2FC", "p_val", "gene") %in% colnames(dge_results))) {
    cat("Required columns missing from DGE results\n")
    return(NULL)
  }
  
  # Create significance categories
  dge_results$significance <- "NS"
  dge_results$significance[dge_results$p_val < 0.05 & abs(dge_results$avg_log2FC) > 0.25] <- "Significant"
  dge_results$significance[dge_results$p_val < 0.001 & abs(dge_results$avg_log2FC) > 0.5] <- "Highly Significant"
  
  # Highlight specific genes if provided
  if (!is.null(highlight_genes)) {
    s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
    trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
    
    dge_results$significance[dge_results$gene %in% s1p_genes] <- "S1P_Genes"
    dge_results$significance[dge_results$gene %in% trafficking_genes] <- "Trafficking_Genes"
  }
  
  # Create volcano plot
  p <- ggplot(dge_results, aes(x = avg_log2FC, y = -log10(p_val), color = significance)) +
    geom_point(alpha = 0.6, size = 1.2) +
    scale_color_manual(values = c(
      "NS" = "grey70",
      "Significant" = "black",
      "Highly Significant" = "darkred",
      "S1P_Genes" = "red",
      "Trafficking_Genes" = "blue"
    )) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", alpha = 0.7) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", alpha = 0.7) +
    labs(
      title = title,
      x = "Average log2 Fold Change",
      y = "-log10(p-value)",
      color = "Significance"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      legend.title = element_text(size = 12),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    )
  
  # Add gene labels for top genes
  if (!is.null(highlight_genes)) {
    top_genes <- dge_results[dge_results$gene %in% highlight_genes, ]
    if (nrow(top_genes) > 0) {
      p <- p + geom_text_repel(
        data = top_genes,
        aes(label = gene),
        size = 3,
        max.overlaps = 20,
        box.padding = 0.3,
        point.padding = 0.3,
        color = "black"
      )
    }
  }
  
  # Save if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, dpi = 300)
    cat(paste("Volcano plot saved to:", output_file, "\n"))
  }
  
  return(p)
}

#' Create volcano plots for all contrasts and methods
#'
#' @param dge_results_list List of DGE results
#' @param output_dir Directory to save plots
#' @return None (saves plots to files)
create_all_volcano_plots <- function(dge_results_list, output_dir = ".") {
  
  cat("Creating volcano plots for all contrasts and methods...\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define genes to highlight
  highlight_genes <- c("S1pr1", "S1pr4", "S1pr5", "Cxcr4", "Sell", "Ccr7", "Cd69", "Klf2")
  
  # Create volcano plots for each method and contrast
  for (method_name in names(dge_results_list)) {
    method_results <- dge_results_list[[method_name]]
    
    if (is.list(method_results)) {
      for (contrast_name in names(method_results)) {
        contrast_data <- method_results[[contrast_name]]
        
        if (is.data.frame(contrast_data) && nrow(contrast_data) > 0) {
          title <- paste("Volcano Plot:", contrast_name, "-", method_name)
          output_file <- file.path(output_dir, paste0("Volcano_", method_name, "_", contrast_name, ".pdf"))
          
          create_enhanced_volcano(contrast_data, title = title, 
                                 output_file = output_file, 
                                 highlight_genes = highlight_genes)
        }
      }
    }
  }
  
  cat("All volcano plots created successfully!\n")
}

# ============================================================================
# NEBULOSA DENSITY PLOT FUNCTIONS
# ============================================================================

#' Create enhanced Nebulosa density plots
#'
#' @param seu_obj Seurat object
#' @param genes_list Vector of genes to plot
#' @param title_prefix Prefix for plot titles
#' @param output_file Path to save the plot
#' @param width Width of the output file
#' @param height Height of the output file
#' @return List of ggplot objects
create_nebulosa_density_plots <- function(seu_obj, genes_list, title_prefix = "", 
                                         output_file = NULL, width = 15, height = 10) {
  
  cat("Creating Nebulosa density plots...\n")
  
  # Filter genes present in the dataset
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in dataset\n")
    return(NULL)
  }
  
  plots <- list()
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }
  
  for (gene in genes_present) {
    tryCatch({
      # Overall density plot
      p1 <- plot_density(seu_obj, features = gene, reduction = "umap", size = 0.8, pal = "viridis") +
        ggtitle(paste(title_prefix, gene, "- Overall Density")) +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
      
      # Density plot by treatment
      p2 <- plot_density(seu_obj, features = gene, reduction = "umap", size = 0.6, pal = "viridis") +
        facet_wrap(~treatment, ncol = 3) +
        ggtitle(paste(title_prefix, gene, "- By Treatment")) +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
      
      # Density plot by compartment if available
      if ("compartment" %in% colnames(seu_obj@meta.data)) {
        p3 <- plot_density(seu_obj, features = gene, reduction = "umap", size = 0.6, pal = "viridis") +
          facet_wrap(~compartment, ncol = 2) +
          ggtitle(paste(title_prefix, gene, "- By Compartment")) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
        
        # Combined plot
        combined_plot <- (p1 | p3) / p2
      } else {
        combined_plot <- p1 / p2
      }
      
      plots[[gene]] <- combined_plot
      
      if (!is.null(output_file)) {
        print(combined_plot)
      }
      
    }, error = function(e) {
      cat(paste("Error creating Nebulosa plot for", gene, ":", e$message, "\n"))
    })
  }
  
  if (!is.null(output_file)) {
    dev.off()
    cat(paste("Nebulosa plots saved to:", output_file, "\n"))
  }
  
  return(plots)
}

#' Create comprehensive Nebulosa plots for all gene sets
#'
#' @param seu_obj Seurat object
#' @param output_dir Directory to save plots
#' @return None (saves plots to files)
create_comprehensive_nebulosa_plots <- function(seu_obj, output_dir = ".") {
  
  cat("Creating comprehensive Nebulosa plots...\n")
  
  # Define gene sets
  s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
  trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
  retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")
  activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
  
  # Create plots for each gene set
  create_nebulosa_density_plots(seu_obj, s1p_genes, "S1P:", 
                               file.path(output_dir, "Nebulosa_S1P_Genes.pdf"))
  
  create_nebulosa_density_plots(seu_obj, trafficking_genes, "Trafficking:", 
                               file.path(output_dir, "Nebulosa_Trafficking_Genes.pdf"))
  
  create_nebulosa_density_plots(seu_obj, retention_genes, "Retention:", 
                               file.path(output_dir, "Nebulosa_Retention_Genes.pdf"))
  
  create_nebulosa_density_plots(seu_obj, activation_genes, "Activation:", 
                               file.path(output_dir, "Nebulosa_Activation_Genes.pdf"))
  
  cat("Comprehensive Nebulosa plots created successfully!\n")
}

# ============================================================================
# VIOLIN AND BOX PLOT FUNCTIONS
# ============================================================================

#' Create enhanced violin plots
#'
#' @param seu_obj Seurat object
#' @param genes_list Vector of genes to plot
#' @param group_by Column to group by
#' @param split_by Column to split by
#' @param title_prefix Prefix for plot titles
#' @param output_file Path to save the plot
#' @param width Width of the output file
#' @param height Height of the output file
#' @return List of ggplot objects
create_enhanced_violin_plots <- function(seu_obj, genes_list, group_by = "treatment", 
                                       split_by = "compartment", title_prefix = "",
                                       output_file = NULL, width = 12, height = 8) {
  
  cat("Creating enhanced violin plots...\n")
  
  # Filter genes present in the dataset
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in dataset\n")
    return(NULL)
  }
  
  plots <- list()
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }
  
  for (gene in genes_present) {
    tryCatch({
      # Basic violin plot by treatment
      p1 <- VlnPlot(seu_obj, features = gene, group.by = group_by, 
                   pt.size = 0.1, combine = FALSE)[[1]] +
        ggtitle(paste(title_prefix, gene, "by", group_by)) +
        theme_minimal() +
        theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
      
      # Violin plot by compartment if available
      if (!is.null(split_by) && split_by %in% colnames(seu_obj@meta.data)) {
        p2 <- VlnPlot(seu_obj, features = gene, group.by = split_by, 
                     pt.size = 0.1, combine = FALSE)[[1]] +
          ggtitle(paste(title_prefix, gene, "by", split_by)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
        
        # Combined violin plot
        p3 <- VlnPlot(seu_obj, features = gene, group.by = group_by, split.by = split_by,
                     pt.size = 0.1, combine = FALSE)[[1]] +
          ggtitle(paste(title_prefix, gene, "by", group_by, "and", split_by)) +
          theme_minimal() +
          theme(plot.title = element_text(size = 12, hjust = 0.5, face = "bold"))
        
        combined_plot <- (p1 | p2) / p3
      } else {
        combined_plot <- p1
      }
      
      plots[[gene]] <- combined_plot
      
      if (!is.null(output_file)) {
        print(combined_plot)
      }
      
    }, error = function(e) {
      cat(paste("Error creating violin plot for", gene, ":", e$message, "\n"))
    })
  }
  
  if (!is.null(output_file)) {
    dev.off()
    cat(paste("Violin plots saved to:", output_file, "\n"))
  }
  
  return(plots)
}

# ============================================================================
# RIDGE PLOT FUNCTIONS
# ============================================================================

#' Create ridge plots for gene expression
#'
#' @param seu_obj Seurat object
#' @param genes_list Vector of genes to plot
#' @param group_by Column to group by
#' @param title_prefix Prefix for plot titles
#' @param output_file Path to save the plot
#' @param width Width of the output file
#' @param height Height of the output file
#' @return ggplot object
create_ridge_plots <- function(seu_obj, genes_list, group_by = "treatment", 
                              title_prefix = "", output_file = NULL, 
                              width = 12, height = 8) {
  
  cat("Creating ridge plots...\n")
  
  # Filter genes present in the dataset
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in dataset\n")
    return(NULL)
  }
  
  plots <- list()
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }
  
  for (gene in genes_present) {
    tryCatch({
      # Get expression data
      expr_data <- GetAssayData(seu_obj, slot = "data")[gene, ]
      
      # Create data frame for plotting
      plot_data <- data.frame(
        expression = expr_data,
        group = seu_obj@meta.data[[group_by]],
        stringsAsFactors = FALSE
      )
      
      # Create ridge plot
      p <- ggplot(plot_data, aes(x = expression, y = group, fill = group)) +
        geom_density_ridges(alpha = 0.7, scale = 0.9) +
        scale_fill_viridis_d() +
        labs(
          title = paste(title_prefix, gene, "Expression Distribution"),
          x = "Expression Level",
          y = group_by,
          fill = group_by
        ) +
        theme_minimal() +
        theme(
          plot.title = element_text(size = 12, hjust = 0.5, face = "bold"),
          legend.position = "none"
        )
      
      plots[[gene]] <- p
      
      if (!is.null(output_file)) {
        print(p)
      }
      
    }, error = function(e) {
      cat(paste("Error creating ridge plot for", gene, ":", e$message, "\n"))
    })
  }
  
  if (!is.null(output_file)) {
    dev.off()
    cat(paste("Ridge plots saved to:", output_file, "\n"))
  }
  
  return(plots)
}

# ============================================================================
# CORRELATION PLOT FUNCTIONS
# ============================================================================

#' Create correlation plots between S1P receptors and trafficking markers
#'
#' @param seu_obj Seurat object
#' @param output_file Path to save the plot
#' @param width Width of the output file
#' @param height Height of the output file
#' @return ggplot object
create_correlation_plots <- function(seu_obj, output_file = NULL, width = 12, height = 10) {
  
  cat("Creating correlation plots...\n")
  
  # Define gene sets
  s1p_genes <- c("S1pr1", "S1pr4", "S1pr5")
  trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Cd69", "Klf2")
  
  # Filter genes present in the dataset
  all_genes <- unique(c(s1p_genes, trafficking_genes))
  genes_present <- all_genes[all_genes %in% rownames(seu_obj)]
  
  if (length(genes_present) < 2) {
    cat("Not enough genes found for correlation analysis\n")
    return(NULL)
  }
  
  # Get expression data
  expr_data <- GetAssayData(seu_obj, slot = "data")[genes_present, ]
  
  # Calculate correlation matrix
  cor_matrix <- cor(t(expr_data), use = "complete.obs")
  
  # Convert to long format for plotting
  cor_long <- expand.grid(Gene1 = rownames(cor_matrix), Gene2 = colnames(cor_matrix))
  cor_long$Correlation <- as.vector(cor_matrix)
  
  # Create correlation heatmap
  p <- ggplot(cor_long, aes(x = Gene1, y = Gene2, fill = Correlation)) +
    geom_tile(color = "white", size = 0.5) +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                        midpoint = 0, limit = c(-1, 1), space = "Lab",
                        name = "Pearson\nCorrelation") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1),
      axis.text.y = element_text(size = 10),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      legend.position = "bottom"
    ) +
    labs(title = "Correlation between S1P Receptors and Trafficking Markers") +
    coord_fixed()
  
  # Save if output file is specified
  if (!is.null(output_file)) {
    ggsave(output_file, p, width = width, height = height, dpi = 300)
    cat(paste("Correlation plot saved to:", output_file, "\n"))
  }
  
  return(p)
}

# ============================================================================
# FEATURE PLOT FUNCTIONS
# ============================================================================

#' Create enhanced feature plots
#'
#' @param seu_obj Seurat object
#' @param genes_list Vector of genes to plot
#' @param output_file Path to save the plot
#' @param width Width of the output file
#' @param height Height of the output file
#' @return List of ggplot objects
create_enhanced_feature_plots <- function(seu_obj, genes_list, output_file = NULL, 
                                         width = 15, height = 10) {
  
  cat("Creating enhanced feature plots...\n")
  
  # Filter genes present in the dataset
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if (length(genes_present) == 0) {
    cat("No genes found in dataset\n")
    return(NULL)
  }
  
  plots <- list()
  
  if (!is.null(output_file)) {
    pdf(output_file, width = width, height = height)
  }
  
  for (gene in genes_present) {
    tryCatch({
      # Standard feature plot
      p1 <- FeaturePlot(seu_obj, features = gene, reduction = "umap", 
                       pt.size = 0.5, order = TRUE) +
        ggtitle(paste("Feature Plot:", gene)) +
        theme_minimal()
      
      # Feature plot split by treatment
      p2 <- FeaturePlot(seu_obj, features = gene, reduction = "umap", 
                       pt.size = 0.5, order = TRUE, split.by = "treatment") +
        ggtitle(paste("Feature Plot:", gene, "by Treatment")) +
        theme_minimal()
      
      # Feature plot split by compartment if available
      if ("compartment" %in% colnames(seu_obj@meta.data)) {
        p3 <- FeaturePlot(seu_obj, features = gene, reduction = "umap", 
                         pt.size = 0.5, order = TRUE, split.by = "compartment") +
          ggtitle(paste("Feature Plot:", gene, "by Compartment")) +
          theme_minimal()
        
        combined_plot <- p1 / p2 / p3
      } else {
        combined_plot <- p1 / p2
      }
      
      plots[[gene]] <- combined_plot
      
      if (!is.null(output_file)) {
        print(combined_plot)
      }
      
    }, error = function(e) {
      cat(paste("Error creating feature plot for", gene, ":", e$message, "\n"))
    })
  }
  
  if (!is.null(output_file)) {
    dev.off()
    cat(paste("Feature plots saved to:", output_file, "\n"))
  }
  
  return(plots)
}

# ============================================================================
# MAIN VISUALIZATION PIPELINE
# ============================================================================

#' Run comprehensive visualization pipeline
#'
#' @param seu_obj Seurat object
#' @param dge_results_list List of DGE results
#' @param output_dir Directory to save visualizations
#' @return None (saves visualizations to files)
run_comprehensive_visualization_pipeline <- function(seu_obj, dge_results_list = NULL, 
                                                    output_dir = "visualizations") {
  
  cat("============================================================================\n")
  cat("STARTING COMPREHENSIVE VISUALIZATION PIPELINE\n")
  cat("============================================================================\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Define gene sets
  s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
  trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
  all_goi <- unique(c(s1p_genes, trafficking_genes))
  
  # Create pathway-specific heatmaps
  cat("\n--- Creating Pathway Heatmaps ---\n")
  create_pathway_heatmaps(seu_obj, output_dir)
  
  # Create volcano plots if DGE results provided
  if (!is.null(dge_results_list)) {
    cat("\n--- Creating Volcano Plots ---\n")
    create_all_volcano_plots(dge_results_list, output_dir)
  }
  
  # Create Nebulosa plots
  cat("\n--- Creating Nebulosa Plots ---\n")
  create_comprehensive_nebulosa_plots(seu_obj, output_dir)
  
  # Create violin plots
  cat("\n--- Creating Violin Plots ---\n")
  create_enhanced_violin_plots(seu_obj, s1p_genes, title_prefix = "S1P:", 
                              output_file = file.path(output_dir, "Violin_S1P_Genes.pdf"))
  
  create_enhanced_violin_plots(seu_obj, trafficking_genes, title_prefix = "Trafficking:", 
                              output_file = file.path(output_dir, "Violin_Trafficking_Genes.pdf"))
  
  # Create ridge plots
  cat("\n--- Creating Ridge Plots ---\n")
  create_ridge_plots(seu_obj, s1p_genes, title_prefix = "S1P:", 
                    output_file = file.path(output_dir, "Ridge_S1P_Genes.pdf"))
  
  # Create correlation plots
  cat("\n--- Creating Correlation Plots ---\n")
  create_correlation_plots(seu_obj, output_file = file.path(output_dir, "Correlation_S1P_Trafficking.pdf"))
  
  # Create feature plots
  cat("\n--- Creating Feature Plots ---\n")
  create_enhanced_feature_plots(seu_obj, s1p_genes, 
                               output_file = file.path(output_dir, "Feature_S1P_Genes.pdf"))
  
  cat("\n============================================================================\n")
  cat("COMPREHENSIVE VISUALIZATION PIPELINE COMPLETE\n")
  cat("============================================================================\n")
  
  cat(paste("All visualizations saved to:", output_dir, "\n"))
}

cat("Enhanced visualization module loaded successfully!\n")
cat("Use run_comprehensive_visualization_pipeline() to execute the complete visualization suite.\n")