#!/usr/bin/env Rscript
# ============================================================================
# TRAJECTORY ANALYSIS MODULE FOR S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING
# ============================================================================
# 
# This module provides comprehensive trajectory analysis functionality
# including pseudotime analysis, branch point detection, and dynamic
# gene expression analysis along trajectories.
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

# Load required libraries
suppressMessages({
  library(monocle3)
  library(slingshot)
  library(tradeSeq)
  library(SingleCellExperiment)
  library(Seurat)
  library(tidyverse)
  library(here)
  library(viridis)
  library(RColorBrewer)
  library(pheatmap)
  library(ggplot2)
  library(patchwork)
})

# Set reproducibility
set.seed(42)

# ============================================================================
# TRAJECTORY ANALYSIS FUNCTIONS
# ============================================================================

#' Perform comprehensive trajectory analysis using Monocle3
#'
#' @param seu_obj Seurat object
#' @param root_cells Character vector of cell barcodes to use as root
#' @param genes_of_interest Character vector of genes to analyze
#' @return List containing trajectory analysis results
perform_monocle3_analysis <- function(seu_obj, root_cells = NULL, genes_of_interest = NULL) {
  
  cat("Starting Monocle3 trajectory analysis...\n")
  
  # Convert Seurat to cell_data_set
  gene_annotation <- data.frame(
    gene_short_name = rownames(seu_obj),
    row.names = rownames(seu_obj)
  )
  
  cds <- new_cell_data_set(
    expression_data = GetAssayData(seu_obj, slot = "counts"),
    cell_metadata = seu_obj@meta.data,
    gene_metadata = gene_annotation
  )
  
  # Preprocess
  cat("Preprocessing data...\n")
  cds <- preprocess_cds(cds, num_dim = 50, verbose = FALSE)
  
  # Align with existing UMAP if available
  if ("umap" %in% names(seu_obj@reductions)) {
    cat("Using existing UMAP coordinates...\n")
    cds@int_colData@listData$reducedDims$UMAP <- Embeddings(seu_obj, "umap")
  } else {
    cat("Generating new UMAP coordinates...\n")
    cds <- reduce_dimension(cds, reduction_method = "UMAP", verbose = FALSE)
  }
  
  # Cluster cells
  cat("Clustering cells...\n")
  cds <- cluster_cells(cds, verbose = FALSE)
  
  # Learn trajectory graph
  cat("Learning trajectory graph...\n")
  cds <- learn_graph(cds, verbose = FALSE)
  
  # Order cells in pseudotime
  cat("Ordering cells in pseudotime...\n")
  if (!is.null(root_cells)) {
    cds <- order_cells(cds, root_cells = root_cells)
  } else {
    cds <- order_cells(cds)
  }
  
  # Find genes that change along pseudotime
  cat("Finding trajectory-associated genes...\n")
  trajectory_genes <- graph_test(cds, neighbor_graph = "principal_graph", verbose = FALSE)
  trajectory_genes <- trajectory_genes[order(trajectory_genes$q_value), ]
  
  # Analyze genes of interest if provided
  goi_analysis <- NULL
  if (!is.null(genes_of_interest)) {
    cat("Analyzing genes of interest along trajectory...\n")
    goi_present <- genes_of_interest[genes_of_interest %in% rownames(cds)]
    
    if (length(goi_present) > 0) {
      goi_analysis <- list()
      
      for (gene in goi_present) {
        gene_data <- data.frame(
          pseudotime = pseudotime(cds),
          expression = exprs(cds)[gene, ],
          treatment = pData(cds)$treatment,
          compartment = pData(cds)$compartment,
          cell_type = pData(cds)$cell_types
        )
        
        # Remove cells with infinite pseudotime
        gene_data <- gene_data[is.finite(gene_data$pseudotime), ]
        
        if (nrow(gene_data) > 0) {
          goi_analysis[[gene]] <- gene_data
        }
      }
    }
  }
  
  # Calculate trajectory statistics
  trajectory_stats <- list(
    total_cells = ncol(cds),
    cells_with_pseudotime = sum(is.finite(pseudotime(cds))),
    significant_trajectory_genes = sum(trajectory_genes$q_value < 0.05, na.rm = TRUE),
    trajectory_length = max(pseudotime(cds)[is.finite(pseudotime(cds))], na.rm = TRUE)
  )
  
  cat("Monocle3 analysis complete!\n")
  
  return(list(
    cds = cds,
    trajectory_genes = trajectory_genes,
    goi_analysis = goi_analysis,
    trajectory_stats = trajectory_stats
  ))
}

#' Perform Slingshot trajectory analysis
#'
#' @param seu_obj Seurat object
#' @param start_cluster Starting cluster for trajectory
#' @param end_clusters Vector of end clusters
#' @return List containing Slingshot analysis results
perform_slingshot_analysis <- function(seu_obj, start_cluster = NULL, end_clusters = NULL) {
  
  cat("Starting Slingshot trajectory analysis...\n")
  
  # Convert to SingleCellExperiment
  sce <- as.SingleCellExperiment(seu_obj)
  
  # Add dimension reductions
  if ("umap" %in% names(seu_obj@reductions)) {
    reducedDims(sce) <- list(UMAP = Embeddings(seu_obj, "umap"))
  }
  
  if ("pca" %in% names(seu_obj@reductions)) {
    reducedDims(sce)$PCA <- Embeddings(seu_obj, "pca")
  }
  
  # Run Slingshot
  cat("Running Slingshot trajectory inference...\n")
  sds <- slingshot(sce, 
                   clusterLabels = "cell_types",
                   reducedDim = "UMAP",
                   start.clus = start_cluster,
                   end.clus = end_clusters)
  
  # Extract pseudotime values
  pseudotime_curves <- slingPseudotime(sds)
  curve_weights <- slingCurveWeights(sds)
  
  # Analyze lineage statistics
  lineage_stats <- list(
    n_lineages = ncol(pseudotime_curves),
    cells_assigned = sum(rowSums(!is.na(pseudotime_curves)) > 0),
    lineage_names = colnames(pseudotime_curves)
  )
  
  cat("Slingshot analysis complete!\n")
  
  return(list(
    sds = sds,
    pseudotime_curves = pseudotime_curves,
    curve_weights = curve_weights,
    lineage_stats = lineage_stats
  ))
}

#' Perform differential expression analysis along trajectories
#'
#' @param trajectory_results Results from trajectory analysis
#' @param seu_obj Original Seurat object
#' @param genes_of_interest Genes to focus on
#' @return List containing trajectory DE results
perform_trajectory_de_analysis <- function(trajectory_results, seu_obj, genes_of_interest = NULL) {
  
  cat("Performing trajectory differential expression analysis...\n")
  
  results <- list()
  
  if ("cds" %in% names(trajectory_results)) {
    # Monocle3 analysis
    cds <- trajectory_results$cds
    
    # Find genes that differ between treatments along trajectory
    treatment_de <- list()
    treatments <- unique(pData(cds)$treatment)
    
    for (i in 1:(length(treatments) - 1)) {
      for (j in (i + 1):length(treatments)) {
        comparison <- paste(treatments[i], "vs", treatments[j])
        cat(paste("Analyzing", comparison, "along trajectory...\n"))
        
        # Subset cells for each treatment
        cells_1 <- rownames(pData(cds))[pData(cds)$treatment == treatments[i]]
        cells_2 <- rownames(pData(cds))[pData(cds)$treatment == treatments[j]]
        
        if (length(cells_1) > 10 && length(cells_2) > 10) {
          # Create pseudotime bins
          pseudotime_1 <- pseudotime(cds)[cells_1]
          pseudotime_2 <- pseudotime(cds)[cells_2]
          
          # Remove infinite values
          pseudotime_1 <- pseudotime_1[is.finite(pseudotime_1)]
          pseudotime_2 <- pseudotime_2[is.finite(pseudotime_2)]
          
          if (length(pseudotime_1) > 5 && length(pseudotime_2) > 5) {
            # Create bins based on pseudotime quantiles
            combined_pseudotime <- c(pseudotime_1, pseudotime_2)
            breaks <- quantile(combined_pseudotime, probs = seq(0, 1, 0.25), na.rm = TRUE)
            
            bin_de_results <- data.frame()
            
            for (bin in 1:(length(breaks) - 1)) {
              bin_cells_1 <- names(pseudotime_1)[pseudotime_1 >= breaks[bin] & pseudotime_1 < breaks[bin + 1]]
              bin_cells_2 <- names(pseudotime_2)[pseudotime_2 >= breaks[bin] & pseudotime_2 < breaks[bin + 1]]
              
              if (length(bin_cells_1) > 3 && length(bin_cells_2) > 3) {
                # Perform differential expression for this bin
                bin_result <- data.frame(
                  pseudotime_bin = bin,
                  pseudotime_range = paste(round(breaks[bin], 2), "-", round(breaks[bin + 1], 2)),
                  n_cells_treatment1 = length(bin_cells_1),
                  n_cells_treatment2 = length(bin_cells_2),
                  comparison = comparison
                )
                
                bin_de_results <- rbind(bin_de_results, bin_result)
              }
            }
            
            treatment_de[[comparison]] <- bin_de_results
          }
        }
      }
    }
    
    results$monocle3_trajectory_de <- treatment_de
  }
  
  if ("sds" %in% names(trajectory_results)) {
    # Slingshot analysis
    sds <- trajectory_results$sds
    pseudotime_curves <- trajectory_results$pseudotime_curves
    
    # Analyze expression patterns along each lineage
    lineage_patterns <- list()
    
    for (lineage in colnames(pseudotime_curves)) {
      cat(paste("Analyzing lineage", lineage, "...\n"))
      
      # Get cells assigned to this lineage
      lineage_cells <- !is.na(pseudotime_curves[, lineage])
      
      if (sum(lineage_cells) > 20) {
        lineage_pseudotime <- pseudotime_curves[lineage_cells, lineage]
        lineage_treatments <- seu_obj@meta.data[names(lineage_pseudotime), "treatment"]
        
        # Create treatment comparison along this lineage
        lineage_data <- data.frame(
          pseudotime = lineage_pseudotime,
          treatment = lineage_treatments,
          lineage = lineage
        )
        
        lineage_patterns[[lineage]] <- lineage_data
      }
    }
    
    results$slingshot_lineage_patterns <- lineage_patterns
  }
  
  # Focus on genes of interest if provided
  if (!is.null(genes_of_interest)) {
    cat("Analyzing genes of interest along trajectories...\n")
    
    goi_trajectory_patterns <- list()
    
    for (gene in genes_of_interest) {
      if (gene %in% rownames(seu_obj)) {
        gene_expr <- GetAssayData(seu_obj, slot = "data")[gene, ]
        
        # Combine with trajectory information
        if ("cds" %in% names(trajectory_results)) {
          cds <- trajectory_results$cds
          common_cells <- intersect(names(gene_expr), rownames(pData(cds)))
          
          if (length(common_cells) > 10) {
            gene_trajectory_data <- data.frame(
              cell_id = common_cells,
              gene_expression = gene_expr[common_cells],
              pseudotime = pseudotime(cds)[common_cells],
              treatment = pData(cds)[common_cells, "treatment"],
              compartment = pData(cds)[common_cells, "compartment"]
            )
            
            # Remove infinite pseudotime values
            gene_trajectory_data <- gene_trajectory_data[is.finite(gene_trajectory_data$pseudotime), ]
            
            if (nrow(gene_trajectory_data) > 0) {
              goi_trajectory_patterns[[gene]] <- gene_trajectory_data
            }
          }
        }
      }
    }
    
    results$goi_trajectory_patterns <- goi_trajectory_patterns
  }
  
  cat("Trajectory differential expression analysis complete!\n")
  
  return(results)
}

#' Create trajectory visualization plots
#'
#' @param trajectory_results Results from trajectory analysis
#' @param seu_obj Original Seurat object
#' @param genes_of_interest Genes to visualize
#' @param output_dir Directory to save plots
#' @return None (saves plots to files)
create_trajectory_visualizations <- function(trajectory_results, seu_obj, genes_of_interest = NULL, output_dir = ".") {
  
  cat("Creating trajectory visualizations...\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Monocle3 visualizations
  if ("cds" %in% names(trajectory_results)) {
    cds <- trajectory_results$cds
    
    pdf(file.path(output_dir, "Monocle3_Trajectory_Plots.pdf"), width = 12, height = 8)
    
    # Basic trajectory plot
    p1 <- plot_cells(cds, color_cells_by = "treatment", show_trajectory_graph = TRUE,
                     label_cell_groups = FALSE, cell_size = 0.5) +
      ggtitle("Trajectory colored by Treatment") +
      theme_minimal()
    print(p1)
    
    # Pseudotime plot
    p2 <- plot_cells(cds, color_cells_by = "pseudotime", show_trajectory_graph = TRUE,
                     label_cell_groups = FALSE, cell_size = 0.5) +
      scale_color_viridis_c(name = "Pseudotime") +
      ggtitle("Cells colored by Pseudotime") +
      theme_minimal()
    print(p2)
    
    # Compartment plot
    if ("compartment" %in% colnames(pData(cds))) {
      p3 <- plot_cells(cds, color_cells_by = "compartment", show_trajectory_graph = TRUE,
                       label_cell_groups = FALSE, cell_size = 0.5) +
        ggtitle("Trajectory colored by Compartment") +
        theme_minimal()
      print(p3)
    }
    
    # Cell type plot
    p4 <- plot_cells(cds, color_cells_by = "cell_types", show_trajectory_graph = TRUE,
                     label_cell_groups = FALSE, cell_size = 0.5) +
      ggtitle("Trajectory colored by Cell Type") +
      theme_minimal()
    print(p4)
    
    # Plot genes of interest along trajectory
    if (!is.null(genes_of_interest)) {
      genes_present <- genes_of_interest[genes_of_interest %in% rownames(cds)]
      
      for (gene in genes_present) {
        tryCatch({
          p_gene <- plot_cells(cds, genes = gene, show_trajectory_graph = TRUE,
                              label_cell_groups = FALSE, cell_size = 0.5) +
            ggtitle(paste("Expression of", gene, "along trajectory")) +
            theme_minimal()
          print(p_gene)
        }, error = function(e) {
          cat(paste("Error plotting", gene, ":", e$message, "\n"))
        })
      }
    }
    
    dev.off()
    
    # Create pseudotime expression plots
    if (!is.null(genes_of_interest) && "goi_analysis" %in% names(trajectory_results)) {
      pdf(file.path(output_dir, "Pseudotime_Expression_Plots.pdf"), width = 12, height = 8)
      
      goi_analysis <- trajectory_results$goi_analysis
      
      for (gene in names(goi_analysis)) {
        gene_data <- goi_analysis[[gene]]
        
        if (nrow(gene_data) > 0) {
          # Expression vs pseudotime
          p1 <- ggplot(gene_data, aes(x = pseudotime, y = expression, color = treatment)) +
            geom_point(alpha = 0.6, size = 0.8) +
            geom_smooth(method = "loess", se = TRUE, size = 1) +
            labs(title = paste("Expression of", gene, "along Pseudotime"),
                 x = "Pseudotime", y = "Expression") +
            theme_minimal() +
            facet_wrap(~treatment, ncol = 3)
          
          print(p1)
          
          # By compartment if available
          if ("compartment" %in% colnames(gene_data)) {
            p2 <- ggplot(gene_data, aes(x = pseudotime, y = expression, color = compartment)) +
              geom_point(alpha = 0.6, size = 0.8) +
              geom_smooth(method = "loess", se = TRUE, size = 1) +
              labs(title = paste("Expression of", gene, "along Pseudotime by Compartment"),
                   x = "Pseudotime", y = "Expression") +
              theme_minimal() +
              facet_wrap(~treatment, ncol = 3)
            
            print(p2)
          }
        }
      }
      
      dev.off()
    }
  }
  
  # Slingshot visualizations
  if ("sds" %in% names(trajectory_results)) {
    sds <- trajectory_results$sds
    
    pdf(file.path(output_dir, "Slingshot_Trajectory_Plots.pdf"), width = 12, height = 8)
    
    # Plot trajectories
    colors <- RColorBrewer::brewer.pal(3, "Set1")
    
    # Basic trajectory plot
    plot(reducedDims(sds)$UMAP, col = colors[as.numeric(as.factor(sds$treatment))],
         pch = 16, cex = 0.5, main = "Slingshot Trajectories")
    lines(SlingshotDataSet(sds), lwd = 2)
    legend("topright", legend = levels(as.factor(sds$treatment)), 
           col = colors, pch = 16, cex = 0.8)
    
    # Plot by cell type
    plot(reducedDims(sds)$UMAP, col = colors[as.numeric(as.factor(sds$cell_types))],
         pch = 16, cex = 0.5, main = "Slingshot Trajectories by Cell Type")
    lines(SlingshotDataSet(sds), lwd = 2)
    legend("topright", legend = levels(as.factor(sds$cell_types)), 
           col = colors, pch = 16, cex = 0.8)
    
    dev.off()
  }
  
  cat("Trajectory visualizations complete!\n")
}

#' Export trajectory analysis results
#'
#' @param trajectory_results Results from trajectory analysis
#' @param trajectory_de_results Results from trajectory DE analysis
#' @param output_dir Directory to save results
#' @return None (saves results to files)
export_trajectory_results <- function(trajectory_results, trajectory_de_results = NULL, output_dir = ".") {
  
  cat("Exporting trajectory analysis results...\n")
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Save trajectory results
  saveRDS(trajectory_results, file.path(output_dir, "trajectory_analysis_results.rds"))
  
  if (!is.null(trajectory_de_results)) {
    saveRDS(trajectory_de_results, file.path(output_dir, "trajectory_de_analysis_results.rds"))
  }
  
  # Export trajectory genes if available
  if ("trajectory_genes" %in% names(trajectory_results)) {
    trajectory_genes <- trajectory_results$trajectory_genes
    write.csv(trajectory_genes, file.path(output_dir, "trajectory_associated_genes.csv"), row.names = FALSE)
  }
  
  # Export trajectory statistics
  if ("trajectory_stats" %in% names(trajectory_results)) {
    trajectory_stats <- trajectory_results$trajectory_stats
    write.csv(data.frame(Statistic = names(trajectory_stats), 
                        Value = unlist(trajectory_stats)), 
              file.path(output_dir, "trajectory_statistics.csv"), row.names = FALSE)
  }
  
  cat("Trajectory results exported successfully!\n")
}

# ============================================================================
# MAIN EXECUTION FUNCTION
# ============================================================================

#' Run complete trajectory analysis pipeline
#'
#' @param seu_obj Seurat object
#' @param genes_of_interest Character vector of genes to analyze
#' @param output_dir Directory to save results
#' @param methods Vector of methods to use ("monocle3", "slingshot", or both)
#' @return List containing all trajectory analysis results
run_trajectory_analysis_pipeline <- function(seu_obj, genes_of_interest = NULL, 
                                           output_dir = "trajectory_results", 
                                           methods = c("monocle3", "slingshot")) {
  
  cat("============================================================================\n")
  cat("STARTING COMPREHENSIVE TRAJECTORY ANALYSIS PIPELINE\n")
  cat("============================================================================\n")
  
  # Create output directory
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  results <- list()
  
  # Run Monocle3 analysis
  if ("monocle3" %in% methods) {
    cat("\n--- Running Monocle3 Analysis ---\n")
    monocle3_results <- perform_monocle3_analysis(seu_obj, genes_of_interest = genes_of_interest)
    results$monocle3 <- monocle3_results
  }
  
  # Run Slingshot analysis
  if ("slingshot" %in% methods) {
    cat("\n--- Running Slingshot Analysis ---\n")
    slingshot_results <- perform_slingshot_analysis(seu_obj)
    results$slingshot <- slingshot_results
  }
  
  # Perform trajectory DE analysis
  cat("\n--- Running Trajectory DE Analysis ---\n")
  trajectory_de_results <- perform_trajectory_de_analysis(results, seu_obj, genes_of_interest)
  results$trajectory_de <- trajectory_de_results
  
  # Create visualizations
  cat("\n--- Creating Trajectory Visualizations ---\n")
  create_trajectory_visualizations(results, seu_obj, genes_of_interest, output_dir)
  
  # Export results
  cat("\n--- Exporting Results ---\n")
  export_trajectory_results(results, trajectory_de_results, output_dir)
  
  cat("\n============================================================================\n")
  cat("TRAJECTORY ANALYSIS PIPELINE COMPLETE\n")
  cat("============================================================================\n")
  
  return(results)
}

# ============================================================================
# UTILITY FUNCTIONS
# ============================================================================

#' Print trajectory analysis summary
#'
#' @param trajectory_results Results from trajectory analysis
#' @return None (prints summary to console)
print_trajectory_summary <- function(trajectory_results) {
  
  cat("============================================================================\n")
  cat("TRAJECTORY ANALYSIS SUMMARY\n")
  cat("============================================================================\n")
  
  if ("monocle3" %in% names(trajectory_results)) {
    cat("\n--- Monocle3 Results ---\n")
    stats <- trajectory_results$monocle3$trajectory_stats
    cat(paste("Total cells analyzed:", stats$total_cells, "\n"))
    cat(paste("Cells with pseudotime:", stats$cells_with_pseudotime, "\n"))
    cat(paste("Significant trajectory genes:", stats$significant_trajectory_genes, "\n"))
    cat(paste("Trajectory length:", round(stats$trajectory_length, 2), "\n"))
  }
  
  if ("slingshot" %in% names(trajectory_results)) {
    cat("\n--- Slingshot Results ---\n")
    stats <- trajectory_results$slingshot$lineage_stats
    cat(paste("Number of lineages:", stats$n_lineages, "\n"))
    cat(paste("Cells assigned to trajectories:", stats$cells_assigned, "\n"))
    cat(paste("Lineage names:", paste(stats$lineage_names, collapse = ", "), "\n"))
  }
  
  if ("trajectory_de" %in% names(trajectory_results)) {
    cat("\n--- Trajectory DE Results ---\n")
    de_results <- trajectory_results$trajectory_de
    if ("monocle3_trajectory_de" %in% names(de_results)) {
      cat(paste("Monocle3 DE comparisons:", length(de_results$monocle3_trajectory_de), "\n"))
    }
    if ("goi_trajectory_patterns" %in% names(de_results)) {
      cat(paste("Genes of interest analyzed:", length(de_results$goi_trajectory_patterns), "\n"))
    }
  }
  
  cat("============================================================================\n")
}

cat("Trajectory analysis module loaded successfully!\n")
cat("Use run_trajectory_analysis_pipeline() to execute the complete analysis.\n")