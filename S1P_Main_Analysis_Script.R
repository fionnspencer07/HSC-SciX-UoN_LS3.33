#!/usr/bin/env Rscript
# ============================================================================
# COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS
# MAIN EXECUTION SCRIPT
# ============================================================================
# 
# This is the main execution script that integrates all analysis modules
# including DGE analysis, trajectory analysis, and comprehensive visualizations
# to address the research question about S1P receptor externalization and
# lymphocyte trafficking in DMG presence.
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

# Set up environment
library(here)
set.seed(42)

# Source required modules
source(here("S1P_Comprehensive_Analysis_Pipeline.R"))
source(here("S1P_Trajectory_Analysis_Module.R"))
source(here("S1P_Enhanced_Visualization_Module.R"))

# Additional packages for enhanced functionality
suppressMessages({
  library(openxlsx)
  library(patchwork)
  library(future)
  library(future.apply)
  library(progressr)
})

# Enable parallel processing
plan(multisession, workers = 4)

# ============================================================================
# MAIN EXECUTION FUNCTION
# ============================================================================

#' Execute the complete S1P receptor and lymphocyte trafficking analysis pipeline
#'
#' @param data_file Path to the RDS file containing the Seurat object
#' @param output_base_dir Base directory for all outputs
#' @param run_trajectory_analysis Logical, whether to run trajectory analysis
#' @param run_visualization Logical, whether to run visualization pipeline
#' @param cleanup_workspace Logical, whether to clean up workspace after analysis
#' @return List containing all analysis results
execute_comprehensive_s1p_analysis <- function(data_file = NULL, 
                                              output_base_dir = "S1P_Comprehensive_Results",
                                              run_trajectory_analysis = TRUE,
                                              run_visualization = TRUE,
                                              cleanup_workspace = FALSE) {
  
  cat("============================================================================\n")
  cat("COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS\n")
  cat("============================================================================\n")
  
  # Start timing
  start_time <- Sys.time()
  
  # ============================================================================
  # 1. DATA LOADING AND PREPARATION
  # ============================================================================
  
  cat("\n=== STEP 1: DATA LOADING AND PREPARATION ===\n")
  
  # Use default data file if not specified
  if (is.null(data_file)) {
    data_file <- here("Data", "seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
  }
  
  if (!file.exists(data_file)) {
    stop("Data file not found. Please ensure the RDS file exists at the specified path.")
  }
  
  # Load the Seurat object
  cat("Loading Seurat object...\n")
  seu <- readRDS(data_file)
  cat(paste("Loaded Seurat object with", ncol(seu), "cells and", nrow(seu), "genes\n"))
  
  # Subset to NK, T, and NKT cells
  seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))
  cat(paste("After cell type filtering:", ncol(seu_NKT), "cells\n"))
  
  # Focus on relevant treatments
  seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "NAIVE", "SHAM"))
  cat(paste("After treatment filtering:", ncol(seu_NKT_focused), "cells\n"))
  
  # Rename treatments for clarity
  seu_NKT_focused$treatment <- factor(seu_NKT_focused$treatment, 
                                     levels = c("NAIVE", "SHAM", "UT"),
                                     labels = c("NAIVE", "SHAM", "DMG"))
  
  # Set default assay
  DefaultAssay(seu_NKT_focused) <- "RNA"
  
  # Print data structure
  cat("Final data structure:\n")
  print(table(seu_NKT_focused$treatment, seu_NKT_focused$compartment))
  
  # ============================================================================
  # 2. GENE SETS AND CONTRASTS DEFINITION
  # ============================================================================
  
  cat("\n=== STEP 2: GENE SETS AND CONTRASTS DEFINITION ===\n")
  
  # Define comprehensive gene sets
  s1p_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", "Sphk1", "Sphk2", "Sgpl1")
  trafficking_genes <- c("Cxcr4", "Sell", "Ccr7", "Itgal", "Cd62l", "Psgl1", "Lfa1")
  retention_genes <- c("Cxcr4", "Vcam1", "Vla4", "Sdf1")
  activation_genes <- c("Cd69", "Cd25", "Cd44", "Cd95", "Klf2")
  chemokine_genes <- c("Ccl5", "Ccl4", "Cxcl10", "Cxcl9")
  cytokine_genes <- c("Ifng", "Il2", "Tnfa", "Il10")
  
  # Combined genes of interest
  all_goi <- unique(c(s1p_genes, trafficking_genes, retention_genes, 
                     activation_genes, chemokine_genes, cytokine_genes))
  
  # Check gene presence
  genes_present <- all_goi[all_goi %in% rownames(seu_NKT_focused)]
  genes_missing <- all_goi[!all_goi %in% rownames(seu_NKT_focused)]
  
  cat(paste("Total genes of interest:", length(all_goi), "\n"))
  cat(paste("Genes present in dataset:", length(genes_present), "\n"))
  cat(paste("Missing genes:", length(genes_missing), "\n"))
  
  if (length(genes_missing) > 0) {
    cat("Missing genes:\n")
    cat(paste(genes_missing, collapse = ", "), "\n")
  }
  
  # Define contrasts
  contrasts <- list(
    "DMG_vs_SHAM" = c("DMG", "SHAM"),
    "DMG_vs_NAIVE" = c("DMG", "NAIVE"),
    "SHAM_vs_NAIVE" = c("SHAM", "NAIVE")
  )
  
  cat(paste("Contrasts defined:", length(contrasts), "\n"))
  
  # ============================================================================
  # 3. CREATE OUTPUT DIRECTORIES
  # ============================================================================
  
  cat("\n=== STEP 3: CREATING OUTPUT DIRECTORIES ===\n")
  
  # Create comprehensive output structure
  if (!dir.exists(output_base_dir)) {
    dir.create(output_base_dir, recursive = TRUE)
  }
  
  # Create subdirectories
  dir.create(file.path(output_base_dir, "Excel_Outputs"), showWarnings = FALSE)
  dir.create(file.path(output_base_dir, "DGE_Results"), showWarnings = FALSE)
  dir.create(file.path(output_base_dir, "Trajectory_Analysis"), showWarnings = FALSE)
  dir.create(file.path(output_base_dir, "Visualizations"), showWarnings = FALSE)
  dir.create(file.path(output_base_dir, "Statistical_Results"), showWarnings = FALSE)
  dir.create(file.path(output_base_dir, "Summary_Reports"), showWarnings = FALSE)
  
  cat(paste("Output directories created in:", output_base_dir, "\n"))
  
  # ============================================================================
  # 4. DIFFERENTIAL GENE EXPRESSION ANALYSIS
  # ============================================================================
  
  cat("\n=== STEP 4: DIFFERENTIAL GENE EXPRESSION ANALYSIS ===\n")
  
  # Initialize results storage
  dge_results <- list()
  
  # 4.1 Wilcoxon Test Analysis
  cat("4.1 Running Wilcoxon test analysis...\n")
  
  with_progress({
    p <- progressor(along = seq_len(length(contrasts)))
    
    wilcox_results <- list()
    
    for (i in seq_along(contrasts)) {
      contrast <- contrasts[[i]]
      contrast_name <- names(contrasts)[i]
      
      p(sprintf("Processing contrast: %s", contrast_name))
      
      # Set identity to treatment
      Idents(seu_NKT_focused) <- "treatment"
      
      tryCatch({
        markers <- FindMarkers(seu_NKT_focused,
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
          
          wilcox_results[[contrast_name]] <- markers
          cat(paste("  Found", nrow(markers), "differential genes for", contrast_name, "\n"))
        }
      }, error = function(e) {
        cat(paste("  Error in", contrast_name, ":", e$message, "\n"))
      })
    }
  })
  
  dge_results$wilcoxon <- wilcox_results
  
  # 4.2 Muscat Analysis
  cat("4.2 Running Muscat pseudobulk analysis...\n")
  
  tryCatch({
    # Convert to SingleCellExperiment
    sce <- as.SingleCellExperiment(seu_NKT_focused)
    
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
    
    # Run differential expression
    muscat_results <- pbDS(pb, 
                          method = c("edgeR", "DESeq2"),
                          filter = "both",
                          min_cells = 10,
                          verbose = FALSE)
    
    # Extract results
    edgeR_results <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cmp = FALSE)
    deseq2_results <- resDS(sce, muscat_results, bind = "row", frq = FALSE, cmp = FALSE, method = "DESeq2")
    
    # Add method information
    if (nrow(edgeR_results) > 0) {
      edgeR_results$method <- "Muscat_edgeR"
    }
    if (nrow(deseq2_results) > 0) {
      deseq2_results$method <- "Muscat_DESeq2"
    }
    
    dge_results$muscat_edgeR <- edgeR_results
    dge_results$muscat_deseq2 <- deseq2_results
    
    cat(paste("  Muscat edgeR results:", nrow(edgeR_results), "genes\n"))
    cat(paste("  Muscat DESeq2 results:", nrow(deseq2_results), "genes\n"))
    
  }, error = function(e) {
    cat(paste("Error in Muscat analysis:", e$message, "\n"))
  })
  
  # Save DGE results
  saveRDS(dge_results, file.path(output_base_dir, "DGE_Results", "comprehensive_dge_results.rds"))
  
  # ============================================================================
  # 5. TRAJECTORY ANALYSIS
  # ============================================================================
  
  trajectory_results <- NULL
  
  if (run_trajectory_analysis) {
    cat("\n=== STEP 5: TRAJECTORY ANALYSIS ===\n")
    
    tryCatch({
      trajectory_results <- run_trajectory_analysis_pipeline(
        seu_NKT_focused,
        genes_of_interest = genes_present,
        output_dir = file.path(output_base_dir, "Trajectory_Analysis"),
        methods = c("monocle3", "slingshot")
      )
      
      cat("Trajectory analysis completed successfully!\n")
      
    }, error = function(e) {
      cat(paste("Error in trajectory analysis:", e$message, "\n"))
      cat("Continuing with remaining analyses...\n")
    })
  }
  
  # ============================================================================
  # 6. COMPREHENSIVE VISUALIZATIONS
  # ============================================================================
  
  if (run_visualization) {
    cat("\n=== STEP 6: COMPREHENSIVE VISUALIZATIONS ===\n")
    
    tryCatch({
      # Run comprehensive visualization pipeline
      run_comprehensive_visualization_pipeline(
        seu_NKT_focused,
        dge_results_list = list(Wilcoxon = wilcox_results),
        output_dir = file.path(output_base_dir, "Visualizations")
      )
      
      cat("Comprehensive visualizations completed successfully!\n")
      
    }, error = function(e) {
      cat(paste("Error in visualization pipeline:", e$message, "\n"))
      cat("Continuing with remaining analyses...\n")
    })
  }
  
  # ============================================================================
  # 7. EXCEL OUTPUT GENERATION
  # ============================================================================
  
  cat("\n=== STEP 7: EXCEL OUTPUT GENERATION ===\n")
  
  # Function to create comprehensive Excel output
  create_comprehensive_excel_output <- function(dge_results, method_name, trajectory_data = NULL) {
    
    cat(paste("Creating Excel output for", method_name, "...\n"))
    
    excel_file <- file.path(output_base_dir, "Excel_Outputs", paste0("DGE_", method_name, "_Complete_Analysis.xlsx"))
    
    wb <- createWorkbook()
    
    # Initialize tracking variables
    sheets_created <- 0
    
    # Get method-specific results
    if (method_name == "Wilcoxon" && "wilcoxon" %in% names(dge_results)) {
      method_results <- dge_results$wilcoxon
    } else if (method_name == "Muscat_DESeq2" && "muscat_deseq2" %in% names(dge_results)) {
      method_results <- dge_results$muscat_deseq2
    } else if (method_name == "Muscat_edgeR" && "muscat_edgeR" %in% names(dge_results)) {
      method_results <- dge_results$muscat_edgeR
    } else {
      cat(paste("No results found for method:", method_name, "\n"))
      return(NULL)
    }
    
    # Sheets 1-3: Complete results for each contrast
    for (contrast_name in names(contrasts)) {
      sheet_name <- paste0("All_Genes_", contrast_name)
      
      if (method_name == "Wilcoxon" && contrast_name %in% names(method_results)) {
        addWorksheet(wb, sheet_name)
        writeData(wb, sheet_name, method_results[[contrast_name]])
        sheets_created <- sheets_created + 1
      } else if (method_name %in% c("Muscat_DESeq2", "Muscat_edgeR") && is.data.frame(method_results)) {
        # For Muscat, filter by contrast
        if ("contrast" %in% colnames(method_results)) {
          contrast_data <- method_results[method_results$contrast == contrast_name, ]
        } else {
          # If no contrast column, use all data for first contrast only
          if (contrast_name == names(contrasts)[1]) {
            contrast_data <- method_results
          } else {
            contrast_data <- data.frame()
          }
        }
        
        if (nrow(contrast_data) > 0) {
          addWorksheet(wb, sheet_name)
          writeData(wb, sheet_name, contrast_data)
          sheets_created <- sheets_created + 1
        }
      }
    }
    
    # Sheet 4: S1P genes across all contrasts
    s1p_sheet_data <- data.frame()
    for (contrast_name in names(contrasts)) {
      if (method_name == "Wilcoxon" && contrast_name %in% names(method_results)) {
        s1p_data <- method_results[[contrast_name]][method_results[[contrast_name]]$gene %in% s1p_genes, ]
        if (nrow(s1p_data) > 0) {
          s1p_sheet_data <- rbind(s1p_sheet_data, s1p_data)
        }
      }
    }
    
    if (nrow(s1p_sheet_data) > 0) {
      addWorksheet(wb, "S1P_Genes_All_Contrasts")
      writeData(wb, "S1P_Genes_All_Contrasts", s1p_sheet_data)
      sheets_created <- sheets_created + 1
    }
    
    # Sheet 5: Trafficking genes across all contrasts
    trafficking_sheet_data <- data.frame()
    for (contrast_name in names(contrasts)) {
      if (method_name == "Wilcoxon" && contrast_name %in% names(method_results)) {
        trafficking_data <- method_results[[contrast_name]][method_results[[contrast_name]]$gene %in% trafficking_genes, ]
        if (nrow(trafficking_data) > 0) {
          trafficking_sheet_data <- rbind(trafficking_sheet_data, trafficking_data)
        }
      }
    }
    
    if (nrow(trafficking_sheet_data) > 0) {
      addWorksheet(wb, "Trafficking_Genes_All_Contrasts")
      writeData(wb, "Trafficking_Genes_All_Contrasts", trafficking_sheet_data)
      sheets_created <- sheets_created + 1
    }
    
    # Sheet 6: Summary statistics
    summary_stats <- data.frame(
      Contrast = names(contrasts),
      Method = method_name,
      Total_Genes_Tested = 0,
      Significant_Genes_p005 = 0,
      Significant_Genes_FDR005 = 0,
      Upregulated_Genes = 0,
      Downregulated_Genes = 0,
      stringsAsFactors = FALSE
    )
    
    for (i in 1:length(contrasts)) {
      contrast_name <- names(contrasts)[i]
      if (method_name == "Wilcoxon" && contrast_name %in% names(method_results)) {
        data <- method_results[[contrast_name]]
        summary_stats$Total_Genes_Tested[i] <- nrow(data)
        summary_stats$Significant_Genes_p005[i] <- sum(data$p_val < 0.05, na.rm = TRUE)
        if ("FDR" %in% colnames(data)) {
          summary_stats$Significant_Genes_FDR005[i] <- sum(data$FDR < 0.05, na.rm = TRUE)
        }
        summary_stats$Upregulated_Genes[i] <- sum(data$avg_log2FC > 0 & data$p_val < 0.05, na.rm = TRUE)
        summary_stats$Downregulated_Genes[i] <- sum(data$avg_log2FC < 0 & data$p_val < 0.05, na.rm = TRUE)
      }
    }
    
    addWorksheet(wb, "Summary_Stats_By_Contrast")
    writeData(wb, "Summary_Stats_By_Contrast", summary_stats)
    sheets_created <- sheets_created + 1
    
    # Sheet 7: Genes of interest by compartment
    if ("compartment" %in% colnames(seu_NKT_focused@meta.data)) {
      compartment_analysis <- data.frame(
        Gene = genes_present,
        WB_Mean_Expression = 0,
        BM_Mean_Expression = 0,
        Compartment_Difference = 0,
        Gene_Category = "Other",
        stringsAsFactors = FALSE
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
        }
      }
      
      addWorksheet(wb, "Genes_of_Interest_by_Compartment")
      writeData(wb, "Genes_of_Interest_by_Compartment", compartment_analysis)
      sheets_created <- sheets_created + 1
    }
    
    # Sheets 8-9: Top up/down regulated genes
    top_up <- data.frame()
    top_down <- data.frame()
    
    for (contrast_name in names(contrasts)) {
      if (method_name == "Wilcoxon" && contrast_name %in% names(method_results)) {
        data <- method_results[[contrast_name]]
        
        # Top upregulated
        up_genes <- data[data$avg_log2FC > 0 & data$p_val < 0.05, ]
        if (nrow(up_genes) > 0) {
          up_genes <- up_genes[order(up_genes$avg_log2FC, decreasing = TRUE), ]
          top_up <- rbind(top_up, head(up_genes, 20))
        }
        
        # Top downregulated
        down_genes <- data[data$avg_log2FC < 0 & data$p_val < 0.05, ]
        if (nrow(down_genes) > 0) {
          down_genes <- down_genes[order(down_genes$avg_log2FC), ]
          top_down <- rbind(top_down, head(down_genes, 20))
        }
      }
    }
    
    if (nrow(top_up) > 0) {
      addWorksheet(wb, "Top_Upregulated_Each_Contrast")
      writeData(wb, "Top_Upregulated_Each_Contrast", top_up)
      sheets_created <- sheets_created + 1
    }
    
    if (nrow(top_down) > 0) {
      addWorksheet(wb, "Top_Downregulated_Each_Contrast")
      writeData(wb, "Top_Downregulated_Each_Contrast", top_down)
      sheets_created <- sheets_created + 1
    }
    
    # Sheet 10: Trajectory-associated genes
    trajectory_genes_data <- data.frame(
      Gene = "Analysis_Pending",
      Pseudotime_Correlation = NA,
      P_Value = NA,
      FDR = NA,
      Note = "Trajectory analysis results to be integrated",
      stringsAsFactors = FALSE
    )
    
    # Add trajectory results if available
    if (!is.null(trajectory_data) && "monocle3" %in% names(trajectory_data)) {
      if ("trajectory_genes" %in% names(trajectory_data$monocle3)) {
        traj_genes <- trajectory_data$monocle3$trajectory_genes
        if (nrow(traj_genes) > 0) {
          trajectory_genes_data <- traj_genes[order(traj_genes$q_value), ]
          trajectory_genes_data <- head(trajectory_genes_data, 100)  # Top 100 genes
        }
      }
    }
    
    addWorksheet(wb, "Trajectory_Associated_Genes")
    writeData(wb, "Trajectory_Associated_Genes", trajectory_genes_data)
    sheets_created <- sheets_created + 1
    
    # Sheet 11: Dynamic S1P and trafficking genes
    dynamic_genes <- data.frame(
      Gene = genes_present,
      S1P_Pathway = genes_present %in% s1p_genes,
      Trafficking_Related = genes_present %in% trafficking_genes,
      Retention_Related = genes_present %in% retention_genes,
      Activation_Related = genes_present %in% activation_genes,
      Trajectory_Pattern = "To_be_analyzed",
      Dynamic_Score = NA,
      stringsAsFactors = FALSE
    )
    
    addWorksheet(wb, "Dynamic_S1P_Trafficking_Genes")
    writeData(wb, "Dynamic_S1P_Trafficking_Genes", dynamic_genes)
    sheets_created <- sheets_created + 1
    
    # Save workbook
    if (sheets_created > 0) {
      saveWorkbook(wb, excel_file, overwrite = TRUE)
      cat(paste("Excel file saved:", excel_file, "(", sheets_created, "sheets)\n"))
    } else {
      cat(paste("No data to save for", method_name, "\n"))
    }
    
    return(excel_file)
  }
  
  # Create Excel outputs for each method
  excel_files <- list()
  
  if ("wilcoxon" %in% names(dge_results)) {
    excel_files[["Wilcoxon"]] <- create_comprehensive_excel_output(dge_results, "Wilcoxon", trajectory_results)
  }
  
  if ("muscat_deseq2" %in% names(dge_results)) {
    excel_files[["Muscat_DESeq2"]] <- create_comprehensive_excel_output(dge_results, "Muscat_DESeq2", trajectory_results)
  }
  
  if ("muscat_edgeR" %in% names(dge_results)) {
    excel_files[["Muscat_edgeR"]] <- create_comprehensive_excel_output(dge_results, "Muscat_edgeR", trajectory_results)
  }
  
  # ============================================================================
  # 8. FINAL SUMMARY AND REPORTING
  # ============================================================================
  
  cat("\n=== STEP 8: FINAL SUMMARY AND REPORTING ===\n")
  
  # Calculate analysis timing
  end_time <- Sys.time()
  analysis_duration <- end_time - start_time
  
  # Generate comprehensive summary
  analysis_summary <- list(
    analysis_info = list(
      start_time = start_time,
      end_time = end_time,
      duration = analysis_duration,
      r_version = R.version.string,
      analysis_date = Sys.Date()
    ),
    data_summary = list(
      total_cells = ncol(seu_NKT_focused),
      total_genes = nrow(seu_NKT_focused),
      cell_types = unique(seu_NKT_focused$cell_types),
      treatments = unique(seu_NKT_focused$treatment),
      compartments = unique(seu_NKT_focused$compartment)
    ),
    gene_sets = list(
      s1p_genes = s1p_genes,
      trafficking_genes = trafficking_genes,
      retention_genes = retention_genes,
      activation_genes = activation_genes,
      chemokine_genes = chemokine_genes,
      cytokine_genes = cytokine_genes,
      total_goi = length(all_goi),
      genes_present = length(genes_present),
      genes_missing = length(genes_missing)
    ),
    analysis_results = list(
      contrasts_analyzed = length(contrasts),
      dge_methods = names(dge_results),
      trajectory_analysis = !is.null(trajectory_results),
      visualization_completed = run_visualization,
      excel_files = excel_files
    ),
    output_structure = list(
      base_directory = output_base_dir,
      excel_outputs = list.files(file.path(output_base_dir, "Excel_Outputs"), pattern = "\\.xlsx$"),
      visualization_files = list.files(file.path(output_base_dir, "Visualizations"), pattern = "\\.pdf$"),
      trajectory_files = if (dir.exists(file.path(output_base_dir, "Trajectory_Analysis"))) {
        list.files(file.path(output_base_dir, "Trajectory_Analysis"), pattern = "\\.pdf$")
      } else NULL
    )
  )
  
  # Save comprehensive summary
  saveRDS(analysis_summary, file.path(output_base_dir, "Summary_Reports", "comprehensive_analysis_summary.rds"))
  
  # Create human-readable summary report
  summary_text <- c(
    "============================================================================",
    "COMPREHENSIVE S1P RECEPTOR AND LYMPHOCYTE TRAFFICKING ANALYSIS SUMMARY",
    "============================================================================",
    "",
    paste("Analysis completed on:", Sys.Date()),
    paste("Total analysis time:", round(as.numeric(analysis_duration), 2), "minutes"),
    "",
    "DATA SUMMARY:",
    paste("- Total cells analyzed:", ncol(seu_NKT_focused)),
    paste("- Total genes:", nrow(seu_NKT_focused)),
    paste("- Cell types:", paste(unique(seu_NKT_focused$cell_types), collapse = ", ")),
    paste("- Treatments:", paste(unique(seu_NKT_focused$treatment), collapse = ", ")),
    paste("- Compartments:", paste(unique(seu_NKT_focused$compartment), collapse = ", ")),
    "",
    "GENE SETS ANALYZED:",
    paste("- S1P pathway genes:", length(s1p_genes)),
    paste("- Trafficking genes:", length(trafficking_genes)),
    paste("- Retention genes:", length(retention_genes)),
    paste("- Activation genes:", length(activation_genes)),
    paste("- Total genes of interest:", length(all_goi)),
    paste("- Genes present in dataset:", length(genes_present)),
    "",
    "ANALYSIS METHODS:",
    paste("- DGE methods:", paste(names(dge_results), collapse = ", ")),
    paste("- Trajectory analysis:", ifelse(!is.null(trajectory_results), "Completed", "Not performed")),
    paste("- Comprehensive visualizations:", ifelse(run_visualization, "Completed", "Not performed")),
    "",
    "OUTPUT FILES:",
    paste("- Excel files:", length(excel_files)),
    paste("- Visualization PDFs:", length(list.files(file.path(output_base_dir, "Visualizations"), pattern = "\\.pdf$"))),
    "",
    "RESULTS DIRECTORY:",
    paste("All results saved to:", output_base_dir),
    "",
    "============================================================================"
  )
  
  # Save summary report
  writeLines(summary_text, file.path(output_base_dir, "Summary_Reports", "analysis_summary.txt"))
  
  # Print summary to console
  cat("\n")
  cat(paste(summary_text, collapse = "\n"))
  cat("\n")
  
  # Clean up workspace if requested
  if (cleanup_workspace) {
    cat("Cleaning up workspace...\n")
    rm(list = setdiff(ls(), c("analysis_summary", "output_base_dir")))
  }
  
  # Return comprehensive results
  return(list(
    seu_object = seu_NKT_focused,
    dge_results = dge_results,
    trajectory_results = trajectory_results,
    analysis_summary = analysis_summary,
    output_directory = output_base_dir
  ))
}

# ============================================================================
# SCRIPT EXECUTION
# ============================================================================

# Main execution
cat("Starting comprehensive S1P receptor and lymphocyte trafficking analysis...\n")

# Execute the complete analysis pipeline
tryCatch({
  results <- execute_comprehensive_s1p_analysis(
    data_file = NULL,  # Uses default data file
    output_base_dir = "S1P_Comprehensive_Results",
    run_trajectory_analysis = TRUE,
    run_visualization = TRUE,
    cleanup_workspace = FALSE
  )
  
  cat("\n============================================================================\n")
  cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
  cat("============================================================================\n")
  
}, error = function(e) {
  cat("\n============================================================================\n")
  cat("ERROR IN ANALYSIS EXECUTION:\n")
  cat("============================================================================\n")
  cat(paste("Error message:", e$message, "\n"))
  cat("Please check the error logs and try again.\n")
  cat("============================================================================\n")
})

cat("Script execution completed.\n")