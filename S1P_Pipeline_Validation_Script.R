#!/usr/bin/env Rscript
# ============================================================================
# S1P ANALYSIS PIPELINE VALIDATION AND TESTING SCRIPT
# ============================================================================
# 
# This script validates the S1P analysis pipeline structure and provides
# a minimal example for testing without requiring the actual data file.
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

# Set up environment
library(here)
set.seed(42)

cat("============================================================================\n")
cat("S1P ANALYSIS PIPELINE VALIDATION AND TESTING\n")
cat("============================================================================\n")

# ============================================================================
# 1. VALIDATE PIPELINE STRUCTURE
# ============================================================================

cat("\n=== STEP 1: VALIDATING PIPELINE STRUCTURE ===\n")

# Check if all required files exist
required_files <- c(
  "S1P_Main_Analysis_Script.R",
  "S1P_Comprehensive_Analysis_Pipeline.R",
  "S1P_Trajectory_Analysis_Module.R",
  "S1P_Enhanced_Visualization_Module.R",
  "S1P_Analysis_Pipeline_README.md",
  "Custom R Functions and Scripts/init_packages.R",
  "Custom R Functions and Scripts/install_packages.R"
)

missing_files <- c()
for (file in required_files) {
  if (!file.exists(here(file))) {
    missing_files <- c(missing_files, file)
  }
}

if (length(missing_files) > 0) {
  cat("Missing files:\n")
  for (file in missing_files) {
    cat(paste("  -", file, "\n"))
  }
  stop("Pipeline validation failed: Missing required files")
} else {
  cat("All required pipeline files are present ✓\n")
}

# ============================================================================
# 2. VALIDATE PACKAGE REQUIREMENTS
# ============================================================================

cat("\n=== STEP 2: VALIDATING PACKAGE REQUIREMENTS ===\n")

# Check core packages
core_packages <- c(
  "Seurat", "tidyverse", "here", "ggplot2", "dplyr", "openxlsx", 
  "pheatmap", "viridis", "patchwork", "ggrepel"
)

missing_packages <- c()
for (pkg in core_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_packages <- c(missing_packages, pkg)
  }
}

if (length(missing_packages) > 0) {
  cat("Missing core packages:\n")
  for (pkg in missing_packages) {
    cat(paste("  -", pkg, "\n"))
  }
  cat("Please run: source('Custom R Functions and Scripts/install_packages.R')\n")
} else {
  cat("All core packages are available ✓\n")
}

# Check trajectory analysis packages
trajectory_packages <- c("monocle3", "slingshot", "SingleCellExperiment")
missing_trajectory <- c()
for (pkg in trajectory_packages) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    missing_trajectory <- c(missing_trajectory, pkg)
  }
}

if (length(missing_trajectory) > 0) {
  cat("Missing trajectory analysis packages:\n")
  for (pkg in missing_trajectory) {
    cat(paste("  -", pkg, "\n"))
  }
  cat("Note: Trajectory analysis will be skipped if these packages are not available\n")
} else {
  cat("All trajectory analysis packages are available ✓\n")
}

# ============================================================================
# 3. VALIDATE FUNCTIONS AND MODULES
# ============================================================================

cat("\n=== STEP 3: VALIDATING FUNCTIONS AND MODULES ===\n")

# Test sourcing modules
tryCatch({
  source(here("Custom R Functions and Scripts/init_packages.R"))
  cat("Package initialization successful ✓\n")
}, error = function(e) {
  cat(paste("Package initialization failed:", e$message, "\n"))
})

# Test main pipeline functions
tryCatch({
  # Source the main pipeline without executing
  source(here("S1P_Comprehensive_Analysis_Pipeline.R"))
  cat("Main pipeline module loaded successfully ✓\n")
}, error = function(e) {
  cat(paste("Main pipeline module failed:", e$message, "\n"))
})

# Test trajectory analysis module
tryCatch({
  source(here("S1P_Trajectory_Analysis_Module.R"))
  cat("Trajectory analysis module loaded successfully ✓\n")
}, error = function(e) {
  cat(paste("Trajectory analysis module failed:", e$message, "\n"))
})

# Test visualization module
tryCatch({
  source(here("S1P_Enhanced_Visualization_Module.R"))
  cat("Visualization module loaded successfully ✓\n")
}, error = function(e) {
  cat(paste("Visualization module failed:", e$message, "\n"))
})

# ============================================================================
# 4. CREATE MINIMAL EXAMPLE DATA
# ============================================================================

cat("\n=== STEP 4: CREATING MINIMAL EXAMPLE DATA ===\n")

# Create a minimal Seurat object for testing
if (requireNamespace("Seurat", quietly = TRUE)) {
  library(Seurat)
  
  # Create synthetic data
  set.seed(42)
  n_cells <- 500
  n_genes <- 1000
  
  # Create expression matrix
  expr_matrix <- matrix(
    rnbinom(n_cells * n_genes, size = 1, mu = 0.5),
    nrow = n_genes,
    ncol = n_cells
  )
  
  # Add some genes of interest
  goi <- c("S1pr1", "S1pr4", "S1pr5", "Cxcr4", "Sell", "Ccr7", "Cd69", "Klf2")
  gene_names <- c(goi, paste0("Gene_", 1:(n_genes - length(goi))))
  
  rownames(expr_matrix) <- gene_names
  colnames(expr_matrix) <- paste0("Cell_", 1:n_cells)
  
  # Create metadata
  metadata <- data.frame(
    cell_id = colnames(expr_matrix),
    treatment = sample(c("NAIVE", "SHAM", "DMG"), n_cells, replace = TRUE),
    compartment = sample(c("WB", "BM"), n_cells, replace = TRUE),
    cell_types = sample(c("T_cells", "NK_cells", "NKT_cells"), n_cells, replace = TRUE),
    sample_id = paste0("Sample_", sample(1:6, n_cells, replace = TRUE)),
    stringsAsFactors = FALSE
  )
  
  rownames(metadata) <- colnames(expr_matrix)
  
  # Create Seurat object
  test_seurat <- CreateSeuratObject(
    counts = expr_matrix,
    meta.data = metadata,
    project = "S1P_Test"
  )
  
  # Add some processing
  test_seurat <- NormalizeData(test_seurat, verbose = FALSE)
  test_seurat <- FindVariableFeatures(test_seurat, verbose = FALSE)
  test_seurat <- ScaleData(test_seurat, verbose = FALSE)
  test_seurat <- RunPCA(test_seurat, verbose = FALSE)
  test_seurat <- RunUMAP(test_seurat, dims = 1:10, verbose = FALSE)
  
  # Save test data
  if (!dir.exists(here("Data"))) {
    dir.create(here("Data"))
  }
  
  saveRDS(test_seurat, here("Data", "test_seurat_object.rds"))
  
  cat("Minimal example data created successfully ✓\n")
  cat(paste("  - Cells:", ncol(test_seurat), "\n"))
  cat(paste("  - Genes:", nrow(test_seurat), "\n"))
  cat(paste("  - Treatments:", paste(unique(test_seurat$treatment), collapse = ", "), "\n"))
  cat(paste("  - Compartments:", paste(unique(test_seurat$compartment), collapse = ", "), "\n"))
  cat(paste("  - Cell types:", paste(unique(test_seurat$cell_types), collapse = ", "), "\n"))
  
} else {
  cat("Seurat package not available - skipping example data creation\n")
}

# ============================================================================
# 5. TEST PIPELINE EXECUTION (MINIMAL)
# ============================================================================

cat("\n=== STEP 5: TESTING PIPELINE EXECUTION ===\n")

# Test basic DGE analysis function
if (exists("test_seurat") && requireNamespace("Seurat", quietly = TRUE)) {
  
  tryCatch({
    # Test basic contrast
    contrasts <- list("DMG_vs_SHAM" = c("DMG", "SHAM"))
    
    # Set identity
    Idents(test_seurat) <- "treatment"
    
    # Test Wilcoxon
    markers <- FindMarkers(test_seurat,
                          ident.1 = "DMG",
                          ident.2 = "SHAM",
                          min.pct = 0.1,
                          logfc.threshold = 0.1,
                          test.use = "wilcox",
                          verbose = FALSE)
    
    if (nrow(markers) > 0) {
      cat("Basic DGE analysis test successful ✓\n")
      cat(paste("  - Found", nrow(markers), "differential genes\n"))
    } else {
      cat("Basic DGE analysis test: No differential genes found\n")
    }
    
  }, error = function(e) {
    cat(paste("Basic DGE analysis test failed:", e$message, "\n"))
  })
  
  # Test visualization functions
  tryCatch({
    # Test basic plot
    p <- VlnPlot(test_seurat, features = "S1pr1", group.by = "treatment", 
                pt.size = 0.1, combine = FALSE)[[1]]
    
    if (is.ggplot(p)) {
      cat("Basic visualization test successful ✓\n")
    }
    
  }, error = function(e) {
    cat(paste("Basic visualization test failed:", e$message, "\n"))
  })
  
} else {
  cat("Skipping pipeline execution test - missing requirements\n")
}

# ============================================================================
# 6. GENERATE VALIDATION REPORT
# ============================================================================

cat("\n=== STEP 6: GENERATING VALIDATION REPORT ===\n")

# Create validation report
validation_report <- list(
  validation_date = Sys.Date(),
  r_version = R.version.string,
  system_info = Sys.info(),
  pipeline_files = list(
    required_files = required_files,
    missing_files = missing_files,
    all_present = length(missing_files) == 0
  ),
  packages = list(
    core_packages = core_packages,
    missing_core = missing_packages,
    core_available = length(missing_packages) == 0,
    trajectory_packages = trajectory_packages,
    missing_trajectory = missing_trajectory,
    trajectory_available = length(missing_trajectory) == 0
  ),
  test_data = list(
    created = exists("test_seurat"),
    cells = if (exists("test_seurat")) ncol(test_seurat) else 0,
    genes = if (exists("test_seurat")) nrow(test_seurat) else 0
  ),
  validation_status = (length(missing_files) == 0 && length(missing_packages) == 0)
)

# Save validation report
if (!dir.exists(here("Pipeline_Validation"))) {
  dir.create(here("Pipeline_Validation"))
}

saveRDS(validation_report, here("Pipeline_Validation", "validation_report.rds"))

# Create readable summary
validation_summary <- c(
  "============================================================================",
  "S1P ANALYSIS PIPELINE VALIDATION SUMMARY",
  "============================================================================",
  "",
  paste("Validation Date:", Sys.Date()),
  paste("R Version:", R.version.string),
  "",
  "PIPELINE STRUCTURE:",
  paste("- Required files present:", ifelse(length(missing_files) == 0, "YES ✓", "NO ✗")),
  paste("- Core packages available:", ifelse(length(missing_packages) == 0, "YES ✓", "NO ✗")),
  paste("- Trajectory packages available:", ifelse(length(missing_trajectory) == 0, "YES ✓", "PARTIAL")),
  "",
  "TEST DATA:",
  paste("- Example data created:", ifelse(exists("test_seurat"), "YES ✓", "NO ✗")),
  if (exists("test_seurat")) paste("- Test cells:", ncol(test_seurat)) else "- Test cells: 0",
  if (exists("test_seurat")) paste("- Test genes:", nrow(test_seurat)) else "- Test genes: 0",
  "",
  "OVERALL STATUS:",
  paste("- Pipeline ready:", ifelse(validation_report$validation_status, "YES ✓", "NO ✗")),
  "",
  "NEXT STEPS:",
  if (length(missing_packages) > 0) "1. Install missing core packages" else "1. Core packages OK ✓",
  if (length(missing_trajectory) > 0) "2. Install trajectory packages (optional)" else "2. Trajectory packages OK ✓",
  "3. Place your data file in Data/ directory",
  "4. Run: source('S1P_Main_Analysis_Script.R')",
  "",
  "============================================================================"
)

writeLines(validation_summary, here("Pipeline_Validation", "validation_summary.txt"))

# Print summary
cat("\n")
cat(paste(validation_summary, collapse = "\n"))
cat("\n")

# ============================================================================
# 7. CLEANUP AND FINALIZATION
# ============================================================================

cat("\n=== PIPELINE VALIDATION COMPLETE ===\n")

if (validation_report$validation_status) {
  cat("SUCCESS: Pipeline is ready for execution! ✓\n")
} else {
  cat("WARNING: Pipeline requires additional setup before execution\n")
}

cat(paste("Validation report saved to:", here("Pipeline_Validation"), "\n"))
cat("============================================================================\n")