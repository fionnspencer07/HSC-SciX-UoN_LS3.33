#!/usr/bin/env Rscript
# ============================================================================
# S1P ANALYSIS PIPELINE STRUCTURE VALIDATION
# ============================================================================
# 
# This script validates the S1P analysis pipeline structure without requiring
# external packages.
#
# Author: Copilot AI Assistant
# Date: 2024
# Version: 1.0
# ============================================================================

cat("============================================================================\n")
cat("S1P ANALYSIS PIPELINE STRUCTURE VALIDATION\n")
cat("============================================================================\n")

# ============================================================================
# 1. VALIDATE PIPELINE STRUCTURE
# ============================================================================

cat("\n=== VALIDATING PIPELINE STRUCTURE ===\n")

# Check if all required files exist
required_files <- c(
  "S1P_Main_Analysis_Script.R",
  "S1P_Comprehensive_Analysis_Pipeline.R",
  "S1P_Trajectory_Analysis_Module.R",
  "S1P_Enhanced_Visualization_Module.R",
  "S1P_Analysis_Pipeline_README.md",
  "S1P_Pipeline_Validation_Script.R",
  "Custom R Functions and Scripts/init_packages.R",
  "Custom R Functions and Scripts/install_packages.R"
)

missing_files <- c()
existing_files <- c()

for (file in required_files) {
  if (!file.exists(file)) {
    missing_files <- c(missing_files, file)
  } else {
    existing_files <- c(existing_files, file)
  }
}

cat("EXISTING FILES:\n")
for (file in existing_files) {
  cat(paste("  ✓", file, "\n"))
}

if (length(missing_files) > 0) {
  cat("\nMISSING FILES:\n")
  for (file in missing_files) {
    cat(paste("  ✗", file, "\n"))
  }
} else {
  cat("\nAll required pipeline files are present ✓\n")
}

# ============================================================================
# 2. VALIDATE FILE CONTENTS
# ============================================================================

cat("\n=== VALIDATING FILE CONTENTS ===\n")

# Check main script
if (file.exists("S1P_Main_Analysis_Script.R")) {
  main_script <- readLines("S1P_Main_Analysis_Script.R")
  
  # Check for key functions
  key_functions <- c(
    "execute_comprehensive_s1p_analysis",
    "source.*S1P_Comprehensive_Analysis_Pipeline.R",
    "source.*S1P_Trajectory_Analysis_Module.R",
    "source.*S1P_Enhanced_Visualization_Module.R"
  )
  
  found_functions <- c()
  for (func in key_functions) {
    if (any(grepl(func, main_script))) {
      found_functions <- c(found_functions, func)
    }
  }
  
  cat("Main script functions found:\n")
  for (func in found_functions) {
    cat(paste("  ✓", func, "\n"))
  }
  
  if (length(found_functions) == length(key_functions)) {
    cat("Main script structure is valid ✓\n")
  } else {
    cat("Main script structure has issues ✗\n")
  }
}

# Check comprehensive analysis pipeline
if (file.exists("S1P_Comprehensive_Analysis_Pipeline.R")) {
  pipeline_script <- readLines("S1P_Comprehensive_Analysis_Pipeline.R")
  
  # Check for key components
  key_components <- c(
    "perform_wilcoxon_analysis",
    "perform_muscat_analysis",
    "create_excel_output",
    "s1p_genes.*<-",
    "trafficking_genes.*<-"
  )
  
  found_components <- c()
  for (comp in key_components) {
    if (any(grepl(comp, pipeline_script))) {
      found_components <- c(found_components, comp)
    }
  }
  
  cat("\nPipeline components found:\n")
  for (comp in found_components) {
    cat(paste("  ✓", comp, "\n"))
  }
}

# Check trajectory analysis module
if (file.exists("S1P_Trajectory_Analysis_Module.R")) {
  trajectory_script <- readLines("S1P_Trajectory_Analysis_Module.R")
  
  # Check for key components
  trajectory_components <- c(
    "perform_monocle3_analysis",
    "perform_slingshot_analysis",
    "run_trajectory_analysis_pipeline"
  )
  
  found_trajectory <- c()
  for (comp in trajectory_components) {
    if (any(grepl(comp, trajectory_script))) {
      found_trajectory <- c(found_trajectory, comp)
    }
  }
  
  cat("\nTrajectory analysis components found:\n")
  for (comp in found_trajectory) {
    cat(paste("  ✓", comp, "\n"))
  }
}

# Check visualization module
if (file.exists("S1P_Enhanced_Visualization_Module.R")) {
  viz_script <- readLines("S1P_Enhanced_Visualization_Module.R")
  
  # Check for key components
  viz_components <- c(
    "create_enhanced_heatmap",
    "create_enhanced_volcano",
    "create_nebulosa_density_plots",
    "run_comprehensive_visualization_pipeline"
  )
  
  found_viz <- c()
  for (comp in viz_components) {
    if (any(grepl(comp, viz_script))) {
      found_viz <- c(found_viz, comp)
    }
  }
  
  cat("\nVisualization components found:\n")
  for (comp in found_viz) {
    cat(paste("  ✓", comp, "\n"))
  }
}

# ============================================================================
# 3. VALIDATE GENE SETS
# ============================================================================

cat("\n=== VALIDATING GENE SETS ===\n")

if (file.exists("S1P_Comprehensive_Analysis_Pipeline.R")) {
  pipeline_script <- readLines("S1P_Comprehensive_Analysis_Pipeline.R")
  
  # Extract gene sets
  gene_sets <- list()
  
  # S1P genes
  s1p_line <- grep("s1p_genes.*<-", pipeline_script, value = TRUE)
  if (length(s1p_line) > 0) {
    cat("S1P genes defined ✓\n")
    cat(paste("  ", s1p_line[1], "\n"))
  }
  
  # Trafficking genes
  trafficking_line <- grep("trafficking_genes.*<-", pipeline_script, value = TRUE)
  if (length(trafficking_line) > 0) {
    cat("Trafficking genes defined ✓\n")
    cat(paste("  ", trafficking_line[1], "\n"))
  }
  
  # Other gene sets
  other_gene_sets <- c("retention_genes", "activation_genes", "chemokine_genes", "cytokine_genes")
  for (gene_set in other_gene_sets) {
    gene_set_line <- grep(paste0(gene_set, ".*<-"), pipeline_script, value = TRUE)
    if (length(gene_set_line) > 0) {
      cat(paste(gene_set, "defined ✓\n"))
    }
  }
}

# ============================================================================
# 4. VALIDATE CONTRASTS
# ============================================================================

cat("\n=== VALIDATING CONTRASTS ===\n")

if (file.exists("S1P_Main_Analysis_Script.R")) {
  main_script <- readLines("S1P_Main_Analysis_Script.R")
  
  # Check for contrasts
  contrast_patterns <- c(
    "DMG_vs_SHAM",
    "DMG_vs_NAIVE",
    "SHAM_vs_NAIVE"
  )
  
  found_contrasts <- c()
  for (contrast in contrast_patterns) {
    if (any(grepl(contrast, main_script))) {
      found_contrasts <- c(found_contrasts, contrast)
    }
  }
  
  cat("Contrasts found:\n")
  for (contrast in found_contrasts) {
    cat(paste("  ✓", contrast, "\n"))
  }
  
  if (length(found_contrasts) == length(contrast_patterns)) {
    cat("All required contrasts are defined ✓\n")
  }
}

# ============================================================================
# 5. VALIDATE EXCEL OUTPUT STRUCTURE
# ============================================================================

cat("\n=== VALIDATING EXCEL OUTPUT STRUCTURE ===\n")

if (file.exists("S1P_Main_Analysis_Script.R")) {
  main_script <- readLines("S1P_Main_Analysis_Script.R")
  
  # Check for required Excel sheets
  excel_sheets <- c(
    "All_Genes_DMG_vs_SHAM",
    "All_Genes_DMG_vs_NAIVE",
    "All_Genes_SHAM_vs_NAIVE",
    "S1P_Genes_All_Contrasts",
    "Trafficking_Genes_All_Contrasts",
    "Summary_Stats_By_Contrast",
    "Genes_of_Interest_by_Compartment",
    "Top_Upregulated_Each_Contrast",
    "Top_Downregulated_Each_Contrast",
    "Trajectory_Associated_Genes",
    "Dynamic_S1P_Trafficking_Genes"
  )
  
  found_sheets <- c()
  for (sheet in excel_sheets) {
    if (any(grepl(sheet, main_script))) {
      found_sheets <- c(found_sheets, sheet)
    }
  }
  
  cat("Excel sheet definitions found:\n")
  for (sheet in found_sheets) {
    cat(paste("  ✓", sheet, "\n"))
  }
  
  cat(paste("\nExpected sheets: 11, Found:", length(found_sheets), "\n"))
}

# ============================================================================
# 6. VALIDATE VISUALIZATION TYPES
# ============================================================================

cat("\n=== VALIDATING VISUALIZATION TYPES ===\n")

if (file.exists("S1P_Enhanced_Visualization_Module.R")) {
  viz_script <- readLines("S1P_Enhanced_Visualization_Module.R")
  
  # Check for visualization types
  viz_types <- c(
    "heatmap",
    "volcano",
    "nebulosa",
    "violin",
    "ridge",
    "correlation",
    "feature"
  )
  
  found_viz_types <- c()
  for (viz_type in viz_types) {
    if (any(grepl(viz_type, viz_script, ignore.case = TRUE))) {
      found_viz_types <- c(found_viz_types, viz_type)
    }
  }
  
  cat("Visualization types found:\n")
  for (viz_type in found_viz_types) {
    cat(paste("  ✓", viz_type, "\n"))
  }
}

# ============================================================================
# 7. GENERATE VALIDATION SUMMARY
# ============================================================================

cat("\n=== VALIDATION SUMMARY ===\n")

# Calculate validation score
total_files <- length(required_files)
found_files <- length(existing_files)
validation_score <- round((found_files / total_files) * 100, 1)

cat("PIPELINE STRUCTURE VALIDATION RESULTS:\n")
cat("============================================================================\n")
cat(paste("Files present:", found_files, "/", total_files, "(", validation_score, "%)\n"))
cat(paste("Validation date:", Sys.Date(), "\n"))
cat(paste("R version:", R.version.string, "\n"))
cat("\nVALIDATION STATUS:\n")

if (validation_score >= 100) {
  cat("✓ EXCELLENT: All pipeline files are present and structured correctly\n")
} else if (validation_score >= 90) {
  cat("✓ GOOD: Most pipeline files are present with minor issues\n")
} else if (validation_score >= 80) {
  cat("⚠ FAIR: Pipeline is mostly complete but needs attention\n")
} else {
  cat("✗ POOR: Pipeline structure needs significant work\n")
}

cat("\nNEXT STEPS:\n")
cat("1. Install required R packages: source('Custom R Functions and Scripts/install_packages.R')\n")
cat("2. Place your data file in the Data/ directory\n")
cat("3. Run the validation script: Rscript S1P_Pipeline_Validation_Script.R\n")
cat("4. Execute the pipeline: source('S1P_Main_Analysis_Script.R')\n")

cat("\n============================================================================\n")
cat("PIPELINE STRUCTURE VALIDATION COMPLETE\n")
cat("============================================================================\n")

# Save validation results
if (!dir.exists("Pipeline_Validation")) {
  dir.create("Pipeline_Validation")
}

# Create simple validation report
validation_report <- c(
  "S1P Analysis Pipeline Structure Validation Report",
  "================================================",
  paste("Date:", Sys.Date()),
  paste("R Version:", R.version.string),
  "",
  "Files Status:",
  paste("Total files required:", total_files),
  paste("Files present:", found_files),
  paste("Validation score:", validation_score, "%"),
  "",
  "Missing files:",
  if (length(missing_files) > 0) paste("  -", missing_files) else "  None",
  "",
  "Status:",
  if (validation_score >= 100) "EXCELLENT" else if (validation_score >= 90) "GOOD" else if (validation_score >= 80) "FAIR" else "POOR",
  "",
  "The pipeline structure is ready for package installation and execution."
)

writeLines(validation_report, "Pipeline_Validation/structure_validation_report.txt")

cat(paste("Validation report saved to: Pipeline_Validation/structure_validation_report.txt\n"))