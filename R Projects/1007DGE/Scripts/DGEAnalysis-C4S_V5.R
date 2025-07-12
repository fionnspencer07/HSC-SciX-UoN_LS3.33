library(here)
#init packages
source(here("Custom R Functions and Scripts", "init_packages.R"))  # will only work if the rproj opened is HSC SciX UoN LS3.33

# Load additional packages for Excel export
library(openxlsx)
library(pheatmap)
library(ComplexHeatmap)

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

# ============================================================================
# COMPREHENSIVE CONDITION-COMPARTMENT ANALYSIS
# All possible comparisons: WB vs BM, treatments within compartments
# ============================================================================

print("=== COMPREHENSIVE CONDITION-COMPARTMENT ANALYSIS ===")

# Create results directory with your specified path
if(!dir.exists(here("R Projects", "1007DGE", "DGEAnalysis-C4S_V4"))) {
  dir.create(here("R Projects", "1007DGE", "DGEAnalysis-C4S_V4"), recursive = TRUE)
}

# Create tissue mapping based on your sample information
create_tissue_mapping <- function(seu_obj) {
  library_codes <- seu_obj$orig.ident
  tissue_mapping <- case_when(
    library_codes %in% c("LIB2A_DUNLAB_TV", "LIB4A_DUNLAB_TV") ~ "WB",
    library_codes %in% c("LIB2E2_DUNLAB_TV", "LIB2F2_DUNLAB_TV", 
                         "LIB2G2_DUNLAB_TV", "LIB2H2_DUNLAB_TV") ~ "BM",
    TRUE ~ "Unknown"
  )
  return(tissue_mapping)
}

# Add tissue information to metadata
seu_NKT_focused$compartment <- create_tissue_mapping(seu_NKT_focused)
seu_NKT_focused$condition <- seu_NKT_focused$treatment

# Create comprehensive grouping variables (FIXED - using sample_name)
seu_NKT_focused$condition_compartment <- paste(seu_NKT_focused$condition, seu_NKT_focused$compartment, sep = "_")
seu_NKT_focused$sample_condition_compartment <- paste(seu_NKT_focused$sample_name, seu_NKT_focused$compartment, sep = "_")

# Create a proper sample_id for muscat compatibility
seu_NKT_focused$sample_id <- paste(seu_NKT_focused$treatment, seu_NKT_focused$orig.ident, sep = "_")

# Verify the setup
print("=== DATA STRUCTURE VERIFICATION ===")
print("Condition (Treatment) distribution:")
print(table(seu_NKT_focused$condition))

print("Compartment (Tissue) distribution:")
print(table(seu_NKT_focused$compartment))

print("Condition-Compartment combinations:")
print(table(seu_NKT_focused$condition_compartment))

print("Sample distribution:")
sample_summary <- seu_NKT_focused@meta.data %>%
  group_by(sample_name, condition, compartment, condition_compartment) %>%  # FIXED - using sample_name
  summarise(n_cells = n(), .groups = "drop")
print(sample_summary)

# Set default assay
DefaultAssay(seu_NKT_focused) <- "RNA"

# ============================================================================
# GENES OF INTEREST DEFINITION
# ============================================================================

# Your genes of interest - defined early for volcano plot highlighting
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")

all_goi <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

# Create GOI categories for better visualization
goi_categories <- data.frame(
  gene = all_goi,
  category = case_when(
    all_goi %in% s1pr_genes ~ "S1P Receptors",
    all_goi %in% egress_genes & !all_goi %in% s1pr_genes ~ "Egress Markers",
    all_goi %in% retention_genes ~ "Retention Markers",
    all_goi %in% trafficking_genes ~ "Trafficking Markers",
    TRUE ~ "Other"
  )
)

print("=== GENES OF INTEREST ===")
print("S1P Receptors:")
print(s1pr_genes)
print("Trafficking genes:")
print(trafficking_genes)
print("Retention genes:")
print(retention_genes)
print("Egress genes:")
print(egress_genes)

# ============================================================================
# ANALYSIS FRAMEWORK - ALL COMPARISON TYPES
# ============================================================================

all_results <- list()

print("=== ANALYSIS PLAN ===")
print("1. WB vs BM within each condition (UT, NAIVE, SHAM)")
print("2. WB vs BM across all conditions combined")  
print("3. Condition comparisons within WB compartment")
print("4. Condition comparisons within BM compartment")
print("5. All pairwise condition comparisons within each compartment")

# ============================================================================
# 1. COMPARTMENT COMPARISONS: WB vs BM
# ============================================================================

print("\n=== 1. COMPARTMENT COMPARISONS: WB vs BM ===")

# 1A. WB vs BM within each condition
for(condition in unique(seu_NKT_focused$condition)) {
  print(paste("Processing WB vs BM within condition:", condition))
  
  # Subset to condition
  seu_condition <- subset(seu_NKT_focused, subset = condition == condition)
  
  # Check compartment availability
  compartments_available <- unique(seu_condition$compartment)
  
  if(length(compartments_available) >= 2 && all(c("WB", "BM") %in% compartments_available)) {
    
    # For each cell type
    for(cell_type in unique(seu_condition$cell_types)) {
      print(paste("  Cell type:", cell_type))
      
      # Subset to cell type
      seu_subset <- subset(seu_condition, subset = cell_types == cell_type)
      
      # Check cell counts
      compartment_counts <- table(seu_subset$compartment)
      print(paste("    Cell counts:", paste(names(compartment_counts), compartment_counts, sep = ":", collapse = ", ")))
      
      if(all(compartment_counts >= 10)) {
        
        # Set identity to compartment
        Idents(seu_subset) <- "compartment"
        
        # Compare WB vs BM
        tryCatch({
          markers <- FindMarkers(seu_subset,
                                 ident.1 = "WB",
                                 ident.2 = "BM",
                                 min.pct = 0.1,
                                 logfc.threshold = 0.1,
                                 test.use = "wilcox",
                                 verbose = FALSE)
          
          if(nrow(markers) > 0) {
            markers$gene <- rownames(markers)
            markers$comparison_type <- "Compartment"
            markers$comparison <- "WB_vs_BM"
            markers$condition <- condition
            markers$cell_type <- cell_type
            markers$method <- "Wilcoxon"
            
            result_name <- paste("Compartment", condition, cell_type, "WB_vs_BM", sep = "_")
            all_results[[result_name]] <- markers
            
            print(paste("    Found", nrow(markers), "differential genes"))
          }
        }, error = function(e) {
          print(paste("    Error:", e$message))
        })
      } else {
        print(paste("    Skipping - insufficient cells"))
      }
    }
  } else {
    print(paste("  Skipping", condition, "- missing WB or BM samples"))
  }
}

# 1B. WB vs BM across all conditions combined
print("\n--- WB vs BM across all conditions combined ---")

for(cell_type in unique(seu_NKT_focused$cell_types)) {
  print(paste("Processing cell type:", cell_type, "(all conditions)"))
  
  # Subset to cell type
  seu_subset <- subset(seu_NKT_focused, subset = cell_types == cell_type)
  
  # Check cell counts
  compartment_counts <- table(seu_subset$compartment)
  print(paste("  Cell counts:", paste(names(compartment_counts), compartment_counts, sep = ":", collapse = ", ")))
  
  if(all(compartment_counts >= 30)) {
    
    # Set identity to compartment
    Idents(seu_subset) <- "compartment"
    
    # Compare WB vs BM
    tryCatch({
      markers <- FindMarkers(seu_subset,
                             ident.1 = "WB",
                             ident.2 = "BM",
                             min.pct = 0.1,
                             logfc.threshold = 0.1,
                             test.use = "wilcox",
                             verbose = FALSE)
      
      if(nrow(markers) > 0) {
        markers$gene <- rownames(markers)
        markers$comparison_type <- "Compartment"
        markers$comparison <- "WB_vs_BM"
        markers$condition <- "All_conditions"
        markers$cell_type <- cell_type
        markers$method <- "Wilcoxon"
        
        result_name <- paste("Compartment", "All_conditions", cell_type, "WB_vs_BM", sep = "_")
        all_results[[result_name]] <- markers
        
        print(paste("  Found", nrow(markers), "differential genes"))
      }
    }, error = function(e) {
      print(paste("  Error:", e$message))
    })
  }
}

# ============================================================================
# 2. CONDITION COMPARISONS WITHIN COMPARTMENTS
# ============================================================================

print("\n=== 2. CONDITION COMPARISONS WITHIN COMPARTMENTS ===")

# Define condition comparisons
condition_comparisons <- list(
  c("NAIVE", "UT"),
  c("SHAM", "UT"),
  c("SHAM", "NAIVE")
)

# 2A. Condition comparisons within WB compartment
print("\n--- Condition comparisons within WB compartment ---")

seu_WB <- subset(seu_NKT_focused, subset = compartment == "WB")

for(comp in condition_comparisons) {
  comparison_name <- paste(comp[1], "vs", comp[2])
  print(paste("Processing", comparison_name, "within WB"))
  
  # Check if both conditions are available in WB
  conditions_available <- unique(seu_WB$condition)
  
  if(all(comp %in% conditions_available)) {
    
    # For each cell type
    for(cell_type in unique(seu_WB$cell_types)) {
      print(paste("  Cell type:", cell_type))
      
      # Subset to cell type
      seu_subset <- subset(seu_WB, subset = cell_types == cell_type)
      
      # Further subset to the two conditions being compared
      seu_subset <- subset(seu_subset, subset = condition %in% comp)
      
      # Check cell counts
      condition_counts <- table(seu_subset$condition)
      print(paste("    Cell counts:", paste(names(condition_counts), condition_counts, sep = ":", collapse = ", ")))
      
      if(all(condition_counts >= 10)) {
        
        # Set identity to condition
        Idents(seu_subset) <- "condition"
        
        # Compare conditions
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
            markers$comparison_type <- "Condition_within_Compartment"
            markers$comparison <- comparison_name
            markers$compartment <- "WB"
            markers$cell_type <- cell_type
            markers$method <- "Wilcoxon"
            
            result_name <- paste("Condition_WB", cell_type, comparison_name, sep = "_")
            all_results[[result_name]] <- markers
            
            print(paste("    Found", nrow(markers), "differential genes"))
          }
        }, error = function(e) {
          print(paste("    Error:", e$message))
        })
      }
    }
  } else {
    print(paste("  Skipping", comparison_name, "- conditions not available in WB"))
  }
}

# 2B. Condition comparisons within BM compartment
print("\n--- Condition comparisons within BM compartment ---")

seu_BM <- subset(seu_NKT_focused, subset = compartment == "BM")

for(comp in condition_comparisons) {
  comparison_name <- paste(comp[1], "vs", comp[2])
  print(paste("Processing", comparison_name, "within BM"))
  
  # Check if both conditions are available in BM
  conditions_available <- unique(seu_BM$condition)
  
  if(all(comp %in% conditions_available)) {
    
    # For each cell type
    for(cell_type in unique(seu_BM$cell_types)) {
      print(paste("  Cell type:", cell_type))
      
      # Subset to cell type
      seu_subset <- subset(seu_BM, subset = cell_types == cell_type)
      
      # Further subset to the two conditions being compared
      seu_subset <- subset(seu_subset, subset = condition %in% comp)
      
      # Check cell counts
      condition_counts <- table(seu_subset$condition)
      print(paste("    Cell counts:", paste(names(condition_counts), condition_counts, sep = ":", collapse = ", ")))
      
      if(all(condition_counts >= 10)) {
        
        # Set identity to condition
        Idents(seu_subset) <- "condition"
        
        # Compare conditions
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
            markers$comparison_type <- "Condition_within_Compartment"
            markers$comparison <- comparison_name
            markers$compartment <- "BM"
            markers$cell_type <- cell_type
            markers$method <- "Wilcoxon"
            
            result_name <- paste("Condition_BM", cell_type, comparison_name, sep = "_")
            all_results[[result_name]] <- markers
            
            print(paste("    Found", nrow(markers), "differential genes"))
          }
        }, error = function(e) {
          print(paste("    Error:", e$message))
        })
      }
    }
  } else {
    print(paste("  Skipping", comparison_name, "- conditions not available in BM"))
  }
}

# ============================================================================
# 3. ALL PAIRWISE COMPARISONS - CONDITION-COMPARTMENT COMBINATIONS
# ============================================================================

print("\n=== 3. ALL PAIRWISE COMPARISONS - CONDITION-COMPARTMENT COMBINATIONS ===")

# Define all pairwise comparisons requested
pairwise_comparisons <- list(
  c("NAIVE_WB", "NAIVE_BM"),      # 1. naive - whole blood vs naive - bone marrow
  c("NAIVE_WB", "SHAM_WB"),       # 2. naive - whole blood vs sham - whole blood
  c("NAIVE_WB", "SHAM_BM"),       # 3. naive - whole blood vs sham - bone marrow
  c("NAIVE_WB", "UT_WB"),         # 4. naive - whole blood vs untreated - whole blood
  c("NAIVE_WB", "UT_BM"),         # 5. naive - whole blood vs untreated - bone marrow
  c("NAIVE_BM", "SHAM_WB"),       # 6. naive - bone marrow vs sham - whole blood
  c("NAIVE_BM", "SHAM_BM"),       # 7. naive - bone marrow vs sham - bone marrow
  c("NAIVE_BM", "UT_WB"),         # 8. naive - bone marrow vs untreated - whole blood
  c("NAIVE_BM", "UT_BM"),         # 9. naive - bone marrow vs untreated - bone marrow
  c("SHAM_WB", "SHAM_BM"),        # 10. sham - whole blood vs sham - bone marrow
  c("SHAM_WB", "UT_WB"),          # 11. sham - whole blood vs untreated - whole blood
  c("SHAM_WB", "UT_BM"),          # 12. sham - whole blood vs untreated - bone marrow
  c("SHAM_BM", "UT_WB"),          # 13. sham - bone marrow vs untreated - whole blood
  c("SHAM_BM", "UT_BM"),          # 14. sham - bone marrow vs untreated - bone marrow
  c("UT_WB", "UT_BM")             # 15. untreated - whole blood vs untreated - bone marrow
)

# Track volcano plots created
volcano_plots_created <- 0

# Perform all pairwise comparisons
for(comp in pairwise_comparisons) {
  comparison_name <- paste(comp[1], "vs", comp[2])
  print(paste("Processing pairwise comparison:", comparison_name))
  
  # Check if both condition-compartment combinations are available
  cc_available <- unique(seu_NKT_focused$condition_compartment)
  
  if(all(comp %in% cc_available)) {
    
    # For each cell type
    for(cell_type in unique(seu_NKT_focused$cell_types)) {
      print(paste("  Cell type:", cell_type))
      
      # Subset to cell type
      seu_subset <- subset(seu_NKT_focused, subset = cell_types == cell_type)
      
      # Further subset to the two condition-compartment combinations being compared
      seu_subset <- subset(seu_subset, subset = condition_compartment %in% comp)
      
      # Check cell counts
      cc_counts <- table(seu_subset$condition_compartment)
      print(paste("    Cell counts:", paste(names(cc_counts), cc_counts, sep = ":", collapse = ", ")))
      
      if(all(cc_counts >= 10)) {
        
        # Set identity to condition_compartment
        Idents(seu_subset) <- "condition_compartment"
        
        # Compare condition-compartment combinations
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
            markers$comparison_type <- "Pairwise_Condition_Compartment"
            markers$comparison <- comparison_name
            markers$group1 <- comp[1]
            markers$group2 <- comp[2]
            markers$cell_type <- cell_type
            markers$method <- "Wilcoxon"
            
            result_name <- paste("Pairwise", gsub("_", "", comp[1]), "vs", gsub("_", "", comp[2]), cell_type, sep = "_")
            all_results[[result_name]] <- markers
            
            print(paste("    Found", nrow(markers), "differential genes"))
          }
        }, error = function(e) {
          print(paste("    Error:", e$message))
        })
      } else {
        print(paste("    Skipping - insufficient cells"))
      }
    }
  } else {
    print(paste("  Skipping", comparison_name, "- condition-compartment combinations not available"))
  }
}

# ============================================================================
# 4. PSEUDOBULK ANALYSIS - BOTH COMPARISON TYPES
# ============================================================================

print("\n=== 4. PSEUDOBULK ANALYSIS ===")

# 4A. Pseudobulk for compartment comparisons (WB vs BM)
print("--- Pseudobulk: WB vs BM ---")

sce_compartment <- as.SingleCellExperiment(seu_NKT_focused)
colData(sce_compartment)$cluster_id <- sce_compartment$cell_types
# FIXED - Use sample_name for original sample ID, then modify for muscat
colData(sce_compartment)$sample_id <- paste(sce_compartment$sample_name, sce_compartment$compartment, sep = "_")
colData(sce_compartment)$group_id <- sce_compartment$compartment
colData(sce_compartment)$condition <- sce_compartment$condition

sce_compartment <- prepSCE(sce_compartment, 
                           kid = "cluster_id",
                           sid = "sample_id", 
                           gid = "group_id",
                           drop = TRUE)

print("Sample summary for compartment comparison:")
print(metadata(sce_compartment)$experiment_info)

pb_compartment <- aggregateData(sce_compartment, assay = "counts", fun = "sum")

# Run muscat for compartment comparison
print("Running muscat edgeR for compartment comparison...")
muscat_compartment_edgeR <- pbDS(pb_compartment, method = "edgeR", filter = "both", min_cells = 10)

print("Running muscat DESeq2 for compartment comparison...")
muscat_compartment_DESeq2 <- pbDS(pb_compartment, method = "DESeq2", filter = "both", min_cells = 10)

# Extract results
muscat_comp_edgeR_res <- resDS(sce_compartment, muscat_compartment_edgeR, bind = "row", frq = FALSE, cmp = FALSE)
muscat_comp_DESeq2_res <- resDS(sce_compartment, muscat_compartment_DESeq2, bind = "row", frq = FALSE, cmp = FALSE)

if(nrow(muscat_comp_edgeR_res) > 0) {
  muscat_comp_edgeR_res$comparison_type <- "Compartment"
  muscat_comp_edgeR_res$method <- "muscat_edgeR"
  all_results[["muscat_compartment_edgeR"]] <- muscat_comp_edgeR_res
  print(paste("Compartment muscat edgeR:", nrow(muscat_comp_edgeR_res), "results"))
}

if(nrow(muscat_comp_DESeq2_res) > 0) {
  muscat_comp_DESeq2_res$comparison_type <- "Compartment"
  muscat_comp_DESeq2_res$method <- "muscat_DESeq2"
  all_results[["muscat_compartment_DESeq2"]] <- muscat_comp_DESeq2_res
  print(paste("Compartment muscat DESeq2:", nrow(muscat_comp_DESeq2_res), "results"))
}

# 4B. Pseudobulk for condition comparisons within each compartment
for(compartment in c("WB", "BM")) {
  print(paste("--- Pseudobulk: Conditions within", compartment, "---"))
  
  # Subset to compartment
  seu_comp <- subset(seu_NKT_focused, subset = compartment == compartment)
  
  # Check if we have multiple conditions
  conditions_available <- unique(seu_comp$condition)
  
  if(length(conditions_available) >= 2) {
    
    sce_condition <- as.SingleCellExperiment(seu_comp)
    colData(sce_condition)$cluster_id <- sce_condition$cell_types
    # FIXED - Use sample_name as the base sample identifier
    colData(sce_condition)$sample_id <- sce_condition$sample_name
    colData(sce_condition)$group_id <- sce_condition$condition
    
    sce_condition <- prepSCE(sce_condition, 
                             kid = "cluster_id",
                             sid = "sample_id", 
                             gid = "group_id",
                             drop = TRUE)
    
    print(paste("Sample summary for", compartment, "condition comparison:"))
    print(metadata(sce_condition)$experiment_info)
    
    pb_condition <- aggregateData(sce_condition, assay = "counts", fun = "sum")
    
    # Run muscat for condition comparison within compartment
    tryCatch({
      print(paste("Running muscat edgeR for conditions within", compartment))
      muscat_cond_edgeR <- pbDS(pb_condition, method = "edgeR", filter = "both", min_cells = 5)
      muscat_cond_edgeR_res <- resDS(sce_condition, muscat_cond_edgeR, bind = "row", frq = FALSE, cmp = FALSE)
      
      if(nrow(muscat_cond_edgeR_res) > 0) {
        muscat_cond_edgeR_res$comparison_type <- "Condition_within_Compartment"
        muscat_cond_edgeR_res$compartment <- compartment
        muscat_cond_edgeR_res$method <- "muscat_edgeR"
        
        result_name <- paste("muscat_condition", compartment, "edgeR", sep = "_")
        all_results[[result_name]] <- muscat_cond_edgeR_res
        
        print(paste("Condition muscat edgeR in", compartment, ":", nrow(muscat_cond_edgeR_res), "results"))
      }
    }, error = function(e) {
      print(paste("Error in muscat for", compartment, ":", e$message))
    })
    
    tryCatch({
      print(paste("Running muscat DESeq2 for conditions within", compartment))
      muscat_cond_DESeq2 <- pbDS(pb_condition, method = "DESeq2", filter = "both", min_cells = 5)
      muscat_cond_DESeq2_res <- resDS(sce_condition, muscat_cond_DESeq2, bind = "row", frq = FALSE, cmp = FALSE)
      
      if(nrow(muscat_cond_DESeq2_res) > 0) {
        muscat_cond_DESeq2_res$comparison_type <- "Condition_within_Compartment"
        muscat_cond_DESeq2_res$compartment <- compartment
        muscat_cond_DESeq2_res$method <- "muscat_DESeq2"
        
        result_name <- paste("muscat_condition", compartment, "DESeq2", sep = "_")
        all_results[[result_name]] <- muscat_cond_DESeq2_res
        
        print(paste("Condition muscat DESeq2 in", compartment, ":", nrow(muscat_cond_DESeq2_res), "results"))
      }
    }, error = function(e) {
      print(paste("Error in muscat DESeq2 for", compartment, ":", e$message))
    })
  }
}

# ============================================================================
# 5. ENHANCED VOLCANO PLOTS WITH GOI EMPHASIS
# ============================================================================

print("\n=== 5. ENHANCED VOLCANO PLOTS WITH GOI EMPHASIS ===")

# Load required libraries
library(ggplot2)
library(ggrepel)
library(dplyr)

# Function to create enhanced volcano plot with GOI emphasis
create_enhanced_volcano_plot <- function(markers_df, comparison_name, cell_type) {
  
  # Handle different p-value columns
  pval_col <- "p_val_adj"
  if(!"p_val_adj" %in% colnames(markers_df)) {
    if("p_adj.loc" %in% colnames(markers_df)) {
      pval_col <- "p_adj.loc"
    } else if("FDR" %in% colnames(markers_df)) {
      pval_col <- "FDR"
    }
  }
  
  # Handle different logFC columns
  logfc_col <- "avg_log2FC"
  if(!"avg_log2FC" %in% colnames(markers_df)) {
    if("logFC" %in% colnames(markers_df)) {
      logfc_col <- "logFC"
    }
  }
  
  # Add volcano plot data
  markers_df$neg_log10_p <- -log10(markers_df[[pval_col]])
  markers_df$significant <- markers_df[[pval_col]] < 0.05 & abs(markers_df[[logfc_col]]) > 0.25
  
  # Identify GOI in the results
  markers_df$is_goi <- markers_df$gene %in% all_goi
  
  # Add GOI category information
  markers_df <- markers_df %>%
    left_join(goi_categories, by = "gene")
  
  # Color coding - prioritize GOI
  markers_df$point_color <- case_when(
    markers_df$is_goi & markers_df[[logfc_col]] > 0.25 & markers_df[[pval_col]] < 0.05 ~ "GOI Up-regulated",
    markers_df$is_goi & markers_df[[logfc_col]] < -0.25 & markers_df[[pval_col]] < 0.05 ~ "GOI Down-regulated",
    markers_df$is_goi ~ "GOI Not Significant",
    markers_df[[logfc_col]] > 0.25 & markers_df[[pval_col]] < 0.05 ~ "Up-regulated",
    markers_df[[logfc_col]] < -0.25 & markers_df[[pval_col]] < 0.05 ~ "Down-regulated",
    TRUE ~ "Not significant"
  )
  
  # Set point sizes - GOI larger
  markers_df$point_size <- ifelse(markers_df$is_goi, 2.5, 1.5)
  
  # Set point alpha - GOI fully opaque
  markers_df$point_alpha <- ifelse(markers_df$is_goi, 1.0, 0.6)
  
  # Color palette
  color_palette <- c(
    "GOI Up-regulated" = "#FF0000",      # Bright red
    "GOI Down-regulated" = "#0000FF",    # Bright blue
    "GOI Not Significant" = "#800080",   # Purple
    "Up-regulated" = "#FF6666",          # Light red
    "Down-regulated" = "#6666FF",        # Light blue
    "Not significant" = "#D3D3D3"        # Light grey
  )
  
  # Create volcano plot
  p <- ggplot(markers_df, aes_string(x = logfc_col, y = "neg_log10_p")) +
    geom_point(aes(color = point_color, size = point_size, alpha = point_alpha)) +
    scale_color_manual(values = color_palette) +
    scale_size_identity() +
    scale_alpha_identity() +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.7) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "black", alpha = 0.7) +
    labs(
      title = paste("Volcano Plot:", comparison_name),
      subtitle = paste("Cell Type:", cell_type, "| GOI:", sum(markers_df$is_goi), 
                      "| Up:", sum(markers_df$point_color == "Up-regulated" | markers_df$point_color == "GOI Up-regulated"),
                      "| Down:", sum(markers_df$point_color == "Down-regulated" | markers_df$point_color == "GOI Down-regulated")),
      x = "Log2 Fold Change",
      y = "-Log10(Adjusted P-value)",
      color = "Gene Type"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 12),
      legend.position = "right",
      legend.title = element_text(face = "bold"),
      panel.grid.minor = element_blank()
    )
  
  # Add labels for ALL GOI that are present
  goi_present <- markers_df %>%
    filter(is_goi) %>%
    arrange(desc(abs(.data[[logfc_col]])))
  
  if(nrow(goi_present) > 0) {
    p <- p + geom_text_repel(
      data = goi_present,
      aes(label = gene),
      size = 3,
      fontface = "bold",
      max.overlaps = 30,
      box.padding = 0.5,
      point.padding = 0.5,
      segment.color = "black",
      segment.size = 0.5,
      min.segment.length = 0.1,
      force = 2,
      seed = 42
    )
    
    print(paste("  -> Labeled", nrow(goi_present), "GOI:", paste(goi_present$gene, collapse = ", ")))
  } else {
    print("  -> No GOI found in this comparison")
  }
  
  return(p)
}

# Create volcano plots for ALL results
print("Creating volcano plots for all comparisons...")

volcano_count <- 0

for(result_name in names(all_results)) {
  result_df <- all_results[[result_name]]
  
  if(nrow(result_df) > 0) {
    # Extract comparison info
    comparison_name <- ifelse("comparison" %in% colnames(result_df), 
                             result_df$comparison[1], 
                             "Unknown_Comparison")
    cell_type <- ifelse("cell_type" %in% colnames(result_df), 
                       result_df$cell_type[1], 
                       ifelse("cluster_id" %in% colnames(result_df),
                              result_df$cluster_id[1],
                              "Unknown_Cell_Type"))
    
    print(paste("Creating volcano plot for:", result_name))
    print(paste("  Comparison:", comparison_name))
    print(paste("  Cell type:", cell_type))
    print(paste("  Total genes:", nrow(result_df)))
    
    # Create volcano plot
    tryCatch({
      p_volcano <- create_enhanced_volcano_plot(result_df, comparison_name, cell_type)
      
      # Create safe filename
      safe_filename <- gsub("[^A-Za-z0-9_]", "_", result_name)
      filename <- paste0("volcano_", safe_filename, ".pdf")
      
      # Save volcano plot
      ggsave(here("R Projects", "1007DGE", "DGEAnalysis-C4S_V4", filename), 
             p_volcano, width = 12, height = 8, dpi = 300)
      
      volcano_count <- volcano_count + 1
      print(paste("  -> Saved:", filename))
      
    }, error = function(e) {
      print(paste("  -> Error creating volcano plot:", e$message))
    })
  }
}

print(paste("=== VOLCANO PLOTS SUMMARY ==="))
print(paste("Total volcano plots created:", volcano_count))

# ============================================================================
# 6. COMPREHENSIVE HEATMAPS FOR ALL METHODS (FIXED)
# ============================================================================

print("\n=== 6. COMPREHENSIVE HEATMAPS FOR ALL METHODS (FIXED) ===")

# Function to create method-specific heatmaps
create_method_heatmap <- function(results_list, method_name, output_dir) {
  
  print(paste("Creating heatmap for method:", method_name))
  
  # Filter results for this method
  method_results <- results_list[grepl(method_name, names(results_list), ignore.case = TRUE)]
  
  if(length(method_results) == 0) {
    print(paste("No results found for method:", method_name))
    return(NULL)
  }
  
  # Collect all significant results for this method
  all_significant <- data.frame()
  
  for(result_name in names(method_results)) {
    result_df <- method_results[[result_name]]
    
    if(nrow(result_df) > 0) {
      # Handle different column names
      pval_col <- "p_val_adj"
      logfc_col <- "avg_log2FC"
      
      if(!"p_val_adj" %in% colnames(result_df)) {
        if("p_adj.loc" %in% colnames(result_df)) {
          pval_col <- "p_adj.loc"
        } else if("FDR" %in% colnames(result_df)) {
          pval_col <- "FDR"
        }
      }
      
      if(!"avg_log2FC" %in% colnames(result_df)) {
        if("logFC" %in% colnames(result_df)) {
          logfc_col <- "logFC"
        }
      }
      
      # Filter for significant genes
      if(pval_col %in% colnames(result_df) && logfc_col %in% colnames(result_df)) {
        sig_genes <- result_df %>%
          filter(.data[[pval_col]] < 0.05 & abs(.data[[logfc_col]]) > 0.25) %>%
          arrange(desc(abs(.data[[logfc_col]]))) %>%
          head(20)  # Top 20 per comparison
        
        if(nrow(sig_genes) > 0) {
          # Ensure we have the columns we need
          sig_genes$result_name <- result_name
          sig_genes$pval_used <- sig_genes[[pval_col]]
          sig_genes$logfc_used <- sig_genes[[logfc_col]]
          
          # Add comparison info if available
          if("comparison_type" %in% colnames(sig_genes)) {
            sig_genes$comparison_type <- sig_genes$comparison_type
          } else {
            sig_genes$comparison_type <- "Unknown"
          }
          
          if("cell_type" %in% colnames(sig_genes)) {
            sig_genes$cell_type <- sig_genes$cell_type
          } else if("cluster_id" %in% colnames(sig_genes)) {
            sig_genes$cell_type <- sig_genes$cluster_id
          } else {
            sig_genes$cell_type <- "Unknown"
          }
          
          all_significant <- rbind(all_significant, sig_genes)
        }
      }
    }
  }
  
  if(nrow(all_significant) > 0) {
    print(paste("Found", nrow(all_significant), "significant results for", method_name))
    
    # Create heatmap data
    heatmap_data <- all_significant %>%
      select(gene, logfc_used, result_name, comparison_type, cell_type) %>%
      # Take top genes overall
      group_by(gene) %>%
      summarise(
        max_abs_fc = max(abs(logfc_used)),
        n_comparisons = n(),
        .groups = "drop"
      ) %>%
      arrange(desc(max_abs_fc)) %>%
      head(50)  # Top 50 genes overall
    
    # Create plot data
    plot_data <- all_significant %>%
      filter(gene %in% heatmap_data$gene) %>%
      mutate(
        comparison_clean = gsub("_", " ", result_name),
        fc_capped = pmax(pmin(logfc_used, 3), -3)  # Cap at ±3 for visualization
      )
    
    # Create heatmap using ggplot2
    tryCatch({
      p_heatmap <- ggplot(plot_data, aes(x = comparison_clean, y = gene, fill = fc_capped)) +
        geom_tile(color = "white", size = 0.1) +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                            midpoint = 0, name = "Log2FC\n(capped ±3)") +
        facet_wrap(~comparison_type, scales = "free_x", ncol = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
              axis.text.y = element_text(size = 8),
              plot.title = element_text(hjust = 0.5),
              strip.text = element_text(size = 10, face = "bold"),
              legend.position = "right") +
        labs(title = paste("Top Differentially Expressed Genes -", method_name),
             x = "Comparison",
             y = "Gene")
      
      # Save heatmap
      filename <- paste0("heatmap_", gsub("[^A-Za-z0-9]", "_", method_name), ".pdf")
      ggsave(file.path(output_dir, filename), p_heatmap, width = 20, height = 16)
      
      print(paste("Created heatmap for", method_name, "->", filename))
      
    }, error = function(e) {
      print(paste("Error creating heatmap for", method_name, ":", e$message))
    })
    
    # Also create a pheatmap version
    tryCatch({
      # Create matrix for pheatmap
      heatmap_matrix <- plot_data %>%
        select(gene, comparison_clean, fc_capped) %>%
        pivot_wider(names_from = comparison_clean, values_from = fc_capped, values_fill = 0) %>%
        column_to_rownames("gene") %>%
        as.matrix()
      
      if(nrow(heatmap_matrix) > 1 && ncol(heatmap_matrix) > 1) {
        # Create pheatmap
        filename_pheatmap <- paste0("pheatmap_", gsub("[^A-Za-z0-9]", "_", method_name), ".pdf")
        pdf(file.path(output_dir, filename_pheatmap), width = 16, height = 12)
        
        pheatmap(heatmap_matrix,
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 breaks = seq(-3, 3, length.out = 101),
                 cluster_rows = TRUE,
                 cluster_cols = TRUE,
                 fontsize_row = 8,
                 fontsize_col = 8,
                 main = paste("Top Differentially Expressed Genes -", method_name),
                 angle_col = 45)
        
        dev.off()
        print(paste("Created pheatmap for", method_name, "->", filename_pheatmap))
      }
      
    }, error = function(e) {
      print(paste("Error creating pheatmap for", method_name, ":", e$message))
    })
    
  } else {
    print(paste("No significant results found for", method_name))
  }
}

# Create heatmaps for each method
output_dir <- here("R Projects", "1007DGE", "DGEAnalysis-C4S_V4")

create_method_heatmap(all_results, "Wilcoxon", output_dir)
create_method_heatmap(all_results, "muscat_edgeR", output_dir)
create_method_heatmap(all_results, "muscat_DESeq2", output_dir)

# ============================================================================
# 7. GENES OF INTEREST ANALYSIS
# ============================================================================

print("\n=== 7. GENES OF INTEREST ANALYSIS ===")

print("Genes of interest:")
print(all_goi)

# Extract GOI from all results with proper handling for different result types
goi_all_results <- list()
goi_summary <- data.frame()

for(result_name in names(all_results)) {
  result_df <- all_results[[result_name]]
  
  print(paste("Processing:", result_name))
  print(paste("  Rows:", nrow(result_df)))
  
  # Handle different result types
  if(grepl("muscat", result_name)) {
    print("  -> This is a muscat result")
    
    # For muscat results, gene names are in the 'gene' column
    if("gene" %in% colnames(result_df)) {
      goi_subset <- result_df[result_df$gene %in% all_goi, ]
      print(paste("  -> Found", nrow(goi_subset), "GOI in muscat results"))
    } else {
      print("  -> No 'gene' column found in muscat results")
      goi_subset <- data.frame() # Empty if no gene column
    }
    
  } else {
    print("  -> This is a Wilcoxon result")
    
    # For Wilcoxon results, gene names are in rownames or 'gene' column
    if("gene" %in% colnames(result_df)) {
      goi_subset <- result_df[result_df$gene %in% all_goi, ]
    } else {
      # Gene names are in rownames
      goi_subset <- result_df[rownames(result_df) %in% all_goi, ]
      if(nrow(goi_subset) > 0) {
        goi_subset$gene <- rownames(goi_subset)
      }
    }
    print(paste("  -> Found", nrow(goi_subset), "GOI in Wilcoxon results"))
  }
  
  # Store if we found any GOI
  if(nrow(goi_subset) > 0) {
    # Add result source info
    goi_subset$result_source <- result_name
    goi_all_results[[result_name]] <- goi_subset
    
    # Add to summary
    goi_summary <- rbind(goi_summary, goi_subset)
    print(paste("  -> Added to summary, total GOI now:", nrow(goi_summary)))
  }
}

# Check what we got
print("\n=== GOI EXTRACTION SUMMARY ===")
print(paste("Total GOI results:", nrow(goi_summary)))

if(nrow(goi_summary) > 0) {
  print("GOI by gene:")
  print(sort(table(goi_summary$gene), decreasing = TRUE))
  
  print("GOI by method:")
  if("method" %in% colnames(goi_summary)) {
    print(table(goi_summary$method))
  }
  
  print("GOI by comparison type:")
  if("comparison_type" %in% colnames(goi_summary)) {
    print(table(goi_summary$comparison_type))
  }
}

# ============================================================================
# 8. COMPREHENSIVE EXCEL EXPORT
# ============================================================================

print("\n=== 8. COMPREHENSIVE EXCEL EXPORT ===")

# Create comprehensive workbook
wb <- createWorkbook()

# Separate results by method
wilcoxon_results <- all_results[!grepl("muscat", names(all_results))]
muscat_edgeR_results <- all_results[grepl("muscat.*edgeR", names(all_results))]
muscat_DESeq2_results <- all_results[grepl("muscat.*DESeq2", names(all_results))]

# Separate by comparison type
compartment_results <- all_results[sapply(all_results, function(x) {
  "comparison_type" %in% colnames(x) && any(x$comparison_type == "Compartment")
})]

condition_results <- all_results[sapply(all_results, function(x) {
  "comparison_type" %in% colnames(x) && any(x$comparison_type == "Condition_within_Compartment")
})]

pairwise_results <- all_results[sapply(all_results, function(x) {
  "comparison_type" %in% colnames(x) && any(x$comparison_type == "Pairwise_Condition_Compartment")
})]

# Summary statistics
summary_stats <- data.frame()

for(result_name in names(all_results)) {
  result_df <- all_results[[result_name]]
  
  if(nrow(result_df) > 0) {
    # Handle different p-value and logFC columns
    pval_col <- "p_val_adj"
    logfc_col <- "avg_log2FC"
    
    if(!"p_val_adj" %in% colnames(result_df)) {
      if("p_adj.loc" %in% colnames(result_df)) {
        pval_col <- "p_adj.loc"
      } else if("FDR" %in% colnames(result_df)) {
        pval_col <- "FDR"
      }
    }
    
    if(!"avg_log2FC" %in% colnames(result_df)) {
      if("logFC" %in% colnames(result_df)) {
        logfc_col <- "logFC"
      }
    }
    
    # Count GOI in results
    goi_in_results <- sum(result_df$gene %in% all_goi)
    goi_significant <- if(pval_col %in% colnames(result_df)) {
      sum(result_df$gene %in% all_goi & result_df[[pval_col]] < 0.05)
    } else {
      0
    }
    
    stats <- data.frame(
      Analysis = result_name,
      Total_Genes = nrow(result_df),
      Significant_Genes = if(pval_col %in% colnames(result_df)) sum(result_df[[pval_col]] < 0.05) else 0,
      Upregulated = if(all(c(pval_col, logfc_col) %in% colnames(result_df))) {
        sum(result_df[[pval_col]] < 0.05 & result_df[[logfc_col]] > 0.25)
      } else 0,
      Downregulated = if(all(c(pval_col, logfc_col) %in% colnames(result_df))) {
        sum(result_df[[pval_col]] < 0.05 & result_df[[logfc_col]] < -0.25)
      } else 0,
      GOI_Present = goi_in_results,
      GOI_Significant = goi_significant,
      Comparison_Type = ifelse("comparison_type" %in% colnames(result_df), result_df$comparison_type[1], "Unknown"),
      Cell_Type = ifelse("cell_type" %in% colnames(result_df), result_df$cell_type[1], 
                        ifelse("cluster_id" %in% colnames(result_df), result_df$cluster_id[1], "Unknown")),
      Method = ifelse("method" %in% colnames(result_df), result_df$method[1], "Unknown")
    )
    
    summary_stats <- rbind(summary_stats, stats)
  }
}

# Add sheets to workbook
addWorksheet(wb, "Summary_Statistics")
writeData(wb, "Summary_Statistics", summary_stats)

# Add method-specific sheets
if(length(wilcoxon_results) > 0) {
  wilcoxon_combined <- do.call(rbind, wilcoxon_results)
  addWorksheet(wb, "Wilcoxon_All_Results")
  writeData(wb, "Wilcoxon_All_Results", wilcoxon_combined)
}

if(length(muscat_edgeR_results) > 0) {
  muscat_edgeR_combined <- do.call(rbind, muscat_edgeR_results)
  addWorksheet(wb, "Muscat_edgeR_All_Results")
  writeData(wb, "Muscat_edgeR_All_Results", muscat_edgeR_combined)
}

if(length(muscat_DESeq2_results) > 0) {
  muscat_DESeq2_combined <- do.call(rbind, muscat_DESeq2_results)
  addWorksheet(wb, "Muscat_DESeq2_All_Results")
  writeData(wb, "Muscat_DESeq2_All_Results", muscat_DESeq2_combined)
}

# Add comparison type sheets
if(length(compartment_results) > 0) {
  compartment_combined <- do.call(rbind, compartment_results)
  addWorksheet(wb, "Compartment_Comparisons")
  writeData(wb, "Compartment_Comparisons", compartment_combined)
}

if(length(condition_results) > 0) {
  condition_combined <- do.call(rbind, condition_results)
  addWorksheet(wb, "Condition_Comparisons")
  writeData(wb, "Condition_Comparisons", condition_combined)
}

if(length(pairwise_results) > 0) {
  pairwise_combined <- do.call(rbind, pairwise_results)
  addWorksheet(wb, "Pairwise_Comparisons")
  writeData(wb, "Pairwise_Comparisons", pairwise_combined)
}

# Add GOI sheets
if(nrow(goi_summary) > 0) {
  addWorksheet(wb, "GOI_All_Results")
  writeData(wb, "GOI_All_Results", goi_summary)
  
  # GOI by method
  goi_wilcoxon <- goi_summary[!grepl("muscat", goi_summary$result_source), ]
  if(nrow(goi_wilcoxon) > 0) {
    addWorksheet(wb, "GOI_Wilcoxon")
    writeData(wb, "GOI_Wilcoxon", goi_wilcoxon)
  }
  
  goi_muscat_edgeR <- goi_summary[grepl("muscat.*edgeR", goi_summary$result_source), ]
  if(nrow(goi_muscat_edgeR) > 0) {
    addWorksheet(wb, "GOI_Muscat_edgeR")
    writeData(wb, "GOI_Muscat_edgeR", goi_muscat_edgeR)
  }
  
  goi_muscat_DESeq2 <- goi_summary[grepl("muscat.*DESeq2", goi_summary$result_source), ]
  if(nrow(goi_muscat_DESeq2) > 0) {
    addWorksheet(wb, "GOI_Muscat_DESeq2")
    writeData(wb, "GOI_Muscat_DESeq2", goi_muscat_DESeq2)
  }
  
  # GOI by comparison type
  goi_compartment <- goi_summary[grepl("Compartment", goi_summary$result_source), ]
  if(nrow(goi_compartment) > 0) {
    addWorksheet(wb, "GOI_Compartment")
    writeData(wb, "GOI_Compartment", goi_compartment)
  }
  
  goi_condition <- goi_summary[grepl("Condition", goi_summary$result_source), ]
  if(nrow(goi_condition) > 0) {
    addWorksheet(wb, "GOI_Condition")
    writeData(wb, "GOI_Condition", goi_condition)
  }
  
  goi_pairwise <- goi_summary[grepl("Pairwise", goi_summary$result_source), ]
  if(nrow(goi_pairwise) > 0) {
    addWorksheet(wb, "GOI_Pairwise")
    writeData(wb, "GOI_Pairwise", goi_pairwise)
  }
}

# Save the comprehensive workbook
saveWorkbook(wb, here("R Projects", "1007DGE", "DGEAnalysis-C4S_V4", "Comprehensive_DGE_Analysis_Results.xlsx"), overwrite = TRUE)

print("=== COMPREHENSIVE EXCEL EXPORT COMPLETE ===")
print("Generated comprehensive Excel file with multiple sheets:")
print("- Summary_Statistics")
print("- Wilcoxon_All_Results")
print("- Muscat_edgeR_All_Results")
print("- Muscat_DESeq2_All_Results")
print("- Compartment_Comparisons")
print("- Condition_Comparisons")
print("- Pairwise_Comparisons")
print("- GOI_All_Results")
print("- GOI_Wilcoxon")
print("- GOI_Muscat_edgeR")
print("- GOI_Muscat_DESeq2")
print("- GOI_Compartment")
print("- GOI_Condition")
print("- GOI_Pairwise")

# ============================================================================
# 9. COMPREHENSIVE SUMMARY REPORT
# ============================================================================

print("\n=== COMPREHENSIVE ANALYSIS COMPLETE ===")
print("Generated files in R Projects/1007DGE/DGEAnalysis-C4S_V4/:")
print("- Comprehensive_DGE_Analysis_Results.xlsx (all results in multiple sheets)")
print("- Individual volcano plots for each comparison")
print("- Method-specific heatmaps (ggplot2 and pheatmap versions)")

print("\n=== FINAL SUMMARY STATISTICS ===")
print(paste("Total analyses performed:", length(all_results)))

# Summary by comparison type
comparison_type_summary <- summary_stats %>%
  group_by(Comparison_Type) %>%
  summarise(
    N_Analyses = n(),
    Total_Genes = sum(Total_Genes),
    Total_Significant = sum(Significant_Genes),
    Total_Upregulated = sum(Upregulated),
    Total_Downregulated = sum(Downregulated),
    Total_GOI_Present = sum(GOI_Present),
    Total_GOI_Significant = sum(GOI_Significant),
    .groups = "drop"
  )

print("Summary by comparison type:")
print(comparison_type_summary)

# Method summary
method_summary <- summary_stats %>%
  group_by(Method) %>%
  summarise(
    N_Analyses = n(),
    Total_Genes = sum(Total_Genes),
    Total_Significant = sum(Significant_Genes),
    Total_GOI_Present = sum(GOI_Present),
    Total_GOI_Significant = sum(GOI_Significant),
    .groups = "drop"
  )

print("Summary by method:")
print(method_summary)

print(paste("Total volcano plots created:", volcano_count))

if(nrow(goi_summary) > 0) {
  print("\nGenes of Interest Summary:")
  print(paste("Total GOI results:", nrow(goi_summary)))
  print("GOI distribution by gene:")
  print(sort(table(goi_summary$gene), decreasing = TRUE))
  
  # Show GOI in pairwise comparisons
  pairwise_goi <- goi_summary[grepl("Pairwise", goi_summary$result_source), ]
  if(nrow(pairwise_goi) > 0) {
    print("\nGOI in pairwise comparisons:")
    print(paste("Total pairwise GOI results:", nrow(pairwise_g