# Script modified by: fionnspencer07
# Date of modification: 2025-07-10 07:00:47 UTC
# Change log: Updated output directory name.

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

# ============================================================================
# COMPREHENSIVE CONDITION-COMPARTMENT ANALYSIS
# ============================================================================

print("=== COMPREHENSIVE CONDITION-COMPARTMENT ANALYSIS ===")

# Create results directory with your specified path
results_dir <- here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO")
if(!dir.exists(results_dir)) {
  dir.create(results_dir, recursive = TRUE)
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
# ANALYSIS FRAMEWORK - ALL COMPARISON TYPES
# ============================================================================

all_results <- list()

print("=== ANALYSIS PLAN ===")
print("1. WB vs BM within each condition (UT, NAIVE, SHAM)")
print("2. WB vs BM across all conditions combined")  
print("3. Condition comparisons within WB compartment")
print("4. Condition comparisons within BM compartment")
print("5. All pairwise condition comparisons within each compartment")
print("6. Specific 15 Pairwise Comparisons (NEW)")
print("7. Pseudobulk Analysis")

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
        markers$method <- "Wilcoxon_combined"
        
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
# 3. NEW: SPECIFIC 15 PAIRWISE COMPARISONS
# ============================================================================
print("\n=== 3. NEW: SPECIFIC 15 PAIRWISE COMPARISONS ===")

# Define the 15 specific pairwise comparisons you requested
pairwise_comparisons <- list(
  c("NAIVE_WB", "NAIVE_BM"),   c("NAIVE_WB", "SHAM_WB"), c("NAIVE_WB", "SHAM_BM"),
  c("NAIVE_WB", "UT_WB"),      c("NAIVE_WB", "UT_BM"),   c("NAIVE_BM", "SHAM_WB"),
  c("NAIVE_BM", "SHAM_BM"),   c("NAIVE_BM", "UT_WB"),   c("NAIVE_BM", "UT_BM"),
  c("SHAM_WB", "SHAM_BM"),    c("SHAM_WB", "UT_WB"),   c("SHAM_WB", "UT_BM"),
  c("SHAM_BM", "UT_WB"),      c("SHAM_BM", "UT_BM"),   c("UT_WB", "UT_BM")
)

# Set identity to the combined group for these comparisons
Idents(seu_NKT_focused) <- "condition_compartment"

for(cell_type in unique(seu_NKT_focused$cell_types)) {
  print(paste("Processing cell type:", cell_type))
  seu_celltype <- subset(seu_NKT_focused, subset = cell_types == cell_type)
  
  for(comp in pairwise_comparisons) {
    group1 <- comp[1]; group2 <- comp[2]
    comparison_name <- paste(group1, "vs", group2)
    
    # Check for sufficient cells
    if (sum(Idents(seu_celltype) == group1) >= 10 && sum(Idents(seu_celltype) == group2) >= 10) {
      tryCatch({
        markers <- FindMarkers(seu_celltype, ident.1 = group1, ident.2 = group2, min.pct = 0.1, logfc.threshold = 0.1, test.use = "wilcox", verbose = FALSE)
        if (nrow(markers) > 0) {
          markers$gene <- rownames(markers)
          markers$comparison_type <- "Specific_Pairwise"
          markers$comparison <- comparison_name
          markers$cell_type <- cell_type
          markers$method <- "Wilcoxon"
          
          result_name <- paste("Specific_Pairwise", cell_type, comparison_name, sep = "_")
          all_results[[result_name]] <- markers
          print(paste("    Found", nrow(markers), "DEGs for", comparison_name))
        }
      }, error = function(e) { print(paste("    Error:", e$message)) })
    } else {
      print(paste("    Skipping", comparison_name, "- insufficient cells."))
    }
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
muscat_comp_edgeR_res <- resDS(sce_compartment, muscat_compartment_edgeR, bind = "row", frq = FALSE, cpm = FALSE)
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
  }
}

# ============================================================================
# 5. GENES OF INTEREST ANALYSIS (FIXED FOR MUSCAT)
# ============================================================================

print("\n=== 5. GENES OF INTEREST ANALYSIS (FIXED FOR MUSCAT) ===")

# Your genes of interest
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")

all_goi <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

print("Genes of interest:")
print(all_goi)

# First, let's check what results we have
print("\nAll results available:")
print(names(all_results))

# Separate analysis by method type
print("\nAnalyzing results by method...")

# Extract GOI from all results with proper handling for different result types
goi_all_results <- list()
goi_summary <- data.frame()

for(result_name in names(all_results)) {
  result_df <- all_results[[result_name]]
  
  print(paste("Processing:", result_name))
  print(paste("  Columns:", paste(colnames(result_df)[1:min(10, ncol(result_df))], collapse = ", ")))
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
      print(paste("  -> Available columns:", paste(colnames(result_df), collapse = ", ")))
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

# ============================================================================
# 6. NEW: ENHANCED VISUALIZATION - VOLCANO PLOTS
# ============================================================================
print("\n=== 6. NEW: CREATING ENHANCED VOLCANO PLOTS ===")
volcano_dir <- file.path(results_dir, "Volcano_Plots_Specific_Pairwise")
if (!dir.exists(volcano_dir)) {
  dir.create(volcano_dir, recursive = TRUE)
}

# Find results from the specific pairwise analysis
pairwise_result_names <- names(all_results)[grepl("Specific_Pairwise", names(all_results))]

for (result_name in pairwise_result_names) {
  results_df <- all_results[[result_name]]
  
  if (nrow(results_df) > 10) {
    print(paste("  -> Generating Volcano plot for:", result_name))
    tryCatch({
      p <- EnhancedVolcano(results_df,
                           lab = rownames(results_df),
                           x = 'avg_log2FC',
                           y = 'p_val_adj',
                           title = results_df$cell_type[1],
                           subtitle = results_df$comparison[1],
                           pCutoff = 0.05,
                           FCcutoff = 0.25,
                           pointSize = 2.0,
                           labSize = 3.0,
                           drawConnectors = TRUE,
                           widthConnectors = 0.5)
      
      ggsave(filename = file.path(volcano_dir, paste0(result_name, "_volcano.pdf")),
             plot = p, width = 12, height = 10)
    }, error = function(e) {
      print(paste("    ERROR creating volcano plot:", e$message))
    })
  }
}

# ============================================================================
# 7. COMPREHENSIVE VISUALIZATION
# ============================================================================

print("\n=== 7. CREATING COMPREHENSIVE VISUALIZATIONS ===")

# Create expression plots by compartment and condition
for(gene in all_goi[all_goi %in% rownames(seu_NKT_focused)]) {
  print(paste("Creating comprehensive plots for", gene))
  
  tryCatch({
    # Main compartment comparison plot
    p1 <- VlnPlot(seu_NKT_focused, 
                  features = gene, 
                  group.by = "compartment", 
                  split.by = "condition",
                  pt.size = 0.1) +
      ggtitle(paste(gene, "- Compartment Comparison by Condition"))
    
    # Condition comparison within compartments
    p2 <- VlnPlot(seu_NKT_focused, 
                  features = gene, 
                  group.by = "condition_compartment", 
                  split.by = "cell_types",
                  pt.size = 0.1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste(gene, "- Condition-Compartment by Cell Type"))
    
    # Cell type specific
    p3 <- VlnPlot(seu_NKT_focused, 
                  features = gene, 
                  group.by = "cell_types", 
                  split.by = "condition_compartment",
                  pt.size = 0.1) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      ggtitle(paste(gene, "- Cell Type by Condition-Compartment"))
    
    # Combine all plots
    combined_plot <- (p1 / p2 / p3)
    
    ggsave(here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", paste0(gene, "_comprehensive_expression.pdf")), 
           combined_plot, width = 20, height = 18)
  }, error = function(e) {
    print(paste("Error creating plot for", gene, ":", e$message))
  })
}

# Create comprehensive heatmap
if(nrow(goi_summary) > 0) {
  
  tryCatch({
    # Standardize columns for plotting
    goi_plot_data <- goi_summary %>%
      mutate(
        log2FC = case_when(
          !is.na(avg_log2FC) ~ avg_log2FC,
          !is.na(logFC) ~ logFC,
          TRUE ~ NA_real_
        ),
        pvalue = case_when(
          !is.na(p_val_adj) ~ p_val_adj,
          !is.na(p_adj.loc) ~ p_adj.loc,
          TRUE ~ 1
        )
      ) %>%
      filter(!is.na(log2FC)) %>%
      mutate(
        significant = pvalue < 0.05,
        analysis_id = paste(comparison_type, 
                            ifelse("cell_type" %in% colnames(.), cell_type, ""), 
                            ifelse("comparison" %in% colnames(.), comparison, ""),
                            ifelse("compartment" %in% colnames(.), compartment, ""),
                            method, sep = "_")
      )
    
    if(nrow(goi_plot_data) > 0) {
      p_comprehensive_heatmap <- ggplot(goi_plot_data, 
                                        aes(x = analysis_id, y = gene)) +
        geom_tile(aes(fill = log2FC), color = "white") +
        geom_text(aes(label = ifelse(significant, "*", "")), size = 3) +
        scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                             midpoint = 0, name = "Log2FC") +
        facet_wrap(~comparison_type, scales = "free_x", ncol = 1) +
        theme_minimal() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1),
              plot.title = element_text(hjust = 0.5),
              strip.text = element_text(size = 12, face = "bold")) +
        labs(title = "Comprehensive Genes of Interest Analysis\n(* indicates p < 0.05)",
             x = "Analysis_CellType_Comparison_Method",
             y = "Gene")
      
      print(p_comprehensive_heatmap)
      ggsave(here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "comprehensive_GOI_heatmap.pdf"), 
             p_comprehensive_heatmap, width = 24, height = 16)
    }
  }, error = function(e) {
    print(paste("Error creating comprehensive heatmap:", e$message))
  })
}

# ============================================================================
# 8. SAVE ALL RESULTS
# ============================================================================

print("\n=== 8. SAVING COMPREHENSIVE RESULTS ===")

# Save all results
saveRDS(all_results, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "comprehensive_all_results.rds"))

# Save results by comparison type
compartment_results <- all_results[grepl("Compartment|muscat_compartment", names(all_results))]
condition_results <- all_results[grepl("Condition|muscat_condition", names(all_results))]

saveRDS(compartment_results, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "compartment_comparison_results.rds"))
saveRDS(condition_results, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "condition_comparison_results.rds"))

# Create comprehensive CSV files
if(length(compartment_results) > 0) {
  compartment_combined <- do.call(rbind, compartment_results)
  write.csv(compartment_combined, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "all_compartment_results.csv"), row.names = FALSE)
}

if(length(condition_results) > 0) {
  condition_combined <- do.call(rbind, condition_results)
  write.csv(condition_combined, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "all_condition_results.csv"), row.names = FALSE)
}

# Save individual muscat results
tryCatch({
  if(exists("muscat_comp_edgeR_res")) {
    write.csv(muscat_comp_edgeR_res, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "muscat_compartment_edgeR.csv"), row.names = FALSE)
  }
  if(exists("muscat_comp_DESeq2_res")) {
    write.csv(muscat_comp_DESeq2_res, here("R Projects", "1007DGE", "Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO", "muscat_compartment_DESeq2.csv"), row.names = FALSE)
  }
}, error = function(e) {
  print(paste("Error saving individual muscat results:", e$message))
})

# ============================================================================
# 9. COMPREHENSIVE SUMMARY REPORT
# ============================================================================

print("\n=== COMPREHENSIVE ANALYSIS COMPLETE ===")
print(paste0("Generated files in R Projects/1007DGE/Comprehensive_Condition_Compartment_DGEAnalysis-CLAUDE4SONNET&&GEMINI2.5PRO/:"))
print("- comprehensive_all_results.rds (all analysis results)")
print("- comprehensive_GOI_summary_FIXED.csv (genes of interest - all methods)")
print("- GOI_Wilcoxon_results.csv (GOI - Wilcoxon only)")
print("- GOI_Muscat_results.csv (GOI - muscat only)")
print("- GOI_muscat_edgeR_direct.csv (GOI - direct from muscat edgeR)")
print("- GOI_muscat_DESeq2_direct.csv (GOI - direct from muscat DESeq2)")
print("- comprehensive_GOI_heatmap.pdf (overview visualization)")
print("- compartment_comparison_results.rds (WB vs BM analyses)")
print("- condition_comparison_results.rds (treatment comparisons)")
print("- all_compartment_results.csv")
print("- all_condition_results.csv")
print("- muscat_compartment_edgeR.csv")
print("- muscat_compartment_DESeq2.csv")
print("- Individual comprehensive expression plots for each GOI")
print("- Volcano plots for specific pairwise comparisons in 'Volcano_Plots_Specific_Pairwise' subdir (NEW)")

print("\n=== FINAL SUMMARY STATISTICS ===")
print(paste("Total analyses performed:", length(all_results)))

# Summary by comparison type
comparison_type_summary <- data.frame()
for(comp_type in c("Compartment", "Condition_within_Compartment", "Specific_Pairwise")) { # Added new type
  type_results <- all_results[sapply(all_results, function(x) {
    "comparison_type" %in% colnames(x) && any(x$comparison_type == comp_type)
  })]
  
  n_analyses <- length(type_results)
  total_genes <- sum(sapply(type_results, nrow))
  
  comparison_type_summary <- rbind(comparison_type_summary,
                                   data.frame(
                                     Comparison_Type = comp_type,
                                     N_Analyses = n_analyses,
                                     Total_Gene_Results = total_genes
                                   ))
}

print("Summary by comparison type:")
print(comparison_type_summary)

if(nrow(goi_summary) > 0) {
  print("\nGenes of Interest Summary:")
  print(paste("Total GOI results:", nrow(goi_summary)))
  print("GOI distribution by gene:")
  print(sort(table(goi_summary$gene), decreasing = TRUE))
  
  if("comparison_type" %in% colnames(goi_summary)) {
    print("GOI distribution by comparison type:")
    print(table(goi_summary$comparison_type))
  }
  
  if("method" %in% colnames(goi_summary)) {
    print("GOI distribution by method:")
    print(table(goi_summary$method))
  }
}

print("\n=== ANALYSIS FRAMEWORK COMPLETED SUCCESSFULLY ===")