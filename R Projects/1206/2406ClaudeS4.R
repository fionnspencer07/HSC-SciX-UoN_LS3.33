# Focused analysis for Question 2: S1PR and lymphocyte trafficking in DMG
source("initialise_packages.R")

# Custom colors (keeping your original)
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

# Load data
seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))

# FOCUSED ANALYSIS: Only Naive, Sham, UT conditions
conditions_focus <- c("Naive", "Sham", "UT")
sample_types <- c("WB", "BM")

# Filter to focused conditions only
focused_samples <- as.vector(outer(conditions_focus, sample_types, paste, sep = "_"))
seu_focused <- subset(seu_NKT, subset = sample_name %in% focused_samples)

print(paste("Focused dataset contains", ncol(seu_focused), "cells"))
print("Sample distribution:")
table(seu_focused$sample_name)

# Define genes of interest for S1PR/trafficking analysis
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")  # Tissue retention markers
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")  # Egress/circulation markers

all_trafficking_genes <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

# Key comparisons for trafficking analysis
trafficking_comparisons <- c(
  "Naive_BM_vs_Naive_WB",   # Normal baseline trafficking
  "Sham_BM_vs_Sham_WB",    # Surgery control trafficking  
  "UT_BM_vs_UT_WB"         # DMG effect on trafficking
)

# Run focused DEG analysis
focused_deg_results <- list()

for (comp in trafficking_comparisons) {
  pair <- strsplit(comp, "_vs_")[[1]]
  
  message("Running trafficking comparison: ", comp)
  
  res <- wilcoxauc(seu_focused, group_by = "sample_name", groups_use = pair) %>%
    filter(group == pair[1]) %>%
    arrange(group, -logFC)
  
  focused_deg_results[[comp]] <- res
}

# Function to extract and summarize trafficking genes
summarize_trafficking_genes <- function(deg_results, genes_of_interest) {
  
  trafficking_summary <- lapply(names(deg_results), function(comp) {
    res <- deg_results[[comp]]
    trafficking_res <- res %>% 
      filter(feature %in% genes_of_interest) %>%
      select(feature, logFC, pval, padj, pct_in, pct_out) %>%
      mutate(comparison = comp,
             significant = padj < 0.05,
             direction = ifelse(logFC > 0, "Higher_in_BM", "Higher_in_WB"))
    
    return(trafficking_res)
  }) %>% bind_rows()
  
  return(trafficking_summary)
}

# Generate trafficking gene summary
trafficking_summary <- summarize_trafficking_genes(focused_deg_results, all_trafficking_genes)

# Key findings table
key_findings <- trafficking_summary %>%
  filter(significant == TRUE) %>%
  arrange(comparison, -abs(logFC)) %>%
  select(comparison, feature, logFC, padj, direction)

print("Significant trafficking gene changes:")
print(key_findings)

# Function to create focused visualizations
create_trafficking_plots <- function(seu_obj, genes_list, title_suffix = "") {
  
  # 1. Average expression heatmap
  avg_exp <- AverageExpression(seu_obj, 
                               features = genes_list,
                               group.by = "sample_name")$RNA
  
  # Prepare annotation
  sample_info <- data.frame(
    sample = colnames(avg_exp),
    condition = gsub("_(BM|WB)$", "", colnames(avg_exp)),
    compartment = ifelse(grepl("_BM$", colnames(avg_exp)), "Bone_Marrow", "Whole_Blood")
  )
  rownames(sample_info) <- sample_info$sample
  sample_info$condition <- factor(sample_info$condition, levels = conditions_focus)
  
  # Create heatmap
  p1 <- pheatmap(avg_exp,
                 annotation_col = sample_info[, c("condition", "compartment")],
                 scale = "row",
                 cluster_cols = FALSE,
                 main = paste("Trafficking Gene Expression", title_suffix),
                 silent = TRUE)
  
  return(p1)
}

# Create S1PR-focused heatmap
# Get proper average expression (not sums)
get_average_expression_correct <- function(seu_obj, genes_list) {
  
  # Check which genes are present
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  
  if(length(genes_present) == 0) {
    warning("No genes found")
    return(NULL)
  }
  
  cat("Analyzing genes:", paste(genes_present, collapse = ", "), "\n")
  
  # Method 1: Use AverageExpression (older function but gives true averages)
  if(exists("AverageExpression")) {
    avg_exp <- AverageExpression(seu_obj, 
                                 features = genes_present,
                                 group.by = "sample_name")$RNA
  } else {
    # Method 2: Calculate manually
    exp_data <- FetchData(seu_obj, vars = c(genes_present, "sample_name"))
    
    avg_exp <- exp_data %>%
      group_by(sample_name) %>%
      summarise(across(all_of(genes_present), mean, na.rm = TRUE), .groups = "drop") %>%
      column_to_rownames("sample_name") %>%
      t()
  }
  
  return(avg_exp)
}

# Get correct S1PR average expression
s1pr_avg_correct <- get_average_expression_correct(seu_focused, s1pr_genes)

if(!is.null(s1pr_avg_correct)) {
  cat("\nCORRECT S1PR Average Expression:\n")
  print(round(s1pr_avg_correct, 3))
  
  # Create a proper heatmap
  library(pheatmap)
  
  # Prepare sample annotation
  sample_info <- data.frame(
    condition = gsub("_(BM|WB)$", "", colnames(s1pr_avg_correct)),
    compartment = ifelse(grepl("_BM$", colnames(s1pr_avg_correct)), "BM", "WB")
  )
  rownames(sample_info) <- colnames(s1pr_avg_correct)
  
  # Create heatmap
  pheatmap(s1pr_avg_correct, 
           annotation_col = sample_info,
           scale = "row",
           cluster_cols = FALSE,
           main = "S1PR Average Expression (Corrected)",
           cellwidth = 50,
           cellheight = 20)
}

# Let's also look at the individual sample distributions
plot_s1pr1_distribution <- function(seu_obj) {
  
  # Get S1PR1 expression data
  s1pr1_data <- FetchData(seu_obj, vars = c("S1pr1", "sample_name")) %>%
    rownames_to_column("cell_id") %>%
    separate(sample_name, into = c("condition", "compartment"), sep = "_", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = c("Naive", "Sham", "UT")),
           compartment = factor(compartment, levels = c("BM", "WB")))
  
  # Calculate means for each group
  means <- s1pr1_data %>%
    group_by(condition, compartment) %>%
    summarise(mean_exp = mean(S1pr1, na.rm = TRUE), .groups = "drop")
  
  cat("S1PR1 mean expression by group:\n")
  print(means)
  
  # Create boxplot
  p1 <- ggplot(s1pr1_data, aes(x = condition, y = S1pr1, fill = compartment)) +
    geom_violin(alpha = 0.05, outlier.size = 0.3) +
    geom_text(data = means, aes(label = round(mean_exp, 2), y = mean_exp + 0.1),
              position = position_dodge(width = 0.75), size = 3) +
    scale_fill_manual(values = c("BM" = "#ff7f00", "WB" = "#1f78b4")) +
    labs(title = "S1PR1 Expression Distribution",
         subtitle = "Numbers show mean expression",
         x = "Condition", 
         y = "S1PR1 Expression",
         fill = "Compartment") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  
  return(list(plot = p1, means = means))
}

# Generate the corrected S1PR1 analysis
s1pr1_analysis <- plot_s1pr1_distribution(seu_focused)
print(s1pr1_analysis$plot)

# Check if the pattern makes biological sense now
cat("\nBiological interpretation check:\n")
cat("Expected pattern: S1PR1 should be HIGHER in WB than BM (promotes egress)\n")
cat("Your results:\n")
print(s1pr1_analysis$means)

# Calculate BM vs WB ratios for each condition
ratios <- s1pr1_analysis$means %>%
  pivot_wider(names_from = compartment, values_from = mean_exp) %>%
  mutate(WB_vs_BM_ratio = WB / BM) %>%
  select(condition, WB_vs_BM_ratio)

cat("\nWB vs BM expression ratios (>1 means higher in blood, as expected):\n")
print(ratios)

# Create comprehensive trafficking heatmap  
trafficking_heatmap <- create_trafficking_plots(seu_focused, all_trafficking_genes, "- All Trafficking Genes")

# Function to plot specific gene comparisons
plot_gene_expression_boxplot <- function(seu_obj, gene, title = NULL) {
  
  if (is.null(title)) title <- paste("Expression of", gene)
  
  plot_data <- FetchData(seu_obj, vars = c(gene, "sample_name")) %>%
    rownames_to_column("cell_id") %>%
    separate(sample_name, into = c("condition", "compartment"), sep = "_", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = conditions_focus),
           compartment = factor(compartment, levels = c("BM", "WB")))
  
  ggplot(plot_data, aes(x = condition, y = .data[[gene]], fill = compartment)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
    scale_fill_manual(values = c("BM" = "#ff7f00", "WB" = "#1f78b4")) +
    labs(title = title,
         x = "Condition", 
         y = "Expression Level",
         fill = "Compartment") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
}

# Create key gene plots
s1pr1_plot <- plot_gene_expression_boxplot(seu_focused, "S1pr1", "S1PR1 Expression Across Conditions")
ccr7_plot <- plot_gene_expression_boxplot(seu_focused, "Ccr7", "CCR7 Expression Across Conditions")
cd69_plot <- plot_gene_expression_boxplot(seu_focused, "Cd69", "CD69 Expression Across Conditions")

# Statistical summary for presentation
create_results_summary <- function(trafficking_summary) {
  
  # Count significant changes by comparison
  sig_counts <- trafficking_summary %>%
    filter(significant == TRUE) %>%
    group_by(comparison, direction) %>%
    summarise(count = n(), .groups = "drop")
  
  # S1PR specific findings
  s1pr_findings <- trafficking_summary %>%
    filter(feature %in% s1pr_genes, significant == TRUE) %>%
    select(comparison, feature, logFC, padj, direction)
  
  # Key trafficking genes
  key_trafficking_findings <- trafficking_summary %>%
    filter(feature %in% c("S1pr1", "Ccr7", "Cd69", "Klf2"), significant == TRUE) %>%
    select(comparison, feature, logFC, padj, direction)
  
  return(list(
    significant_counts = sig_counts,
    s1pr_findings = s1pr_findings,
    key_findings = key_trafficking_findings
  ))
}

results_summary <- create_results_summary(trafficking_summary)

# Print summary for supervisors
cat("\n=== TRAFFICKING ANALYSIS SUMMARY ===\n")
cat("Research Question: Is S1PR dysregulation in DMG direct or indirect?\n\n")

cat("Significant gene changes by comparison:\n")
print(results_summary$significant_counts)

cat("\nS1PR family findings:\n")
print(results_summary$s1pr_findings)

cat("\nKey trafficking genes:\n")
print(results_summary$key_findings)

# Save key results for presentation
saveRDS(focused_deg_results, "focused_trafficking_deg_results.rds")
saveRDS(trafficking_summary, "trafficking_gene_summary.rds")
saveRDS(results_summary, "results_summary_for_supervisors.rds")

# Display plots
print(s1pr1_plot)
print(ccr7_plot) 
print(cd69_plot)