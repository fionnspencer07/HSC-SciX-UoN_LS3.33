# Focused analysis for Question 2: S1PR and lymphocyte trafficking in DMG
# Using ComplexHeatmap with improved layout and colors
# Author: fionnspencer07
# Date: 2025-06-25

source("initialise_packages.R")

# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

# Load data
seu <- readRDS("seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))

# FOCUSED ANALYSIS: Only Naive, Sham, UT conditions
conditions_focus <- c("Naive", "Sham", "UT")
sample_types <- c("BM", "WB")  # Reorder so BM comes first

# Filter to focused conditions only
focused_samples <- as.vector(outer(conditions_focus, sample_types, paste, sep = "_"))
seu_focused <- subset(seu_NKT, subset = sample_name %in% focused_samples)

# Define genes of interest for S1PR/trafficking analysis
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
trafficking_genes <- c("Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")
retention_genes <- c("Cd69", "Itgae", "Cxcr4")
egress_genes <- c("S1pr1", "Klf2", "Sell", "Ccr7")

all_trafficking_genes <- unique(c(s1pr_genes, trafficking_genes, retention_genes, egress_genes))

# Check gene availability
genes_present <- all_trafficking_genes[all_trafficking_genes %in% rownames(seu_focused)]
genes_missing <- all_trafficking_genes[!all_trafficking_genes %in% rownames(seu_focused)]

if(length(genes_missing) > 0) {
  warning("Missing genes: ", paste(genes_missing, collapse = ", "))
}

# Update gene lists to only include present genes
s1pr_genes_present <- s1pr_genes[s1pr_genes %in% rownames(seu_focused)]
all_trafficking_genes_present <- genes_present

# Key comparisons for trafficking analysis
trafficking_comparisons <- c(
  "Naive_BM_vs_Naive_WB",
  "Sham_BM_vs_Sham_WB",
  "UT_BM_vs_UT_WB"
)

# Run DEG analysis
for (comp in trafficking_comparisons) {
  pair <- strsplit(comp, "_vs_")[[1]]
  
  res <- suppressMessages(
    wilcoxauc(seu_focused, group_by = "sample_name", groups_use = pair) %>%
      filter(group == pair[1]) %>%
      arrange(group, -logFC)
  )
  
  assign(paste0("deg_", gsub("[^A-Za-z0-9]", "_", comp)), res)
}

# Extract trafficking gene summary
get_trafficking_summary <- function() {
  deg_objects <- ls(pattern = "^deg_", envir = .GlobalEnv)
  
  all_results <- list()
  for(obj_name in deg_objects) {
    comp_name <- gsub("^deg_", "", obj_name)
    
    res <- get(obj_name, envir = .GlobalEnv)
    trafficking_res <- res %>% 
      filter(feature %in% all_trafficking_genes_present) %>%
      select(feature, logFC, pval, padj, pct_in, pct_out) %>%
      mutate(comparison = trafficking_comparisons[which(grepl(gsub("_", ".*", strsplit(comp_name, " vs ")[[1]][1]), trafficking_comparisons))],
             significant = padj < 0.05,
             direction = ifelse(logFC > 0, "Higher_in_BM", "Higher_in_WB"))
    
    all_results[[comp_name]] <- trafficking_res
  }
  
  rm(list = deg_objects, envir = .GlobalEnv)
  return(bind_rows(all_results))
}

trafficking_summary <- get_trafficking_summary()

# Key findings
key_findings <- trafficking_summary %>%
  filter(significant == TRUE) %>%
  arrange(comparison, -abs(logFC)) %>%
  select(comparison, feature, logFC, padj, direction)

if(nrow(key_findings) > 0) {
  key_findings
}

# IMPROVED: ComplexHeatmap with better layout and colors
create_complexheatmap <- function(seu_obj, genes_list, title_suffix = "", 
                                  cluster_genes = TRUE, show_gene_names = TRUE) {
  
  genes_present <- genes_list[genes_list %in% rownames(seu_obj)]
  if(length(genes_present) == 0) return(invisible(NULL))
  
  avg_exp <- AverageExpression(seu_obj, features = genes_present, group.by = "sample_name")$RNA
  avg_exp_scaled <- t(scale(t(avg_exp)))
  avg_exp_scaled[is.nan(avg_exp_scaled)] <- 0
  
  # Get sample names and parse them
  sample_names <- colnames(avg_exp_scaled)
  sample_parts <- strsplit(sample_names, "_")
  
  sample_info <- data.frame(
    sample = sample_names,
    condition = sapply(sample_parts, function(x) {
      if(length(x) >= 2) {
        paste(x[1:(length(x)-1)], collapse = "_")
      } else {
        x[1]
      }
    }),
    compartment = sapply(sample_parts, function(x) {
      x[length(x)]
    }),
    stringsAsFactors = FALSE
  )
  
  # Set rownames and create factors
  rownames(sample_info) <- sample_info$sample
  sample_info$condition <- factor(sample_info$condition, levels = conditions_focus)
  sample_info$compartment <- factor(sample_info$compartment, levels = c("BM", "WB"))
  
  # Define colors - different scheme for expression
  condition_colors <- c(
    "Naive" = "#2E8B57",    # Sea Green
    "Sham" = "#FFD700",     # Gold  
    "UT" = "#DC143C"        # Crimson
  )
  
  compartment_colors <- c(
    "BM" = "#FF7F00",       # Orange
    "WB" = "#1F78B4"        # Blue
  )
  
  # NEW: Different color scale for expression values (viridis-like)
  expression_colors <- colorRamp2(
    c(min(avg_exp_scaled, na.rm = TRUE), 0, max(avg_exp_scaled, na.rm = TRUE)), 
    c("#440154", "#21908C", "#FDE725")  # Purple -> Teal -> Yellow (viridis)
  )
  
  # Order samples: first by condition, then by compartment within each condition
  sample_order <- sample_info %>%
    arrange(condition, compartment) %>%
    pull(sample)
  
  avg_exp_scaled_ordered <- avg_exp_scaled[, sample_order]
  sample_info_ordered <- sample_info[sample_order, ]
  
  # Create column annotation
  col_anno <- HeatmapAnnotation(
    Condition = sample_info_ordered$condition,
    Compartment = sample_info_ordered$compartment,
    col = list(
      Condition = condition_colors,
      Compartment = compartment_colors
    ),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 12, fontface = "bold"),
    annotation_name_rot = 0,
    simple_anno_size = unit(0.6, "cm"),
    gap = unit(2, "mm"),
    annotation_legend_param = list(
      Condition = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                       labels_gp = gpar(fontsize = 10)),
      Compartment = list(title_gp = gpar(fontsize = 12, fontface = "bold"),
                         labels_gp = gpar(fontsize = 10))
    )
  )
  
  # Create custom column labels (condition names with compartment subdivisions)
  condition_labels <- sample_info_ordered$condition
  compartment_labels <- sample_info_ordered$compartment
  
  # Create compound labels showing both condition and compartment subtly
  custom_labels <- paste0(condition_labels, "\n(", compartment_labels, ")")
  
  # Actually, let's just show condition and let the annotation bars show compartment
  simple_labels <- as.character(condition_labels)
  
  # Create the heatmap
  ht <- Heatmap(
    avg_exp_scaled_ordered,
    name = "Expression\nZ-score",
    col = expression_colors,
    
    # Column settings
    top_annotation = col_anno,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_labels = simple_labels,  # Just show conditions
    column_names_gp = gpar(fontsize = 11, fontface = "bold"),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = paste("Trafficking Gene Expression", title_suffix),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    # Row settings
    cluster_rows = cluster_genes,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_side = "right",
    
    # Heatmap body
    rect_gp = gpar(col = "white", lwd = 1),
    
    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(5, "cm"),
      legend_width = unit(1.5, "cm"),
      at = c(-2, -1, 0, 1, 2),
      labels = c("-2", "-1", "0", "1", "2")
    ),
    
    # Border
    border = TRUE
  )
  
  # Add column split to visually separate conditions
  column_split <- sample_info_ordered$condition
  
  ht_with_split <- Heatmap(
    avg_exp_scaled_ordered,
    name = "Expression\nZ-score",
    col = expression_colors,
    
    # Column settings with split
    top_annotation = col_anno,
    cluster_columns = FALSE,
    column_split = column_split,
    column_gap = unit(3, "mm"),
    show_column_names = TRUE,
    column_labels = compartment_labels,  # Show BM/WB within each condition
    column_names_gp = gpar(fontsize = 10, fontface = "bold"),
    column_names_rot = 0,
    column_names_centered = TRUE,
    column_title = paste("Trafficking Gene Expression", title_suffix),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    
    # Row settings
    cluster_rows = cluster_genes,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 10, fontface = "bold"),
    row_names_side = "right",
    
    # Heatmap body
    rect_gp = gpar(col = "white", lwd = 1),
    
    # Legend
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 12, fontface = "bold"),
      labels_gp = gpar(fontsize = 10),
      legend_height = unit(5, "cm"),
      legend_width = unit(1.5, "cm")
    ),
    
    # Border
    border = TRUE
  )
  
  # Draw with proper spacing
  draw(ht_with_split, 
       annotation_legend_side = "right",
       heatmap_legend_side = "right",
       gap = unit(5, "mm"))
  
  return(invisible(NULL))
}

# Enhanced boxplot function with matching colors
create_enhanced_boxplot <- function(seu_obj, gene, title = NULL) {
  
  if(!gene %in% rownames(seu_obj)) return(invisible(NULL))
  if (is.null(title)) title <- paste("Expression of", gene)
  
  plot_data <- FetchData(seu_obj, vars = c(gene, "sample_name")) %>%
    rownames_to_column("cell_id") %>%
    separate(sample_name, into = c("condition", "compartment"), sep = "_", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = conditions_focus),
           compartment = factor(compartment, levels = c("BM", "WB")))
  
  stats <- plot_data %>%
    group_by(condition, compartment) %>%
    summarise(mean_exp = mean(.data[[gene]], na.rm = TRUE), 
              n_cells = n(), .groups = "drop")
  
  # Use the same compartment colors as heatmap
  compartment_colors <- c("BM" = "#FF7F00", "WB" = "#1F78B4")
  
  p <- ggplot(plot_data, aes(x = condition, y = .data[[gene]], fill = compartment)) +
    geom_boxplot(alpha = 0.8, outlier.size = 0.3, width = 0.6) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.2), 
               alpha = 0.1, size = 0.1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                 position = position_dodge(width = 0.6), 
                 color = "black", fill = "white") +
    scale_fill_manual(values = compartment_colors) +
    labs(title = title,
         subtitle = paste("Mean expression shown as diamonds | n =", sum(stats$n_cells), "cells"),
         x = "Condition", y = "Expression Level", fill = "Compartment") +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 0, hjust = 0.5, size = 11),
      axis.text.y = element_text(size = 11),
      axis.title = element_text(size = 12, face = "bold"),
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      legend.position = "bottom",
      legend.title = element_text(size = 12, face = "bold"),
      legend.text = element_text(size = 11),
      panel.grid.major = element_line(alpha = 0.3),
      panel.grid.minor = element_line(alpha = 0.1)
    )
  
  print(p)
  return(invisible(NULL))
}

# Create visualizations
create_complexheatmap(seu_focused, all_trafficking_genes_present, "- All Trafficking Genes")

if(length(s1pr_genes_present) > 0) {
  create_complexheatmap(seu_focused, s1pr_genes_present, "- S1PR Family", cluster_genes = FALSE)
}

egress_genes_present <- egress_genes[egress_genes %in% rownames(seu_focused)]
if(length(egress_genes_present) > 0) {
  create_complexheatmap(seu_focused, egress_genes_present, "- Egress Markers")
}

retention_genes_present <- retention_genes[retention_genes %in% rownames(seu_focused)]
if(length(retention_genes_present) > 0) {
  create_complexheatmap(seu_focused, retention_genes_present, "- Retention Markers")
}

# Create boxplots for key genes
key_genes_for_plots <- c("S1pr1", "Ccr7", "Cd69", "Klf2")
key_genes_available <- key_genes_for_plots[key_genes_for_plots %in% rownames(seu_focused)]

for(gene in key_genes_available) {
  create_enhanced_boxplot(seu_focused, gene, paste(gene, "Expression Across Conditions"))
}

# S1PR1 biological interpretation
if("S1pr1" %in% rownames(seu_focused)) {
  s1pr1_data <- FetchData(seu_focused, vars = c("S1pr1", "sample_name")) %>%
    rownames_to_column("cell_id") %>%
    separate(sample_name, into = c("condition", "compartment"), sep = "_", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = c("Naive", "Sham", "UT")),
           compartment = factor(compartment, levels = c("BM", "WB")))
  
  s1pr1_means <- s1pr1_data %>%
    group_by(condition, compartment) %>%
    summarise(mean_exp = mean(S1pr1, na.rm = TRUE), .groups = "drop")
  
  s1pr1_ratios <- s1pr1_means %>%
    select(condition, compartment, mean_exp) %>%
    pivot_wider(names_from = compartment, values_from = mean_exp) %>%
    mutate(WB_vs_BM_ratio = WB / BM,
           biological_expectation = ifelse(WB_vs_BM_ratio > 1, "Expected", "Unexpected")) %>%
    select(condition, WB_vs_BM_ratio, biological_expectation)
  
  s1pr1_ratios
}