# Focused analysis for Question 2: S1PR and lymphocyte trafficking in DMG

source("/workspaces/HSC-SciX-UoN_LS3.33/Custom R Functions and Scripts/init_packages.R")

source("/workspaces/HSC-SciX-UoN_LS3.33/Custom R Functions and Scripts/install_packages.R")
# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs", "plasma_cells", "progenitor_cells", "NKT_cells")

# Load data
seu <- readRDS("/workspaces/HSC-SciX-UoN_LS3.33/Data/seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds")
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))

# FOCUSED ANALYSIS: Only Naive, Sham, UT conditions
conditions_focus <- c("Naive", "Sham", "UT")
sample_types <- c("WB", "BM")

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

# Only show missing genes if there are any
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

# Key findings - only show if there are results
key_findings <- trafficking_summary %>%
  filter(significant == TRUE) %>%
  arrange(comparison, -abs(logFC)) %>%
  select(comparison, feature, logFC, padj, direction)

# Only print results if there are significant findings
if(nrow(key_findings) > 0) {
  key_findings
}

# FIXED: ComplexHeatmap function with correct annotation parsing
create_complexheatmap <- function(seu_focused, genes_list, title_suffix = "", 
                                  cluster_genes = TRUE, show_gene_names = TRUE) {
  
  genes_present <- genes_list[genes_list %in% rownames(seu_focused)]
  if(length(genes_present) == 0) return(invisible(NULL))
  
  avg_exp <- AverageExpression(seu_focused, features = genes_present, group.by = "sample_name")$RNA
  avg_exp_scaled <- t(scale(t(avg_exp)))
  avg_exp_scaled[is.nan(avg_exp_scaled)] <- 0
  
  # FIXED: Proper parsing of sample names
  sample_names <- colnames(avg_exp_scaled)
  
  # Parse condition and compartment correctly
  sample_parts <- strsplit(sample_names, "_")
  
  sample_info <- data.frame(
    sample = sample_names,
    condition = sapply(sample_parts, function(x) {
      if(length(x) >= 2) {
        paste(x[1:(length(x)-1)], collapse = "_")  # Everything except the last part
      } else {
        x[1]
      }
    }),
    compartment = sapply(sample_parts, function(x) {
      x[length(x)]  # Last part is compartment
    }),
    stringsAsFactors = FALSE
  )
  
  # Set rownames for annotation
  rownames(sample_info) <- sample_info$sample
  
  # Ensure factors have correct levels
  sample_info$condition <- factor(sample_info$condition, levels = conditions_focus)
  sample_info$compartment <- factor(sample_info$compartment, levels = c("BM", "WB"))
  
  # Create colors
  unique_conditions <- levels(sample_info$condition)
  condition_colors <- brewer.pal(min(length(unique_conditions), 11), "Set3")[1:length(unique_conditions)]
  names(condition_colors) <- unique_conditions
  compartment_colors <- c("BM" = "#ff7f00", "WB" = "#1f78b4")
  
  
  # Create column annotation
  col_anno <- HeatmapAnnotation(
    Condition = sample_info$condition,
    Compartment = sample_info$compartment,
    col = list(Condition = condition_colors, Compartment = compartment_colors),
    annotation_name_side = "left",
    annotation_name_gp = gpar(fontsize = 10),
    simple_anno_size = unit(0.5, "cm")
  )
  
  # Create color scale
  col_fun <- colorRamp2(
    c(min(avg_exp_scaled, na.rm = TRUE), 0, max(avg_exp_scaled, na.rm = TRUE)), 
    c("blue", "white", "red")
  )
  
  # Create heatmap
  ht <- Heatmap(
    avg_exp_scaled,
    name = "Scaled\nExpression",
    col = col_fun,
    top_annotation = col_anno,
    cluster_columns = FALSE,
    show_column_names = TRUE,
    column_names_gp = gpar(fontsize = 9),
    column_names_rot = 45,
    column_title = paste("Trafficking Gene Expression", title_suffix),
    column_title_gp = gpar(fontsize = 12, fontface = "bold"),
    cluster_rows = cluster_genes,
    show_row_names = show_gene_names,
    row_names_gp = gpar(fontsize = 9),
    rect_gp = gpar(col = "white", lwd = 0.5),
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 10, fontface = "bold"),
      labels_gp = gpar(fontsize = 8),
      legend_height = unit(4, "cm")
    )
  )
  
  draw(ht)
  return(invisible(NULL))
}

# Enhanced boxplot function
create_enhanced_boxplot <- function(seu_focused, gene, title = NULL) {
  
  if(!gene %in% rownames(seu_focused)) return(invisible(NULL))
  if (is.null(title)) title <- paste("Expression of", gene)
  
  plot_data <- FetchData(seu_focused, vars = c(gene, "sample_name")) %>%
    rownames_to_column("cell_id") %>%
    separate(sample_name, into = c("condition", "compartment"), sep = "_", remove = FALSE) %>%
    mutate(condition = factor(condition, levels = conditions_focus),
           compartment = factor(compartment, levels = c("BM", "WB")))
  
  stats <- plot_data %>%
    group_by(condition, compartment) %>%
    summarise(mean_exp = mean(.data[[gene]], na.rm = TRUE), 
              n_cells = n(), .groups = "drop")
  
  p <- ggplot(plot_data, aes(x = condition, y = .data[[gene]], fill = compartment)) +
    geom_boxplot(alpha = 0.7, outlier.size = 0.3, width = 0.6) +
    geom_point(position = position_jitterdodge(dodge.width = 0.6, jitter.width = 0.2), 
               alpha = 0.1, size = 0.1) +
    stat_summary(fun = mean, geom = "point", shape = 23, size = 3, 
                 position = position_dodge(width = 0.6), 
                 color = "black", fill = "white") +
    scale_fill_manual(values = c("BM" = "#ff7f00", "WB" = "#1f78b4")) +
    labs(title = title,
         subtitle = paste("Mean expression shown as diamonds | n =", sum(stats$n_cells), "cells"),
         x = "Condition", y = "Expression Level", fill = "Compartment") +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          plot.title = element_text(size = 14, face = "bold"),
          plot.subtitle = element_text(size = 11),
          legend.position = "bottom")
  
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

# S1PR1 biological interpretation - only show results
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
  
  # Show the results
  s1pr1_ratios
}

