# Make sure you have these packages loaded

source("initialise_packages.r")
source("focused_deg-Gemini2.5PRO.R")

# Your Seurat object should be 'seu_focused' from the previous script

# The gene you want to plot
gene_to_plot <- "S1pr1" # Use the correct gene symbol for mouse

# Create the Violin Plot
VlnPlot(seu_focused, 
        features = gene_to_plot, 
        group.by = "sample_name",
        # To ensure the order is logical for comparison
        pt.size = 0.1) + # Adjust point size to make it less busy
  labs(title = "S1pr1 Expression and T-cell Trafficking",
       y = "S1pr1 Expression Level",
       x = "Sample Condition") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  # Add a statistical comparison
  stat_compare_means(comparisons = list(c("UT_BM", "UT_WB"), c("Naive_BM", "Naive_WB")), label = "p.signif")


# Genes of interest for lymphocyte trafficking and activation
genes_for_dotplot <- c("S1pr1", "Cxcr4", "Cd69", "Pdcd1", "Ctla4")

# Create the Dot Plot
DotPlot(seu_focused, 
        features = genes_for_dotplot, 
        group.by = "sample_name") +
  labs(title = "Trafficking and Activation Marker Expression",
       y = "Sample Condition",
       x = "Genes") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))


# Extract the results for your key comparison from the list you created
key_comparison_results <- deg_results[["UT_BM_vs_UT_WB"]]

# Add a column to indicate significance
key_comparison_results$significant <- ifelse(key_comparison_results$padj < 0.05 & abs(key_comparison_results$logFC) > 0.25, "Yes", "No")

# Create the plot
ggplot(key_comparison_results, aes(x = logFC, y = -log10(padj))) +
  geom_point(aes(color = significant), alpha = 0.6) +
  scale_color_manual(values = c("Yes" = "red", "No" = "grey")) +
  labs(title = "DEGs in Bone Marrow vs. Whole Blood (Untreated)",
       x = "Log2 Fold Change",
       y = "-log10 Adjusted P-value") +
  theme_minimal() +
  # Add a label for your key gene
  geom_text_repel(data = subset(key_comparison_results, gene == "S1pr1"),
                  aes(label = gene),
                  box.padding = 0.5,
                  point.padding = 0.5)