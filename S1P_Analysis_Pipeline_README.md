# Comprehensive S1P Receptor and Lymphocyte Trafficking Analysis Pipeline

## Overview

This comprehensive analysis pipeline addresses the research question: **Is the loss of Sphingosine-1-Phosphate (S1P) receptor externalization on T-cells and impaired lymphocyte trafficking from bone marrow a direct or indirect consequence of DMG presence?**

The pipeline includes:
- Multiple differential gene expression methods (Wilcoxon, Muscat+DESeq2, Muscat+edgeR)
- Trajectory analysis with Monocle3 and Slingshot
- Comprehensive visualizations (heatmaps, volcano plots, Nebulosa density plots, etc.)
- Excel outputs with 11 sheets per method
- Statistical analysis and reporting

## Pipeline Components

### 1. Main Analysis Script
- **`S1P_Main_Analysis_Script.R`**: Main execution script that runs the complete pipeline

### 2. Core Analysis Modules
- **`S1P_Comprehensive_Analysis_Pipeline.R`**: Core DGE analysis and Excel output generation
- **`S1P_Trajectory_Analysis_Module.R`**: Trajectory analysis using Monocle3 and Slingshot
- **`S1P_Enhanced_Visualization_Module.R`**: Comprehensive visualization functions

### 3. Supporting Files
- **`Custom R Functions and Scripts/init_packages.R`**: Package initialization
- **`Custom R Functions and Scripts/install_packages.R`**: Package installation

## Required Input Data

The pipeline expects the following input file:
- `Data/seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds`

### Data Requirements
- **Cell types**: NK/T cells from mouse blood/bone marrow
- **Conditions**: NAIVE, SHAM, DMG (untreated DIPG)
- **Compartments**: WB (whole blood) vs BM (bone marrow)
- **Required metadata columns**: 
  - `cell_types`
  - `treatment`
  - `compartment`
  - `sample_id`

## Gene Sets Analyzed

### S1P Pathway Genes
- S1pr1, S1pr2, S1pr3, S1pr4, S1pr5
- Sphk1, Sphk2, Sgpl1

### Trafficking Genes
- Cxcr4, Sell, Ccr7, Itgal, Cd62l, Psgl1, Lfa1

### Bone Marrow Retention Genes
- Cxcr4, Vcam1, Vla4, Sdf1

### T-cell Activation Genes
- Cd69, Cd25, Cd44, Cd95, Klf2

### Chemokine Genes
- Ccl5, Ccl4, Cxcl10, Cxcl9

### Cytokine Genes
- Ifng, Il2, Tnfa, Il10

## Installation and Setup

### 1. Install Required Packages

```r
# Run the package installation script
source("Custom R Functions and Scripts/install_packages.R")
```

### 2. Additional Packages for Trajectory Analysis

```r
# Install trajectory analysis packages
if (!require("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

BiocManager::install(c("monocle3", "slingshot", "tradeSeq", "SingleCellExperiment"))
```

## Usage

### Quick Start

```r
# Run the complete analysis pipeline
source("S1P_Main_Analysis_Script.R")
```

### Customized Execution

```r
# Load the main script functions
source("S1P_Main_Analysis_Script.R")

# Run with custom parameters
results <- execute_comprehensive_s1p_analysis(
  data_file = "path/to/your/data.rds",
  output_base_dir = "Custom_Results",
  run_trajectory_analysis = TRUE,
  run_visualization = TRUE,
  cleanup_workspace = FALSE
)
```

### Running Individual Components

```r
# DGE analysis only
source("S1P_Comprehensive_Analysis_Pipeline.R")

# Trajectory analysis only
source("S1P_Trajectory_Analysis_Module.R")
results <- run_trajectory_analysis_pipeline(seu_obj, genes_of_interest)

# Visualizations only
source("S1P_Enhanced_Visualization_Module.R")
run_comprehensive_visualization_pipeline(seu_obj, dge_results)
```

## Output Structure

The pipeline generates the following output structure:

```
S1P_Comprehensive_Results/
├── Excel_Outputs/
│   ├── DGE_Wilcoxon_Complete_Analysis.xlsx
│   ├── DGE_Muscat_DESeq2_Complete_Analysis.xlsx
│   └── DGE_Muscat_edgeR_Complete_Analysis.xlsx
├── DGE_Results/
│   └── comprehensive_dge_results.rds
├── Trajectory_Analysis/
│   ├── Monocle3_Trajectory_Plots.pdf
│   ├── Slingshot_Trajectory_Plots.pdf
│   ├── Pseudotime_Expression_Plots.pdf
│   └── trajectory_analysis_results.rds
├── Visualizations/
│   ├── S1P_Pathway_Heatmap.pdf
│   ├── Trafficking_Genes_Heatmap.pdf
│   ├── Volcano_Plots_Wilcoxon.pdf
│   ├── Nebulosa_S1P_Genes.pdf
│   ├── Violin_S1P_Genes.pdf
│   └── Correlation_S1P_Trafficking.pdf
├── Statistical_Results/
│   └── comprehensive_statistical_results.rds
└── Summary_Reports/
    ├── comprehensive_analysis_summary.rds
    └── analysis_summary.txt
```

## Excel Output Structure

Each DGE method generates an Excel file with 11 sheets:

1. **All_Genes_DMG_vs_SHAM**: Complete DGE results
2. **All_Genes_DMG_vs_NAIVE**: Complete DGE results
3. **All_Genes_SHAM_vs_NAIVE**: Complete DGE results
4. **S1P_Genes_All_Contrasts**: S1P pathway genes highlighted
5. **Trafficking_Genes_All_Contrasts**: Trafficking genes highlighted
6. **Summary_Stats_By_Contrast**: Statistical summaries
7. **Genes_of_Interest_by_Compartment**: WB vs BM analysis
8. **Top_Upregulated_Each_Contrast**: Top upregulated genes
9. **Top_Downregulated_Each_Contrast**: Top downregulated genes
10. **Trajectory_Associated_Genes**: Pseudotime-correlated genes
11. **Dynamic_S1P_Trafficking_Genes**: Trajectory-specific patterns

## Trajectory Analysis Features

### Monocle3 Analysis
- Pseudotime trajectory inference
- S1P receptor expression dynamics along trajectories
- Trajectory-associated gene identification
- Condition-specific trajectory differences

### Slingshot Analysis
- Multiple lineage trajectory inference
- Branch point analysis
- Lineage-specific expression patterns
- Trajectory weight calculations

### Dynamic Gene Expression
- Time-series-like analysis of S1P pathway genes
- Trafficking marker expression dynamics
- Trajectory-specific gene signatures
- Condition-dependent trajectory disruptions

## Visualization Types

### Heatmaps (Enhanced Implementation)
- Pathway-specific heatmaps
- Clustered by sample and condition
- Split by compartment (WB vs BM)
- Gene category annotations

### Volcano Plots
- Per DGE method and contrast
- S1P genes highlighted (red)
- Trafficking genes highlighted (blue)
- Top genes labeled

### Nebulosa Density Plots
- Individual gene density plots
- Split by treatment and compartment
- Combined S1P signature density

### Violin and Box Plots
- Expression by condition and compartment
- Statistical comparisons
- T-cell subset analysis

### Trajectory Visualizations
- Pseudotime trajectory plots
- Expression along trajectories
- Branch point analysis
- Condition-specific patterns

### Additional Visualizations
- Ridge plots
- Feature plots
- Correlation heatmaps
- Pathway analysis plots

## Statistical Features

- Multiple testing correction (FDR, Bonferroni)
- Effect size calculations (Cohen's d)
- Confidence intervals
- Trajectory statistical testing
- Pseudotime correlation analysis

## Performance Considerations

The pipeline is designed to handle large datasets efficiently:
- Parallel processing support
- Progress indicators
- Memory-efficient operations
- Modular architecture for selective execution

## Troubleshooting

### Common Issues

1. **Missing packages**: Run the installation script first
2. **Data file not found**: Check the path to your RDS file
3. **Memory issues**: Reduce the number of genes or use subsampling
4. **Long execution time**: Run components separately or use parallel processing

### Error Handling

The pipeline includes comprehensive error handling:
- Graceful degradation if components fail
- Detailed error messages
- Continuation of analysis after non-critical errors

## Citation

If you use this pipeline in your research, please cite:

```
S1P Receptor and Lymphocyte Trafficking Analysis Pipeline
Comprehensive analysis of Sphingosine-1-Phosphate receptor externalization 
and lymphocyte trafficking patterns in DMG presence
```

## Support

For questions or issues, please:
1. Check the error logs in the output directory
2. Review the analysis summary report
3. Verify input data format and requirements

## Version Information

- Version: 1.0
- Date: 2024
- R Version: 4.3.3 or higher
- Platform: Cross-platform (Linux, macOS, Windows)

## License

This pipeline is provided as-is for research purposes. Please ensure compliance with all relevant software licenses for the underlying packages used.