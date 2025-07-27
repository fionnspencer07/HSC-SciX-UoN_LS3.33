---
output_dir: R Projects/2107claude/output/rmd_outs
output:
  html_document:
    toc: true
    toc_depth: 6
    toc_float: true
    code_folding: hide
    keep_md: true
    number_sections: true
  pdf_document:
    toc: true
    toc_depth: '6'
---


``` r
library(here) 
```

```
## here() starts at C:/Users/fionn/Newcastle Grammar School/NGS SciX-2025 UoN DMG Cancer Signalling - Bioinformatics and Proteomics/Fionn Git/HSC-SciX-UoN_LS3.33
```

``` r
# Set global knitr chunk options
knitr::opts_chunk$set(
    error     = TRUE,        # Ensures that error messages are being shown in the output
    fig.align = "center",    # Ensures that all figures are aligned to the centre in the output
    message   = TRUE,       # <--- Recommend FALSE for cleaner reports, set to TRUE for development
    warning   = TRUE,       # <--- Recommend FALSE for cleaner reports, set to TRUE for development
    autodep   = TRUE,        # Automatically detects dependencies for caching
    cache     = TRUE,        # <--- IMPORTANT: Enable caching globally here!
    results   = "markup",    # Displays the output as normal formatted text
    echo      = TRUE,        # Ensures that the R code itself is displayed in the document
    dev       = "svg",       # Graphics output will be in Scalable Vector Graphics (SVG) format
    out.width = "90%",       # Sets the width of the output image or plot to 90% of the available width 
    out.extra = "keepaspectratio=true" # Adds extra options to the output, preserving aspect ratio
)

# Set global option for future package (important for large parallel computations like Seurat/MAST)
mib <- 9000*1024^2 # 9GB in bytes
options(future.globals.maxSize = mib) 
```


``` r
library(here)
```

# Set Up Working Environment

## Set Seed


``` r
seed <- "210825"
set.seed(seed)
message("seed has been set to 210825 ")
```

```
## seed has been set to 210825
```

## Load Packages & Functions


``` r
source(here("Custom R Functions and Scripts", "init_packages.R"))
```

```
## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
## ✔ dplyr     1.1.4     ✔ readr     2.1.5
## ✔ forcats   1.0.0     ✔ stringr   1.5.1
## ✔ ggplot2   3.5.2     ✔ tibble    3.3.0
## ✔ lubridate 1.9.4     ✔ tidyr     1.3.1
## ✔ purrr     1.1.0     
## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
## ✖ dplyr::filter() masks stats::filter()
## ✖ dplyr::lag()    masks stats::lag()
## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors
## 
## Attaching package: 'scales'
## 
## 
## The following object is masked from 'package:purrr':
## 
##     discard
## 
## 
## The following object is masked from 'package:readr':
## 
##     col_factor
## 
## 
## 
## Attaching package: 'kableExtra'
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     group_rows
## 
## 
## Loading required package: SeuratObject
## 
## Loading required package: sp
## 
## 'SeuratObject' was built under R 4.5.0 but the current version is
## 4.5.1; it is recomended that you reinstall 'SeuratObject' as the ABI
## for R may have changed
## 
## 
## Attaching package: 'SeuratObject'
## 
## 
## The following object is masked from 'package:DT':
## 
##     JS
## 
## 
## The following objects are masked from 'package:base':
## 
##     intersect, t
## 
## 
## 
## Attaching package: 'Seurat'
## 
## 
## The following object is masked from 'package:DT':
## 
##     JS
## 
## 
## Loading required package: SummarizedExperiment
## 
## Loading required package: MatrixGenerics
## 
## Loading required package: matrixStats
## 
## 
## Attaching package: 'matrixStats'
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     count
## 
## 
## 
## Attaching package: 'MatrixGenerics'
## 
## 
## The following objects are masked from 'package:matrixStats':
## 
##     colAlls, colAnyNAs, colAnys, colAvgsPerRowSet, colCollapse,
##     colCounts, colCummaxs, colCummins, colCumprods, colCumsums,
##     colDiffs, colIQRDiffs, colIQRs, colLogSumExps, colMadDiffs,
##     colMads, colMaxs, colMeans2, colMedians, colMins, colOrderStats,
##     colProds, colQuantiles, colRanges, colRanks, colSdDiffs, colSds,
##     colSums2, colTabulates, colVarDiffs, colVars, colWeightedMads,
##     colWeightedMeans, colWeightedMedians, colWeightedSds,
##     colWeightedVars, rowAlls, rowAnyNAs, rowAnys, rowAvgsPerColSet,
##     rowCollapse, rowCounts, rowCummaxs, rowCummins, rowCumprods,
##     rowCumsums, rowDiffs, rowIQRDiffs, rowIQRs, rowLogSumExps,
##     rowMadDiffs, rowMads, rowMaxs, rowMeans2, rowMedians, rowMins,
##     rowOrderStats, rowProds, rowQuantiles, rowRanges, rowRanks,
##     rowSdDiffs, rowSds, rowSums2, rowTabulates, rowVarDiffs, rowVars,
##     rowWeightedMads, rowWeightedMeans, rowWeightedMedians,
##     rowWeightedSds, rowWeightedVars
## 
## 
## Loading required package: GenomicRanges
## 
## Loading required package: stats4
## 
## Loading required package: BiocGenerics
## 
## Loading required package: generics
## 
## 
## Attaching package: 'generics'
## 
## 
## The following object is masked from 'package:lubridate':
## 
##     as.difftime
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     explain
## 
## 
## The following objects are masked from 'package:base':
## 
##     as.difftime, as.factor, as.ordered, intersect, is.element, setdiff,
##     setequal, union
## 
## 
## 
## Attaching package: 'BiocGenerics'
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     combine
## 
## 
## The following objects are masked from 'package:stats':
## 
##     IQR, mad, sd, var, xtabs
## 
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, aperm, append, as.data.frame, basename, cbind,
##     colnames, dirname, do.call, duplicated, eval, evalq, Filter, Find,
##     get, grep, grepl, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rownames, sapply, saveRDS, table, tapply, unique,
##     unsplit, which.max, which.min
## 
## 
## Loading required package: S4Vectors
## 
## 
## Attaching package: 'S4Vectors'
## 
## 
## The following objects are masked from 'package:lubridate':
## 
##     second, second<-
## 
## 
## The following objects are masked from 'package:dplyr':
## 
##     first, rename
## 
## 
## The following object is masked from 'package:tidyr':
## 
##     expand
## 
## 
## The following object is masked from 'package:utils':
## 
##     findMatches
## 
## 
## The following objects are masked from 'package:base':
## 
##     expand.grid, I, unname
## 
## 
## Loading required package: IRanges
## 
## 
## Attaching package: 'IRanges'
## 
## 
## The following object is masked from 'package:sp':
## 
##     %over%
## 
## 
## The following object is masked from 'package:glue':
## 
##     trim
## 
## 
## The following object is masked from 'package:lubridate':
## 
##     %within%
## 
## 
## The following objects are masked from 'package:dplyr':
## 
##     collapse, desc, slice
## 
## 
## The following object is masked from 'package:purrr':
## 
##     reduce
## 
## 
## The following object is masked from 'package:grDevices':
## 
##     windows
## 
## 
## Loading required package: GenomeInfoDb
## 
## Loading required package: Biobase
## 
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## 
## 
## Attaching package: 'Biobase'
## 
## 
## The following object is masked from 'package:MatrixGenerics':
## 
##     rowMedians
## 
## 
## The following objects are masked from 'package:matrixStats':
## 
##     anyMissing, rowMedians
## 
## 
## 
## Attaching package: 'SummarizedExperiment'
## 
## 
## The following object is masked from 'package:Seurat':
## 
##     Assays
## 
## 
## The following object is masked from 'package:SeuratObject':
## 
##     Assays
## 
## 
## Loading required package: SingleCellExperiment
## 
## Loading required package: scuttle
## 
## 
## Attaching package: 'monocle3'
## 
## 
## The following objects are masked from 'package:Biobase':
## 
##     exprs, fData, fData<-, pData, pData<-
## 
## 
## Loading required package: princurve
## 
## Loading required package: TrajectoryUtils
## 
## 
## Attaching package: 'TrajectoryUtils'
## 
## 
## The following object is masked from 'package:scran':
## 
##     createClusterMST
## 
## 
## Loading required package: Rcpp
## 
## Loading required package: data.table
## 
## 
## Attaching package: 'data.table'
## 
## 
## The following object is masked from 'package:SummarizedExperiment':
## 
##     shift
## 
## 
## The following object is masked from 'package:GenomicRanges':
## 
##     shift
## 
## 
## The following object is masked from 'package:IRanges':
## 
##     shift
## 
## 
## The following objects are masked from 'package:S4Vectors':
## 
##     first, second
## 
## 
## The following objects are masked from 'package:lubridate':
## 
##     hour, isoweek, mday, minute, month, quarter, second, wday, week,
##     yday, year
## 
## 
## The following objects are masked from 'package:dplyr':
## 
##     between, first, last
## 
## 
## The following object is masked from 'package:purrr':
## 
##     transpose
## 
## 
## 
## Attaching package: 'presto'
## 
## 
## The following object is masked from 'package:monocle3':
## 
##     top_markers
## 
## 
## 
## Attaching package: 'glmGamPoi'
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     vars
## 
## 
## The following object is masked from 'package:ggplot2':
## 
##     vars
## 
## 
## 
## Attaching package: 'edge'
## 
## 
## The following object is masked from 'package:monocle3':
## 
##     fit_models
## 
## 
## 
## Attaching package: 'scrapper'
## 
## 
## The following object is masked from 'package:scran':
## 
##     scoreMarkers
## 
## 
## The following objects are masked from 'package:scater':
## 
##     aggregateAcrossCells, normalizeCounts
## 
## 
## The following objects are masked from 'package:scuttle':
## 
##     aggregateAcrossCells, normalizeCounts
## 
## 
## 
## Attaching package: 'DESeq2'
## 
## 
## The following object is masked from 'package:scater':
## 
##     fpkm
## 
## 
## Loading required package: limma
## 
## 
## Attaching package: 'limma'
## 
## 
## The following object is masked from 'package:DESeq2':
## 
##     plotMA
## 
## 
## The following object is masked from 'package:scater':
## 
##     plotMDS
## 
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
## 
## 
## 
## Attaching package: 'edgeR'
## 
## 
## The following object is masked from 'package:SingleCellExperiment':
## 
##     cpm
## 
## 
## 
## 
## clusterProfiler v4.16.0 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
## 
## Please cite:
## 
## Guangchuang Yu, Li-Gen Wang, Yanyan Han and Qing-Yu He.
## clusterProfiler: an R package for comparing biological themes among
## gene clusters. OMICS: A Journal of Integrative Biology. 2012,
## 16(5):284-287
## 
## 
## Attaching package: 'clusterProfiler'
## 
## 
## The following object is masked from 'package:IRanges':
## 
##     slice
## 
## 
## The following object is masked from 'package:S4Vectors':
## 
##     rename
## 
## 
## The following object is masked from 'package:purrr':
## 
##     simplify
## 
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## 
## Loading required package: AnnotationDbi
## 
## 
## Attaching package: 'AnnotationDbi'
## 
## 
## The following object is masked from 'package:clusterProfiler':
## 
##     select
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     select
## 
## 
## 
## 
## 
## Attaching package: 'celldex'
## 
## 
## The following objects are masked from 'package:SingleR':
## 
##     BlueprintEncodeData, DatabaseImmuneCellExpressionData,
##     HumanPrimaryCellAtlasData, ImmGenData, MonacoImmuneData,
##     MouseRNAseqData, NovershternHematopoieticData
## 
## 
## 
## Attaching package: 'cowplot'
## 
## 
## The following object is masked from 'package:patchwork':
## 
##     align_plots
## 
## 
## The following object is masked from 'package:ggpubr':
## 
##     get_legend
## 
## 
## The following object is masked from 'package:lubridate':
## 
##     stamp
## 
## 
## 
## Attaching package: 'gridExtra'
## 
## 
## The following object is masked from 'package:Biobase':
## 
##     combine
## 
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     combine
## 
## 
## The following object is masked from 'package:dplyr':
## 
##     combine
## 
## 
## Loading required package: viridisLite
## 
## 
## Attaching package: 'viridis'
## 
## 
## The following object is masked from 'package:scales':
## 
##     viridis_pal
## 
## 
## Loading required package: grid
## 
## ========================================
## ComplexHeatmap version 2.25.2
## Bioconductor page: http://bioconductor.org/packages/ComplexHeatmap/
## Github page: https://github.com/jokergoo/ComplexHeatmap
## Documentation: http://jokergoo.github.io/ComplexHeatmap-reference
## 
## If you use it in published research, please cite either one:
## - Gu, Z. Complex Heatmap Visualization. iMeta 2022.
## - Gu, Z. Complex heatmaps reveal patterns and correlations in multidimensional 
##     genomic data. Bioinformatics 2016.
## 
## 
## The new InteractiveComplexHeatmap package can directly export static 
## complex heatmaps into an interactive Shiny app with zero effort. Have a try!
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(ComplexHeatmap))
## ========================================
## 
## 
## 
## Attaching package: 'pheatmap'
## 
## 
## The following object is masked from 'package:ComplexHeatmap':
## 
##     pheatmap
## 
## 
## ========================================
## circlize version 0.4.16
## CRAN page: https://cran.r-project.org/package=circlize
## Github page: https://github.com/jokergoo/circlize
## Documentation: https://jokergoo.github.io/circlize_book/book/
## 
## If you use it in published research, please cite:
## Gu, Z. circlize implements and enhances circular visualization
##   in R. Bioinformatics 2014.
## 
## This message can be suppressed by:
##   suppressPackageStartupMessages(library(circlize))
## ========================================
## 
## 
## Loading required package: futile.logger
## 
## 
## Attaching package: 'VennDiagram'
## 
## 
## The following object is masked from 'package:ggpubr':
## 
##     rotate
## 
## 
## 
## Attaching package: 'hdf5r'
## 
## 
## The following object is masked from 'package:SummarizedExperiment':
## 
##     values
## 
## 
## The following object is masked from 'package:GenomicRanges':
## 
##     values
## 
## 
## The following object is masked from 'package:S4Vectors':
## 
##     values
## 
## 
## The following object is masked from 'package:purrr':
## 
##     flatten_df
```

# Configuration

## Custom Colours


``` r
# Custom colors
man_cols <- c("#86b0cc",  "#f3e65d", "#d5c1e7", "#eeb84c", "#82c39e", "#525252","#4d9f6b", "#b3939e", "#e76031", "#e9944b")
names(man_cols) <- c("B_cells", "NK_cells", "monocytes", "T_cells", "neutrophils", "megakaryocytes", "pDCs",
                     "plasma_cells", "progenitor_cells", "NKT_cells")
```

## Output Folders

### Primary


``` r
dir.create(here("R Projects", "2107claude", "output"), recursive = TRUE, showWarnings = FALSE)
```

### DGE


``` r
dir.create(here("R Projects", "2107claude", "output", "de_analysis", "Find Markers"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "de_analysis", "Presto"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "de_analysis", "t.test"), recursive = TRUE, showWarnings = FALSE)
```

### Visualisation


``` r
dir.create(here("R Projects", "2107claude", "output", "visualisation", "box plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "visualisation", "complex heatmaps"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "visualisation", "nebulosa"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "visualisation", "violin plots"), recursive = TRUE, showWarnings = FALSE)
dir.create(here("R Projects", "2107claude", "output", "visualisation", "volcano"), recursive = TRUE, showWarnings = FALSE)
```

# Import & Pre-processing of Data

## Load RDS file

### Load Original RDS

Importing PhD Clara Savary's RDS & Initialising the object as a Seurat Object


``` r
seu <- readRDS(here("Data",          "seu_NKT_blood_mouse_PPK_ONC_harmony_integrated_filtered_v2_20231_20236_annot.rds"))
```

### Subset

Initial Subset to isolate T Cells, Natural Killer T Cells and Natural Killer Cells and label


``` r
seu_NKT <- subset(seu, subset = cell_types %in% c("T_cells", "NK_cells", "NKT_cells"))
```

### Subset Dataset Again

Second Subset -\> subsetting to remove treatment conditions involving theraputics as this is outside the bounds of the scope of the investigation.

Including:

-   DMG Mice (represented as untreated or UT)

-   Control: No Cancer (Treatment and Disease Naive Group)

-   Surgery Control (SHAM - Injected with Saline - to make sure that the stress of surgery isn't a factor and can be compared between Naive, UT and Sham)

Excluding:

-   ONC201 - Dordaviprone -\> DRD2 antagonist and ClpP agonist

    -   <https://doi.org/10.1158/0008-5472.CAN-23-0186>

    -   <https://doi.org/10.1158/2159-8290.CD-23-0131>

-   Dex - Dexmethasone -\> Corticosteroid

-   Dex & ONC201 Combo


``` r
seu_NKT_focused <- subset(seu_NKT, subset = treatment %in% c("UT", "NAIVE", "SHAM"))
```

For later scripting - load seu_NKT_focused

#### Export

Use either base or BioCgenerics

##### BioCgenerics


``` r
BiocGenerics::saveRDS(seu_NKT_focused, file = here("Data", "seu_NKT_focused.rds"))
```

OR

##### Base


``` r
base::saveRDS(seu_NKT_focused, file = here("Data", "seu_NKT_focused.rds") )
```

#### Load rds


``` r
seu_NKT_focused <- readRDS(here("Data", "seu_NKT_focused.rds"))
```

# Differential Gene Expression Analysis

## T Cells


``` r
seu_T <- subset(seu_NKT_focused, subset = T_clusters %in% "T_cells")
```

### Presto

#### Wilcoxon signed-rank test


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- wilcoxauc(seu_T, group_by = "sample_name", groups_use = pairwise_comparisons[[i]])
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(group == pairwise_comparisons[[i]][[1]]) %>%
    filter(feature %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(group, -logFC)
  
  writeData(wb, sheet = i, x = deg_T)
}
```

```
## Warning: The `slot` argument of `GetAssayData()` is deprecated as of SeuratObject 5.0.0.
## ℹ Please use the `layer` argument instead.
## ℹ The deprecated feature was likely used in the presto package.
##   Please report the issue to the authors.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Presto", "DE_presto-wilcoxauc_filtered.xlsx"), overwrite = TRUE)
```

### Find Markers

#### Wilcox Presto


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(seu_T, group.by = "sample_name", ident.1 = pairwise_comparisons[[i]][1], ident.2 = pairwise_comparisons[[i]][2], test.use = "wilcox")
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(-avg_log2FC)
  
  # Add gene names as a column
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  writeData(wb, sheet = i, x = deg_T)
}
```

```
## Warning: `PackageCheck()` was deprecated in SeuratObject 5.0.0.
## ℹ Please use `rlang::check_installed()` instead.
## ℹ The deprecated feature was likely used in the Seurat package.
##   Please report the issue at <https://github.com/satijalab/seurat/issues>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

``` r
saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_wilcox_filtered.xlsx"), overwrite = TRUE)
```

#### Wilcoxon Limma


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(seu_T, group.by = "sample_name", ident.1 = pairwise_comparisons[[i]][1], ident.2 = pairwise_comparisons[[i]][2], test.use = "wilcox_limma")
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(-avg_log2FC)
  
  # Add gene names as a column
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  writeData(wb, sheet = i, x = deg_T)
}

saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_wilcox-limma_filtered.xlsx"), overwrite = TRUE)
```

#### deseq2


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(seu_T, group.by = "sample_name", ident.1 = pairwise_comparisons[[i]][1], ident.2 = pairwise_comparisons[[i]][2], test.use = "DESeq2")
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(-avg_log2FC)
  
  # Add gene names as a column
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  writeData(wb, sheet = i, x = deg_T)
}

saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_deseq2_filtered.xlsx"), overwrite = TRUE)
```

#### edgeR - Empirical Methods for DE (Negative Binomial Distribution)


``` r
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(seu_T, group.by = "sample_name", ident.1 = pairwise_comparisons[[i]][1], ident.2 = pairwise_comparisons[[i]][2], test.use = "edgeR")
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(-avg_log2FC)
  
  # Add gene names as a column
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  writeData(wb, sheet = i, x = deg_T)
}

saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_edgeR_filtered.xlsx"), overwrite = TRUE)
```

#### MAST


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(seu_T, group.by = "sample_name", ident.1 = pairwise_comparisons[[i]][1], ident.2 = pairwise_comparisons[[i]][2], test.use = "MAST")
  
  addWorksheet(wb, sheetName = i)
  
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%   # <- Filter only your selected genes
    arrange(-avg_log2FC)
  
  # Add gene names as a column
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  writeData(wb, sheet = i, x = deg_T)
}
```

```
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
## 
## Done!
```

```
## Combining coefficients and standard errors
```

```
## Calculating log-fold changes
```

```
## Calculating likelihood ratio tests
```

```
## Refitting on reduced model...
```

```
## 
## Done!
```

``` r
saveWorkbook(wb, here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_MAST_filtered.xlsx"), overwrite = TRUE)
```

#### Student's T Test


``` r
# Define the gene list you want to focus on
selected_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5", 
                    "Ccr7", "Sell", "Cd69", "Klf2", "Cxcr4", "Cd44", "Itgae", "Itgal")

pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

# Create an Excel workbook
wb <- createWorkbook()

for (i in names(pairwise_comparisons)) {
  
  deg_T <- FindMarkers(
    seu_T,
    group.by = "sample_name",
    ident.1 = pairwise_comparisons[[i]][1],
    ident.2 = pairwise_comparisons[[i]][2],
    test.use = "t"
  )
  
  addWorksheet(wb, sheetName = i)
  
  # Filter for selected genes and sort by log2FC
  deg_T <- deg_T %>%
    filter(rownames(.) %in% selected_genes) %>%
    arrange(-avg_log2FC)
  
  # Add gene names as a column and reorder
  deg_T$gene <- rownames(deg_T)
  deg_T <- deg_T %>% dplyr::select(gene, dplyr::everything())
  
  # Write to Excel sheet
  writeData(wb, sheet = i, x = deg_T)
}

# Save the workbook
saveWorkbook(
  wb,
  here("R Projects", "2107claude", "output", "de_analysis", "Find Markers", "DE_t_filtered.xlsx"),
  overwrite = TRUE
)
```

### Welch Two Sample t-test

#### S1pr1


``` r
# Define the gene of interest
gene_of_interest <- "S1pr1"

# Define pairwise comparisons
pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

# Initialize list to store results for FDR calculation
all_pvalues <- numeric(length(pairwise_comparisons))
results_list <- list()

# Create workbook
wb <- createWorkbook()

# Loop through each pairwise comparison
for (i in seq_along(pairwise_comparisons)) {
  comparison_name <- names(pairwise_comparisons)[i]
  groups <- pairwise_comparisons[[i]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  cat("Processing comparison:", comparison_name, "\n")
  
  # Extract expression data for the gene of interest
  gene_expression <- FetchData(seu_T, vars = gene_of_interest, layer = "data")
  
  # Get sample names from metadata
  sample_names <- seu_T@meta.data$sample_name
  
  # Extract expression values for each group
  group1_cells <- which(sample_names == group1)
  group2_cells <- which(sample_names == group2)
  
  group1_expression <- gene_expression[group1_cells, gene_of_interest]
  group2_expression <- gene_expression[group2_cells, gene_of_interest]
  
  # Check if both groups have data
  if (length(group1_expression) == 0 || length(group2_expression) == 0) {
    cat("Warning: No cells found for one or both groups in", comparison_name, "\n")
    next
  }
  
  # Perform Welch's t-test
  t_test_result <- t.test(group1_expression, group2_expression, var.equal = FALSE)
  
  # Calculate means
  mean_group1 <- mean(group1_expression, na.rm = TRUE)
  mean_group2 <- mean(group2_expression, na.rm = TRUE)
  
  # Calculate log2 fold change (group1 - group2)
  log2_fc <- mean_group1 - mean_group2
  
  # Extract results
  t_statistic <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  conf_int_lower <- t_test_result$conf.int[1]
  conf_int_upper <- t_test_result$conf.int[2]
  
  # Store p-value for FDR calculation
  all_pvalues[i] <- p_value
  
  # Create results data frame
  results_df <- data.frame(
    Gene = gene_of_interest,
    Comparison = comparison_name,
    Group1 = group1,
    Group2 = group2,
    N_Group1 = length(group1_expression),
    N_Group2 = length(group2_expression),
    Mean_Group1 = mean_group1,
    Mean_Group2 = mean_group2,
    Log2_Fold_Change = log2_fc,
    t_Statistic = as.numeric(t_statistic),
    Degrees_of_Freedom = as.numeric(df),
    P_Value = p_value,
    Conf_Int_Lower = conf_int_lower,
    Conf_Int_Upper = conf_int_upper,
    stringsAsFactors = FALSE
  )
  
  # Store results for later FDR adjustment
  results_list[[i]] <- results_df
  
  # Add worksheet to workbook
  addWorksheet(wb, sheetName = comparison_name)
  writeData(wb, sheet = comparison_name, x = results_df)
}
```

```
## Processing comparison: NAIVE_WB_vs_NAIVE_BM 
## Processing comparison: NAIVE_WB_vs_SHAM_WB 
## Processing comparison: NAIVE_WB_vs_SHAM_BM 
## Processing comparison: NAIVE_WB_vs_UT_WB 
## Processing comparison: NAIVE_WB_vs_UT_BM 
## Processing comparison: NAIVE_BM_vs_SHAM_WB 
## Processing comparison: NAIVE_BM_vs_SHAM_BM 
## Processing comparison: NAIVE_BM_vs_UT_WB 
## Processing comparison: NAIVE_BM_vs_UT_BM 
## Processing comparison: SHAM_WB_vs_SHAM_BM 
## Processing comparison: SHAM_WB_vs_UT_WB 
## Processing comparison: SHAM_WB_vs_UT_BM 
## Processing comparison: SHAM_BM_vs_UT_WB 
## Processing comparison: SHAM_BM_vs_UT_BM 
## Processing comparison: UT_WB_vs_UT_BM
```

``` r
# Calculate FDR-adjusted p-values
fdr_adjusted <- p.adjust(all_pvalues, method = "fdr")

# Add FDR-adjusted p-values to each worksheet
for (i in seq_along(results_list)) {
  if (!is.null(results_list[[i]])) {
    comparison_name <- names(pairwise_comparisons)[i]
    
    # Add FDR column to existing data
    results_list[[i]]$FDR_Adjusted_P_Value <- fdr_adjusted[i]
    
    # Clear and rewrite the worksheet with updated data
    removeWorksheet(wb, comparison_name)
    addWorksheet(wb, sheetName = comparison_name)
    writeData(wb, sheet = comparison_name, x = results_list[[i]])
    
    # Print summary for this comparison
    cat("\n", comparison_name, ":\n")
    cat("  t =", round(results_list[[i]]$t_Statistic, 4), "\n")
    cat("  df =", round(results_list[[i]]$Degrees_of_Freedom, 2), "\n")
    cat("  p-value =", format(results_list[[i]]$P_Value, scientific = TRUE), "\n")
    cat("  FDR p-value =", format(results_list[[i]]$FDR_Adjusted_P_Value, scientific = TRUE), "\n")
    cat("  95% CI: [", round(results_list[[i]]$Conf_Int_Lower, 4), ", ", 
        round(results_list[[i]]$Conf_Int_Upper, 4), "]\n")
    cat("  Log2 FC =", round(results_list[[i]]$Log2_Fold_Change, 4), "\n")
  }
}
```

```
## 
##  NAIVE_WB_vs_NAIVE_BM :
##   t = 3.7211 
##   df = 594.76 
##   p-value = 2.171856e-04 
##   FDR p-value = 2.71482e-04 
##   95% CI: [ 0.1107 ,  0.3582 ]
##   Log2 FC = 0.2345 
## 
##  NAIVE_WB_vs_SHAM_WB :
##   t = -2.4876 
##   df = 955.65 
##   p-value = 1.303036e-02 
##   FDR p-value = 1.503503e-02 
##   95% CI: [ -0.2254 ,  -0.0266 ]
##   Log2 FC = -0.126 
## 
##  NAIVE_WB_vs_SHAM_BM :
##   t = 4.9144 
##   df = 371.12 
##   p-value = 1.338628e-06 
##   FDR p-value = 2.231046e-06 
##   95% CI: [ 0.1958 ,  0.4571 ]
##   Log2 FC = 0.3265 
## 
##  NAIVE_WB_vs_UT_WB :
##   t = 0.3244 
##   df = 626.64 
##   p-value = 7.457086e-01 
##   FDR p-value = 7.457086e-01 
##   95% CI: [ -0.0726 ,  0.1014 ]
##   Log2 FC = 0.0144 
## 
##  NAIVE_WB_vs_UT_BM :
##   t = -8.5054 
##   df = 523.95 
##   p-value = 1.903605e-16 
##   FDR p-value = 7.138518e-16 
##   95% CI: [ -0.4431 ,  -0.2768 ]
##   Log2 FC = -0.36 
## 
##  NAIVE_BM_vs_SHAM_WB :
##   t = -6.3817 
##   df = 481.39 
##   p-value = 4.121104e-10 
##   FDR p-value = 8.830937e-10 
##   95% CI: [ -0.4714 ,  -0.2495 ]
##   Log2 FC = -0.3605 
## 
##  NAIVE_BM_vs_SHAM_BM :
##   t = 1.2963 
##   df = 380.3 
##   p-value = 1.95663e-01 
##   FDR p-value = 2.096389e-01 
##   95% CI: [ -0.0475 ,  0.2316 ]
##   Log2 FC = 0.092 
## 
##  NAIVE_BM_vs_UT_WB :
##   t = -4.3262 
##   df = 328.2 
##   p-value = 2.015908e-05 
##   FDR p-value = 3.023863e-05 
##   95% CI: [ -0.3202 ,  -0.12 ]
##   Log2 FC = -0.2201 
## 
##  NAIVE_BM_vs_UT_BM :
##   t = -12.0936 
##   df = 286.41 
##   p-value = 1.771736e-27 
##   FDR p-value = 1.328802e-26 
##   95% CI: [ -0.6912 ,  -0.4977 ]
##   Log2 FC = -0.5944 
## 
##  SHAM_WB_vs_SHAM_BM :
##   t = 7.5065 
##   df = 280.16 
##   p-value = 8.12291e-13 
##   FDR p-value = 2.030728e-12 
##   95% CI: [ 0.3338 ,  0.5711 ]
##   Log2 FC = 0.4525 
## 
##  SHAM_WB_vs_UT_WB :
##   t = 4.0803 
##   df = 1762.42 
##   p-value = 4.697728e-05 
##   FDR p-value = 6.405993e-05 
##   95% CI: [ 0.0729 ,  0.2079 ]
##   Log2 FC = 0.1404 
## 
##  SHAM_WB_vs_UT_BM :
##   t = -7.3562 
##   df = 1332.51 
##   p-value = 3.294636e-13 
##   FDR p-value = 9.883908e-13 
##   95% CI: [ -0.2964 ,  -0.1716 ]
##   Log2 FC = -0.234 
## 
##  SHAM_BM_vs_UT_WB :
##   t = -5.6689 
##   df = 198.06 
##   p-value = 5.02713e-08 
##   FDR p-value = 9.425868e-08 
##   95% CI: [ -0.4207 ,  -0.2035 ]
##   Log2 FC = -0.3121 
## 
##  SHAM_BM_vs_UT_BM :
##   t = -12.8383 
##   df = 176.32 
##   p-value = 4.612805e-27 
##   FDR p-value = 2.306402e-26 
##   95% CI: [ -0.792 ,  -0.5809 ]
##   Log2 FC = -0.6864 
## 
##  UT_WB_vs_UT_BM :
##   t = -18.5054 
##   df = 5436.21 
##   p-value = 3.408491e-74 
##   FDR p-value = 5.112737e-73 
##   95% CI: [ -0.414 ,  -0.3347 ]
##   Log2 FC = -0.3744
```

``` r
# Save the workbook
output_path <- here("R Projects", "2107claude", "output", "de_analysis", "t.test", paste0(gene_of_interest, ".xlsx"))

# Create directory if it doesn't exist
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Save workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)
```

#### S1pr2


``` r
# Define the gene of interest
gene_of_interest <- "S1pr2"

# Define pairwise comparisons
pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

# Initialize list to store results for FDR calculation
all_pvalues <- numeric(length(pairwise_comparisons))
results_list <- list()

# Create workbook
wb <- createWorkbook()

# Loop through each pairwise comparison
for (i in seq_along(pairwise_comparisons)) {
  comparison_name <- names(pairwise_comparisons)[i]
  groups <- pairwise_comparisons[[i]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  cat("Processing comparison:", comparison_name, "\n")
  
  # Extract expression data for the gene of interest
  gene_expression <- FetchData(seu_T, vars = gene_of_interest, layer = "data")
  
  # Get sample names from metadata
  sample_names <- seu_T@meta.data$sample_name
  
  # Extract expression values for each group
  group1_cells <- which(sample_names == group1)
  group2_cells <- which(sample_names == group2)
  
  group1_expression <- gene_expression[group1_cells, gene_of_interest]
  group2_expression <- gene_expression[group2_cells, gene_of_interest]
  
  # Check if both groups have data
  if (length(group1_expression) == 0 || length(group2_expression) == 0) {
    cat("Warning: No cells found for one or both groups in", comparison_name, "\n")
    next
  }
  
  # Perform Welch's t-test
  t_test_result <- t.test(group1_expression, group2_expression, var.equal = FALSE)
  
  # Calculate means
  mean_group1 <- mean(group1_expression, na.rm = TRUE)
  mean_group2 <- mean(group2_expression, na.rm = TRUE)
  
  # Calculate log2 fold change (group1 - group2)
  log2_fc <- mean_group1 - mean_group2
  
  # Extract results
  t_statistic <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  conf_int_lower <- t_test_result$conf.int[1]
  conf_int_upper <- t_test_result$conf.int[2]
  
  # Store p-value for FDR calculation
  all_pvalues[i] <- p_value
  
  # Create results data frame
  results_df <- data.frame(
    Gene = gene_of_interest,
    Comparison = comparison_name,
    Group1 = group1,
    Group2 = group2,
    N_Group1 = length(group1_expression),
    N_Group2 = length(group2_expression),
    Mean_Group1 = mean_group1,
    Mean_Group2 = mean_group2,
    Log2_Fold_Change = log2_fc,
    t_Statistic = as.numeric(t_statistic),
    Degrees_of_Freedom = as.numeric(df),
    P_Value = p_value,
    Conf_Int_Lower = conf_int_lower,
    Conf_Int_Upper = conf_int_upper,
    stringsAsFactors = FALSE
  )
  
  # Store results for later FDR adjustment
  results_list[[i]] <- results_df
  
  # Add worksheet to workbook
  addWorksheet(wb, sheetName = comparison_name)
  writeData(wb, sheet = comparison_name, x = results_df)
}
```

```
## Processing comparison: NAIVE_WB_vs_NAIVE_BM 
## Processing comparison: NAIVE_WB_vs_SHAM_WB 
## Processing comparison: NAIVE_WB_vs_SHAM_BM 
## Processing comparison: NAIVE_WB_vs_UT_WB 
## Processing comparison: NAIVE_WB_vs_UT_BM 
## Processing comparison: NAIVE_BM_vs_SHAM_WB 
## Processing comparison: NAIVE_BM_vs_SHAM_BM 
## Processing comparison: NAIVE_BM_vs_UT_WB 
## Processing comparison: NAIVE_BM_vs_UT_BM 
## Processing comparison: SHAM_WB_vs_SHAM_BM 
## Processing comparison: SHAM_WB_vs_UT_WB 
## Processing comparison: SHAM_WB_vs_UT_BM 
## Processing comparison: SHAM_BM_vs_UT_WB 
## Processing comparison: SHAM_BM_vs_UT_BM 
## Processing comparison: UT_WB_vs_UT_BM
```

``` r
# Calculate FDR-adjusted p-values
fdr_adjusted <- p.adjust(all_pvalues, method = "fdr")

# Add FDR-adjusted p-values to each worksheet
for (i in seq_along(results_list)) {
  if (!is.null(results_list[[i]])) {
    comparison_name <- names(pairwise_comparisons)[i]
    
    # Add FDR column to existing data
    results_list[[i]]$FDR_Adjusted_P_Value <- fdr_adjusted[i]
    
    # Clear and rewrite the worksheet with updated data
    removeWorksheet(wb, comparison_name)
    addWorksheet(wb, sheetName = comparison_name)
    writeData(wb, sheet = comparison_name, x = results_list[[i]])
    
    # Print summary for this comparison
    cat("\n", comparison_name, ":\n")
    cat("  t =", round(results_list[[i]]$t_Statistic, 4), "\n")
    cat("  df =", round(results_list[[i]]$Degrees_of_Freedom, 2), "\n")
    cat("  p-value =", format(results_list[[i]]$P_Value, scientific = TRUE), "\n")
    cat("  FDR p-value =", format(results_list[[i]]$FDR_Adjusted_P_Value, scientific = TRUE), "\n")
    cat("  95% CI: [", round(results_list[[i]]$Conf_Int_Lower, 4), ", ", 
        round(results_list[[i]]$Conf_Int_Upper, 4), "]\n")
    cat("  Log2 FC =", round(results_list[[i]]$Log2_Fold_Change, 4), "\n")
  }
}
```

```
## 
##  NAIVE_WB_vs_NAIVE_BM :
##   t = -1.6017 
##   df = 329.24 
##   p-value = 1.101748e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0538 ,  0.0055 ]
##   Log2 FC = -0.0242 
## 
##  NAIVE_WB_vs_SHAM_WB :
##   t = -1.1969 
##   df = 1173.01 
##   p-value = 2.315805e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0221 ,  0.0053 ]
##   Log2 FC = -0.0084 
## 
##  NAIVE_WB_vs_SHAM_BM :
##   t = -1.0005 
##   df = 231.12 
##   p-value = 3.181368e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0374 ,  0.0122 ]
##   Log2 FC = -0.0126 
## 
##  NAIVE_WB_vs_UT_WB :
##   t = -0.5922 
##   df = 666.06 
##   p-value = 5.539232e-01 
##   FDR p-value = 6.391422e-01 
##   95% CI: [ -0.0145 ,  0.0078 ]
##   Log2 FC = -0.0034 
## 
##  NAIVE_WB_vs_UT_BM :
##   t = -0.2869 
##   df = 518.41 
##   p-value = 7.743287e-01 
##   FDR p-value = 7.743287e-01 
##   95% CI: [ -0.012 ,  0.009 ]
##   Log2 FC = -0.0015 
## 
##  NAIVE_BM_vs_SHAM_WB :
##   t = 1.0581 
##   df = 318.18 
##   p-value = 2.908202e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0136 ,  0.0452 ]
##   Log2 FC = 0.0158 
## 
##  NAIVE_BM_vs_SHAM_BM :
##   t = 0.6348 
##   df = 420.78 
##   p-value = 5.259343e-01 
##   FDR p-value = 6.391422e-01 
##   95% CI: [ -0.0243 ,  0.0474 ]
##   Log2 FC = 0.0116 
## 
##  NAIVE_BM_vs_UT_WB :
##   t = 1.447 
##   df = 273.77 
##   p-value = 1.490491e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0075 ,  0.0491 ]
##   Log2 FC = 0.0208 
## 
##  NAIVE_BM_vs_UT_BM :
##   t = 1.5901 
##   df = 263.52 
##   p-value = 1.130122e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0054 ,  0.0507 ]
##   Log2 FC = 0.0226 
## 
##  SHAM_WB_vs_SHAM_BM :
##   t = -0.3408 
##   df = 220.01 
##   p-value = 7.336062e-01 
##   FDR p-value = 7.743287e-01 
##   95% CI: [ -0.0287 ,  0.0202 ]
##   Log2 FC = -0.0042 
## 
##  SHAM_WB_vs_UT_WB :
##   t = 0.948 
##   df = 1604.34 
##   p-value = 3.432844e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0053 ,  0.0153 ]
##   Log2 FC = 0.005 
## 
##  SHAM_WB_vs_UT_BM :
##   t = 1.3994 
##   df = 1210.62 
##   p-value = 1.619471e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0027 ,  0.0164 ]
##   Log2 FC = 0.0068 
## 
##  SHAM_BM_vs_UT_WB :
##   t = 0.7869 
##   df = 176.2 
##   p-value = 4.324087e-01 
##   FDR p-value = 6.391422e-01 
##   95% CI: [ -0.0139 ,  0.0323 ]
##   Log2 FC = 0.0092 
## 
##  SHAM_BM_vs_UT_BM :
##   t = 0.9576 
##   df = 166.32 
##   p-value = 3.396496e-01 
##   FDR p-value = 5.721406e-01 
##   95% CI: [ -0.0117 ,  0.0339 ]
##   Log2 FC = 0.0111 
## 
##  UT_WB_vs_UT_BM :
##   t = 0.6763 
##   df = 4825.68 
##   p-value = 4.988886e-01 
##   FDR p-value = 6.391422e-01 
##   95% CI: [ -0.0035 ,  0.0072 ]
##   Log2 FC = 0.0018
```

``` r
# Save the workbook
output_path <- here("R Projects", "2107claude", "output", "de_analysis", "t.test", paste0(gene_of_interest, ".xlsx"))

# Create directory if it doesn't exist
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Save workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)
```

#### S1pr3


``` r
# Define the gene of interest
gene_of_interest <- "S1pr3"

# Define pairwise comparisons
pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

# Initialize list to store results for FDR calculation
all_pvalues <- numeric(length(pairwise_comparisons))
results_list <- list()

# Create workbook
wb <- createWorkbook()

# Loop through each pairwise comparison
for (i in seq_along(pairwise_comparisons)) {
  comparison_name <- names(pairwise_comparisons)[i]
  groups <- pairwise_comparisons[[i]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  cat("Processing comparison:", comparison_name, "\n")
  
  # Extract expression data for the gene of interest
  gene_expression <- FetchData(seu_T, vars = gene_of_interest, layer = "data")
  
  # Get sample names from metadata
  sample_names <- seu_T@meta.data$sample_name
  
  # Extract expression values for each group
  group1_cells <- which(sample_names == group1)
  group2_cells <- which(sample_names == group2)
  
  group1_expression <- gene_expression[group1_cells, gene_of_interest]
  group2_expression <- gene_expression[group2_cells, gene_of_interest]
  
  # Check if both groups have data
  if (length(group1_expression) == 0 || length(group2_expression) == 0) {
    cat("Warning: No cells found for one or both groups in", comparison_name, "\n")
    next
  }
  
  # Perform Welch's t-test
  t_test_result <- t.test(group1_expression, group2_expression, var.equal = FALSE)
  
  # Calculate means
  mean_group1 <- mean(group1_expression, na.rm = TRUE)
  mean_group2 <- mean(group2_expression, na.rm = TRUE)
  
  # Calculate log2 fold change (group1 - group2)
  log2_fc <- mean_group1 - mean_group2
  
  # Extract results
  t_statistic <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  conf_int_lower <- t_test_result$conf.int[1]
  conf_int_upper <- t_test_result$conf.int[2]
  
  # Store p-value for FDR calculation
  all_pvalues[i] <- p_value
  
  # Create results data frame
  results_df <- data.frame(
    Gene = gene_of_interest,
    Comparison = comparison_name,
    Group1 = group1,
    Group2 = group2,
    N_Group1 = length(group1_expression),
    N_Group2 = length(group2_expression),
    Mean_Group1 = mean_group1,
    Mean_Group2 = mean_group2,
    Log2_Fold_Change = log2_fc,
    t_Statistic = as.numeric(t_statistic),
    Degrees_of_Freedom = as.numeric(df),
    P_Value = p_value,
    Conf_Int_Lower = conf_int_lower,
    Conf_Int_Upper = conf_int_upper,
    stringsAsFactors = FALSE
  )
  
  # Store results for later FDR adjustment
  results_list[[i]] <- results_df
  
  # Add worksheet to workbook
  addWorksheet(wb, sheetName = comparison_name)
  writeData(wb, sheet = comparison_name, x = results_df)
}
```

```
## Processing comparison: NAIVE_WB_vs_NAIVE_BM 
## Processing comparison: NAIVE_WB_vs_SHAM_WB 
## Processing comparison: NAIVE_WB_vs_SHAM_BM 
## Processing comparison: NAIVE_WB_vs_UT_WB 
## Processing comparison: NAIVE_WB_vs_UT_BM 
## Processing comparison: NAIVE_BM_vs_SHAM_WB 
## Processing comparison: NAIVE_BM_vs_SHAM_BM 
## Processing comparison: NAIVE_BM_vs_UT_WB 
## Processing comparison: NAIVE_BM_vs_UT_BM 
## Processing comparison: SHAM_WB_vs_SHAM_BM 
## Processing comparison: SHAM_WB_vs_UT_WB 
## Processing comparison: SHAM_WB_vs_UT_BM 
## Processing comparison: SHAM_BM_vs_UT_WB 
## Processing comparison: SHAM_BM_vs_UT_BM 
## Processing comparison: UT_WB_vs_UT_BM
```

``` r
# Calculate FDR-adjusted p-values
fdr_adjusted <- p.adjust(all_pvalues, method = "fdr")

# Add FDR-adjusted p-values to each worksheet
for (i in seq_along(results_list)) {
  if (!is.null(results_list[[i]])) {
    comparison_name <- names(pairwise_comparisons)[i]
    
    # Add FDR column to existing data
    results_list[[i]]$FDR_Adjusted_P_Value <- fdr_adjusted[i]
    
    # Clear and rewrite the worksheet with updated data
    removeWorksheet(wb, comparison_name)
    addWorksheet(wb, sheetName = comparison_name)
    writeData(wb, sheet = comparison_name, x = results_list[[i]])
    
    # Print summary for this comparison
    cat("\n", comparison_name, ":\n")
    cat("  t =", round(results_list[[i]]$t_Statistic, 4), "\n")
    cat("  df =", round(results_list[[i]]$Degrees_of_Freedom, 2), "\n")
    cat("  p-value =", format(results_list[[i]]$P_Value, scientific = TRUE), "\n")
    cat("  FDR p-value =", format(results_list[[i]]$FDR_Adjusted_P_Value, scientific = TRUE), "\n")
    cat("  95% CI: [", round(results_list[[i]]$Conf_Int_Lower, 4), ", ", 
        round(results_list[[i]]$Conf_Int_Upper, 4), "]\n")
    cat("  Log2 FC =", round(results_list[[i]]$Log2_Fold_Change, 4), "\n")
  }
}
```

```
## 
##  NAIVE_WB_vs_NAIVE_BM :
##   t = NaN 
##   df = NaN 
##   p-value = NaN 
##   FDR p-value = NaN 
##   95% CI: [ NaN ,  NaN ]
##   Log2 FC = 0 
## 
##  NAIVE_WB_vs_SHAM_WB :
##   t = -1 
##   df = 1041 
##   p-value = 3.175429e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0042 ,  0.0013 ]
##   Log2 FC = -0.0014 
## 
##  NAIVE_WB_vs_SHAM_BM :
##   t = -1 
##   df = 162 
##   p-value = 3.188018e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0222 ,  0.0073 ]
##   Log2 FC = -0.0075 
## 
##  NAIVE_WB_vs_UT_WB :
##   t = NaN 
##   df = NaN 
##   p-value = NaN 
##   FDR p-value = NaN 
##   95% CI: [ NaN ,  NaN ]
##   Log2 FC = 0 
## 
##  NAIVE_WB_vs_UT_BM :
##   t = -1 
##   df = 7129 
##   p-value = 3.173444e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0005 ,  0.0002 ]
##   Log2 FC = -0.0002 
## 
##  NAIVE_BM_vs_SHAM_WB :
##   t = -1 
##   df = 1041 
##   p-value = 3.175429e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0042 ,  0.0013 ]
##   Log2 FC = -0.0014 
## 
##  NAIVE_BM_vs_SHAM_BM :
##   t = -1 
##   df = 162 
##   p-value = 3.188018e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0222 ,  0.0073 ]
##   Log2 FC = -0.0075 
## 
##  NAIVE_BM_vs_UT_WB :
##   t = NaN 
##   df = NaN 
##   p-value = NaN 
##   FDR p-value = NaN 
##   95% CI: [ NaN ,  NaN ]
##   Log2 FC = 0 
## 
##  NAIVE_BM_vs_UT_BM :
##   t = -1 
##   df = 7129 
##   p-value = 3.173444e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0005 ,  0.0002 ]
##   Log2 FC = -0.0002 
## 
##  SHAM_WB_vs_SHAM_BM :
##   t = -0.7984 
##   df = 173.58 
##   p-value = 4.257429e-01 
##   FDR p-value = 4.257429e-01 
##   95% CI: [ -0.0211 ,  0.0089 ]
##   Log2 FC = -0.0061 
## 
##  SHAM_WB_vs_UT_WB :
##   t = 1 
##   df = 1041 
##   p-value = 3.175429e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0013 ,  0.0042 ]
##   Log2 FC = 0.0014 
## 
##  SHAM_WB_vs_UT_BM :
##   t = 0.8738 
##   df = 1071.14 
##   p-value = 3.824271e-01 
##   FDR p-value = 4.171932e-01 
##   95% CI: [ -0.0015 ,  0.004 ]
##   Log2 FC = 0.0012 
## 
##  SHAM_BM_vs_UT_WB :
##   t = 1 
##   df = 162 
##   p-value = 3.188018e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0073 ,  0.0222 ]
##   Log2 FC = 0.0075 
## 
##  SHAM_BM_vs_UT_BM :
##   t = 0.9772 
##   df = 162.16 
##   p-value = 3.299055e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0075 ,  0.0221 ]
##   Log2 FC = 0.0073 
## 
##  UT_WB_vs_UT_BM :
##   t = -1 
##   df = 7129 
##   p-value = 3.173444e-01 
##   FDR p-value = 3.958866e-01 
##   95% CI: [ -0.0005 ,  0.0002 ]
##   Log2 FC = -0.0002
```

``` r
# Save the workbook
output_path <- here("R Projects", "2107claude", "output", "de_analysis", "t.test", paste0(gene_of_interest, ".xlsx"))

# Create directory if it doesn't exist
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Save workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)
```

#### S1pr4


``` r
# Define the gene of interest
gene_of_interest <- "S1pr4"

# Define pairwise comparisons
pairwise_comparisons <- list(
  "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
  "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
  "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
  "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
  "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
  "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
  "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
  "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
  "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
  "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
  "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
  "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
  "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
  "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
  "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
)

# Initialize list to store results for FDR calculation
all_pvalues <- numeric(length(pairwise_comparisons))
results_list <- list()

# Create workbook
wb <- createWorkbook()

# Loop through each pairwise comparison
for (i in seq_along(pairwise_comparisons)) {
  comparison_name <- names(pairwise_comparisons)[i]
  groups <- pairwise_comparisons[[i]]
  group1 <- groups[1]
  group2 <- groups[2]
  
  cat("Processing comparison:", comparison_name, "\n")
  
  # Extract expression data for the gene of interest
  gene_expression <- FetchData(seu_T, vars = gene_of_interest, layer = "data")
  
  # Get sample names from metadata
  sample_names <- seu_T@meta.data$sample_name
  
  # Extract expression values for each group
  group1_cells <- which(sample_names == group1)
  group2_cells <- which(sample_names == group2)
  
  group1_expression <- gene_expression[group1_cells, gene_of_interest]
  group2_expression <- gene_expression[group2_cells, gene_of_interest]
  
  # Check if both groups have data
  if (length(group1_expression) == 0 || length(group2_expression) == 0) {
    cat("Warning: No cells found for one or both groups in", comparison_name, "\n")
    next
  }
  
  # Perform Welch's t-test
  t_test_result <- t.test(group1_expression, group2_expression, var.equal = FALSE)
  
  # Calculate means
  mean_group1 <- mean(group1_expression, na.rm = TRUE)
  mean_group2 <- mean(group2_expression, na.rm = TRUE)
  
  # Calculate log2 fold change (group1 - group2)
  log2_fc <- mean_group1 - mean_group2
  
  # Extract results
  t_statistic <- t_test_result$statistic
  df <- t_test_result$parameter
  p_value <- t_test_result$p.value
  conf_int_lower <- t_test_result$conf.int[1]
  conf_int_upper <- t_test_result$conf.int[2]
  
  # Store p-value for FDR calculation
  all_pvalues[i] <- p_value
  
  # Create results data frame
  results_df <- data.frame(
    Gene = gene_of_interest,
    Comparison = comparison_name,
    Group1 = group1,
    Group2 = group2,
    N_Group1 = length(group1_expression),
    N_Group2 = length(group2_expression),
    Mean_Group1 = mean_group1,
    Mean_Group2 = mean_group2,
    Log2_Fold_Change = log2_fc,
    t_Statistic = as.numeric(t_statistic),
    Degrees_of_Freedom = as.numeric(df),
    P_Value = p_value,
    Conf_Int_Lower = conf_int_lower,
    Conf_Int_Upper = conf_int_upper,
    stringsAsFactors = FALSE
  )
  
  # Store results for later FDR adjustment
  results_list[[i]] <- results_df
  
  # Add worksheet to workbook
  addWorksheet(wb, sheetName = comparison_name)
  writeData(wb, sheet = comparison_name, x = results_df)
}
```

```
## Processing comparison: NAIVE_WB_vs_NAIVE_BM 
## Processing comparison: NAIVE_WB_vs_SHAM_WB 
## Processing comparison: NAIVE_WB_vs_SHAM_BM 
## Processing comparison: NAIVE_WB_vs_UT_WB 
## Processing comparison: NAIVE_WB_vs_UT_BM 
## Processing comparison: NAIVE_BM_vs_SHAM_WB 
## Processing comparison: NAIVE_BM_vs_SHAM_BM 
## Processing comparison: NAIVE_BM_vs_UT_WB 
## Processing comparison: NAIVE_BM_vs_UT_BM 
## Processing comparison: SHAM_WB_vs_SHAM_BM 
## Processing comparison: SHAM_WB_vs_UT_WB 
## Processing comparison: SHAM_WB_vs_UT_BM 
## Processing comparison: SHAM_BM_vs_UT_WB 
## Processing comparison: SHAM_BM_vs_UT_BM 
## Processing comparison: UT_WB_vs_UT_BM
```

``` r
# Calculate FDR-adjusted p-values
fdr_adjusted <- p.adjust(all_pvalues, method = "fdr")

# Add FDR-adjusted p-values to each worksheet
for (i in seq_along(results_list)) {
  if (!is.null(results_list[[i]])) {
    comparison_name <- names(pairwise_comparisons)[i]
    
    # Add FDR column to existing data
    results_list[[i]]$FDR_Adjusted_P_Value <- fdr_adjusted[i]
    
    # Clear and rewrite the worksheet with updated data
    removeWorksheet(wb, comparison_name)
    addWorksheet(wb, sheetName = comparison_name)
    writeData(wb, sheet = comparison_name, x = results_list[[i]])
    
    # Print summary for this comparison
    cat("\n", comparison_name, ":\n")
    cat("  t =", round(results_list[[i]]$t_Statistic, 4), "\n")
    cat("  df =", round(results_list[[i]]$Degrees_of_Freedom, 2), "\n")
    cat("  p-value =", format(results_list[[i]]$P_Value, scientific = TRUE), "\n")
    cat("  FDR p-value =", format(results_list[[i]]$FDR_Adjusted_P_Value, scientific = TRUE), "\n")
    cat("  95% CI: [", round(results_list[[i]]$Conf_Int_Lower, 4), ", ", 
        round(results_list[[i]]$Conf_Int_Upper, 4), "]\n")
    cat("  Log2 FC =", round(results_list[[i]]$Log2_Fold_Change, 4), "\n")
  }
}
```

```
## 
##  NAIVE_WB_vs_NAIVE_BM :
##   t = 0.5509 
##   df = 515.72 
##   p-value = 5.819313e-01 
##   FDR p-value = 7.274141e-01 
##   95% CI: [ -0.0644 ,  0.1146 ]
##   Log2 FC = 0.0251 
## 
##  NAIVE_WB_vs_SHAM_WB :
##   t = -0.1805 
##   df = 991.06 
##   p-value = 8.568284e-01 
##   FDR p-value = 9.370495e-01 
##   95% CI: [ -0.0713 ,  0.0593 ]
##   Log2 FC = -0.006 
## 
##  NAIVE_WB_vs_SHAM_BM :
##   t = -1.0785 
##   df = 272.33 
##   p-value = 2.817669e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.1641 ,  0.0479 ]
##   Log2 FC = -0.0581 
## 
##  NAIVE_WB_vs_UT_WB :
##   t = -0.8283 
##   df = 640.53 
##   p-value = 4.078026e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.0806 ,  0.0328 ]
##   Log2 FC = -0.0239 
## 
##  NAIVE_WB_vs_UT_BM :
##   t = -0.1579 
##   df = 521.38 
##   p-value = 8.745795e-01 
##   FDR p-value = 9.370495e-01 
##   95% CI: [ -0.0582 ,  0.0495 ]
##   Log2 FC = -0.0043 
## 
##  NAIVE_BM_vs_SHAM_WB :
##   t = -0.7378 
##   df = 425.36 
##   p-value = 4.610393e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.1139 ,  0.0517 ]
##   Log2 FC = -0.0311 
## 
##  NAIVE_BM_vs_SHAM_BM :
##   t = -1.3922 
##   df = 343.6 
##   p-value = 1.64751e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.2007 ,  0.0343 ]
##   Log2 FC = -0.0832 
## 
##  NAIVE_BM_vs_UT_WB :
##   t = -1.2644 
##   df = 310.8 
##   p-value = 2.070486e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.1253 ,  0.0273 ]
##   Log2 FC = -0.049 
## 
##  NAIVE_BM_vs_UT_BM :
##   t = -0.7808 
##   df = 277.43 
##   p-value = 4.35567e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.1036 ,  0.0447 ]
##   Log2 FC = -0.0294 
## 
##  SHAM_WB_vs_SHAM_BM :
##   t = -1.0208 
##   df = 226.04 
##   p-value = 3.084501e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.1526 ,  0.0484 ]
##   Log2 FC = -0.0521 
## 
##  SHAM_WB_vs_UT_WB :
##   t = -0.7736 
##   df = 1760.75 
##   p-value = 4.392892e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.0633 ,  0.0275 ]
##   Log2 FC = -0.0179 
## 
##  SHAM_WB_vs_UT_BM :
##   t = 0.0789 
##   df = 1299.61 
##   p-value = 9.371243e-01 
##   FDR p-value = 9.371243e-01 
##   95% CI: [ -0.0401 ,  0.0434 ]
##   Log2 FC = 0.0017 
## 
##  SHAM_BM_vs_UT_WB :
##   t = 0.7078 
##   df = 181.94 
##   p-value = 4.800052e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.0611 ,  0.1294 ]
##   Log2 FC = 0.0342 
## 
##  SHAM_BM_vs_UT_BM :
##   t = 1.1344 
##   df = 169.15 
##   p-value = 2.582346e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.0398 ,  0.1473 ]
##   Log2 FC = 0.0538 
## 
##  UT_WB_vs_UT_BM :
##   t = 1.4628 
##   df = 5173.52 
##   p-value = 1.435819e-01 
##   FDR p-value = 6.545525e-01 
##   95% CI: [ -0.0067 ,  0.0459 ]
##   Log2 FC = 0.0196
```

``` r
# Save the workbook
output_path <- here("R Projects", "2107", "output", "de_analysis", "t.test", paste0(gene_of_interest, ".xlsx"))

# Create directory if it doesn't exist
dir.create(dirname(output_path), recursive = TRUE, showWarnings = FALSE)

# Save workbook
saveWorkbook(wb, file = output_path, overwrite = TRUE)
```

# Visualisation

## T Cells

### Volcano Plots


``` r
# Function to create volcano plots for S1PR genes
create_s1pr_volcano_plots <- function(seu_T) {
  # Define S1PR genes
  s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
  
  # Define all 15 pairwise comparisons
  pairwise_comparisons <- list(
    "NAIVE_WB_vs_NAIVE_BM"     = c("Naive_WB", "Naive_BM"),
    "NAIVE_WB_vs_SHAM_WB"      = c("Naive_WB", "Sham_WB"),
    "NAIVE_WB_vs_SHAM_BM"      = c("Naive_WB", "Sham_BM"),
    "NAIVE_WB_vs_UT_WB"        = c("Naive_WB", "UT_WB"),
    "NAIVE_WB_vs_UT_BM"        = c("Naive_WB", "UT_BM"),
    "NAIVE_BM_vs_SHAM_WB"      = c("Naive_BM", "Sham_WB"),
    "NAIVE_BM_vs_SHAM_BM"      = c("Naive_BM", "Sham_BM"),
    "NAIVE_BM_vs_UT_WB"        = c("Naive_BM", "UT_WB"),
    "NAIVE_BM_vs_UT_BM"        = c("Naive_BM", "UT_BM"),
    "SHAM_WB_vs_SHAM_BM"       = c("Sham_WB", "Sham_BM"),
    "SHAM_WB_vs_UT_WB"         = c("Sham_WB", "UT_WB"),
    "SHAM_WB_vs_UT_BM"         = c("Sham_WB", "UT_BM"),
    "SHAM_BM_vs_UT_WB"         = c("Sham_BM", "UT_WB"),
    "SHAM_BM_vs_UT_BM"         = c("Sham_BM", "UT_BM"),
    "UT_WB_vs_UT_BM"           = c("UT_WB", "UT_BM")
  )
  
  volcano_plots <- list()
  
  for (comparison_name in names(pairwise_comparisons)) {
    groups <- pairwise_comparisons[[comparison_name]]
    
    cat("Creating volcano plot for:", comparison_name, "\n")
    
    # Run differential expression analysis
    de_results <- FindMarkers(
      seu_T,
      group.by = "sample_name",
      ident.1 = groups[1],
      ident.2 = groups[2],
      test.use = "wilcox",
      logfc.threshold = 0,
      min.pct = 0
    )
    
    # Add gene names
    de_results$gene <- rownames(de_results)
    
    # Fix zero p-values to avoid volcano plot warnings
    # Replace 0 p-values with a very small number
    min_nonzero_p <- min(de_results$p_val[de_results$p_val > 0], na.rm = TRUE)
    min_nonzero_p_adj <- min(de_results$p_val_adj[de_results$p_val_adj > 0], na.rm = TRUE)
    
    de_results$p_val[de_results$p_val == 0] <- min_nonzero_p * 0.1
    de_results$p_val_adj[de_results$p_val_adj == 0] <- min_nonzero_p_adj * 0.1
    
    # Create volcano plot
    volcano_plot <- EnhancedVolcano(
      de_results,
      lab = de_results$gene,
      x = 'avg_log2FC',
      y = 'p_val_adj',
      title = paste('Volcano Plot:', comparison_name),
      subtitle = paste('T Cells:', groups[1], 'vs', groups[2]),
      pCutoff = 0.05,
      FCcutoff = 0.25,
      pointSize = 1.5,
      labSize = 3,
      # Highlight S1PR genes
      selectLab = s1pr_genes,
      boxedLabels = TRUE,
      colAlpha = 0.7,
      legendPosition = 'right',
      legendLabSize = 10,
      legendIconSize = 3.0,
      drawConnectors = TRUE,
      widthConnectors = 0.3,
      colConnectors = 'black'
    )
    
    volcano_plots[[comparison_name]] <- volcano_plot
    print(volcano_plot)
    
    # Save individual volcano plot as SVG
    ggsave(
      filename = here("R Projects", "2107claude", "output", "visualisation", "volcano", 
                      paste0("volcano_", comparison_name, ".svg")),
      plot = volcano_plot,
      device = "svg",
      width = 10, height = 8
    )
  }
  
  cat("\nTotal volcano plots created:", length(volcano_plots), "\n")
  
  return(volcano_plots)
}

# Create volcano plots
volcano_plots <- create_s1pr_volcano_plots(seu_T)
```

```
## Creating volcano plot for: NAIVE_WB_vs_NAIVE_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-1.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_WB_vs_SHAM_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-2.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_WB_vs_SHAM_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-3.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_WB_vs_UT_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-4.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_WB_vs_UT_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-5.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_BM_vs_SHAM_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-6.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_BM_vs_SHAM_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-7.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_BM_vs_UT_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-8.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: NAIVE_BM_vs_UT_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-9.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: SHAM_WB_vs_SHAM_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-10.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: SHAM_WB_vs_UT_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-11.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: SHAM_WB_vs_UT_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-12.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: SHAM_BM_vs_UT_WB
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-13.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: SHAM_BM_vs_UT_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-14.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating volcano plot for: UT_WB_vs_UT_BM
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-27-15.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## 
## Total volcano plots created: 15
```

### Box Plots


``` r
# Function to create individual box plots for each S1PR gene
create_s1pr_box_plots <- function(seu_T) {
  # Define S1PR genes
  s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
  
  # Check which genes are available
  available_genes <- s1pr_genes[s1pr_genes %in% rownames(seu_T)]
  
  if (length(available_genes) == 0) {
    cat("No S1PR genes found in dataset\n")
    return(NULL)
  }
  
  cat("Creating individual box plots for genes:", paste(available_genes, collapse = ", "), "\n")
  
  # Extract expression data
  expr_data <- GetAssayData(seu_T, slot = "data", layer = "data")[available_genes, , drop = FALSE]
  if (inherits(expr_data, "dgCMatrix")) {
    expr_data <- as.matrix(expr_data)
  }
  
  # Create long format data
  plot_data <- expr_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(seu_T@meta.data %>% tibble::rownames_to_column("cell"), by = "cell")
  
  # Define colors for your actual data values
  treatment_colors <- c(
    "NAIVE" = "#2166AC",
    "SHAM" = "#762A83", 
    "UT" = "#C51B7D"
  )
  
  tissue_colors <- c(
    "Whole_blood" = "#D73027",
    "Bone_marrow" = "#4575B4"
  )
  
  all_plots <- list()
  
  # Create individual plots for each gene
  for (gene in available_genes) {
    gene_data <- plot_data %>% filter(gene == !!gene)
    
    cat("Creating box plots for", gene, "\n")
    
    # Box plot by treatment for this gene
    p_treatment <- ggplot(gene_data, aes(x = treatment, y = expression, fill = treatment)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.2, size = 0.3) +
      scale_fill_manual(values = treatment_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "none",
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Treatment",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression by Treatment"),
        subtitle = "T Cells"
      )
    
    # Box plot by tissue for this gene
    p_tissue <- ggplot(gene_data, aes(x = tissue, y = expression, fill = tissue)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      geom_jitter(width = 0.2, alpha = 0.2, size = 0.3) +
      scale_fill_manual(values = tissue_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "none",
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Tissue",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression by Tissue"),
        subtitle = "T Cells"
      )
    
    # Box plot by sample for this gene
    p_sample <- ggplot(gene_data, aes(x = sample_name, y = expression, fill = treatment)) +
      geom_boxplot(alpha = 0.7, outlier.size = 0.5) +
      scale_fill_manual(values = treatment_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Sample",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression by Sample"),
        subtitle = "T Cells",
        fill = "Treatment"
      )
    
    # Display plots
    print(p_treatment)
    print(p_tissue)
    print(p_sample)
    
    # Save individual plots as SVG
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "box plots", 
                paste0(gene, "_boxplot_treatment.svg")),
           p_treatment, device = "svg", width = 8, height = 6)
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "box plots", 
                paste0(gene, "_boxplot_tissue.svg")),
           p_tissue, device = "svg", width = 8, height = 6)
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "box plots", 
                paste0(gene, "_boxplot_sample.svg")),
           p_sample, device = "svg", width = 10, height = 6)
    
    # Store plots
    all_plots[[paste0(gene, "_treatment")]] <- p_treatment
    all_plots[[paste0(gene, "_tissue")]] <- p_tissue
    all_plots[[paste0(gene, "_sample")]] <- p_sample
  }
  
  return(all_plots)
}

# Create box plots
box_plots <- create_s1pr_box_plots(seu_T)
```

```
## Creating individual box plots for genes: S1pr1, S1pr2, S1pr3, S1pr4, S1pr5 
## Creating box plots for S1pr1
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-1.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-2.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-3.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating box plots for S1pr2
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-4.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-5.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-6.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating box plots for S1pr3
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-7.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-8.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-9.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating box plots for S1pr4
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-10.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-11.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-12.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating box plots for S1pr5
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-13.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-14.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-28-15.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

### Violin Plots


``` r
# Function to create individual violin plots for each S1PR gene
create_s1pr_violin_plots <- function(seu_T) {
  # Define S1PR genes
  s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")
  
  # Check which genes are available
  available_genes <- s1pr_genes[s1pr_genes %in% rownames(seu_T)]
  
  if (length(available_genes) == 0) {
    cat("No S1PR genes found in dataset\n")
    return(NULL)
  }
  
  cat("Creating individual violin plots for genes:", paste(available_genes, collapse = ", "), "\n")
  
  # Extract expression data
  expr_data <- GetAssayData(seu_T, slot = "data", layer = "data")[available_genes, , drop = FALSE]
  if (inherits(expr_data, "dgCMatrix")) {
    expr_data <- as.matrix(expr_data)
  }
  
  # Create long format data
  plot_data <- expr_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    pivot_longer(-gene, names_to = "cell", values_to = "expression") %>%
    left_join(seu_T@meta.data %>% tibble::rownames_to_column("cell"), by = "cell")
  
  # Define colors for your actual data values
  treatment_colors <- c(
    "NAIVE" = "#2166AC",
    "SHAM" = "#762A83", 
    "UT" = "#C51B7D"
  )
  
  tissue_colors <- c(
    "Whole_blood" = "#D73027",
    "Bone_marrow" = "#4575B4"
  )
  
  all_plots <- list()
  
  # Create individual plots for each gene
  for (gene in available_genes) {
    gene_data <- plot_data %>% filter(gene == !!gene)
    
    cat("Creating violin plots for", gene, "\n")
    
    # Violin plot by treatment for this gene
    p_violin_treatment <- ggplot(gene_data, aes(x = treatment, y = expression, fill = treatment)) +
      geom_violin(scale = "width", alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8, fill = "white") +
      stat_summary(fun = median, geom = "point", size = 3, color = "black") +
      scale_fill_manual(values = treatment_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.position = "none",
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Treatment",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression Distribution by Treatment"),
        subtitle = "T Cells - Violin Plot"
      )
    
    # Violin plot by sample for this gene
    p_violin_sample <- ggplot(gene_data, aes(x = sample_name, y = expression, fill = treatment)) +
      geom_violin(scale = "width", alpha = 0.7, trim = FALSE) +
      geom_boxplot(width = 0.1, outlier.shape = NA, alpha = 0.8, fill = "white") +
      scale_fill_manual(values = treatment_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 11),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Sample",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression Distribution by Sample"),
        subtitle = "T Cells - Violin Plot",
        fill = "Treatment"
      )
    
    # Split violin plot comparing tissues within each treatment for this gene
    p_violin_split <- ggplot(gene_data, aes(x = treatment, y = expression, fill = tissue)) +
      geom_violin(position = position_dodge(width = 0.8), alpha = 0.7) +
      geom_boxplot(position = position_dodge(width = 0.8), width = 0.1, 
                   outlier.shape = NA, alpha = 0.8, fill = "white") +
      scale_fill_manual(values = tissue_colors) +
      theme_minimal() +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
        axis.text.y = element_text(size = 11),
        axis.title = element_text(size = 14, face = "bold"),
        plot.title = element_text(size = 16, face = "bold"),
        plot.subtitle = element_text(size = 12),
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = "Treatment",
        y = "Expression (Log-normalized)",
        title = paste(gene, "Expression: Treatment vs Tissue"),
        subtitle = "T Cells - Split Violin Plot",
        fill = "Tissue"
      )
    
    # Display plots
    print(p_violin_treatment)
    print(p_violin_sample)
    print(p_violin_split)
    
    # Save individual plots as SVG
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "violin plots", 
                paste0(gene, "_violin_treatment.svg")),
           p_violin_treatment, device = "svg", width = 8, height = 6)
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "violin plots", 
                paste0(gene, "_violin_sample.svg")),
           p_violin_sample, device = "svg", width = 10, height = 6)
    ggsave(here("R Projects", "2107claude", "output", "visualisation", "violin plots", 
                paste0(gene, "_violin_split.svg")),
           p_violin_split, device = "svg", width = 8, height = 6)
    
    # Store plots
    all_plots[[paste0(gene, "_treatment")]] <- p_violin_treatment
    all_plots[[paste0(gene, "_sample")]] <- p_violin_sample
    all_plots[[paste0(gene, "_split")]] <- p_violin_split
  }
  
  return(all_plots)
}

# Create violin plots
violin_plots <- create_s1pr_violin_plots(seu_T)
```

```
## Creating individual violin plots for genes: S1pr1, S1pr2, S1pr3, S1pr4, S1pr5 
## Creating violin plots for S1pr1
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-1.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-2.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-3.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating violin plots for S1pr2
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-4.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-5.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-6.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating violin plots for S1pr3
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-7.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-8.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-9.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating violin plots for S1pr4
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-10.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-11.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-12.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

```
## Creating violin plots for S1pr5
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-13.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-14.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" /><img src="2107notebook_claude_files/figure-html/unnamed-chunk-29-15.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

### Nebulosa


``` r
# S1pr1
Nebulosa::plot_density(seu_T, "S1pr1") + 
    facet_wrap(.~seu_T$sample_name, ncol = 3) + theme_void()
```

```
## Warning: The `slot` argument of `FetchData()` is deprecated as of SeuratObject 5.0.0.
## ℹ Please use the `layer` argument instead.
## ℹ The deprecated feature was likely used in the Nebulosa package.
##   Please report the issue at
##   <https://github.com/powellgenomicslab/Nebulosa/issues>.
## This warning is displayed once every 8 hours.
## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
## generated.
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-30-1.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# S1pr2
Nebulosa::plot_density(seu_T, "S1pr2") + 
    facet_wrap(.~seu_T$sample_name, ncol = 3) + theme_void()
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-30-2.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# S1pr3
Nebulosa::plot_density(seu_T, "S1pr3") + 
    facet_wrap(.~seu_T$sample_name, ncol = 3) + theme_void()
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-30-3.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# S1pr4
Nebulosa::plot_density(seu_T, "S1pr4") + 
    facet_wrap(.~seu_T$sample_name, ncol = 3) + theme_void()
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-30-4.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# S1pr5
Nebulosa::plot_density(seu_T, "S1pr5") + 
    facet_wrap(.~seu_T$sample_name, ncol = 3) + theme_void()
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-30-5.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />



### Heatmaps


``` r
# Define S1PR gene family
s1pr_genes <- c("S1pr1", "S1pr2", "S1pr3", "S1pr4", "S1pr5")

# Function to check gene availability and get expression data
extract_s1pr_data <- function(seurat_obj, genes = s1pr_genes) {
  # Check which S1PR genes are present in the dataset
  available_genes <- genes[genes %in% rownames(seurat_obj)]
  missing_genes <- genes[!genes %in% rownames(seurat_obj)]
  
  if (length(missing_genes) > 0) {
    cat("Warning: Missing genes:", paste(missing_genes, collapse = ", "), "\n")
  }
  
  if (length(available_genes) == 0) {
    stop("No S1PR genes found in the dataset!")
  }
  
  cat("Found S1PR genes:", paste(available_genes, collapse = ", "), "\n")
  
  # Extract expression data (using layer parameter instead of deprecated slot)
  expr_data <- GetAssayData(seurat_obj, layer = "data")[available_genes, , drop = FALSE]
  
  # Convert sparse matrix to dense if needed
  if (inherits(expr_data, "dgCMatrix")) {
    expr_data <- as.matrix(expr_data)
  }
  
  # Get metadata
  metadata <- seurat_obj@meta.data
  
  return(list(
    expression = expr_data,
    metadata = metadata,
    available_genes = available_genes
  ))
}

# Function to create sample annotation colors
create_annotation_colors <- function(metadata) {
  # Define color palettes for your specific conditions
  treatment_colors <- c(
    "NAIVE" = "#2166AC",
    "SHAM" = "#762A83", 
    "UT" = "#C51B7D"
  )
  
  tissue_colors <- c(
    "Whole_blood" = "#D73027",
    "Bone_marrow" = "#4575B4"
  )
  
  return(list(
    treatment = treatment_colors,
    tissue = tissue_colors
  ))
}

# Function to calculate summary statistics by treatment and tissue
calculate_treatment_tissue_stats <- function(expr_data, metadata) {
  # Convert to long format for easier manipulation
  expr_long <- expr_data %>%
    as.data.frame() %>%
    tibble::rownames_to_column("gene") %>%
    tidyr::pivot_longer(-gene, names_to = "cell", values_to = "expression")
  
  # Add metadata
  expr_long <- expr_long %>%
    dplyr::left_join(metadata %>% tibble::rownames_to_column("cell"), by = "cell")
  
  # Calculate summary statistics by treatment, tissue, and gene
  summary_stats <- expr_long %>%
    dplyr::group_by(gene, treatment, tissue) %>%
    dplyr::summarise(
      mean_expr = mean(expression),
      median_expr = median(expression),
      sd_expr = sd(expression),
      n_cells = n(),
      pct_positive = sum(expression > 0) / n() * 100,
      .groups = "drop"
    ) %>%
    # Create combined treatment_tissue column for better visualization
    dplyr::mutate(treatment_tissue = paste(treatment, tissue, sep = "_"))
  
  return(summary_stats)
}

# Function to create S1PR expression heatmap by treatment and tissue
create_s1pr_treatment_tissue_heatmap <- function(seurat_obj, metric = "mean") {
  # Extract data
  s1pr_data <- extract_s1pr_data(seurat_obj)
  expr_data <- s1pr_data$expression
  metadata <- s1pr_data$metadata
  
  # Calculate summary statistics
  summary_data <- calculate_treatment_tissue_stats(expr_data, metadata)
  
  cat("Summary data structure:\n")
  print(head(summary_data))
  cat("Unique treatment_tissue combinations:", paste(unique(summary_data$treatment_tissue), collapse = ", "), "\n")
  
  # Create summary matrix using explicit dplyr functions
  if (metric == "mean") {
    summary_matrix <- summary_data %>%
      dplyr::select(gene, treatment_tissue, mean_expr) %>%
      tidyr::pivot_wider(names_from = treatment_tissue, values_from = mean_expr) %>%
      tibble::column_to_rownames("gene") %>%
      as.matrix()
    title_suffix <- "Mean Expression"
    legend_title <- "Mean Expression"
  } else if (metric == "pct_positive") {
    summary_matrix <- summary_data %>%
      dplyr::select(gene, treatment_tissue, pct_positive) %>%
      tidyr::pivot_wider(names_from = treatment_tissue, values_from = pct_positive) %>%
      tibble::column_to_rownames("gene") %>%
      as.matrix()
    title_suffix <- "% Positive Cells"
    legend_title <- "% Positive"
  }
  
  cat("Summary matrix dimensions:", dim(summary_matrix), "\n")
  cat("Column names:", colnames(summary_matrix), "\n")
  
  # Create annotation for treatment-tissue combinations
  col_info <- summary_data %>%
    dplyr::select(treatment_tissue, treatment, tissue) %>%
    dplyr::distinct() %>%
    dplyr::arrange(treatment_tissue) %>%
    tibble::column_to_rownames("treatment_tissue")
  
  # Ensure column order matches
  summary_matrix <- summary_matrix[, rownames(col_info), drop = FALSE]
  
  anno_colors <- create_annotation_colors(metadata)
  
  col_anno <- HeatmapAnnotation(
    Treatment = col_info$treatment,
    Tissue = col_info$tissue,
    col = list(
      Treatment = anno_colors$treatment,
      Tissue = anno_colors$tissue
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  # Create color function
  if (metric == "mean") {
    max_val <- max(summary_matrix, na.rm = TRUE)
    col_fun <- colorRamp2(c(0, max_val/2, max_val), c("white", "pink", "red"))
  } else {
    col_fun <- colorRamp2(c(0, 50, 100), c("white", "orange", "red"))
  }
  
  # Create heatmap
  ht <- Heatmap(
    summary_matrix,
    name = legend_title,
    col = col_fun,
    
    # Appearance
    column_title = paste("S1PR", title_suffix, "by Treatment and Tissue"),
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    row_title = "S1PR Genes", 
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = gpar(fontsize = 11, fontface = "italic"),
    column_names_gp = gpar(fontsize = 10),
    
    # Clustering
    cluster_rows = nrow(summary_matrix) > 1,
    cluster_columns = TRUE,
    
    # Annotations
    top_annotation = col_anno,
    
    # Cell text
    cell_fun = function(j, i, x, y, width, height, fill) {
      grid.text(sprintf("%.2f", summary_matrix[i, j]), 
                x, y, gp = gpar(fontsize = 9, col = "black"))
    },
    
    border = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  return(list(
    heatmap = ht,
    summary_matrix = summary_matrix,
    summary_data = summary_data
  ))
}

# Function for detailed single-cell heatmap with subsampling by treatment and tissue
create_s1pr_detailed_treatment_tissue_heatmap <- function(seurat_obj, 
                                                          subsample_cells = 100,
                                                          scale_data = TRUE) {
  
  # Extract data
  s1pr_data <- extract_s1pr_data(seurat_obj)
  expr_data <- s1pr_data$expression
  metadata <- s1pr_data$metadata
  available_genes <- s1pr_data$available_genes
  
  # Subsample cells for visualization - stratified by treatment and tissue
  cat("Subsampling to", subsample_cells, "cells per treatment-tissue combination...\n")
  
  # Add cell names to metadata for easier manipulation
  metadata$cell_name <- rownames(metadata)
  
  # Stratified sampling by treatment and tissue - fixed the pull issue
  sampled_cells <- metadata %>%
    dplyr::mutate(treatment_tissue = paste(treatment, tissue, sep = "_")) %>%
    dplyr::group_by(treatment_tissue) %>%
    dplyr::sample_n(min(n(), subsample_cells)) %>%
    dplyr::ungroup() %>%
    dplyr::pull(cell_name)  # Now pulling a specific column
  
  expr_data <- expr_data[, sampled_cells]
  metadata <- metadata[sampled_cells, ]
  
  # Scale data if requested
  if (scale_data) {
    expr_data <- t(scale(t(expr_data)))
    # Handle any NaN values (genes with zero variance)
    expr_data[is.nan(expr_data)] <- 0
  }
  
  # Create annotation colors
  anno_colors <- create_annotation_colors(metadata)
  
  # Create column annotations
  col_anno <- HeatmapAnnotation(
    Treatment = metadata$treatment,
    Tissue = metadata$tissue,
    col = list(
      Treatment = anno_colors$treatment,
      Tissue = anno_colors$tissue
    ),
    annotation_name_gp = gpar(fontsize = 10, fontface = "bold")
  )
  
  # Order cells by treatment and tissue for better visualization
  metadata <- metadata %>%
    dplyr::mutate(treatment_tissue = paste(treatment, tissue, sep = "_"))
  cell_order <- order(metadata$treatment_tissue)
  expr_data_ordered <- expr_data[, cell_order]
  
  # Create color function for expression
  if (scale_data) {
    col_fun <- colorRamp2(c(-2, 0, 2), c("blue", "white", "red"))
    heatmap_title <- "S1PR Expression (Z-score) by Treatment and Tissue"
    legend_title <- "Z-score"
  } else {
    max_expr <- max(expr_data, na.rm = TRUE)
    col_fun <- colorRamp2(c(0, max_expr/2, max_expr), c("white", "pink", "red"))
    heatmap_title <- "S1PR Expression (Log-normalized) by Treatment and Tissue"
    legend_title <- "Expression"
  }
  
  # Create main heatmap
  ht <- Heatmap(
    expr_data_ordered,
    name = legend_title,
    col = col_fun,
    
    # Row parameters
    row_title = "S1PR Genes",
    row_title_gp = gpar(fontsize = 12, fontface = "bold"),
    row_names_gp = gpar(fontsize = 11, fontface = "italic"),
    cluster_rows = length(available_genes) > 1,
    
    # Column parameters  
    column_title = heatmap_title,
    column_title_gp = gpar(fontsize = 14, fontface = "bold"),
    cluster_columns = FALSE,
    show_column_names = FALSE,
    
    # Annotations
    top_annotation = col_anno,
    
    # Appearance
    border = TRUE,
    heatmap_legend_param = list(
      title_gp = gpar(fontsize = 11, fontface = "bold"),
      labels_gp = gpar(fontsize = 9)
    )
  )
  
  return(ht)
}

# Create directories if they don't exist
dir.create(here("R Projects", "2107claude", "output", "visualisation", "complex heatmaps"), 
           recursive = TRUE, showWarnings = FALSE)

# Create mean expression heatmap by treatment and tissue
cat("Creating S1PR mean expression heatmap by treatment and tissue...\n")
```

```
## Creating S1PR mean expression heatmap by treatment and tissue...
```

``` r
mean_heatmap_result <- create_s1pr_treatment_tissue_heatmap(seu_T, metric = "mean")
```

```
## Found S1PR genes: S1pr1, S1pr2, S1pr3, S1pr4, S1pr5 
## Summary data structure:
## # A tibble: 6 × 9
##   gene  treatment tissue      mean_expr median_expr sd_expr n_cells pct_positive
##   <chr> <chr>     <chr>           <dbl>       <dbl>   <dbl>   <int>        <dbl>
## 1 S1pr1 NAIVE     Bone_marrow     0.484        0      0.773     260         30.4
## 2 S1pr1 NAIVE     Whole_blood     0.718        0      0.875     458         43.7
## 3 S1pr1 SHAM      Bone_marrow     0.392        0      0.668     163         28.8
## 4 S1pr1 SHAM      Whole_blood     0.844        0      0.965    1042         45.7
## 5 S1pr1 UT        Bone_marrow     1.08         1.30   0.920    7130         62.3
## 6 S1pr1 UT        Whole_blood     0.704        0      0.923    2930         40.6
## # ℹ 1 more variable: treatment_tissue <chr>
## Unique treatment_tissue combinations: NAIVE_Bone_marrow, NAIVE_Whole_blood, SHAM_Bone_marrow, SHAM_Whole_blood, UT_Bone_marrow, UT_Whole_blood 
## Summary matrix dimensions: 5 6 
## Column names: NAIVE_Bone_marrow NAIVE_Whole_blood SHAM_Bone_marrow SHAM_Whole_blood UT_Bone_marrow UT_Whole_blood
```

``` r
draw(mean_heatmap_result$heatmap)
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-32-1.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# Save mean expression heatmap as SVG using a different approach
cat("Saving mean expression heatmap...\n")
```

```
## Saving mean expression heatmap...
```

``` r
heatmap_file <- here("R Projects", "2107claude", "output", "visualisation", "complex heatmaps", "s1pr_mean_expression_treatment_tissue_heatmap.svg")
svg(heatmap_file, width = 10, height = 6)
draw(mean_heatmap_result$heatmap)
dev.off()
```

```
## svg 
##   2
```

``` r
cat("Mean expression heatmap saved to:", heatmap_file, "\n")
```

```
## Mean expression heatmap saved to: C:/Users/fionn/Newcastle Grammar School/NGS SciX-2025 UoN DMG Cancer Signalling - Bioinformatics and Proteomics/Fionn Git/HSC-SciX-UoN_LS3.33/R Projects/2107claude/output/visualisation/complex heatmaps/s1pr_mean_expression_treatment_tissue_heatmap.svg
```

``` r
# Create percentage positive heatmap by treatment and tissue
cat("Creating S1PR percentage positive heatmap by treatment and tissue...\n")
```

```
## Creating S1PR percentage positive heatmap by treatment and tissue...
```

``` r
pct_heatmap_result <- create_s1pr_treatment_tissue_heatmap(seu_T, metric = "pct_positive")
```

```
## Found S1PR genes: S1pr1, S1pr2, S1pr3, S1pr4, S1pr5 
## Summary data structure:
## # A tibble: 6 × 9
##   gene  treatment tissue      mean_expr median_expr sd_expr n_cells pct_positive
##   <chr> <chr>     <chr>           <dbl>       <dbl>   <dbl>   <int>        <dbl>
## 1 S1pr1 NAIVE     Bone_marrow     0.484        0      0.773     260         30.4
## 2 S1pr1 NAIVE     Whole_blood     0.718        0      0.875     458         43.7
## 3 S1pr1 SHAM      Bone_marrow     0.392        0      0.668     163         28.8
## 4 S1pr1 SHAM      Whole_blood     0.844        0      0.965    1042         45.7
## 5 S1pr1 UT        Bone_marrow     1.08         1.30   0.920    7130         62.3
## 6 S1pr1 UT        Whole_blood     0.704        0      0.923    2930         40.6
## # ℹ 1 more variable: treatment_tissue <chr>
## Unique treatment_tissue combinations: NAIVE_Bone_marrow, NAIVE_Whole_blood, SHAM_Bone_marrow, SHAM_Whole_blood, UT_Bone_marrow, UT_Whole_blood 
## Summary matrix dimensions: 5 6 
## Column names: NAIVE_Bone_marrow NAIVE_Whole_blood SHAM_Bone_marrow SHAM_Whole_blood UT_Bone_marrow UT_Whole_blood
```

``` r
draw(pct_heatmap_result$heatmap)
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-32-2.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# Save percentage positive heatmap as SVG
cat("Saving percentage positive heatmap...\n")
```

```
## Saving percentage positive heatmap...
```

``` r
pct_file <- here("R Projects", "2107claude", "output", "visualisation", "complex heatmaps", "s1pr_percent_positive_treatment_tissue_heatmap.svg")
svg(pct_file, width = 10, height = 6)
draw(pct_heatmap_result$heatmap)
dev.off()
```

```
## svg 
##   2
```

``` r
cat("Percentage positive heatmap saved to:", pct_file, "\n")
```

```
## Percentage positive heatmap saved to: C:/Users/fionn/Newcastle Grammar School/NGS SciX-2025 UoN DMG Cancer Signalling - Bioinformatics and Proteomics/Fionn Git/HSC-SciX-UoN_LS3.33/R Projects/2107claude/output/visualisation/complex heatmaps/s1pr_percent_positive_treatment_tissue_heatmap.svg
```

``` r
# Create detailed single-cell heatmap by treatment and tissue
cat("Creating detailed single-cell heatmap by treatment and tissue...\n")
```

```
## Creating detailed single-cell heatmap by treatment and tissue...
```

``` r
detailed_heatmap <- create_s1pr_detailed_treatment_tissue_heatmap(seu_T, subsample_cells = 100, scale_data = TRUE)
```

```
## Found S1PR genes: S1pr1, S1pr2, S1pr3, S1pr4, S1pr5 
## Subsampling to 100 cells per treatment-tissue combination...
```

``` r
draw(detailed_heatmap)
```

<img src="2107notebook_claude_files/figure-html/unnamed-chunk-32-3.svg" width="90%" keepaspectratio=true style="display: block; margin: auto;" />

``` r
# Save detailed heatmap as SVG
cat("Saving detailed single-cell heatmap...\n")
```

```
## Saving detailed single-cell heatmap...
```

``` r
detailed_file <- here("R Projects", "2107claude", "output", "visualisation", "complex heatmaps", "s1pr_detailed_treatment_tissue_heatmap.svg")
svg(detailed_file, width = 14, height = 6)
draw(detailed_heatmap)
dev.off()
```

```
## svg 
##   2
```

``` r
cat("Detailed heatmap saved to:", detailed_file, "\n")
```

```
## Detailed heatmap saved to: C:/Users/fionn/Newcastle Grammar School/NGS SciX-2025 UoN DMG Cancer Signalling - Bioinformatics and Proteomics/Fionn Git/HSC-SciX-UoN_LS3.33/R Projects/2107claude/output/visualisation/complex heatmaps/s1pr_detailed_treatment_tissue_heatmap.svg
```

``` r
cat("All S1PR treatment-tissue heatmaps created and saved as SVG!\n")
```

```
## All S1PR treatment-tissue heatmaps created and saved as SVG!
```
