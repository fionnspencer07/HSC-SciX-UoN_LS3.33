source("initialise_packages.r")

source("Creating_Comparisons.r")

saveRDS(deg_results, file = "deg_results.rds")
save

# load results

degresults <- readRDS("deg_results.RDS")

"Naive_vs_Sham_WB" <- deg_results$Naive_WB_vs_Sham_WB

"Naive_vs_UT_WB" <- deg_results$Naive_WB_vs_UT_WB

"Naive_vs_ONC_WB" <- deg_results$Naive_WB_vs_ONC_WB

