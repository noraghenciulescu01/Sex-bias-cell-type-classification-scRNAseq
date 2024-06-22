# load packages
library(edgeR)
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(Matrix)
library(scMerge)
print("successfully loaded packages!")
flush.console()

# load data
pb <- readRDS("aggregated_covid_sce.rds")
print("loaded data!")
flush.console()



lst <- list()
# iterate over all cell types, individually perform DS analysis, save results to lst
for (assay_name in assayNames(pb)){
  pb_subset <- SingleCellExperiment(assays = setNames(list(assay(pb, assay_name)), assay_name),
                                    colData = colData(pb),
                                    rowData = rowData(pb),
                                    metadata = metadata(pb))
  
  # we need to edit the internal column metadata too
  n_cells <- lapply(int_colData(pb)$n_cells, function(x) {
    x[assay_name, drop = FALSE]
  })
  int_colData(pb_subset)$n_cells <- n_cells
  
  # preprocess
  # remove undetected genes
  pb_subset_fil <- pb_subset[rowSums(assay(pb_subset, assay_name)) > 0, ]
  # remove donors with no samples of this cell type
  pb_subset_fil <- pb_subset_fil[, !(colSums(assay(pb_subset_fil, assay_name)) == 0)]  
  
  # DS analysis
  df <- tryCatch({
    res <- pbDS(pb_subset_fil, verbose = FALSE)
    # save res to file
    saveRDS(res, file = paste0("individualDS_covid/result_", assay_name, ".rds"))
    
    
    # Handle result to get the number and percentage of DS genes
    tbl <- res$table[[1]]
    tbl_fil <- lapply(tbl, function(u) {
      u <- dplyr::filter(u, p_adj.loc < 0.05)
      dplyr::arrange(u, p_adj.loc)
    })
    n_de <- vapply(tbl_fil, nrow, numeric(1))
    p_de <- format(n_de / nrow(pb_subset_fil) * 100, digits = 3)
    data.frame("#DS" = n_de, "%DS" = p_de, check.names = FALSE)
  }, error = function(e) {
    # If an error occurs, create a data frame with NaN values
    message("Error in pbDS: ", e$message, " for cell type: ", assay_name)
    data.frame("#DS" = NaN, "%DS" = NaN, check.names = FALSE)
  })
  rownames(df) <- assay_name
  
  # save df to lst
  lst[[assay_name]] <- df
  
}
print("performed DS for all, now stitching together...")
flush.console()


# stitch all dfs together
final_results <- do.call(rbind, lst)
final_results <- final_results[order(as.numeric(final_results$"%DS"), decreasing=TRUE),]

print(final_results)
flush.console()

# save results
saveRDS(final_results, "final_DS_results_covid.rds")
