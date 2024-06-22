# Load packages
library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
print("loaded packages!")
flush.console()

sce = readRDS("/winhome/nghenciulescu/Desktop/covid_sce.rds")
print("loaded object!")
print(sce)
flush.console()

# Pre-processing:
print("removing undetected genes:")
sce <- sce[rowSums(counts(sce) > 0) > 0, ] # remove undetected genes
print(dim(sce))
flush.console()

print("removing outliers:")
qc <- perCellQCMetrics(sce)
ol <- isOutlier(metric = qc$detected, nmads = 2, log = TRUE)
sce <- sce[, !ol] # remove cells with outliers (in terms of number of detected genes)
print(dim(sce))
flush.console()

print("removing lowly expressed genes:")
sce <- sce[rowSums(counts(sce) > 1) >= 10, ] # remove lowly expressed genes
print(dim(sce))
flush.console()


# prepare for muscat expected format
(sce <- prepSCE(sce, 
                kid = "ann_finest_level", # sub-population assignments
                gid = "sex",  # group IDs (male/female)
                sid = "donor_id",   # sample IDs (donor ids)
                drop = TRUE))  # drop all other colData columns


print(sce)
print("prepped object!")
flush.console()

# Aggregate to pseudo-bulk
pb <- aggregateData(sce,
                    assay = "counts", fun = "sum",
                    by = c("cluster_id", "sample_id"))

print(assayNames(pb)) # check: one sheet per subpopulation
print(t(head(assay(pb)))) # check: 1st subpopulation
print("aggregated!")
flush.console()

# save aggregated object
saveRDS(pb, "aggregated_covid_sce.rds")
print("saved object!")
flush.console()

