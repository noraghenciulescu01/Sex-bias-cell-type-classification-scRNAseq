# Note: the covid_finest can be saved from the pre-processing notebook;
# it contains the object pre-processed at the finest level (since we will be running DE analysis at the finest level too)

library(SingleCellExperiment)
library(dplyr)
library(ggplot2)
library(limma)
library(muscat)
library(purrr)
library(scater)
library(zellkonverter)

# Convert COVID-19 .h5ad object to .rds:
sce <- readH5AD("covid_finest.h5ad", reader = "R")
assays(sce)$counts <- assays(sce)$X
assays(sce)$X <- NULL 
saveRDS(sce, "covid_sce.rds")