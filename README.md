This repository is divided into 4 sub-sections of experiments, numbered in the order that they are being referred to in the paper; notebooks or scripts within each folder are also given in order.  Most experiments are performed  in Jupyter notebooks; a few use Python or R scripts. 

The full code for the classification experiment is not given - instead, a sample script is included for every experiment, along with instructions on how to generate subsequent scripts. The classification results are given in full.

This research uses data from the Human Lung Cell Atlas (Sikkema et al., 2023); the data can be downloaded from the CELLxGENE database: https://cellxgene.cziscience.com/collections/6f6d381a-7701-4781-935c-db10d30de293. The data is used in both the .h5ad format (for the analyses in Python) and in the Seurat format (for those in R). Specifically:
- analyses in folders 1., 2., 3., 5., 6., 7., 9. only use the .h5ad data
- analyses in folders 4. and 8. use the Seurat object (or both the Seurat and the .h5ad object)
