## Classification result analysis
The results of the classification experiment, performed as described in the last section, are added in the folder Classification Results. It contains a sub-folder for every seed, which in turn contains a sub-folder for each annotation level, which in turn contains a sub-folder for the KNN and RF classifiers. The results are given as pickle files containing a dictionary of metrics. There is one pickle file for every proportion of female cells in the training set (0, 0.1, ..., 1). <br>

Determining of cell type proportions in the training and test set (performed in notebook 1.) is needed for the classification pattern analysis. After running the notebook, a proportion_dictionary pickle file is saved, which contains the training and test set proportions for each cell type, at all annotation levels and proportions of female cells in the training set. <br>

The general classification results are analyzed in notebook 2; the accuracy and F1 scores plots are saved in Analysis_results/Performance_plots. Individual classification trends are analyzed in notebook 3; the resulting dataframe is saved in Analysis_results. We make use of custom dictionaries of cell names and color mappings to in notebooks 1 and 3 for ease of plotting. <br>