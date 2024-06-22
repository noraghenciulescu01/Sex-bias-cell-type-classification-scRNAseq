import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from scipy.sparse import vstack

import os
import pickle

from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_score
import seaborn as sns

print('successfully imported packages!', flush=True)

path = '/winhome/nghenciulescu/Desktop/ild.h5ad'
adata = anndata.read_h5ad(path)

print('successfully read file!', flush=True)

# remove problematic cell types:
to_remove = ['EC capillary', 'Innate lymphoid cell NK', 'Lymphoid', 'SM activated stress response', 'Smooth muscle FAM83D+', 'Lymphatic EC differentiating', 'Rare']
adata = adata[~adata.obs["ann_level_3"].isin(to_remove)]


embedding = adata.obsm['X_scanvi_emb']
sex_labels = adata.obs['sex']
cell_type_labels = adata.obs['ann_level_3'].astype(str)
cell_type_labels = cell_type_labels.where((cell_type_labels != 'nan') & (cell_type_labels != 'Unknown'), adata.obs['ann_level_2'].astype(str))
classes = sorted(set(cell_type_labels))

from helper_functions import final_train_clf_and_predict, final_evaluate_clf, plot_confusion, fixed_select_indices_by_proportion

# ----------
# Run on data increasing female proportion:
    # (metrics dictionary is pickled)
    # (results are also printed and plots are saved)

prop = 0

print(f"PROPORTION OF FEMALE CELLS: {prop}", flush=True)
print('Training and testing...', flush=True)
male_pred, male_true_labels, female_pred, female_true_labels = final_train_clf_and_predict(embedding, cell_type_labels, sex_labels, prop)
print('Evaluating...', flush=True)
male_metrics = final_evaluate_clf(male_pred, male_true_labels, classes, prop, 'male')
female_metrics = final_evaluate_clf(female_pred, female_true_labels, classes, prop, 'female')

# File path to save the dictionary
male_file_path = f"{''.join(str(prop).split('.'))}_male_metrics.pickle"
female_file_path = f"{''.join(str(prop).split('.'))}_female_metrics.pickle"

# Save the dictionary to a file
with open(male_file_path, 'wb') as file:
    pickle.dump(male_metrics, file)

with open(female_file_path, 'wb') as file:
    pickle.dump(female_metrics, file)
