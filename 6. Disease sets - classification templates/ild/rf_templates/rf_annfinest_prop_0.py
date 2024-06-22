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

mask = (adata.obs['ann_finest_level'] == 'Unknown')
adata.obs['ann_finest_level'] = adata.obs['ann_finest_level'].astype('object')
adata.obs.loc[mask, 'ann_finest_level'] = adata.obs.loc[mask, 'ann_level_4']
adata.obs.loc[mask, 'ann_finest_level'] = adata.obs.loc[mask, 'ann_level_3']
adata.obs.loc[mask, 'ann_finest_level'] = adata.obs.loc[mask, 'ann_level_2']
adata.obs['ann_finest_level'] = adata.obs['ann_finest_level'].astype('category')

# remove problematic cell types:
to_remove = ['AT2 proliferating', 'NK cells', 'T cells proliferating', 'Neuroendocrine', 'Migratory DCs', 'Plasma cells', 'Alveolar Mph CCL3+', 'Hillock-like', 'Deuterosomal', 'Alveolar Mph MT-positive', 'Lymphoid', 'SM activated stress response', 'Smooth muscle FAM83D+', 'Goblet (bronchial)', 'DC1', 'Lymphatic EC differentiating', 'EC general capillary', 'AT0', 'Interstitial Mph perivascular','Multiciliated (nasal)']
adata = adata[~adata.obs["ann_finest_level"].isin(to_remove)]


counts = adata.X
sex_labels = adata.obs['sex']
cell_type_labels = adata.obs['ann_finest_level'].astype(str)
classes = sorted(set(cell_type_labels))


from helper_functions import final_train_clf_and_predict, final_evaluate_clf, plot_confusion, fixed_select_indices_by_proportion

# ----------
# Run on data increasing female proportion:
    # (metrics dictionary is pickled)
    # (results are also printed and plots are saved)

prop = 0

print(f"PROPORTION OF FEMALE CELLS: {prop}", flush=True)
print('Training and testing...', flush=True)
male_pred, male_true_labels, female_pred, female_true_labels = final_train_clf_and_predict(counts, cell_type_labels, sex_labels, prop, classifier = 'rf')
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
