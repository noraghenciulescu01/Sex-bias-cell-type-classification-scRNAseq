import numpy as np
import pandas as pd
import scanpy as sc
import anndata
import matplotlib.pyplot as plt
from scipy.sparse import vstack

import os

from sklearn.model_selection import train_test_split
from sklearn.neighbors import KNeighborsClassifier
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import accuracy_score, f1_score, confusion_matrix, precision_score
import seaborn as sns


# Train and evaluate functions:

# with equalized test sets!
def final_train_clf_and_predict(X, y, sex_labels, proportion_female, classifier='knn', k = 30):
    '''
    Parameters:
    X = expression matrix; matrix of shape n_obs x n_vars
    y = (cell type) labels; array/list of shape n_obs
    sex_labels = 'male' or 'female' label for each entry; array/list of shape n_obs
    proportion_female = desired proportion of female cells; float between 0 and 1
    classifier = 'knn' or 'rf'
    k = value if using knn; default 30
    ----------
    
    Returns:

    male_pred, female_pred = array of shape n_obs; contains the prediction on the male 
                                or female test set
    y_male_test, y_female_test = array of shape n_obs; contains the true labels of the 
                                male or female test set
    
    '''
    
    np.random.seed(42)
    
    print('entered training function!', flush = True)
    
    male_indices = np.where(sex_labels == 'male')[0]
    female_indices = np.where(sex_labels == 'female')[0]

    X_male = X[male_indices]
    y_male = y[male_indices]
    X_female = X[female_indices]
    y_female = y[female_indices]
    
    
    X_female_train, X_female_test, y_female_train, y_female_test = train_test_split(
        X_female, y_female, test_size=0.2, stratify=y_female, random_state=42)
    
    # compute what to pass to test_size to get equal test set size to the female set
    male_proportion = X_female_test.shape[0] / X_male.shape[0]

    X_male_train, X_male_test, y_male_train, y_male_test = train_test_split(
        X_male, y_male, test_size=male_proportion, stratify=y_male, random_state=42)
    

    # merge training sets back together
    # X_train = np.concatenate([X_male_train, X_female_train])
    X_train = vstack([X_male_train, X_female_train])
    y_train = np.concatenate([y_male_train, y_female_train])
    sex_labels_train = ['male'] * X_male_train.shape[0] + ['female'] * X_female_train.shape[0]


    # Select female cells based on proportion_female
    selected_indices = fixed_select_indices_by_proportion(sex_labels_train, proportion_female)
    X_selected = X_train.tocsr()[selected_indices]
    y_selected = y_train[selected_indices]

    print('selected training!', flush = True)

    # Initialize classifier
    if classifier == 'knn':
        clf = KNeighborsClassifier(n_neighbors=k)
    elif classifier == 'rf':
        clf = RandomForestClassifier(n_jobs=-1)
        

    print('initialized classif!', flush = True)
    # Train
    clf.fit(X_selected, y_selected)
        # "calling fit() more than once will overwrite what was learned by any previous fit()"

    print('done training!', flush = True)
    
    # Predict
    male_pred = clf.predict(X_male_test)
    female_pred = clf.predict(X_female_test)
     
    return male_pred, y_male_test, female_pred, y_female_test


def final_evaluate_clf(predictions, true_labels, classes, prop, sex):
    '''
    This is meant to be run separately on the male and female results of the train function.
    ---------
    
    Parameters:

    predictions = predictions on a test set
    true_labels = true labels on that test set
    classes = sorted list of classes
    prop = proportion of female cells
    sex = 'male' or 'female'; dneotes what test set we are working with
        # prop and sex are only used for naming the confusion matrix plots
    
    ----------
    
    Returns:

    accuracy = accuracy score
    f1_per_class = array of shape n_classes; each entry is the f1 score of that class
    median_f1 = median of the f1_per_class array
    precision_per_class = array of shape n_classes; each entry is the precision score of that class
    median_precision = median of the precision_per_class array
    cm = confusion matrix
    cm_normalized = normalized confusion matrix
        (the function saves plots of each confusion matrix)
    
    '''
    n_classes = len(classes)

    # Accuracy scores
    accuracy = accuracy_score(true_labels, predictions)

    # F1 per class
    f1_per_class = robust_f1_score(true_labels, predictions, classes)

    # Median F1
    median_f1 = np.median(f1_per_class)

    # Precision per class
    precision_per_class = precision_score(true_labels, predictions, average=None)

    # Median precision
    median_precision = np.median(precision_per_class)

    # Confusion matrix
    cm = confusion_matrix(true_labels, predictions, labels=classes)
    # Normalized confusion matrix:
    eps = 1e-9  # to avoid division by 0
    cm_normalized = cm.astype('float') / (cm.sum(axis=1)[:, np.newaxis] + eps)
    
    # Create dictionary
    metrics = {
        'accuracy': accuracy,
        'f1_scores': f1_per_class,
        'median_f1': median_f1,
        'precision_scores': precision_per_class,
        'median_precision': median_precision,
        'aggregated_confusion_matrix': cm,
        'normalized_aggregated_confusion_matrix': cm_normalized
    }
        
        
    plot_confusion(cm, classes, f'Prop {prop}, {sex} test set Confusion Matrix', False)
    plot_confusion(cm_normalized, classes, f'Prop {prop}, {sex} test set Normalized Confusion Matrix', True)
    
    
    return metrics



# Helper functions:

def robust_f1_score(y_true, y_pred, labels):
    f1_scores = []
    for label in labels:
        if (y_pred == label).sum() == 0:
            # no predictions for this class -> set F1 to 0
            f1_scores.append(0.0)
        else:
            f1 = f1_score(y_true == label, y_pred == label, zero_division=0)
            f1_scores.append(f1)
    return np.array(f1_scores)

def plot_confusion(confusion_matrix, classes, title, normalize = False):
    '''
    Plot and save the current confusion matrix.
    Make sure to pass a title and the normalize parameter (if the matrix is normalized).
    '''
    if not os.path.exists('cms'):
        os.makedirs('cms')
    if not os.path.exists('norm_cms'):
        os.makedirs('norm_cms')
    
    plt.clf()
    plt.figure(figsize=(10, 8))
    if normalize:
        sns.heatmap(confusion_matrix, annot=True, fmt=".3f", cmap="Blues", xticklabels=classes, yticklabels=classes)
    else:
        sns.heatmap(confusion_matrix, annot=True, fmt="g", cmap="Blues", xticklabels=classes, yticklabels=classes)
    plt.title('Confusion Matrix')
    plt.xlabel('Predicted Label')
    plt.ylabel('True Label')
    plt.title(title)
    if normalize:
        plt.savefig(f'norm_cms/{title}.png', bbox_inches='tight')
    else:
        plt.savefig(f'cms/{title}.png', bbox_inches='tight')
    plt.close()


def fixed_select_indices_by_proportion(sex_labels, proportion_female):
    np.random.seed(42)
    sex_labels_series = pd.Series( (el for el in sex_labels) )
    
    female_indices = np.where(sex_labels_series == 'female')[0]
    male_indices = np.where(sex_labels_series == 'male')[0]
    
    fixed_size = min(len(female_indices), len(male_indices))
    
    np.random.shuffle(female_indices)
    np.random.shuffle(male_indices)

    num_female_cells = int(fixed_size * proportion_female)
    num_male_cells = fixed_size - num_female_cells
        # total will always be fixed_size
        # this works for cases with prop 0% or 100% --> no need to handle them separately
    
    # adjust in case of rounding errors
    num_female_cells = min(num_female_cells, len(female_indices))
    num_male_cells = min(num_male_cells, len(male_indices))

    selected_female_indices = female_indices[:num_female_cells]
    selected_male_indices = male_indices[:num_male_cells]

    return np.concatenate([selected_female_indices, selected_male_indices])