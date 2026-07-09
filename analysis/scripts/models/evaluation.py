import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import warnings
import os

from sklearn.ensemble import RandomForestClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.svm import SVC
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.inspection import permutation_importance
import xgboost as xgb

from sklearn.model_selection import GridSearchCV, StratifiedKFold
from sklearn.metrics import (roc_auc_score, auc, roc_curve, accuracy_score, precision_score, 
                             recall_score, f1_score, confusion_matrix, ConfusionMatrixDisplay, 
                             precision_recall_curve, average_precision_score)
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline

warnings.filterwarnings("ignore", category=FutureWarning)
warnings.filterwarnings("ignore", category=UserWarning)

os.makedirs('results/stats', exist_ok=True)
os.makedirs('results/figures', exist_ok=True)

# 1. Load data and preprocess
data = pd.read_csv("../../data/genepy_matrix_v2_less.csv", index_col=0)

X = data.drop('target', axis=1)
y = data['target']

ratio = float(sum(y == 0)) / sum(y == 1)

# 2. Declare classifiers to evaluate
classifiers = {
    'RandomForest': {
        'model': RandomForestClassifier(random_state=42, class_weight='balanced'),
        'params': {
            'n_estimators': [50, 100, 200],
            'max_depth': [None, 5, 10, 20],
            'min_samples_split': [2, 5, 10],
            'min_samples_leaf': [1, 2, 4],
            'max_features': ['sqrt', 'log2']
        }
    },
    'Logistic Regression': {
        'model': LogisticRegression(penalty='l1', solver='liblinear', random_state=42, class_weight='balanced'),
        'params': {
            'C': np.logspace(-4, 4, 10)
        }
    },
    'SVM': {
        'model': SVC(probability=True, random_state=42, class_weight='balanced'),
        'params': {
            'C': [0.1, 1, 10, 100],
            'kernel': ['linear', 'rbf'],
            'gamma': ['scale', 'auto', 0.001, 0.01, 0.1, 1]
        }
    },
    'KNN': {
        'model': KNeighborsClassifier(),
        'params': {
            'n_neighbors': [3, 5, 7, 9],
            # 'weights': ['uniform', 'distance'],
            # 'metric': ['euclidean', 'manhattan']
        }
    },
    'Naive Bayes': {
        'model': GaussianNB(),
        'params': {
            'var_smoothing': np.logspace(-10, -7, 5)
        }
    },
    'XGBoost': {
        'model': xgb.XGBClassifier(scale_pos_weight=ratio, eval_metric='logloss', random_state=42),
        'params': {
            'n_estimators': [50, 100, 200, 500],
            'max_depth': [3, 5, 7],
            'learning_rate': [0.01, 0.05, 0.1, 0.2],
            'subsample': [0.8, 1.0],
        }
    }
}

# 3. Nested cross-validation configuration

outer_cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
inner_cv = StratifiedKFold(n_splits=3, shuffle=True, random_state=42)

results = []
conf_matrices = {}
roc_plotting_data = {}

print("Starting nested cross-validation for model evaluation...")

# 4. Main loop (Outer loop)
for name, config in classifiers.items():
    print(f'--- Evaluating {name} ---')

    # Initialize lists to store scores and predictions
    best_params_per_fold = []
    outer_scores = {'auc': [], 'accuracy': [], 'precision': [], 'recall': [], 'f1': [], }
    combined_y_true = []
    combined_y_pred = []
    combined_y_prob = []

    for i, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        print(f' Fold {i + 1}/5')

        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]

        pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', config['model'])
    ])
        current_params = {f'classifier__{key}': value for key, value in config['params'].items()}

        # Inner loop for hyperparameter tuning
        grid_search = GridSearchCV(estimator=pipeline, param_grid= current_params, cv=inner_cv, scoring='roc_auc')
        grid_search.fit(X_train, y_train)

        # Store best parameters for this fold
        best_params_per_fold.append({
            'fold': i + 1,
            'params': grid_search.best_params_
        })

        # Best model from inner loop
        best_model = grid_search.best_estimator_
        y_pred = best_model.predict(X_test)
        y_prob = best_model.predict_proba(X_test)[:, 1]

        # Store true and predicted labels for confusion matrix
        outer_scores['auc'].append(roc_auc_score(y_test, y_prob))
        outer_scores['accuracy'].append(accuracy_score(y_test, y_pred))
        outer_scores['precision'].append(precision_score(y_test, y_pred))
        outer_scores['recall'].append(recall_score(y_test, y_pred))
        outer_scores['f1'].append(f1_score(y_test, y_pred))

        combined_y_true.extend(y_test)
        combined_y_pred.extend(y_pred)
        combined_y_prob.extend(y_prob)

    # results aggregation
    results.append({
        'Model': name,
        'Mean AUC': np.mean(outer_scores['auc']),
        'Std AUC': np.std(outer_scores['auc']),
        'Mean Accuracy': np.mean(outer_scores['accuracy']),
        'Mean Recall': np.mean(outer_scores['recall']),
        'Mean Precision': np.mean(outer_scores['precision']),
        'Mean F1': np.mean(outer_scores['f1'])
    })

    # Store best parameters for this model
    with open(f'results/stats/best_params_{name.replace(" ", "_")}.txt', 'w') as f:
        f.write(f"Best hyperparameters for {name} found in each outer fold:\n")
        for entry in best_params_per_fold:
            f.write(f"Fold {entry['fold']}: {entry['params']}\n")

    roc_plotting_data[name] = {
        'true': combined_y_true,
        'prob': combined_y_prob
    }

    conf_matrices[name] = confusion_matrix(combined_y_true, combined_y_pred)



# 5. Display results
results_df = pd.DataFrame(results).sort_values(by='Mean AUC', ascending=False)
print("\n--- MODEL PERFORMANCE COMPARISON ---")
print(results_df)

# 6. Plot confusion matrices
n_models = len(conf_matrices)
fix, axes = plt.subplots(2, 3, figsize=(18, 10))
axes_flat = axes.flatten()

for i, (model_name, matrix) in enumerate(conf_matrices.items()):
    disp = ConfusionMatrixDisplay(confusion_matrix=matrix, display_labels=['Control', 'UC'])
    disp.plot(ax=axes_flat[i], cmap='Blues', colorbar=False)
    axes_flat[i].set_title(f'{model_name} Confusion Matrix')

plt.tight_layout()
plt.savefig('results/figures/model_comparison_cm_v3.png')
plt.show()

print("\n--- GENEROWANIE KONSENSUSU WAŻNOŚCI (PERMUTATION IMPORTANCE) ---")
os.makedirs('results/figures/importance', exist_ok=True)

consensus_data = []

for name, config in classifiers.items():
    print(f"Przetwarzanie: {name}...")
    
    # 1. Trenujemy finalny model na całym zbiorze
    pipeline = Pipeline([
        ('scaler', StandardScaler()),
        ('classifier', config['model'])
    ])
    pipeline.fit(X, y)
    
    # 2. Obliczamy Permutation Importance (scoring='roc_auc')
    r = permutation_importance(pipeline, X, y, n_repeats=10, random_state=42, n_jobs=-1, scoring='roc_auc')
    
    # 3. Zbieramy Top 10 genów dla tego modelu
    for i in r.importances_mean.argsort()[::-1][:10]:
        consensus_data.append({
            'Model': name,
            'Gene': X.columns[i],
            'Importance': r.importances_mean[i]
        })

# 4. Tworzymy tabelę zbiorczą
consensus_df = pd.DataFrame(consensus_data)

# 5. Sprawdzamy, które geny powtarzają się najczęściej w Top 10 wszystkich modeli
gene_votes = consensus_df.groupby('Gene').size().sort_values(ascending=False).reset_index()
gene_votes.columns = ['Gene', 'Model_Votes']

print("\n--- GENY Z NAJWIĘKSZYM KONSENSUSEM (Top 10 we wszystkich modelach) ---")
print(gene_votes.head(10))

# 6. Wizualizacja konsensusu
plt.figure(figsize=(12, 6))
sns.barplot(x='Model_Votes', y='Gene', data=gene_votes.head(15), palette='mako')
plt.title("Geny najczęściej wskazywane jako ważne przez różne modele")
plt.xlabel("Liczba modeli, które uznały gen za istotny (Top 10)")
plt.savefig('results/figures/importance/gene_consensus_ranking_v3.png')
plt.show()

# 7. Generating ROC curves for all models
print("\nGenerating combined ROC curves...")

plt.figure(figsize=(12, 9))

# Color configuration - using a professional palette
colors = sns.color_palette("Set1", len(roc_plotting_data))

# Plotting the "random guess" line
plt.plot([0, 1], [0, 1], linestyle='--', lw=2, color='#7f8c8d', label='Random Guess (AUC = 0.50)', alpha=0.8)

# Loop through each model's results
for i, (name, data) in enumerate(roc_plotting_data.items()):
    # Calculating FPR and TPR based on aggregated data from all folds
    fpr, tpr, _ = roc_curve(data['true'], data['prob'])
    roc_auc = auc(fpr, tpr)
    
    # Plotting the curve for each specific model
    plt.plot(fpr, tpr, lw=3, color=colors[i],
             label=f'{name} (AUC = {roc_auc:.2f})')

# Plot aesthetics and English labels
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.05])
plt.xlabel('False Positive Rate (1 - Specificity)', fontsize=14, labelpad=10)
plt.ylabel('True Positive Rate (Sensitivity)', fontsize=14, labelpad=10)
plt.title('ROC Curve Comparison for UC Classification\n(Nested Cross-Validation Results)', 
          fontsize=16, fontweight='bold', pad=20)
plt.legend(loc="lower right", fontsize=12, frameon=True, shadow=True)
plt.grid(alpha=0.3, linestyle='--')

# Ensure labels and layout are optimized
plt.tight_layout()

# Saving in vector formats for high-quality publication
plt.savefig('results/figures/multi_model_roc_comparison_v3.svg', format='svg', bbox_inches='tight')
plt.savefig('results/figures/multi_model_roc_comparison_v3.pdf', format='pdf', bbox_inches='tight')

plt.show()

print(f"Success! ROC plot has been saved in the results/figures/ folder.")

# -- 8. Generating Precision-Recall curves (English Version) ---
print("\nGenerating collective Precision-Recall plot...")

# Set publication style parameters
plt.figure(figsize=(12, 9))
sns.set_style("whitegrid", {'axes.grid': True, 'grid.linestyle': '--'})
colors = sns.color_palette("Set1", len(roc_plotting_data))

# Calculate baseline (no-skill model)
# Baseline in PR curve is the ratio of positive samples (UC) to total samples
baseline_ratio = sum(y == 1) / len(y)
plt.plot([0, 1], [baseline_ratio, baseline_ratio], linestyle='--', lw=2, color='#7f8c8d', 
         label=f'Baseline (AP = {baseline_ratio:.2f})', alpha=0.8)

# Iterate through each model's combined data
for i, (name, data) in enumerate(roc_plotting_data.items()):
    # Compute Precision-Recall curve pairs
    precision, recall, _ = precision_recall_curve(data['true'], data['prob'])
    # Compute Average Precision (AP) score
    ap_score = average_precision_score(data['true'], data['prob'])
    
    # Plot curve for specific model
    plt.plot(recall, precision, lw=3, color=colors[i],
             label=f'{name} (AP = {ap_score:.2f})')

# Aesthetics and formatting (English labels)
plt.xlim([-0.01, 1.01])
plt.ylim([-0.01, 1.05])
plt.xlabel('Recall (Sensitivity)', fontsize=14, labelpad=10)
plt.ylabel('Precision (Positive Predictive Value)', fontsize=14, labelpad=10)
plt.title('Precision-Recall Curve Comparison for UC Classification\n(Nested Cross-Validation Results)', 
          fontsize=16, fontweight='bold', pad=20)
plt.legend(loc="upper right", fontsize=12, frameon=True, shadow=True)
plt.grid(alpha=0.3, linestyle='--')

# Ensure labels are not cut off
plt.tight_layout()

# Save in vector formats for high-quality printing
plt.savefig('results/figures/multi_model_pr_comparison_v3.svg', format='svg', bbox_inches='tight')
plt.savefig('results/figures/multi_model_pr_comparison_v3.pdf', format='pdf', bbox_inches='tight')

plt.show()

print(f"Success! Precision-Recall plot saved to results/figures/")