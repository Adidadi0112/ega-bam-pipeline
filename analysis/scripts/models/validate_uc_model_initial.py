import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import StratifiedKFold, cross_val_score
from sklearn.metrics import roc_auc_score, roc_curve, auc

# 1. Load data
df = pd.read_csv("/home/adam/projects/ega-bam-pipeline/analysis/genepy_matrix.csv", index_col=0)
X = df.drop('target', axis=1)
y = df['target']

# 2. Configure cross-validation (5 folds)
cv = StratifiedKFold(n_splits=5, shuffle=True, random_state=42)
model = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)

# 3. Calculate stable AUC
auc_scores = cross_val_score(model, X, y, cv=cv, scoring='roc_auc')

print("--- 5-FOLD CROSS-VALIDATION RESULTS ---")
print(f"AUC for every fold: {auc_scores}")
print(f"Mean ROC-AUC: {np.mean(auc_scores):.3f} (+/- {np.std(auc_scores):.3f})")

# 4. Calculate stable feature importance across folds
all_importances = []

tprs = []
base_fpr = np.linspace(0, 1, 101)
plt.figure(figsize=(15, 6))

# Plot 1: ROC curves for each fold
plt.subplot(1, 2, 1)
for i, (train, test) in enumerate(cv.split(X, y)):
    model.fit(X.iloc[train], y.iloc[train])
    y_prob = model.predict_proba(X.iloc[test])[:, 1]
    fpr, tpr, _ = roc_curve(y.iloc[test], y_prob)
    roc_auc = auc(fpr, tpr)
    plt.plot(fpr, tpr, lw=1, alpha=0.3, label=f'ROC fold {i} (AUC = {roc_auc:.2f})')
    
    # For interpolation of TPR at fixed FPR points
    tpr_interp = np.interp(base_fpr, fpr, tpr)
    tpr_interp[0] = 0.0
    tprs.append(tpr_interp)
    
    # Collect feature importances for stability analysis
    all_importances.append(model.feature_importances_)

# Plot mean ROC curve
mean_tprs = np.mean(tprs, axis=0)
plt.plot(base_fpr, mean_tprs, 'b', label=f'Mean ROC (AUC = {np.mean(auc_scores):.2f})', lw=2)
plt.plot([0, 1], [0, 1], 'r--')
plt.title('ROC Curves (5-Fold CV)')
plt.legend()

# 5. Stability analysis of feature importance
mean_importance = np.mean(all_importances, axis=0)
std_importance = np.std(all_importances, axis=0)

importance_df = pd.DataFrame({
    'Gene': X.columns,
    'Mean_Importance': mean_importance,
    'Std': std_importance
}).sort_values(by='Mean_Importance', ascending=False)

# Plot 2: Top 15 most stable important genes
plt.subplot(1, 2, 2)
sns.barplot(x='Mean_Importance', y='Gene', data=importance_df.head(15), palette='viridis')
plt.errorbar(importance_df.head(15)['Mean_Importance'], range(15), 
             xerr=importance_df.head(15)['Std'], fmt='none', c='black', capsize=3)
plt.title('Top 15 Stable Predictors (Mean Importance)')

plt.tight_layout()
plt.savefig('validation_results.png')
print("\nTop 10 most stable predictors:")
print(importance_df[['Gene', 'Mean_Importance']].head(10))