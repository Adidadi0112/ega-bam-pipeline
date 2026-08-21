import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import os
import matplotlib.gridspec as gridspec
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.neighbors import KNeighborsClassifier

# 1. Przygotowanie folderów
os.makedirs('results/figures', exist_ok=True)

# 2. Ładowanie danych
df = pd.read_csv("analysis/data/genepy_matrix_v2.csv", index_col=0)
X = df.drop('target', axis=1)
y = df['target']

target_names = {0: 'Control', 1: 'UC'}
df_plot = df.copy()
df_plot['Group'] = df['target'].map(target_names)

# 3. Trening finalnego modelu (k-NN)
# Wybieramy k=5 (lub Twoją wartość z optymalizacji)
model = KNeighborsClassifier(n_neighbors=3)

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', model)
])

pipeline.fit(X, y)

# 4. Obliczanie wartości SHAP (KernelExplainer)
X_scaled = pipeline.named_steps['scaler'].transform(X)
X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

# KernelExplainer jest wolniejszy, więc używamy zredukowanego tła (np. 50 próbek)
# Pozwala to modelowi zrozumieć "średnie" zachowanie danych
background = shap.kmeans(X_scaled, 10) 

# Definiujemy funkcję predykcji prawdopodobieństwa dla klasy 1 (UC)
predict_fn = lambda x: pipeline.named_steps['classifier'].predict_proba(x)[:, 1]

explainer = shap.KernelExplainer(predict_fn, background)

# Obliczamy SHAP dla wszystkich pacjentów (może to zająć chwilę)
print("Obliczanie wartości SHAP dla k-NN (KernelExplainer)...")
shap_values = explainer.shap_values(X_scaled)

# W KernelExplainer dla jednej funkcji wyjściowej shap_values jest macierzą (n_samples, n_features)
shap_v = shap_values

# 5. Identyfikacja TOP 10 genów wg SHAP
top_idx = np.argsort(np.abs(shap_v).mean(0))[::-1][:10]
top_genes = X.columns[top_idx].tolist()

pd.DataFrame({
    'Gene': X.columns,
    'Mean_Absolute_SHAP': np.abs(shap_v).mean(0)
}).sort_values('Mean_Absolute_SHAP', ascending=False).to_csv(
    'results/stats/knn_shap_feature_ranking_v2.csv', index=False
)

# --- KOMPOZYCJA WYKRESU (PUBLICATION STYLE) ---
fig = plt.figure(figsize=(22, 10)) 
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.4])

# PANEL A: Violin Plots
ax1 = plt.subplot(gs[0])
df_melted = df_plot.melt(id_vars=['Group'], value_vars=top_genes, var_name='Gene', value_name='GenePy Score')

sns.violinplot(
    data=df_melted, x='Gene', y='GenePy Score', hue='Group', 
    split=True, inner="point", palette='muted', ax=ax1, 
    cut=0, bw_method=0.2
)

ax1.set_title("Distributions of GenePy scores (Top 10 k-NN Predictors)", fontsize=16, fontweight='bold', loc='left', pad=20)
ax1.set_xticklabels(top_genes, rotation=45, ha='right', fontsize=12)
ax1.set_ylabel("Raw GenePy Score", fontsize=13)
ax1.grid(axis='y', linestyle='--', alpha=0.3)

# PANEL B: SHAP Summary Plot
ax2 = plt.subplot(gs[1])
plt.sca(ax2)

shap.summary_plot(
    shap_v, X_scaled_df, feature_names=X.columns, 
    show=False, plot_size=None, alpha=0.6
)

plt.title("k-NN SHAP values: Similarity-based contribution to UC", fontsize=16, fontweight='bold', loc='left', pad=20)
ax2.tick_params(labelsize=12)
plt.xlabel("SHAP value (impact on UC probability)", fontsize=13)

plt.tight_layout()
plt.savefig('results/figures/final_analysis_knn_v2.svg', format='svg', dpi=300, bbox_inches='tight')
plt.show()

print("\nAnaliza SHAP dla k-NN zakończona.")
