import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import shap
import os
import matplotlib.gridspec as gridspec
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
import xgboost as xgb

# 1. Przygotowanie folderów
os.makedirs('results/figures', exist_ok=True)

# 2. Ładowanie danych
# Ścieżka zaktualizowana zgodnie z Twoim środowiskiem
df = pd.read_csv("../..data/genepy_matrix.csv", index_col=0)
X = df.drop('target', axis=1)
y = df['target']

# Mapowanie etykiet dla czytelności wykresów
target_names = {0: 'Control', 1: 'UC'}
df_plot = df.copy()
df_plot['Group'] = df['target'].map(target_names)

# 3. Trening finalnego modelu (XGBoost)
ratio = float(sum(y == 0)) / sum(y == 1)

# Parametry na podstawie Twojej hiperparametryzacji i rygoru z publikacji
model = xgb.XGBClassifier(
    n_estimators=500, 
    max_depth=5, 
    learning_rate=0.05, 
    scale_pos_weight=ratio,
    random_state=42,
    eval_metric='logloss'
)

pipeline = Pipeline([
    ('scaler', StandardScaler()),
    ('classifier', model)
])

pipeline.fit(X, y)

# 4. Obliczanie wartości SHAP 
best_clf = pipeline.named_steps['classifier']
X_scaled = pipeline.named_steps['scaler'].transform(X)
X_scaled_df = pd.DataFrame(X_scaled, columns=X.columns)

explainer = shap.TreeExplainer(best_clf)
shap_values = explainer.shap_values(X_scaled)

# Wybieramy wartości dla klasy pozytywnej (UC)
if isinstance(shap_values, list):
    shap_v = shap_values[1]
else:
    shap_v = shap_values

# 5. Identyfikacja TOP 10 genów wg SHAP do Panelu A
top_idx = np.argsort(np.abs(shap_v).mean(0))[::-1][:10]
top_genes = X.columns[top_idx].tolist()

# --- KOMPOZYCJA WYKRESU (PUBLICATION STYLE) ---
# Zwiększona szerokość figury, aby uniknąć ścisku
fig = plt.figure(figsize=(20, 10)) 

# Ustawiamy proporcje: Panel B (SHAP) potrzebuje więcej miejsca na nazwy genów (1 : 1.6)
gs = gridspec.GridSpec(1, 2, width_ratios=[1, 1.6])

# PANEL A: Violin Plots (Surowe dane GenePy)
ax1 = plt.subplot(gs[0])

df_melted = df_plot.melt(id_vars=['Group'], value_vars=top_genes, var_name='Gene', value_name='GenePy Score')

# Naprawa Warningu: Jawne ustawienie lokalizacji ticków
ax1.set_xticks(range(len(top_genes)))

# inner="point" i bw_method=0.2 pomagają zwizualizować rzadkie dane (sparsity)
sns.violinplot(
    data=df_melted, x='Gene', y='GenePy Score', hue='Group', 
    split=True, inner="point", palette='muted', ax=ax1, 
    cut=0, bw_method=0.2
)

ax1.set_title("A] Distributions of GenePy scores (Top 10 Genes)", fontsize=16, fontweight='bold', loc='left', pad=20)
ax1.set_xticklabels(top_genes, rotation=45, ha='right', fontsize=12)
ax1.set_ylabel("Raw GenePy Score", fontsize=13)
ax1.set_xlabel("Gene", fontsize=13)
ax1.grid(axis='y', linestyle='--', alpha=0.3)

# PANEL B: SHAP Summary Plot
ax2 = plt.subplot(gs[1])
plt.sca(ax2) # Ustawienie kontekstu dla biblioteki SHAP

# shap.summary_plot przejmuje kontrolę nad osiami, GridSpec trzyma go w ryzach
shap.summary_plot(
    shap_v, X_scaled_df, feature_names=X.columns, 
    show=False, plot_size=None, alpha=0.6
)

plt.title("B] SHAP values: GenePy scores contribution to UC", fontsize=16, fontweight='bold', loc='left', pad=20)
ax2.tick_params(labelsize=12)
plt.xlabel("SHAP value (impact on model output)", fontsize=13)

# Finalne dopasowanie i zapis
plt.tight_layout()
plt.savefig('results/figures/final_analysis_publication.png', dpi=300, bbox_inches='tight')
plt.show()

print("\nWykres został wygenerowany pomyślnie.")
print("Lokalizacja: results/figures/final_analysis_publication.png")