import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.model_selection import train_test_split
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix, roc_auc_score, roc_curve

# 1. Load data
df = pd.read_csv("genepy_matrix.csv", index_col=0)
X = df.drop('target', axis=1)
y = df['target']

# 2. Split into training and testing sets (stratified to maintain group proportions)
X_train, X_test, y_train, y_test = train_test_split(
    X, y, test_size=0.2, random_state=42, stratify=y
)

# 3. Initialize and train the model
model = RandomForestClassifier(n_estimators=100, class_weight='balanced', random_state=42)
model.fit(X_train, y_train)

# 4. Predictions and evaluation
y_pred = model.predict(X_test)
y_prob = model.predict_proba(X_test)[:, 1]

print("--- CLASSIFICATION REPORT ---")
print(classification_report(y_test, y_pred))
print(f"ROC-AUC Score: {roc_auc_score(y_test, y_prob):.2f}")

# 5. Extract most important genes (Feature Importance)
importances = pd.DataFrame({
    'Gene': X.columns,
    'Importance': model.feature_importances_
}).sort_values(by='Importance', ascending=False)

# --- VISUALIZATION ---
plt.figure(figsize=(15, 5))

# Chart 1: Confusion Matrix
plt.subplot(1, 2, 1)
cm = confusion_matrix(y_test, y_pred)
sns.heatmap(cm, annot=True, fmt='d', cmap='Blues')
plt.title('Confusion Matrix')
plt.xlabel('Predicted')
plt.ylabel('Actual')

# Chart 2: Top 15 most important genes
plt.subplot(1, 2, 2)
sns.barplot(x='Importance', y='Gene', data=importances.head(15), palette='magma')
plt.title('Top 15 Predictor Genes (Feature Importance)')

plt.tight_layout()
plt.savefig('model_results.png')
print("\nTop 10 most important genes:")
print(importances.head(10))