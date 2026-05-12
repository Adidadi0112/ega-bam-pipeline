import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

vep_path = '../../05_annotation/vep_results_vqsr.txt'

# Loading with VEP header support
header_line = 0
with open(vep_path, 'r') as f:
    for i, line in enumerate(f):
        if line.startswith('#Uploaded'):
            header_line = i
            break

df = pd.read_csv(vep_path, sep='\t', skiprows=header_line)
df.columns = [c.lstrip('#') for c in df.columns]

# Extracting CADD if it's in the Extra column
if 'CADD_PHRED' not in df.columns and 'Extra' in df.columns:
    df['CADD_PHRED'] = df['Extra'].str.extract(r'CADD_PHRED=([^;]+)')[0]
df['CADD_PHRED'] = pd.to_numeric(df['CADD_PHRED'], errors='coerce').fillna(0)

# --- VISUALIZATION ---
plt.figure(figsize=(16, 6))

# 1. Variant types (Consequences)
plt.subplot(1, 3, 1)
df['Consequence'].value_counts().head(10).plot(kind='bar', color='skyblue')
plt.title('Top 10 Variant Consequences')
plt.xticks(rotation=45, ha='right')

# 2. Distribution of CADD scores
plt.subplot(1, 3, 2)
sns.histplot(df['CADD_PHRED'], bins=30, kde=True, color='salmon')
plt.axvline(15, color='red', linestyle='--', label='Likely Deleterious (>15)')
plt.title('Distribution of CADD PHRED Scores')
plt.legend()

# 3. Genes with the highest number of variants
plt.subplot(1, 3, 3)
df['SYMBOL'].value_counts().head(15).plot(kind='bar', color='lightgreen')
plt.title('Top 15 Most Variable Genes')

plt.tight_layout()
plt.savefig('annotation_exploration_vqsr.png')
print("CADD PHRED Statistics:")
print(df['CADD_PHRED'].describe())
print("\nTop 5 Most Common Consequences:")
print(df['Consequence'].value_counts().head(5))