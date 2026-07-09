import pandas as pd
import os
import gzip

# 1. Konfiguracja ścieżek na podstawie dokumentacji
VEP_FILE = "../../05_annotation/vep_results_vqsr.txt"
VCF_FILE = "../../02_targeting/uc_vqsr_cleaned.vcf.gz" # Finalny, czysty VCF
MATRIX_FILE = "../../data/genepy_matrix.csv"
OUTPUT_DIR = "results/bio_verification"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Geny, które SHAP wskazał jako najważniejsze
target_genes = ['GALC', 'SCN7A', 'IP6K3', 'GJA3']

# 2. Pobranie i mapowanie pacjentów
print("Krok 0: Mapowanie pacjentów...")
df_matrix = pd.read_csv(MATRIX_FILE, index_col=0)
uc_patients_raw = df_matrix[df_matrix['target'] == 1].index.tolist()

# Mapujemy długie ID (macierz) na krótkie ID (VCF)
# Np. NA18939.mapped... -> NA18939
patient_mapping = {p.split('.')[0]: p for p in uc_patients_raw}
uc_core_ids = set(patient_mapping.keys())

# 3. Wczytywanie VEP (Lokalizacje)
print("Krok 1: Wczytywanie adnotacji VEP...")
def extract_from_extra(extra_str, key):
    if pd.isna(extra_str): return None
    for part in extra_str.split(';'):
        if '=' in part:
            k, v = part.split('=')
            if k == key: return v
    return None

with open(VEP_FILE, 'r') as f:
    skip_count = 0
    for i, line in enumerate(f):
        if line.startswith('#Uploaded_variation'):
            skip_count = i
            break

df_vep = pd.read_csv(VEP_FILE, sep='\t', skiprows=skip_count)
df_vep.rename(columns={'#Uploaded_variation': 'ID'}, inplace=True)
df_vep['SYMBOL'] = df_vep['Extra'].apply(lambda x: extract_from_extra(x, 'SYMBOL'))
df_vep['IMPACT'] = df_vep['Extra'].apply(lambda x: extract_from_extra(x, 'IMPACT'))

interesting = df_vep[df_vep['SYMBOL'].isin(target_genes)].copy()
# Zachowujemy tylko najbardziej szkodliwą adnotację dla każdej pozycji
interesting['imp_val'] = interesting['IMPACT'].map({'HIGH': 3, 'MODERATE': 2, 'LOW': 1}).fillna(0)
interesting = interesting.sort_values('imp_val', ascending=False).drop_duplicates('Location')
variant_map = interesting.set_index('Location')[['SYMBOL', 'Consequence', 'IMPACT']].to_dict('index')

print(f"Znaleziono {len(variant_map)} pozycji wariantów dla wybranych genów.")

# 4. Przeszukiwanie skompresowanego VCF
print(f"Krok 2: Przeszukiwanie VCF ({VCF_FILE})...")
results = []

with gzip.open(VCF_FILE, 'rt') as vcf:
    samples = []
    for line in vcf:
        if line.startswith('##'): continue
        if line.startswith('#CHROM'):
            raw_vcf_samples = line.strip().split('\t')[9:]
            # Normalizujemy nazwy z VCF do "rdzenia"
            samples = [s.split('.')[0] for s in raw_vcf_samples]
            continue
        
        cols = line.strip().split('\t')
        # VEP używa formatu 5:123, VCF może mieć chr5:123
        loc = f"{cols[0].replace('chr', '')}:{cols[1]}"
        
        if loc in variant_map:
            genotypes = cols[9:]
            for i, gt in enumerate(genotypes):
                vcf_core_id = samples[i]
                # Czy ten pacjent jest w naszej grupie UC?
                if vcf_core_id in uc_core_ids:
                    # Czy ma mutację (0/1 lub 1/1)?
                    gt_call = gt.split(':')[0]
                    if gt_call not in ['0/0', './.', '0|0']:
                        res = variant_map[loc].copy()
                        res['Sample_ID'] = patient_mapping[vcf_core_id]
                        res['Location'] = loc
                        res['Genotype'] = gt_call
                        results.append(res)

# 5. Eksport wyników
if results:
    final_df = pd.DataFrame(results)
    final_df.to_csv(os.path.join(OUTPUT_DIR, "final_clinical_evidence.csv"), index=False)
    print(f"\n--- ANALIZA ZAKOŃCZONA ---")
    print(f"Zidentyfikowano {len(final_df)} wariantów u pacjentów UC.")
    print(f"Top 5 wariantów o najwyższym IMPACT:")
    print(final_df.sort_values('IMPACT').head())
else:
    print("\nNadal brak dopasowań. Sprawdź ręcznie: zgrep '#CHROM' " + VCF_FILE + " | cut -f 10")