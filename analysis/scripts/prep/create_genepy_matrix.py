from collections import defaultdict
import math

import pandas as pd
import gzip
import os

# --- PATHS CONFIGURATION ---
VEP_FILE = "../../05_annotation/vep_results_vqsr.txt"
CADD_FILE = "../../05_annotation/cadd_scores.tsv"
VCF_FILE = "../../02_targeting/uc_vqsr_rare_included.vcf.gz"
OUTPUT_FILE = "genepy_matrix.csv"

def get_cadd_df(path):
    skip = 0
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#Chrom'):
                skip = i
                break
    df = pd.read_csv(path, sep='\t', skiprows=skip)
    df.columns = [c.replace('#', '').strip() for c in df.columns]
    return df

def get_vep_df(path):
    skip = 0
    header_line = []
    with open(path, 'r') as f:
        for i, line in enumerate(f):
            if line.startswith('#Uploaded_variation'):
                skip = i
                header_line = line.strip().lstrip('#').split('\t')
                break
    df = pd.read_csv(path, sep='\t', skiprows=skip+1, header=None, names=header_line)
    
    if 'Extra' in df.columns:
        df['SYMBOL_EXT'] = df['Extra'].str.extract(r'SYMBOL=([^;]+)')
        df['Gene_Name'] = df['SYMBOL_EXT'].fillna(df['Gene'])

        # WYCIĄGANIE MAF Z KOLUMNY EXTRA (Priorytet dla East Asian, fallback na Global)
        df['MAF_EAS'] = df['Extra'].str.extract(r'gnomADe_EAS_AF=([^;]+)').astype(float)
        df['MAF_GLOBAL'] = df['Extra'].str.extract(r'gnomADe_AF=([^;]+)').astype(float)
        
        # Łączenie MAF i zabezpieczenie dla wariantów nieznanych (ustawienie na ultra-rzadkie 1e-5)
        df['MAF'] = df['MAF_EAS'].fillna(df['MAF_GLOBAL'])
        df['MAF'] = df['MAF'].fillna(1e-5).replace(0.0, 1e-5)
    else:
        df['Gene_Name'] = df.get('SYMBOL', df['Gene'])
        df['MAF'] = 1e-5
    return df

print("--- PATH 1: Loading CADD ---")
cadd = get_cadd_df(CADD_FILE)
cadd = cadd[['Chrom', 'Pos', 'Ref', 'Alt', 'PHRED']].copy()
cadd['Variant'] = cadd['Chrom'].astype(str) + "_" + cadd['Pos'].astype(str) + "_" + cadd['Ref'] + "/" + cadd['Alt']
print(f"Loaded {len(cadd)} CADD results.")

print("\n--- PATH 2: Loading VEP ---")
vep = get_vep_df(VEP_FILE)
vep['Variant'] = vep['Uploaded_variation'].str.replace(':', '_')
vep = vep[['Variant', 'Gene_Name', 'MAF']].copy()
print(f"Loaded {len(vep)} VEP annotations.")

print("\n--- PATH 3: Merging Data ---")
mapping = vep.merge(cadd[['Variant', 'PHRED']], on='Variant')
mapping_unique = mapping.sort_values('PHRED', ascending=False).drop_duplicates(subset=['Variant', 'Gene_Name'])
variant_to_gene_impact = mapping_unique.set_index(['Variant', 'Gene_Name'])[['PHRED', 'MAF']].to_dict('index')
print(f"Zmapowano {len(variant_to_gene_impact)} unikalnych par wariant-gen.")



print("\n--- PATH 4: Processing Genotypes from VCF ---")
samples = []
with gzip.open(VCF_FILE, 'rt') as f:
    for line in f:
        if line.startswith("#CHROM"):
            samples = line.strip().split('\t')[9:]
            break

unique_genes = list(set([g for (v, g) in variant_to_gene_impact.keys()]))
matrix = pd.DataFrame(0.0, index=samples, columns=unique_genes)

var_to_genes = defaultdict(list)
for v, g in variant_to_gene_impact.keys():
    var_to_genes[v].append(g)

with gzip.open(VCF_FILE, 'rt') as f:
    for line in f:
        if line.startswith("#"): continue
        cols = line.strip().split('\t')
        var_id = f"{cols[0]}_{cols[1]}_{cols[3]}/{cols[4]}"
                
        if var_id in var_to_genes:
            for gene in var_to_genes[var_id]:
                impact_data = variant_to_gene_impact[(var_id, gene)]
                phred = impact_data['PHRED']
                maf = impact_data['MAF']
                
                for i, sample_data in enumerate(cols[9:]):
                    gt = sample_data.split(':')[0]
                    
                    if gt in ['0/1', '1/0']:
                        # HETEROZYGOTA: Wzór GenePy ( Fi = MAF )
                        score = phred * math.log10(1.0 / maf)
                        matrix.at[samples[i], gene] += score
                        
                    elif gt == '1/1':
                        # HOMOZYGOTA: Wzór GenePy ( Fi = MAF^2 )
                        score = phred * math.log10(1.0 / (maf ** 2))
                        matrix.at[samples[i], gene] += score

print("\n--- PATH 5: Adding target column ---")
# Control group contains the phrase '.mapped.ILLUMINA.bwa.JPT.exome'
is_control = matrix.index.str.contains('.mapped.ILLUMINA.bwa.JPT.exome')

# By default, all are UC (1)
matrix['target'] = 1 
# Samples meeting the control condition get 0
matrix.loc[is_control, 'target'] = 0

print(f"Group statistics: UC={sum(matrix['target']==1)}, Control={sum(matrix['target']==0)}")

print(f"\n--- PATH 6: Finalizing ---")
# We remove genes without any mutations (columns with only zeros), but keep the target column
genes_only = matrix.drop('target', axis=1)
active_genes = genes_only.columns[(genes_only != 0).any(axis=0)]
matrix = matrix[list(active_genes) + ['target']]

print(f"Writing matrix with dimensions {matrix.shape}...")
matrix.to_csv(OUTPUT_FILE)
print(f"SUCCESS! File {OUTPUT_FILE} is ready.")