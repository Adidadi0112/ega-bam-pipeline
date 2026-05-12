import pandas as pd
import gseapy as gp
import matplotlib.pyplot as plt

# 1. Load list of top genes from previous analysis (TODO: replace with actual gene names from model feature importance)

top_genes = [
    'MST1', 'TSPAN4', 'LPP', 'HDAC9', 'CSMD1', 'CD80', 'DNMT3B', 
    'CCDC93', 'ASAP2', 'NCOR2', 'IL1R2', 'VWF', 'MYO16', 'ITGA4', 
    'DOCK5', 'PCDH15', 'MYO3B', 'MAP4K3', 'TRPM3', 'COL28A1'
]

print(f"Starting functional analysis for {len(top_genes)} genes...")

# 2. Running Enrichment in GO and KEGG databases
try:
    enr = gp.enrichr(gene_list=top_genes,
                     gene_sets=['GO_Biological_Process_2023', 'KEGG_2021_Human'],
                     organism='human',
                     outdir=None)
    
    # 3. Fetching results and filtering for significance
    results = enr.results
    
    sig_results = results[results['Adjusted P-value'] < 0.05]
    
    if sig_results.empty:
        print("No significant results (P-value < 0.05). Showing top terms (non-significant):")
        sig_results = results.head(10)
    else:
        print(f"Found {len(sig_results)} significant pathways!")

    # 4. Visualization of Top 10 pathways
    from gseapy import dotplot
    
    ax = dotplot(enr.results,
                 column="Adjusted P-value",
                 x='Gene_set', # Grouping by database (GO vs KEGG)
                 size=10,
                 top_term=10,
                 figsize=(10, 8),
                 title="GO & KEGG Enrichment Analysis",
                 xticklabels_rot=45,
                 show_ring=True,
                 marker='o')
    
    plt.savefig('pathway_enrichment.png', bbox_inches='tight')
    print("Plot saved as pathway_enrichment.png")
    
    # Saving the table to CSV
    sig_results.to_csv('enrichment_results.csv', index=False)
    print("Full results table saved in enrichment_results.csv")

except Exception as e:
    print(f"Error while connecting to Enrichr database: {e}")
    print("Make sure you have an internet connection.")