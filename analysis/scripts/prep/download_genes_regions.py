import requests
import time

def create_bed_from_symbols(input_file, output_file):
    server = "https://grch37.rest.ensembl.org"
    headers = {"Content-Type": "application/json"}
    
    # Load gene symbols from the input file
    with open(input_file, 'r') as f:
        # Clean up the content to get a list of gene symbols
        content = f.read().replace(',', ' ')
        genes = [g.strip() for g in content.split() if g.strip()]

    print(f"Starting download for {len(genes)} genes...")

    with open(output_file, 'w') as bed_file:
        for symbol in genes:
            # Endpoint for gene lookup by symbol
            ext = f"/lookup/symbol/homo_sapiens/{symbol}?"
            
            try:
                r = requests.get(server + ext, headers=headers)
                
                # Handle rate limiting
                if r.status_code == 429:
                    retry_after = float(r.headers.get("Retry-After", 1))
                    time.sleep(retry_after)
                    r = requests.get(server + ext, headers=headers)

                if not r.ok:
                    print(f"Skipped: {symbol} (not found in GRCh37)")
                    continue

                data = r.json()
                
                # Fetch: chrom, start, end
                # BED format: chrom (without 'chr'), start-1, end
                chrom = data['seq_region_name']
                start = data['start'] - 1 # BED is 0-based at the start
                end = data['end']
                
                bed_file.write(f"{chrom}\t{start}\t{end}\t{symbol}\n")
                print(f"OK: {symbol} -> {chrom}:{start}-{end}")

            except Exception as e:
                print(f"Error with gene {symbol}: {e}")
            
            # Small delay to avoid overwhelming the API
            time.sleep(0.1)

    print(f"\nSUCCESS! File saved as: {output_file}")

if __name__ == "__main__":
    create_bed_from_symbols("genes_uc.txt", "uc_genes.bed")