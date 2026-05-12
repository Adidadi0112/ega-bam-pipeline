# 1. Znajdź wszystkie pliki .vcf.gz w obu folderach i zapisz do listy
find results control_results -name "*.vcf.gz" > all_vcfs.list

# 2. Połącz wszystko w jeden plik
# Flaga -m none pozwala na łączenie próbek o różnych nazwach (nawet jeśli miałyby te same ID wewnątrz VCF)
bcftools merge -l all_vcfs.list -Oz -o all_samples_merged_raw.vcf.gz

echo "Łączenie zakończone. Plik wynikowy: all_samples_merged_raw.vcf.gz"
