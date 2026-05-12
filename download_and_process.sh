#!/bin/bash

# Plik z linkami (np. wygenerowany wcześniej jpt_mapped_only.txt)
LINKS_FILE="jpt_mapped_only.txt"

while read -r url; do
    filename=$(basename "$url")
    
    echo "################################################"
    echo "Pobieranie: $filename"
    echo "################################################"
    
    # Pobieranie z maksymalną prędkością (aria2c)
    aria2c -x 16 -s 16 -k 1M "$url"
    
    # Uruchomienie analizy
    bash variant_calling.sh -f "$filename"
    
    # Usunięcie surowego BAMa po przetworzeniu
    echo "Usuwanie pliku źródłowego: $filename"
    rm "$filename"
    
done < "$LINKS_FILE"
