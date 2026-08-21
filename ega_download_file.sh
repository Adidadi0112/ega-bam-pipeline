#!/bin/bash

INPUT_FILE="/home/adam/projects/EGAD00001005237/to_download.txt"

if [[ ! -f "$INPUT_FILE" ]]; then
	echo "Failure, file "$INPUT_FILE" don't exist"
	exit 1
fi

echo "Starting downloading files..."

while IFS= read -r id || [[ -n "$id" ]]; do
	clean_id=$(echo "$id" | tr -d ',' | xargs)

	if [[ -n "$clean_id" ]]; then
		echo "---------------------------------------------------"
		echo "Downloading: $clean_id"
		pyega3 -c 8 -cf ~/programs/pyega3_credentials.json fetch "$clean_id" --output-dir ~/projects/EGAD00001005237
	fi
done < "$INPUT_FILE"

echo "---------------------------------------------------"
echo "Download succesful"
