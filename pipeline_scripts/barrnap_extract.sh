#!/usr/bin/env bash
#barrnap_extract.sh 
#
#
# extracts 16S from each .gff and .fasta pair from barrnap output
#
# Usage:
#   bash barrnap_extract.sh \
#       --genus-dir genus_genomes \

for dir in */; do
    echo "Starting $dir"
    for gff in "$dir"*.gff; do

        fasta="${gff%barrnap.gff}rRNA.fasta"
        output="${gff%barrnap.gff}16S.fna"

        if [[ ! -f "$fasta" ]]; then
            echo "[WARN] No matching fasta for $gff — skipping."
            continue
        fi

        # pull 16S coordinate strings directly from GFF headers
        # format matches the barrnap --outseq fasta headers exactly
        grep "16S_rRNA" "$gff" \
            | awk '{print $1":"$4-1"-"$5}' \
            > "${gff}.ids"

        if [[ ! -s "${gff}.ids" ]]; then
            echo "[WARN] No 16S rRNA found in $gff — skipping."
            rm -f "${gff}.ids"
            continue
        fi

        # use awk to extract matching records from the rRNA fasta
        # no samtools or bedtools needed
        awk 'NR==FNR { ids[$1]=1; next }
             /^>/ { keep=0; for (id in ids) if (index($0, id)) keep=1 }
             keep { print }' \
            "${gff}.ids" "$fasta" \
            > "$output"

        if [[ ! -s "$output" ]]; then
            echo "[WARN] No 16S extracted for $gff"
        fi

        rm -f "${gff}.ids"

    done
done

for dir in */;
do 
    echo "Moving 16S from $dir" 
    mkdir "${dir%/}_16S"
    for fna in "$dir"*_16S.fna; do

        if [[ ! -f "$fna" ]]; then
            echo "[WARN] No 16S .fna files found in $dir — skipping."
            continue
        fi

        mv "$fna" "${dir%/}_16S/"

    done
done

