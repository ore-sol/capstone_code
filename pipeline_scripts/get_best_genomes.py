#!/usr/bin/env python3
# get_best_genomes.py
#
# 1. For each genus, finds the single best reference .fna across all ASV comparisons.
# 2. Strips the _genomic_16S.fna suffix to get the genome accession.
# 3. Looks up the taxon ID for that genus from the CSV table.
# 4. Searches genome_dir/taxon_id/ncbi_dataset/data/accession/ for the matching .fna
# 5. Copies the matched file to an output directory.
#
# Usage:
#   python get_best_genomes.py \
#       --input      results/mash_results.tsv \
#       --output     results/mash_best_hits.tsv \
#       --taxon-csv  taxon_table.csv \
#       --genome-dir ncbi_cli/unzipped \
#       --copy-to    results/best_genomes

import argparse
import csv
import shutil
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input",      required=True, help="Full mash results TSV")
    parser.add_argument("--output",     required=True, help="Filtered best hits TSV")
    parser.add_argument("--taxon-csv",  required=True, help="CSV with Genus and TaxonID columns")
    parser.add_argument("--genome-dir", required=True, help="Parent unzipped directory containing taxon ID subfolders")
    parser.add_argument("--copy-to",    required=True, help="Directory to copy matched genome files into")
    return parser.parse_args()


def load_taxon_map(taxon_csv: str) -> dict:
    """
    Reads the CSV and returns a dict of {genus: taxon_id}.
    CSV format: Genus, TaxonID, <metadata1>, <metadata2>
    One genus per taxon ID.
    """
    taxon_map = {}
    with open(taxon_csv, "r") as f:
        reader = csv.reader(f)
        next(reader)                   # skip header
        for row in reader:
            if len(row) < 2:
                continue
            genus, taxon_id = row[0].strip(), row[1].strip()
            taxon_map[genus] = taxon_id
    return taxon_map


def find_best_per_genus(input_tsv: str) -> dict:
    """
    For each genus, finds the row with the lowest Mash distance,
    then lowest p-value, then first occurrence.
    Returns dict of {genus: row}
    """
    best = {}
    with open(input_tsv, "r") as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            genus    = row["Genus"]
            distance = float(row["Mash_Distance"])
            pvalue   = float(row["P_Value"])

            if genus not in best:
                best[genus] = row
            else:
                current_dist = float(best[genus]["Mash_Distance"])
                current_pval = float(best[genus]["P_Value"])

                if (distance, pvalue) < (current_dist, current_pval):
                    best[genus] = row
    return best


def strip_suffix(reference_id: str) -> str:
    
    """
    Strips everything after the bare accession number.
    e.g. GCF_004124235.1_ASM412423v1_genomic_16S.fna
       → GCF_004124235.1   (folder name)
    keeping _ASM*_genomic.fna as the glob pattern for the file inside.
    """
    # accession is always the first two underscore-separated parts: GCF_XXXXXXXXX.X
    parts = reference_id.split("_")
    return f"{parts[0]}_{parts[1]}"
    

def find_genome_file(accession: str, taxon_id: str, genome_dir: str) -> Path | None:
    """
    Looks inside:
      genome_dir/taxon_id/ncbi_dataset/data/accession/
    for a file matching accession*_genomic.fna.
    There is exactly one match per accession folder.
    """
    accession_folder = Path(genome_dir) / taxon_id / "ncbi_dataset" / "data" / accession

    if not accession_folder.is_dir():
        print(f"[WARN] Accession folder not found: {accession_folder}")
        return None

    matches = list(accession_folder.glob(f"{accession}*_genomic.fna"))

    if len(matches) == 0:
        print(f"[WARN] No .fna file found in: {accession_folder}")
        return None
    if len(matches) > 1:
        print(f"[WARN] Multiple .fna files found in {accession_folder} — using first: {matches[0].name}")

    return matches[0]


def main():
    args = parse_args()

    copy_to = Path(args.copy_to)
    copy_to.mkdir(parents=True, exist_ok=True)

    # ── STEP 1: find best reference per genus ────────────────
    best_per_genus = find_best_per_genus(args.input)

    # ── STEP 2: load genus → taxon ID map ────────────────────
    taxon_map = load_taxon_map(args.taxon_csv)

    # ── STEP 3: write best hits TSV + find and copy files ────
    with open(args.output, "w", newline="") as f:
        fieldnames = ["Genus", "TaxonID", "Accession", "Reference_ID", "ASV_ID",
                      "Mash_Distance", "P_Value", "Matching_Hashes", "File_Found"]
        writer = csv.DictWriter(f, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()

        for genus, row in best_per_genus.items():

            accession = strip_suffix(row["Reference_ID"])

            # look up taxon ID for this genus
            taxon_id = taxon_map.get(genus)
            if not taxon_id:
                print(f"[WARN] No taxon ID found for genus '{genus}' in {args.taxon_csv}")
                taxon_id = "unknown"

            # search genome_dir/taxon_id/ncbi_dataset/data/accession/ for .fna
            matched_file = find_genome_file(accession, taxon_id, args.genome_dir)

            if matched_file:
                dest = copy_to / matched_file.name
                shutil.copy2(matched_file, dest)
                print(f"[OK]   {genus} → {matched_file.name}")
                file_found = matched_file.name
            else:
                file_found = "NOT_FOUND"

            writer.writerow({
                "Genus":           genus,
                "TaxonID":         taxon_id,
                "Accession":       accession,
                "Reference_ID":    row["Reference_ID"],
                "ASV_ID":          row["ASV_ID"],
                "Mash_Distance":   row["Mash_Distance"],
                "P_Value":         row["P_Value"],
                "Matching_Hashes": row["Matching_Hashes"],
                "File_Found":      file_found
            })

    print(f"\n[DONE] Best hits written to: {args.output}")
    print(f"[INFO] Genera processed: {len(best_per_genus)}")
    print(f"[INFO] Files copied to:  {copy_to}")


if __name__ == "__main__":
    main()
