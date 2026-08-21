#!/usr/bin/env bash
# rename_genus_folders.sh
#
# Renames folders from taxon ID to genus name using a CSV lookup table.
# CSV format: Genus, TaxonID, <metadata1>, <metadata2>
#
# Usage:
#   bash rename_genus_folders.sh \
#       --genus-dir genus_genomes \
#       --table     genus_taxids.csv

set -euo pipefail

GENUS_DIR=""
TABLE=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genus-dir) GENUS_DIR="$2"; shift 2 ;;
        --table)     TABLE="$2";     shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

if [[ -z "$GENUS_DIR" || -z "$TABLE" ]]; then
    echo "[ERROR] --genus-dir and --table are required."
    exit 1
fi

if [[ ! -d "$GENUS_DIR" ]]; then
    echo "[ERROR] Directory not found: $GENUS_DIR"
    exit 1
fi

if [[ ! -f "$TABLE" ]]; then
    echo "[ERROR] Table not found: $TABLE"
    exit 1
fi

# ─────────────────────────────────────────────
# READ CSV AND RENAME FOLDERS
# IFS=','      sets comma as delimiter
# $1 = genus name (col 1)
# $2 = taxon ID  (col 2)
# tail -n +2   skips the header line
# ─────────────────────────────────────────────

renamed=0
warned=0
skipped=0

while IFS=',' read -r genus_name taxon_id remainder; do

    # Skip invalid/placeholder genus names
    if [[ "$genus_name" == "." || -z "$genus_name" ]]; then
        echo "[SKIP] Invalid genus name '$genus_name' for taxon $taxon_id — skipping."
        (( skipped++ )) || true
        continue
    fi

    src="$GENUS_DIR/${taxon_id}_16S"   # folders are named <taxonID>_16S
    dst="$GENUS_DIR/$genus_name"

    if [[ ! -d "$src" ]]; then
        echo "[WARN] No folder found for taxon ID $taxon_id ($genus_name) — skipping."
        (( warned++ )) || true
        continue
    fi

    if [[ -d "$dst" ]]; then
        echo "[WARN] Destination already exists: $dst — skipping."
        (( warned++ )) || true
        continue
    fi

    mv "$src" "$dst"
    echo "[OK] ${taxon_id}_16S → $genus_name"
    (( renamed++ )) || true

done < <(tail -n +2 "$TABLE")

# --- summmary prints -------------
echo "[DONE] Folder renaming complete."
echo "       Renamed:  $renamed"
echo "       Warned:   $warned"
echo "       Skipped:  $skipped"
