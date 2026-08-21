
#!/usr/bin/env bash
# run_mash_pipeline.sh
#
# Usage:
#   bash run_mash_pipeline.sh \
#       --genus-dir genus_genomes \
#       --asv-fasta genus_ASVs.fasta \
#       --output    results/mash_results.tsv

set -euo pipefail

# --------- checking passed in arguments ------------

GENUS_DIR=""
ASV_FASTA=""
OUTPUT=""

while [[ $# -gt 0 ]]; do
    case "$1" in
        --genus-dir)  GENUS_DIR="$2";  shift 2 ;;
        --asv-fasta)  ASV_FASTA="$2";  shift 2 ;;
        --output)     OUTPUT="$2";     shift 2 ;;
        *) echo "[ERROR] Unknown argument: $1"; exit 1 ;;
    esac
done

# ------- error catching --------------------

if [[ -z "$GENUS_DIR" || -z "$ASV_FASTA" || -z "$OUTPUT" ]]; then
    echo "[ERROR] --genus-dir, --asv-fasta, and --output are all required."
    exit 1
fi

if ! command -v mash &>/dev/null; then
    echo "[ERROR] mash not found in PATH."
    exit 1
fi

if [[ ! -d "$GENUS_DIR" ]]; then
    echo "[ERROR] Genus directory not found: $GENUS_DIR"
    exit 1
fi

if [[ ! -f "$ASV_FASTA" ]]; then
    echo "[ERROR] ASV fasta file not found: $ASV_FASTA"
    exit 1
fi

# --------- set up output -----------------------

THREADS=${NCPUS:-${PBS_NP:-$(nproc)}}

mkdir -p "$(dirname "$OUTPUT")"

echo -e "Genus\tASV_ID\tReference_ID\tMash_Distance\tP_Value\tMatching_Hashes" > "$OUTPUT"

TEMP_ASV=$(mktemp /tmp/asvs_XXXXXX.fasta)
TEMP_DIST=$(mktemp /tmp/mash_dist_XXXXXX.txt)
trap 'rm -f "$TEMP_ASV" "$TEMP_DIST"' EXIT

# ----- looping through genus folders to do mash dist calc -----

for genus_folder in "$GENUS_DIR"/*/; do

    [[ -d "$genus_folder" ]] || continue
    genus_name=$(basename "$genus_folder")

    # ── STEP 1: Sketch all non-empty reference 16S for this genus ──
    if ! ls "$genus_folder"*_16S.fna &>/dev/null; then
        echo "[WARN] No *_16S.fna files in $genus_name — skipping."
        continue
    fi

    if [[ ! -f "$genus_folder/reference.msh" ]]; then

        non_empty_fna=()
        for fna in "$genus_folder"*_16S.fna; do
            if [[ -s "$fna" ]]; then
                non_empty_fna+=("$(basename "$fna")")
            else
                echo "[WARN] Empty file skipped from sketch: $fna"
            fi
        done

        if [[ "${#non_empty_fna[@]}" -eq 0 ]]; then
            echo "[WARN] All *_16S.fna files are empty in $genus_name — skipping."
            continue
        fi

        echo "[CMD] cd $genus_folder && mash sketch -o reference ${non_empty_fna[*]}"
        (
            cd "$genus_folder"
            mash sketch -o reference "${non_empty_fna[@]}"
        )
    else
        echo "[INFO] Sketch already exists for $genus_name — skipping."
    fi

    # ── STEP 2: Extract ASVs matching this genus ──
    awk -v genus="$genus_name" '
        /^>/ { keep = (tolower($0) ~ tolower(genus)) }
        keep  { print }
    ' "$ASV_FASTA" > "$TEMP_ASV"

    asv_count=$(grep -c '^>' "$TEMP_ASV" || true)
    echo "[DEBUG] ASVs matched for $genus_name: $asv_count"

    if [[ ! -s "$TEMP_ASV" ]]; then
        echo "[WARN] No ASVs matched genus '$genus_name' in $ASV_FASTA — skipping."
        continue
    fi

    # ── STEP 3: mash dist — one ASV at a time so Query_ID is the ASV header ──
    while IFS= read -r header; do
        asv_id="${header#>}"

        # extract just this one ASV sequence into a single record temp file
        awk -v id="$asv_id" '
            /^>/ { keep = (substr($0,2) == id) }
            keep { print }
        ' "$TEMP_ASV" > "$TEMP_DIST"

        if [[ ! -s "$TEMP_DIST" ]]; then
            echo "[WARN] Could not extract sequence for $asv_id — skipping."
            continue
        fi

        echo "[CMD] mash dist -p $THREADS $genus_folder/reference.msh $TEMP_DIST"
        mash dist \
            -p "$THREADS" \
            "$genus_folder/reference.msh" \
            "$TEMP_DIST" \
        | while IFS=$'\t' read -r ref_id query_id distance pvalue hashes; do
            printf "%s\t%s\t%s\t%s\t%s\t%s\n" \
                "$genus_name" \
                "$asv_id" \
                "$ref_id" \
                "$distance" \
                "$pvalue" \
                "$hashes"
          done >> "$OUTPUT"

    done < <(grep '^>' "$TEMP_ASV")

done

echo "[DONE] Results written to: $OUTPUT"
