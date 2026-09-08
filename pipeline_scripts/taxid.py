from Bio import Entrez
import time
import pandas as pd

# ── Setup ─────────────────────────────────────────────────────────────────────
Entrez.email = "ore.solanke25@imperial.ac.uk"
# Entrez.api_key = "YOUR_API_KEY_HERE"  # Recommended for large lists

# ── Your DADA2 genera list ────────────────────────────────────────────────────
genera = []
with open('genera_list.txt', 'r') as f:
    for row in f:
        genera.append(row.strip())

def get_bacterial_taxid(genus_name):
    """
    Query NCBI Taxonomy for the TaxID of a genus, restricted to Bacteria only.
    Returns a dict with TaxID, status, and hit count for transparency.
    """
    try:
        # 🔑 Key: AND Bacteria[subtree] restricts results to bacterial taxa only
        query = f"{genus_name}[Genus] AND Bacteria[subtree]"
        handle = Entrez.esearch(db="taxonomy", term=query)
        record = Entrez.read(handle)
        handle.close()

        hit_count = int(record["Count"])

        if hit_count == 0:
            #print(f"  [NOT FOUND]  {genus_name} — no bacterial match")
            return {"Genus": genus_name, "TaxID": None, "Status": "not_found", "Hits": 0}

        elif hit_count == 1:
            #print(f"  [OK]         {genus_name} → TaxID: {record['IdList'][0]}")
            return {"Genus": genus_name, "TaxID": record["IdList"][0], "Status": "ok", "Hits": 1}

        else:
            # Multiple hits even within Bacteria — flag for manual review
            #print(f"  [MULTI-HIT]  {genus_name} → {hit_count} hits: {record['IdList']}")
            return {"Genus": genus_name, "TaxID": record["IdList"][0], "Status": "multi_hit", "Hits": hit_count}

    except Exception as e:
        #print(f"  [ERROR]      {genus_name}: {e}")
        return {"Genus": genus_name, "TaxID": None, "Status": "error", "Hits": None}


# ── Main loop ─────────────────────────────────────────────────────────────────
results = []

for genus in genera:
    result = get_bacterial_taxid(genus)
    results.append(result)
    time.sleep(0.4)  # Respect NCBI rate limits

# ── Output ────────────────────────────────────────────────────────────────────
df = pd.DataFrame(results)

ids = df[
	(df["Status"]=="ok") &
	(df["Hits"]==1) &
	(df["TaxID"].notna())][["TaxID"]].copy()


# Save outputs
df.to_csv("genus_taxids.csv", index=False)

ids.to_csv("genera_ncbi_ids.txt", header=False, index=False)
