# ============================================================
# GapSeq SBML .xml parser
# Outputs:
#   {model}_reactions.csv
#   {model}_metabolites.csv
#   {model}_genes.csv
#   {model}_pathways.csv
#   {model}_summary.csv
#   {model}_transporters.csv
#   all_models_summary.csv
#   all_transporters_summary.csv 
# ============================================================

import xml.etree.ElementTree as ET
import pandas as pd
import os
import sys

NS = {
    "sbml"  : "http://www.sbml.org/sbml/level3/version2/core",
    "fbc"   : "http://www.sbml.org/sbml/level3/version1/fbc/version2",
    "groups": "http://www.sbml.org/sbml/level3/version1/groups/version1",
    "rdf"   : "http://www.w3.org/1999/02/22-rdf-syntax-ns#",
    "bqbiol": "http://biomodels.net/biology-qualifiers/",
}

FBC = "http://www.sbml.org/sbml/level3/version1/fbc/version2"
GRP = "http://www.sbml.org/sbml/level3/version1/groups/version1"

# ============================================================
# HELPERS
# ============================================================

def get_annotations(element):
    annotations = {}
    for qualifier in ["bqbiol:is", "bqbiol:isHomologTo"]:
        bag = element.find(f".//{qualifier}/rdf:Bag", NS)
        if bag is not None:
            for li in bag:
                uri = li.attrib.get(
                    "{http://www.w3.org/1999/02/22-rdf-syntax-ns#}resource", ""
                )
                if "identifiers.org/" not in uri:
                    continue
                clean = uri.replace("http://identifiers.org/", "")
                if "/" in clean:
                    db, entry = clean.split("/", 1)
                elif ":" in clean:
                    db, entry = clean.split(":", 1)
                else:
                    continue
                annotations.setdefault(db, []).append(entry)
    return annotations


def fmt(annotation_dict, key):
    return ", ".join(annotation_dict.get(key, []))


def get_gpr(rxn_el):
    gpa = rxn_el.find("fbc:geneProductAssociation", NS)
    if gpa is None:
        return "No genes", []

    def extract_refs(parent):
        return [
            ref.attrib.get(f"{{{FBC}}}geneProduct", "").replace("G_", "")
            for ref in parent.findall("fbc:geneProductRef", NS)
        ]

    single = gpa.find("fbc:geneProductRef", NS)
    if single is not None:
        gene = single.attrib.get(f"{{{FBC}}}geneProduct", "").replace("G_", "")
        return gene, [gene]

    or_el = gpa.find("fbc:or", NS)
    if or_el is not None:
        genes = extract_refs(or_el)
        return " OR ".join(genes), genes

    and_el = gpa.find("fbc:and", NS)
    if and_el is not None:
        genes = extract_refs(and_el)
        return " AND ".join(genes), genes

    return "No genes", []


def resolve_bound(param_id, param_map):
    return param_map.get(param_id, param_id)


def load_reaction_scores(tbl_file):
    if not os.path.exists(tbl_file):
        return {}

    df = pd.read_csv(tbl_file, sep="\t", skiprows=3)
    df = df[df["bitscore"].notna()]

    best = (
        df.sort_values("bitscore", ascending=False)
          .groupby("rxn", as_index=False)
          .first()
    )

    score_map = {}
    for _, row in best.iterrows():
        score_map[row["rxn"]] = {
            "pident"  : row.get("pident",   ""),
            "evalue"  : row.get("evalue",   ""),
            "bitscore": row.get("bitscore", ""),
            "qcovs"   : row.get("qcovs",    ""),
            "status"  : row.get("status",   ""),
        }

    return score_map


# ============================================================
# PARSERS
# ============================================================

def parse_parameters(model):
    param_map = {}
    for param in model.findall(".//sbml:listOfParameters/sbml:parameter", NS):
        pid = param.attrib.get("id", "")
        val = param.attrib.get("value", "")
        try:
            param_map[pid] = float(val)
        except ValueError:
            param_map[pid] = val
    return param_map

def parse_transporter_tbl(tbl_file, model_id):
    if not os.path.exists(tbl_file):
        return pd.DataFrame()

    df = pd.read_csv(tbl_file, sep="\t", skiprows=3)
    df.insert(0, "Model ID", model_id)

    return df


def parse_reactions(model, param_map, model_id, score_map):
    rows = []
    for rxn in model.findall(".//sbml:listOfReactions/sbml:reaction", NS):

        rxn_id             = rxn.attrib.get("id", "")
        ann                = get_annotations(rxn)
        gpr_str, gene_list = get_gpr(rxn)
        lb_ref             = rxn.attrib.get(f"{{{FBC}}}lowerFluxBound", "")
        ub_ref             = rxn.attrib.get(f"{{{FBC}}}upperFluxBound", "")

        reactant_refs = rxn.findall(".//sbml:listOfReactants/sbml:speciesReference", NS)
        reactants     = ", ".join(sr.attrib.get("species",      "") for sr in reactant_refs)
        reactant_stoi = ", ".join(sr.attrib.get("stoichiometry", "1") for sr in reactant_refs)

        product_refs  = rxn.findall(".//sbml:listOfProducts/sbml:speciesReference", NS)
        products      = ", ".join(sr.attrib.get("species",      "") for sr in product_refs)
        product_stoi  = ", ".join(sr.attrib.get("stoichiometry", "1") for sr in product_refs)

        seed_ids = ann.get("seed.reaction", [])
        scores   = {}
        for sid in seed_ids:
            if sid in score_map:
                scores = score_map[sid]
                break

        rows.append({
            "Model ID"         : model_id,
            "Reaction ID"      : rxn_id,
            "Reaction Name"    : rxn.attrib.get("name",       ""),
            "SBO Term"         : rxn.attrib.get("sboTerm",    ""),
            "Meta ID"          : rxn.attrib.get("metaid",     ""),
            "Reversible"       : rxn.attrib.get("reversible", ""),
            "Lower Bound"      : resolve_bound(lb_ref, param_map),
            "Upper Bound"      : resolve_bound(ub_ref, param_map),
            "Reactants"        : reactants,
            "Reactant Stoich"  : reactant_stoi,
            "Products"         : products,
            "Product Stoich"   : product_stoi,
            "GPR"              : gpr_str,
            "Num Genes"        : len(gene_list),
            "EC"               : fmt(ann, "ec-code"),
            "KEGG Reaction"    : fmt(ann, "kegg.reaction"),
            "SEED Reaction"    : fmt(ann, "seed.reaction"),
            "BiGG Reaction"    : fmt(ann, "bigg.reaction"),
            "MetaNetX Reaction": fmt(ann, "metanetx.reaction"),
            "BioCyc"           : fmt(ann, "biocyc"),
            "pident"           : scores.get("pident",   ""),
            "evalue"           : scores.get("evalue",   ""),
            "bitscore"         : scores.get("bitscore", ""),
            "qcovs"            : scores.get("qcovs",    ""),
            "blast_status"     : scores.get("status",   ""),
        })

    return pd.DataFrame(rows)

def parse_metabolites(model, model_id):
    rows = []
    for sp in model.findall(".//sbml:listOfSpecies/sbml:species", NS):
        ann = get_annotations(sp)
        rows.append({
            "Model ID"         : model_id,
            "Metabolite ID"    : sp.attrib.get("id",   ""),
            "Name"             : sp.attrib.get("name", ""),
            "Meta ID"          : sp.attrib.get("metaid",  ""),
            "SBO Term"         : sp.attrib.get("sboTerm", ""),
            "Compartment"      : sp.attrib.get("compartment", ""),
            "Charge"           : sp.attrib.get(f"{{{FBC}}}charge",          ""),
            "Formula"          : sp.attrib.get(f"{{{FBC}}}chemicalFormula", ""),
            "Boundary"         : sp.attrib.get("boundaryCondition", ""),
            "Constant"         : sp.attrib.get("constant",          ""),
            "KEGG Compound"    : fmt(ann, "kegg.compound"),
            "SEED Compound"    : fmt(ann, "seed.compound"),
            "ChEBI"            : fmt(ann, "CHEBI"),
            "BiGG Metabolite"  : fmt(ann, "bigg.metabolite"),
            "MetaNetX Chemical": fmt(ann, "metanetx.chemical"),
            "InChIKey"         : fmt(ann, "inchikey"),
            "HMDB"             : fmt(ann, "hmdb"),
            "Reactome"         : fmt(ann, "reactome"),
            "BioCyc"           : fmt(ann, "biocyc"),
        })

    return pd.DataFrame(rows)


def parse_genes(model, model_id):
    rows = []
    for gp in model.findall(".//fbc:listOfGeneProducts/fbc:geneProduct", NS):
        ann = get_annotations(gp)
        rows.append({
            "Model ID" : model_id,
            "Gene ID"  : gp.attrib.get(f"{{{FBC}}}id",    "").replace("G_", ""),
            "Label"    : gp.attrib.get(f"{{{FBC}}}label", ""),
            "Name"     : gp.attrib.get(f"{{{FBC}}}name",  ""),
            "Meta ID"  : gp.attrib.get("metaid",  ""),
            "SBO Term" : gp.attrib.get("sboTerm", ""),
            "UniRef"   : fmt(ann, "uniref"),
        })

    return pd.DataFrame(rows)


def parse_pathways(model, model_id):
    rows           = []
    list_of_groups = model.find("groups:listOfGroups", NS)

    if list_of_groups is None:
        return pd.DataFrame()

    for grp in list_of_groups:
        grp_id   = grp.attrib.get(f"{{{GRP}}}id",   "")
        grp_name = grp.attrib.get(f"{{{GRP}}}name", "")
        members  = grp.findall("groups:listOfMembers/groups:member", NS)

        for m in members:
            rows.append({
                "Model ID"    : model_id,
                "Pathway ID"  : grp_id,
                "Pathway Name": grp_name,
                "Reaction ID" : m.attrib.get(f"{{{GRP}}}idRef", ""),
                "Pathway Size": len(members),
            })

    return pd.DataFrame(rows)


def build_summary(rxn_df, met_df, gene_df, pathway_df):

    met_lookup  = met_df.set_index("Metabolite ID")["Name"].to_dict()
    gene_lookup = gene_df.set_index("Gene ID")["UniRef"].to_dict()

    def resolve_met_names(stoich_str):
        if not stoich_str:
            return stoich_str
        parts = []
        for item in stoich_str.split(", "):
            met_id = item.split(" (")[0]
            stoich = item.split("(")[-1].rstrip(")")
            name   = met_lookup.get(met_id, met_id)
            parts.append(f"{name} ({stoich})")
        return ", ".join(parts)

    def resolve_gene_uniref(gpr_str):
        if gpr_str in ("No genes", ""):
            return ""
        genes  = [g.strip() for g in gpr_str.replace(" OR ", " ").replace(" AND ", " ").split()]
        uniref = [gene_lookup.get(g, "") for g in genes if gene_lookup.get(g, "")]
        return "; ".join(uniref)

    if not pathway_df.empty:
        pwy_grouped = pathway_df.groupby("Reaction ID").agg({
            "Pathway ID"  : "; ".join,
            "Pathway Name": "; ".join,
            "Pathway Size": lambda x: "; ".join(str(v) for v in x),
        }).reset_index()
    else:
        pwy_grouped = pd.DataFrame(
            columns=["Reaction ID", "Pathway ID", "Pathway Name", "Pathway Size"]
        )

    summary                = rxn_df.copy()
    summary["Reactants"]   = summary["Reactants"].apply(resolve_met_names)
    summary["Products"]    = summary["Products"].apply(resolve_met_names)
    summary["Gene UniRef"] = summary["GPR"].apply(resolve_gene_uniref)
    summary                = summary.merge(pwy_grouped, on="Reaction ID", how="left")

    return summary


# ============================================================
# MAIN
# ============================================================

def summarize_model(xml_file):

    if not os.path.exists(xml_file):
        return None, None

    try:
        tree  = ET.parse(xml_file)
        root  = tree.getroot()
        model = root.find("sbml:model", NS)
        if model is None:
            return None, None
    except ET.ParseError:
        return None, None

    model_id  = model.attrib.get("id", os.path.basename(xml_file).replace("-draft.xml", ""))
    param_map = parse_parameters(model)

    tbl_file        = xml_file.replace("-draft.xml", "-all-Reactions.tbl")
    score_map       = load_reaction_scores(tbl_file)

    rxn_df     = parse_reactions(model,   param_map, model_id, score_map)
    met_df     = parse_metabolites(model, model_id)
    gene_df    = parse_genes(model,       model_id)
    pathway_df = parse_pathways(model,    model_id)
    summary_df = build_summary(rxn_df, met_df, gene_df, pathway_df)

    transporter_tbl = xml_file.replace("-draft.xml", "-Transporter.tbl")
    transporter_df  = parse_transporter_tbl(transporter_tbl, model_id)

    base = xml_file.replace("-draft.xml", "")
    rxn_df.to_csv(    f"{base}_reactions.csv",   index=False)
    met_df.to_csv(    f"{base}_metabolites.csv", index=False)
    gene_df.to_csv(   f"{base}_genes.csv",       index=False)
    pathway_df.to_csv(f"{base}_pathways.csv",    index=False)
    summary_df.to_csv(f"{base}_summary.csv",     index=False)
    if not transporter_df.empty:
        transporter_df.to_csv(f"{base}_transporters.csv", index=False)

    return summary_df, transporter_df


if __name__ == "__main__":

    model_list = sys.argv[1] if len(sys.argv) > 1 else "draft_list.txt"

    if not os.path.exists(model_list):
        sys.exit(1)

    with open(model_list, "r") as f:
        files = [row.strip() for row in f if row.strip()]

    all_summaries    = []
    all_transporters = []

    for xml_file in files:
        summary_df, transporter_df = summarize_model(xml_file)
        if summary_df is not None:
            all_summaries.append(summary_df)
        if transporter_df is not None and not transporter_df.empty:
            all_transporters.append(transporter_df)

    if all_summaries:
        pd.concat(all_summaries, ignore_index=True).to_csv(
            "all_models_summary.csv", index=False
        )
    if all_transporters:
        pd.concat(all_transporters, ignore_index=True).to_csv(
            "all_transporters_summary.csv", index=False
        )
