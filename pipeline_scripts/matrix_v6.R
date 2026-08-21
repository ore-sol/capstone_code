# ==================== matrix_tensor_construction.R ====================
# Builds three metabolic interaction matrices from gapseq model outputs:
#   1. Uptake matrix    [species x metabolite] binary
#   2. Export matrix    [species x metabolite] binary  
#   3. Leakage tensor   [species x reactant x product] binary
#
# All three require explicit transporter evidence and bound-confirmed
# directionality. See inline comments for logic at each step.


#set up env
rm(list = ls())  # clear all variables

# load packages ----------------------------------------------

library(dplyr)
library(tidyr)
library(readr)
library(tools)

#make output dir 
if(!file.exists("matrices_usethis")){dir.create("matrices_usethis")}


#set up to read in arguments --------------------------------
args     <- commandArgs(trailingOnly = TRUE)
basename <- tools::file_path_sans_ext(args[1])  # strips .csv if passed with extension
basename <- sub("_genera$", "", basename)        # strips trailing _genera if present

cat(sprintf("[INFO] Running for: %s\n", basename))

cat(sprintf("[INFO] args[1]   = %s\n", args[1]))
cat(sprintf("[INFO] basename  = %s\n", basename))
cat(sprintf("[INFO] CSV path  = %s\n", paste0(basename, "_genera.csv")))
cat(sprintf("[INFO] rxn_files = %s\n", paste(list.files("summaries", pattern="_reactions\\.csv$"), collapse=", ")))
cat(sprintf("[INFO] comm_genera = %s\n", paste(read.csv(paste0(basename, "_genera.csv"))[[1]], collapse=", ")))

# run this command: in folder that has ls *_genera.csv | xargs -I{} Rscript matrix_v6.R {}

# load files ----------------------------------------------
rxn_files  <- list.files("../summaries", pattern = "_reactions\\.csv$",    full.names = TRUE)
trp_files  <- list.files("../summaries", pattern = "_transporters\\.csv$", full.names = TRUE)
metb_files  <- list.files("../summaries", pattern = "_metabolites\\.csv$",  full.names = TRUE)

# Extract genome accession from filename — strip path and suffix
get_accession <- function(f, suffix) {
  gsub(suffix, "", basename(f))
}

rxn_files <- get_accession(rxn_files, "_reactions.csv")
trp_files <- get_accession(trp_files, "_transporters.csv")
metb_files <- get_accession(metb_files, "_metabolites.csv")

# Filter to community genera list
comm_genera <- read.csv(paste0(basename, "_genera.csv"))[[1]]

rxn_genera  <- rxn_files[rxn_files %in% comm_genera]
trp_genera  <- trp_files[trp_files %in% comm_genera]
metb_genera <- metb_files[metb_files %in% comm_genera]

# Write genus list for downstream reference 
writeLines(rxn_genera, file.path("matrices_usethis", paste0(basename, "_genera_wIDS.txt")))

genera <- rxn_genera
n_species <- length(rxn_genera)
cat(sprintf("[INFO] Processing %d species\n", n_species))

# metabolite parsing function -------------------------------------------
# Removes any stoichiometric coefficients written as (n), splits on
# comma, trims whitespace, drops empty strings.

parse_mets <- function(x) {
  if (is.na(x) || x == "") return(character(0))
  x <- gsub("\\s*\\(\\d+\\.?\\d*\\)", "", x)
  trimws(unlist(strsplit(x, ",\\s*")))
}

# 2. Per-genome processing function --------------------------------------
#    Returns a list with:
#      $uptake_mets  : character vector of base metabolite IDs
#      $export_mets  : character vector of base metabolite IDs
#      $tensor_pairs : data.frame with columns reactant, product


EXCLUDE_SBO <- c("SBO:0000629", "SBO:0000627")
process_genome <- function(genus) {
  
  cat(sprintf("[PROCESSING] %s\n", genus))
  
  # File paths fully determined by genus name
  rxn <- read.csv(
    file.path("../summaries", paste0(genus, "_reactions.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  trp <- read.csv(
    file.path("../summaries", paste0(genus, "_transporters.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  met <- read.csv(
    file.path("../summaries", paste0(genus, "_metabolites.csv")),
    stringsAsFactors = FALSE, check.names = FALSE
  )
  
  colnames(rxn) <- trimws(colnames(rxn))
  colnames(trp) <- trimws(colnames(trp))
  colnames(met) <- trimws(colnames(met))
  
  # -- Step A: Identify metabolites with confirmed transporters ----------
  trp_exids    <- unique(trp$exid)
  ex_rxn_ids   <- paste0("R_", trp_exids)
  base_from_ex <- sub("^EX_", "", sub("_e0$", "", trp_exids))
  
  ex_map <- data.frame(
    base_met  = base_from_ex,
    ex_rxn_id = ex_rxn_ids,
    stringsAsFactors = FALSE
  )
  
  # -- Step B: Exchange reaction bounds ----------------------------------
  ex_rxns <- rxn[rxn$`Reaction ID` %in% ex_rxn_ids,
                 c("Reaction ID", "Lower Bound", "Upper Bound")]
  colnames(ex_rxns) <- c("ex_rxn_id", "lb_ex", "ub_ex")
  
  ex_map <- merge(ex_map, ex_rxns, by = "ex_rxn_id", all.x = TRUE)
  
  ex_map$uptake_gate1 <- !is.na(ex_map$lb_ex) & ex_map$lb_ex < 0
  ex_map$export_gate1 <- !is.na(ex_map$ub_ex) & ex_map$ub_ex > 0
  
  # -- Step C: Transport reactions via rea -> SEED Reaction join ---------
  trp_expanded <- trp %>%
    select(exid, rea) %>%
    distinct() %>%
    mutate(seed_id = strsplit(rea, ",\\s*")) %>%
    unnest(seed_id) %>%
    mutate(
      seed_id  = trimws(seed_id),
      base_met = sub("^EX_", "", sub("_e0$", "", exid))
    ) %>%
    select(base_met, seed_id) %>%
    distinct()
  
  rxn_expanded <- rxn %>%
    select(`Reaction ID`, `Lower Bound`, `Upper Bound`,
           `SEED Reaction`, Reactants, Products, `SBO Term`) %>%
    filter(!`SBO Term` %in% EXCLUDE_SBO) %>%
    mutate(seed_id = strsplit(`SEED Reaction`, ",\\s*")) %>%
    unnest(seed_id) %>%
    mutate(seed_id = trimws(seed_id)) %>%
    filter(seed_id != "")
  
  transport_rxns <- inner_join(trp_expanded, rxn_expanded, by = "seed_id") %>%
    select(base_met, `Reaction ID`, `Lower Bound`, `Upper Bound`,
           Reactants, Products) %>%
    distinct()
  # Gate 1
  cat(sprintf("  Gate 1 — transporter metabolites found : %d\n", nrow(ex_map)))
  cat(sprintf("  Gate 1 — uptake_gate1 TRUE             : %d\n", sum(ex_map$uptake_gate1)))
  cat(sprintf("  Gate 1 — export_gate1 TRUE             : %d\n", sum(ex_map$export_gate1)))
  
  # -- Step D: Gate 2 — per-metabolite directionality from bounds --------
  transport_rxns$import_gate2 <- FALSE
  transport_rxns$export_gate2 <- FALSE
  
  for (k in seq_len(nrow(transport_rxns))) {
    bm  <- transport_rxns$base_met[k]
    lb  <- transport_rxns$`Lower Bound`[k]
    ub  <- transport_rxns$`Upper Bound`[k]
    rct <- parse_mets(transport_rxns$Reactants[k])
    prd <- parse_mets(transport_rxns$Products[k])
    
    e0_form <- paste0("M_", bm, "_e0")
    c0_form <- paste0("M_", bm, "_c0")
    
    e0_in_reactants <- e0_form %in% rct
    c0_in_products  <- c0_form %in% prd
    c0_in_reactants <- c0_form %in% rct
    e0_in_products  <- e0_form %in% prd
    
    # Forward: e0 → c0 = import
    if (e0_in_reactants && c0_in_products && ub > 0)
      transport_rxns$import_gate2[k] <- TRUE
    # Reverse of e0 → c0 written reaction = export
    if (e0_in_reactants && c0_in_products && lb < 0)
      transport_rxns$export_gate2[k] <- TRUE
    # Forward: c0 → e0 = export
    if (c0_in_reactants && e0_in_products && ub > 0)
      transport_rxns$export_gate2[k] <- TRUE
    # Reverse of c0 → e0 written reaction = import
    if (c0_in_reactants && e0_in_products && lb < 0)
      transport_rxns$import_gate2[k] <- TRUE
  }
  
  gate2_summary <- transport_rxns %>%
    group_by(base_met) %>%
    summarise(
      import_gate2 = any(import_gate2),
      export_gate2 = any(export_gate2),
      .groups = "drop"
    )
  # Gate 2
  cat(sprintf("  Gate 2 — transport reactions found     : %d\n", nrow(transport_rxns)))
  cat(sprintf("  Gate 2 — import_gate2 TRUE             : %d\n", sum(gate2_summary$import_gate2)))
  cat(sprintf("  Gate 2 — export_gate2 TRUE             : %d\n", sum(gate2_summary$export_gate2)))
  # Gate 2 join diagnostic — did the SEED ID join find anything?
  cat(sprintf("  Gate 2 — trp_expanded rows             : %d\n", nrow(trp_expanded)))
  cat(sprintf("  Gate 2 — rxn_expanded rows             : %d\n", nrow(rxn_expanded)))
  cat(sprintf("  Gate 2 — post-join transport_rxns rows : %d\n", nrow(transport_rxns)))
  
  # -- Step E: Gate 3 — intracellular reactions --------------------------
  intracellular <- rxn %>%
    filter(
      !`SBO Term` %in% EXCLUDE_SBO,
      `Upper Bound` > 0
    ) %>%
    rowwise() %>%
    filter({
      all_mets <- c(parse_mets(Reactants), parse_mets(Products))
      length(all_mets) > 0 && !any(grepl("_e0$", all_mets))
    }) %>%
    ungroup()
  
  consumed_c0   <- unique(unlist(lapply(intracellular$Reactants, parse_mets)))
  consumed_c0   <- consumed_c0[grepl("_c0$", consumed_c0)]
  consumed_base <- sub("^M_", "", sub("_c0$", "", consumed_c0))
  
  produced_c0   <- unique(unlist(lapply(intracellular$Products, parse_mets)))
  produced_c0   <- produced_c0[grepl("_c0$", produced_c0)]
  produced_base <- sub("^M_", "", sub("_c0$", "", produced_c0))
  
  # Gate 3
  cat(sprintf("  Gate 3 — intracellular reactions       : %d\n", nrow(intracellular)))
  cat(sprintf("  Gate 3 — unique consumed base mets     : %d\n", length(consumed_base)))
  cat(sprintf("  Gate 3 — unique produced base mets     : %d\n", length(produced_base)))
  
  # -- Step F: Combine all three gates -----------------------------------
  result <- ex_map %>%
    left_join(gate2_summary, by = "base_met") %>%
    mutate(
      import_gate2     = replace_na(import_gate2, FALSE),
      export_gate2     = replace_na(export_gate2, FALSE),
      gate3_uptake     = base_met %in% consumed_base,
      gate3_export     = base_met %in% produced_base,
      
      # Uptake: Gate 1 removed — exchange LB is always 0 in gapseq default 
      # state regardless of true uptake capability. Gate 2 (transport reaction
      # directionality from bounds) is the correct and sufficient evidence
      # that inward flux is physically possible.
      uptake_confirmed = import_gate2 & gate3_uptake,
      
      # Export: Gate 1 retained — UB = 1000 on exchange reactions is
      # meaningful and confirms the boundary is open for secretion.
      export_confirmed = export_gate1 & export_gate2 & gate3_export
    )
  uptake_mets <- unique(result$base_met[result$uptake_confirmed])
  export_mets <- unique(result$base_met[result$export_confirmed])
  
  # Combined
  cat(sprintf("  Final  — uptake confirmed              : %d\n", sum(result$uptake_confirmed)))
  cat(sprintf("  Final  — export confirmed              : %d\n", sum(result$export_confirmed)))
  
  # -- Step G: Leakage tensor pairs --------------------------------------
  exportable_base <- unique(result$base_met[result$export_gate1 & result$export_gate2])
  
  tensor_pairs <- data.frame(
    reactant = character(0),
    product  = character(0),
    stringsAsFactors = FALSE
  )
  for (k in seq_len(nrow(intracellular))) {
    rcts <- parse_mets(intracellular$Reactants[k])
    prds <- parse_mets(intracellular$Products[k])
    
    rcts_c0         <- rcts[grepl("_c0$", rcts)]
    prds_c0         <- prds[grepl("_c0$", prds)]
    prds_base       <- sub("^M_", "", sub("_c0$", "", prds_c0))
    prds_exportable <- prds_base[prds_base %in% exportable_base]  # already base IDs
    rcts_base       <- sub("^M_", "", sub("_c0$", "", rcts_c0))   # strip to base IDs
    
    if (length(rcts_base) == 0 || length(prds_exportable) == 0) next
    
    pairs <- expand.grid(
      reactant = rcts_base,       # "cpd00001"
      product  = prds_exportable, # "cpd00001"
      stringsAsFactors = FALSE
    )
    tensor_pairs <- rbind(tensor_pairs, pairs)
  }
  
  return(list(
    uptake_mets  = uptake_mets,
    export_mets  = export_mets,
    tensor_pairs = distinct(tensor_pairs),
    gate_result  = result
  ))
}

# -----------------------------------------------------------------------
# 3. Run across all genera
# -----------------------------------------------------------------------
results <- vector("list", n_species)
names(results) <- genera

for (i in seq_len(n_species)) {
  g <- genera[i]
  results[[g]] <- tryCatch(
    process_genome(g),
    error = function(e) {
      cat(sprintf("[ERROR] %s: %s\n", g, conditionMessage(e)))
      list(
        uptake_mets  = character(0),
        export_mets  = character(0),
        tensor_pairs = data.frame(reactant = character(0),
                                  product  = character(0),
                                  stringsAsFactors = FALSE)
      )
    }
  )
}


# -----------------------------------------------------------------------
# 4. Build uptake matrix
# -----------------------------------------------------------------------
all_uptake_mets <- unique(unlist(lapply(results, `[[`, "uptake_mets")))

uptake_matrix <- matrix(
  0L,
  nrow = n_species,
  ncol = length(all_uptake_mets),
  dimnames = list(genera, all_uptake_mets)
)

for (g in genera) {
  mets <- results[[g]]$uptake_mets
  if (length(mets) > 0)
    uptake_matrix[g, intersect(mets, all_uptake_mets)] <- 1L
}

# -----------------------------------------------------------------------
# 5. Build export matrix
# -----------------------------------------------------------------------
all_export_mets <- unique(unlist(lapply(results, `[[`, "export_mets")))

export_matrix <- matrix(
  0L,
  nrow = n_species,
  ncol = length(all_export_mets),
  dimnames = list(genera, all_export_mets)
)

for (g in genera) {
  mets <- results[[g]]$export_mets
  if (length(mets) > 0)
    export_matrix[g, intersect(mets, all_export_mets)] <- 1L
}

# -----------------------------------------------------------------------
# 6. Build leakage tensor [species x reactant x product]
# -----------------------------------------------------------------------
all_reactants <- unique(unlist(lapply(results, function(r) r$tensor_pairs$reactant)))
all_products  <- unique(unlist(lapply(results, function(r) r$tensor_pairs$product)))

tensor <- array(
  0L,
  dim      = c(n_species, length(all_reactants), length(all_products)),
  dimnames = list(
    species  = genera,
    reactant = all_reactants,
    product  = all_products
  )
)

for (i in seq_len(n_species)) {
  pairs <- results[[genera[i]]]$tensor_pairs
  if (nrow(pairs) == 0) next
  r_idx <- match(pairs$reactant, all_reactants)
  p_idx <- match(pairs$product,  all_products)
  valid <- !is.na(r_idx) & !is.na(p_idx)
  tensor[cbind(i, r_idx[valid], p_idx[valid])] <- 1L
}

# Double check / Verify species alignment across all three outputs -------------
stopifnot(identical(rownames(uptake_matrix), genera))
stopifnot(identical(rownames(export_matrix), genera))
stopifnot(identical(dimnames(tensor)[[1]],   genera))

cat(sprintf(
  "[DONE] Uptake matrix: %d x %d | Export matrix: %d x %d | Tensor: %d x %d x %d\n",
  nrow(uptake_matrix), ncol(uptake_matrix),
  nrow(export_matrix), ncol(export_matrix),
  dim(tensor)[1], dim(tensor)[2], dim(tensor)[3]
))

# Save outputs --------------------------------------------------------

write_rds(uptake_matrix, file.path("matrices_usethis", paste0("uptake_matrix_", basename, ".rds")))
write_rds(export_matrix, file.path("matrices_usethis", paste0("export_matrix_", basename, ".rds")))
write_rds(tensor, file.path("matrices_usethis", paste0("leakage_tensor_", basename, ".rds")))

cat("[DONE] All matrices saved.\n")

