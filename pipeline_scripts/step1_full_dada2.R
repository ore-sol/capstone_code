# ============================================================
# DADA2 Full Pipeline: QC → ASV → Chimera Removal → Taxonomy
# Features:
#   - Explicit thread count passed to all multithread args
#   - RDS checkpoints after each major step
#   - Live log writing throughout
# ============================================================

# ---- Setup ------------------------------------------------
THREADS <- 32

# Live logging to file
log_file <- file("pipeline.log", open = "wt")
sink(log_file, type = "output")
sink(log_file, type = "message")

cat("============================================================\n")
cat("Pipeline started:", format(Sys.time()), "\n")
cat("Cores requested:", THREADS, "\n")
cat("Cores available:", parallel::detectCores(), "\n")
cat("============================================================\n\n")

# ---- Load packages ----------------------------------------
cat("Loading packages:", format(Sys.time()), "\n")
packages <- list('dada2', 'ggplot2', 'DECIPHER', 'dplyr', 'readr',
                 'BiocManager', 'data.table', 'reshape2')
lapply(packages, library, character.only = TRUE)
cat("Packages loaded\n\n")

# ---- Get sequence files -----------------------------------
cat("Reading sequence files:", format(Sys.time()), "\n")
seq_path <- '../fastqfiles_1'
list.files(seq_path)

fnFs <- sort(list.files(seq_path, pattern = "_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(seq_path, pattern = "_R2_001.fastq.gz", full.names = TRUE))
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)

cat("Forward files found:", length(fnFs), "\n")
cat("Reverse files found:", length(fnRs), "\n")
cat("Sample names:\n")
cat(paste(" ", sample.names, collapse = "\n"), "\n\n")


# -------- plot quality scores to decide where to trim reads --> need to plot to determine these cutoffs !
# the higher the score the better
f.plots <- plotQualityProfile(fnFs[6:7])
r.plots <- plotQualityProfile(fnRs[6:7])

ggsave(f.plots, filename='f-plot.png')
ggsave(r.plots, filename='r-plot.png')

# ---- Filter and Trim --------------------------------------
# maxEE:    Maximum number of expected errors allowed in a read
# rm.phix:  Discard reads matching the phiX genome
# truncQ:   Truncate reads at first instance of quality score <= truncQ
# trimLeft: Removes fixed number of nucleotides from start of read

filtFs <- file.path(seq_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(seq_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

if (!file.exists("filtered_reads.rds")) {
    cat("Starting filterAndTrim:", format(Sys.time()), "\n")
    t0 <- system.time({
        filtered <- filterAndTrim(
            fnFs, filtFs, fnRs, filtRs,
            truncLen = c(250, 250),   # forward, reverse trunc lengths
            maxN     = 0,
            maxEE    = c(2, 2),
            truncQ   = 2,
            rm.phix  = TRUE,
            compress = TRUE,
            multithread = THREADS     # <-- THREADS passed here
        )
    })
    saveRDS(filtered, "filtered_reads.rds")
    cat("filterAndTrim done:", format(Sys.time()), "\n")
    cat("Time taken:", t0["elapsed"], "seconds\n")
    cat("Reads in vs reads out:\n")
    print(head(filtered))
    cat("\n")
} else {
    cat("Checkpoint found — loading filtered_reads.rds\n\n")
    filtered <- readRDS("filtered_reads.rds")
}

# ---- Learn Error Rates ------------------------------------
if (!file.exists("errF.rds") || !file.exists("errR.rds")) {
    cat("Learning error rates:", format(Sys.time()), "\n")
    t1 <- system.time({
        errF <- learnErrors(filtFs, multithread = THREADS)   # <-- THREADS passed here
        errR <- learnErrors(filtRs, multithread = THREADS)   # <-- THREADS passed here
    })
    saveRDS(errF, "errF.rds")
    saveRDS(errR, "errR.rds")
    cat("Error learning done:", format(Sys.time()), "\n")
    cat("Time taken:", t1["elapsed"], "seconds\n\n")
} else {
    cat("Checkpoint found — loading errF.rds and errR.rds\n\n")
    errF <- readRDS("errF.rds")
    errR <- readRDS("errR.rds")
}

# Save error plots to file (avoids display issues on HPC)
cat("Saving error rate plots:", format(Sys.time()), "\n")
png("errF_plot.png", width = 1200, height = 800)
print(plotErrors(errF, nominalQ = TRUE))
dev.off()
png("errR_plot.png", width = 1200, height = 800)
print(plotErrors(errR, nominalQ = TRUE))
dev.off()
cat("Error plots saved\n\n")

# ---- Denoise (DADA2 core) ---------------------------------
if (!file.exists("dadaFs.rds") || !file.exists("dadaRs.rds")) {
    cat("Starting DADA2 denoising:", format(Sys.time()), "\n")
    t2 <- system.time({
        dadaFs <- dada(filtFs, err = errF, multithread = THREADS)  # <-- THREADS passed here
        dadaRs <- dada(filtRs, err = errR, multithread = THREADS)  # <-- THREADS passed here
    })
    saveRDS(dadaFs, "dadaFs.rds")
    saveRDS(dadaRs, "dadaRs.rds")
    cat("Denoising done:", format(Sys.time()), "\n")
    cat("Time taken:", t2["elapsed"], "seconds\n\n")
} else {
    cat("Checkpoint found — loading dadaFs.rds and dadaRs.rds\n\n")
    dadaFs <- readRDS("dadaFs.rds")
    dadaRs <- readRDS("dadaRs.rds")
}

# ---- Merge Paired Reads -----------------------------------
if (!file.exists("mergers.rds")) {
    cat("Merging paired reads:", format(Sys.time()), "\n")
    t3 <- system.time({
        mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE)
    })
    saveRDS(mergers, "mergers.rds")
    cat("Merging done:", format(Sys.time()), "\n")
    cat("Time taken:", t3["elapsed"], "seconds\n\n")
} else {
    cat("Checkpoint found — loading mergers.rds\n\n")
    mergers <- readRDS("mergers.rds")
}

# ---- Make Sequence Table ----------------------------------
if (!file.exists("seqtab_wchimeras.rds")) {
    cat("Making sequence table:", format(Sys.time()), "\n")
    seqtab <- makeSequenceTable(mergers)
    cat("Sequence table dimensions:", dim(seqtab), "\n")
    cat("Sequence length distribution:\n")
    print(table(nchar(getSequences(seqtab))))
    saveRDS(seqtab, "seqtab_wchimeras.rds")
    write.csv(seqtab, "seqtab_wchimeras.csv")
    cat("Sequence table saved\n\n")
} else {
    cat("Checkpoint found — loading seqtab_wchimeras.rds\n\n")
    seqtab <- readRDS("seqtab_wchimeras.rds")
}

# ---- Load seqtab (if resuming from CSV) -------------------
# Uncomment the block below if resuming pipeline from a saved CSV
# cat("Loading seqtab from CSV:", format(Sys.time()), "\n")
# seqtab <- as.matrix(
#   data.frame(fread("seqtab_wchimeras.csv", header = TRUE), row.names = 1)
# )
# cat("seqtab dimensions:", dim(seqtab), "\n")
# cat("  Samples:", nrow(seqtab), "\n")
# cat("  ASVs:",    ncol(seqtab), "\n\n")

# ---- Chimera Removal --------------------------------------
if (!file.exists("seqtab_nochim.rds")) {
    cat("Starting chimera removal:", format(Sys.time()), "\n")
    t4 <- system.time({
        seqtab.nochim <- removeBimeraDenovo(
            seqtab,
            method      = "consensus",
            multithread = THREADS,     # <-- THREADS passed here
            verbose     = TRUE
        )
    })
    saveRDS(seqtab.nochim, "seqtab_nochim.rds")
    cat("Chimera removal done:", format(Sys.time()), "\n")
    cat("Time taken:", t4["elapsed"], "seconds\n")
    cat("Proportion of reads retained:", sum(seqtab.nochim) / sum(seqtab), "\n\n")
} else {
    cat("Checkpoint found — loading seqtab_nochim.rds\n\n")
    seqtab.nochim <- readRDS("seqtab_nochim.rds")
}

# ---- Taxonomy Assignment ----------------------------------
if (!file.exists("taxa_raw.rds")) {
    cat("Starting assignTaxonomy:", format(Sys.time()), "\n")
    cat("This is the slowest step — estimated 2-4hrs with", THREADS, "cores\n")
    t5 <- system.time({
        taxa <- assignTaxonomy(
            seqtab.nochim,
            "silva_nr99_v138.2_toGenus_trainset.fa.gz",
            multithread = THREADS      # <-- THREADS passed here
        )
    })
    saveRDS(taxa, "taxa_raw.rds")
    cat("assignTaxonomy done:", format(Sys.time()), "\n")
    cat("Time taken:", t5["elapsed"], "seconds\n\n")
} else {
    cat("Checkpoint found — loading taxa_raw.rds\n\n")
    taxa <- readRDS("taxa_raw.rds")
}

# ---- Add Species ------------------------------------------
if (!file.exists("taxa_species.rds")) {
    cat("Starting addSpecies:", format(Sys.time()), "\n")
    t6 <- system.time({
        taxa <- addSpecies(taxa, "silva_v138.2_assignSpecies.fa.gz")
    })
    saveRDS(taxa, "taxa_species.rds")
    cat("addSpecies done:", format(Sys.time()), "\n")
    cat("Time taken:", t6["elapsed"], "seconds\n\n")
} else {
    cat("Checkpoint found — loading taxa_species.rds\n\n")
    taxa <- readRDS("taxa_species.rds")
}

# ---- Filter Taxonomy Table --------------------------------
cat("Filtering taxonomy table:", format(Sys.time()), "\n")
family_id <- taxa[!is.na(taxa[, 5]), ] %>% as.data.frame()
genus_id  <- taxa[!is.na(taxa[, 6]), ] %>% as.data.frame()
cat("Family-level IDs:", nrow(family_id), "\n")
cat("Genus-level IDs:",  nrow(genus_id),  "\n\n")

# ---- Process Family Level ---------------------------------
cat("Processing family-level IDs:", format(Sys.time()), "\n")
family_id$ASV    <- rownames(family_id)
rownames(family_id) <- 1:nrow(family_id)
family_counts    <- table(family_id$Family)
write.csv(family_counts, 'family_counts.csv')
write.csv(family_id,     'family-lvl-ids.csv')
cat("Family CSVs written\n\n")

# ---- Process Genus Level ----------------------------------
cat("Processing genus-level IDs:", format(Sys.time()), "\n")
genus_id$ASV    <- rownames(genus_id)
rownames(genus_id) <- 1:nrow(genus_id)
genus_counts    <- table(genus_id$Genus)
write.csv(taxa,         'all_asv_ids.csv')
write.csv(genus_counts, 'genus_counts.csv')
write.csv(genus_id,     'genus-lvl-ids.csv')
cat("Genus CSVs written\n\n")

# ---- Melt ASV Table ---------------------------------------
cat("Melting ASV table:", format(Sys.time()), "\n")
asv_counts          <- seqtab.nochim %>% as.data.frame()
asv_counts$Run      <- row.names(asv_counts)
row.names(asv_counts) <- 1:nrow(asv_counts)
asv_counts          <- melt(asv_counts)
colnames(asv_counts)[2:3] <- c("ASV", "Count")
cat("ASV table melted\n\n")

# ---- Join and Write Final Tables --------------------------
cat("Writing final output tables:", format(Sys.time()), "\n")
all.info.df <- left_join(genus_id, asv_counts, by = 'ASV')
write.csv(all.info.df, 'genus_sample_asv_seqtable.csv')
left_join(family_id, asv_counts, by = 'ASV') %>% write.csv('family_sample_asv_seqtable.csv')
cat("Main tables written\n\n")

# ---- Genera List ------------------------------------------
genera <- all.info.df$Genus %>% unique() %>% as.data.frame()
write_delim(genera, "genera_list.txt")
cat("Genera list written\n\n")

# ---- All ASV fasta ---------------------------------------
#get / format table
library(tidyr)
taxa_table <- readRDS("taxa_raw.rds") %>% as.data.frame()
genus.table <- taxa_table[!is.na(taxa_table$Genus),] 
genus.table <- genus.table[genus.table$Kingdom == "Bacteria" ,]
genus.table$Row_num_id <- 1:nrow(genus.table)

library(readr)
write_rds(genus.table, "filtered_genus_taxaraw.rds")


#format and write fasta
asv_seqs <- row.names(genus.table)  # ASV sequences
asv_info <- paste0(genus.table$Row_num_id,"_",genus.table$Genus)
asv_headers <- paste0(">ASV", "_row_", asv_info)  # ASV identifiers using row num and genus 

fasta_file <- "genus_ASVs.fasta"
writeLines(c(rbind(asv_headers, asv_seqs)), fasta_file)

# ---- Sample Trackings ------------------------------------

f.plots <- plotQualityProfile(fnFs[6:7])
r.plots <- plotQualityProfile(fnRs[6:7])

ggsave(f.plots, filename='f-plot.png')
ggsave(r.plots, filename='r-plot.png')

# -------- tracking # of reads lost at each step (just a measure of good habits to check)

getN <- function(x) sum(getUniques(x))
track <- cbind(filtered, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
print(head(track))

write.csv(track, "sample_tracking.csv")

# ---- Done -------------------------------------------------
cat("============================================================\n")
cat("Pipeline completed:", format(Sys.time()), "\n")
cat("============================================================\n")

sink(type = "output")
sink(type = "message")
close(log_file)

