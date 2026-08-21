#put code description at top like in class tutorials 
# you will still need a mapping file 

#see Notion to see which BioProjects have been downloaded 

# -------- load packages and metadata 
packages <- list('dada2', 'phyloseq', 'ggplot2', 'DECIPHER', 'phangorn','dplyr','readr',
                 'BiocManager', 'phyloseq','reshape2')
lapply(packages, library, character.only=TRUE)

setwd('/media/ore-s/0E10-169F/capstone26/sra_2')

# -------- load / edit seq metadata file -- this step probably wont be in pipeline since metadata varies so much ...
metadata <- read.csv('SraRunTable.csv')
seq_nums <- read_tsv('SRR_Acc_List.txt',col_names=FALSE) #file that has list of downloaded seq 
colnames(seq_nums)[1] = 'Run'

#meta2 <- metadata[c('Run','Sample_ID','Sample.Name')] %>% inner_join(seq_nums, by="Run")

# --------- load seq data
seq_path <- 'data' #file where the seq are located
list.files(seq_path)

#sort forward and reverse reads 
# forward and reverse fastq filenames have format: SAMPLENAME_R1_001.fastq and SAMPLENAME_R2_001.fastq
fnFs <- sort(list.files(seq_path, pattern="_R1_001.fastq.gz", full.names = TRUE))
fnRs <- sort(list.files(seq_path, pattern="_R2_001.fastq.gz", full.names = TRUE))
# Extract sample names, assuming filenames have format: SAMPLENAME_XXX.fastq
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)


# -------- plot quality scores to decide where to trim reads --> need to plot to determine these cutoffs !
# the higher the score the better 
f.plots <- plotQualityProfile(fnFs[6:7])
r.plots <- plotQualityProfile(fnRs[6:7])

ggsave(f.plots, filename='f-plot.png')
ggsave(r.plots, filename='r-plot.png')

# -------- do the filtering : explanations below 
# maxEE: parameter sets the maximum number of “expected errors” allowed in a read
# rm.phix: Default TRUE. If TRUE, discard reads that match against the phiX genome, as determined by isPhiX (phiX is a well-characterized bacteriophage genome)
# trunQ: Default 2. Truncate reads at the first instance of a quality score less than or equal to truncQ
# trimLeft: This removes a fixed number of nucleotides from the start of the read

# Place filtered files in filtered/ subdirectory

filtFs <- file.path(seq_path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(seq_path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

filtered <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=c(240,240), #240 is for forward seq, 240 is for reverse seq
                          maxN=0, maxEE=c(2,2), truncQ=2, rm.phix=TRUE,
                          compress=TRUE, multithread=TRUE)

head(filtered)
# -------- visualize the error rates of the reads 
#(will be used to make a plot produced for each type of seq error- ie A to C artifact; G to T...)
# these error rates will be used in the making of ASVs !
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)
plotErrors(errF, nominalQ=TRUE)

#generate ASVs using filter rate
dadaFs <- dada(filtFs, err=errF, multithread=TRUE)
dadaRs <- dada(filtRs, err=errR, multithread=TRUE)

# -------- merge forward and reverse for full denoised seq 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)

seqtab <- makeSequenceTable(mergers) #making a seq table

#inspect table
dim(seqtab)
table(nchar(getSequences(seqtab)))

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)

#find proportion of chimeras to original ASV table
sum(seqtab.nochim)/sum(seqtab)

# -------- tracking # of reads lost at each step (just a measure of good habits to check)

getN <- function(x) sum(getUniques(x))
track <- cbind(filtered, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
# If processing a single sample, remove the sapply calls: e.g. replace sapply(dadaFs, getN) with getN(dadaFs)
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)


# -------- assign taxonomy based on SILVA database

taxa <- assignTaxonomy(seqtab.nochim, "~/Documents/sra-test-data/silva_nr99_v138.2_toGenus_trainset.fa.gz", multithread=TRUE)

taxa <- addSpecies(taxa, "silva_nr99_v138.2_toSpecies_trainset.fa.gz")

#some filtering to remove unknown IDs
family_id <- taxa[!is.na(taxa[, 5]),] %>% as.data.frame()

genus_id <- taxa[!is.na(taxa[, 6]),] %>% as.data.frame()

#format family level ID table & get counts
family_id$ASV <- rownames(family_id)
rownames(family_id) <- 1:nrow(family_id)

family_counts <- table(family_id$family)

write.csv(family_counts, 'family_counts.csv')
write.csv(family_id, 'family-lvl-ids.csv')

#format genus level ID table & get counts
genus_id$ASV <- rownames(genus_id)
rownames(genus_id) <- 1:nrow(genus_id)

genus_counts <- table(genus_id$Genus)

#write tables to csv
write.csv(taxa, 'all_asv_ids.csv')
write.csv(genus_counts, 'genus_counts.csv')
write.csv(genus_id, 'genus-lvl-ids.csv')

# -------- melt asv table 

asv_counts <- seqtab.nochim %>% as.data.frame()
asv_counts$Run <- row.names(asv_counts)
row.names(asv_counts) <- 1:nrow(asv_counts)

asv_counts <- melt(asv_counts)
colnames(asv_counts)[2:3] <- c("ASV", "Count")

asv_info <- inner_join(asv_counts, meta2)


# -------------- csv writing 
# join to first asv table and write to csv
left_join(genus_id, asv_info, by='ASV') %>% write.csv('genus_sample_asv_seqtable.csv')
left_join(family_id, asv_info, by='ASV') %>% write.csv('family_sample_asv_seqtable.csv')

write.csv(asv_fasta, 'sra_1_inputpicrust2.fna')

# --------------- formatting for picrust2
# Extract ASV sequences from the column names of your ASV table (seqtab.nochim)
asv_seqs <- colnames(seqtab.nochim)

# Create custom headers (e.g., ASV1, ASV2, etc.)
asv_headers <- paste0(">ASV", seq_along(asv_seqs))

# Combine headers and sequences into a single vector
asv_fasta <- c(rbind(asv_headers, asv_seqs))

# Write to a FASTA file
writeLines(asv_fasta, "asv_sequences_fasta.fasta")


# Transpose your ASV table so samples are rows and ASVs are columns
seqtab_transposed <- t(seqtab.nochim)

# Rename rows from sequences to simple ASV IDs
rownames(seqtab_transposed) <- paste0("ASV", seq_along(asv_seqs))

# Add a column for the ASV ID (required by PICRUSt2)
seqtab_final <- cbind(Sequence_ID = rownames(seqtab_transposed), seqtab_transposed)

write.table(seqtab_final, file="asv_biome_table.tsv", sep="\t", row.names=FALSE, quote=FALSE)
