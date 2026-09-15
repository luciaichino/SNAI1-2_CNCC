#!/usr/bin/env Rscript

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 4: background for the motif occurrences.
#
# Unbound copies of the identical 9-mer word. For every word among the bound
# occurrences, the same number of unbound distal occurrences of that exact word are
# sampled without replacement (all of them if fewer exist). Repeated N_DRAWS_KMER
# times with fixed seeds (900 + draw number).
#
# Outputs, written to $INT:
#   bg_kmer.bed            chrom start end name, with name = km{draw}_{i}
#   kmer_availability.tsv  per word: N bound occurrences, number of unbound copies

suppressMessages(library(data.table))
setDTthreads(1)

INT   <- Sys.getenv("INT")
NDRAW <- as.integer(Sys.getenv("N_DRAWS_KMER", "100"))

BED_COLS <- c("chrom", "start", "end", "id", "score", "strand", "pval", "word")


# =============================================================================
# 1. Bound and unbound occurrences
# =============================================================================

bound   <- fread(file.path(INT, "motifs_bound.bed"),          col.names = BED_COLS)
unbound <- fread(file.path(INT, "motifs_unbound_distal.bed"), col.names = BED_COLS)

# How many copies of each word are needed per draw.
need <- bound[, .(N = .N), by = word]

setkey(unbound, word)


# =============================================================================
# 2. One word-matched draw
# =============================================================================

draw_one <- function(k) {
    set.seed(900L + k)

    # Unbound occurrences of the words present among the bound occurrences, each
    # row carrying the number N required for its word.
    candidates <- unbound[need, on = "word", nomatch = 0L]

    # Sample N rows per word (or all rows if fewer than N exist).
    picked <- candidates[, .SD[sample.int(.N, min(.N, N[1]))], by = word]

    picked[, .(chrom, start, end, name = sprintf("km%d_%d", k, .I))]
}


# =============================================================================
# 3. All draws
# =============================================================================

draws <- vector("list", NDRAW)
for (k in seq_len(NDRAW)) {
    draws[[k]] <- draw_one(k)
}

bg <- rbindlist(draws)
setorder(bg, chrom, start, end)

fwrite(bg, file.path(INT, "bg_kmer.bed"), sep = "\t", col.names = FALSE)


# =============================================================================
# 4. Availability table: how many unbound copies exist per word
# =============================================================================

available <- unbound[, .(available = .N), by = word]
avail     <- merge(need, available, all.x = TRUE)

fwrite(avail, file.path(INT, "kmer_availability.tsv"), sep = "\t")

cat(sprintf("word-matched background: %d draws, %d words, %d occurrences per draw\n",
            NDRAW, nrow(need), nrow(bg) / NDRAW))
