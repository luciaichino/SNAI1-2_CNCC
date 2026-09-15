#!/usr/bin/env Rscript

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 7: statistics behind the eight bars and the four nested brackets.
#
# Statistic per element set (elements with >= 50% of bases aligned):
#   phastCons470way  fraction of aligned bases with score >= 0.8, pooled over elements
#   LINSIGHT         mean of the per-element mean score
#
# Background tests (draw z-test): the observed value is compared with the mean and
# SD of the same statistic across background draws; p = one-sided normal tail.
#   all peaks           vs the 200 open-chromatin draws
#   K27ac peak subset   vs the same draws restricted to the regions paired to those peaks
#   all motifs          vs the 100 word-matched draws
#   K27ac motif subset  vs a word-matched null taken from the same draws at the
#                       subset's word composition (one sample per draw)
#
# Nested tests: the subset against N_DRAWS_NESTED composition-matched samples of
# the elements of its parent set that are not in the subset
#   peaks   matched on width band x GC tertile
#   motifs  matched on 9-mer word
#
# Outputs, written to $RES:
#   final_comparison_values.tsv   panel group n_fg obs bg_mean bg_sd n_draws z p
#   final_comparison_nested.tsv   panel comparison obs bg_mean bg_sd n_draws z p
#   covariate_balance.tsv         comparison covariate focal_mean bg_mean smd

suppressMessages(library(data.table))
setDTthreads(1)

INT    <- Sys.getenv("INT")
RES    <- Sys.getenv("RES")
NNEST  <- as.integer(Sys.getenv("N_DRAWS_NESTED", "1000"))
COVMIN <- 0.5                       # minimum aligned fraction per element

MOTIF_COLS <- c("chrom", "start", "end", "id", "score", "strand", "pval", "word")


# =============================================================================
# 1. Helpers
# =============================================================================

# Score tables written by step 06, restricted to sufficiently aligned elements.
read_cons <- function(nm) {
    d <- fread(file.path(INT, "scores", paste0(nm, ".phastCons470way.tsv")),
               col.names = c("name", "size", "covered", "mean", "bp_hi"))
    d[covered / size >= COVMIN]
}

read_lin <- function(nm) {
    d <- fread(file.path(INT, "scores", paste0(nm, ".linsight.tsv")),
               col.names = c("name", "size", "covered", "mean"))
    d[covered / size >= COVMIN]
}

# The two summary statistics.
frac_hi <- function(d) sum(d$bp_hi) / sum(d$covered)
mean_sc <- function(d) mean(d$mean)

# Background names are "<prefix><draw>_<index>", e.g. atm17_2371 or km3_508.
draw_id <- function(name) as.integer(sub("^[a-z]+([0-9]+)_.*$",     "\\1", name))
pair_id <- function(name) as.integer(sub("^[a-z]+[0-9]+_([0-9]+)$", "\\1", name))

# One-sided draw z-test of an observed value against a vector of null values.
ztest <- function(obs, draws) {
    draws <- draws[is.finite(draws)]
    z <- (obs - mean(draws)) / sd(draws)

    list(obs     = obs,
         bg_mean = mean(draws),
         bg_sd   = sd(draws),
         n_draws = length(draws),
         z       = z,
         p       = pnorm(z, lower.tail = FALSE))
}

# Join `want` (which carries the required count N per group) onto `pool` by the
# column(s) `on`, then sample within each group defined by `by` as many rows as N
# asks for (all of them if fewer exist). Returns the ids of the sampled rows.
sample_by_group <- function(pool, want, on, by = on) {
    candidates <- pool[want, on = on, nomatch = 0L]
    picked     <- candidates[, .SD[sample.int(.N, min(.N, N[1]))], by = by]
    picked$id
}


# =============================================================================
# 2. Element sets
# =============================================================================

# -- peaks and the K27ac subset --------------------------------------------------
pidx    <- fread(file.path(INT, "peak_index.tsv"))
k27_ids <- readLines(file.path(INT, "peaks_k27ac.ids"))
k27_idx <- pidx[id %in% k27_ids, idx]

peaks_bed <- fread(file.path(INT, "peaks.bed"), col.names = c("chrom", "start", "end", "peak"))
k27_bed   <- peaks_bed[peak %in% k27_ids]
setkey(k27_bed, chrom, start, end)

# -- motifs and the K27ac subset (motifs lying inside a K27ac peak) --------------
motifs <- fread(file.path(INT, "motifs_bound.bed"), col.names = MOTIF_COLS)

hits      <- foverlaps(motifs, k27_bed, by.x = c("chrom", "start", "end"), nomatch = 0L)
k27_motif <- unique(hits$id)

cat(sprintf("peaks %d (K27ac subset %d); motifs %d (in K27ac peaks %d)\n",
            nrow(pidx), length(k27_idx), nrow(motifs), length(k27_motif)))

# -- word of every element of the word-matched background -----------------------
# bg_kmer.bed carries no word; it is recovered by matching coordinates against the
# unbound occurrences it was sampled from.
unbound <- fread(file.path(INT, "motifs_unbound_distal.bed"), select = c(1:3, 8),
                 col.names = c("chrom", "start", "end", "word"))

bg_kmer <- fread(file.path(INT, "bg_kmer.bed"), col.names = c("chrom", "start", "end", "name"))
bg_kmer <- merge(bg_kmer, unbound, by = c("chrom", "start", "end"))
bg_kmer[, draw := draw_id(name)]
setnames(bg_kmer, "name", "id")


# =============================================================================
# 3. Null samples, drawn once and reused for both tracks
# =============================================================================

# -- (a) word-matched null for the K27ac motif subset ---------------------------
# Per draw and per word, as many background elements as the subset has motifs of
# that word. One null sample per draw.
need_words <- motifs[id %in% k27_motif, .(N = .N), by = word]

set.seed(2026)
k27_null_names <- sample_by_group(bg_kmer, need_words, on = "word", by = c("draw", "word"))

# -- (b) nested null for the peaks: rest of the peaks, matched on width band x GC tertile
pidx[, gct     := cut(gc, quantile(gc, 0:3 / 3), include.lowest = TRUE, labels = FALSE)]
pidx[, stratum := paste(band, gct)]

want_strata <- pidx[idx %in% k27_idx, .N, by = stratum]
rest_peaks  <- pidx[!idx %in% k27_idx]

set.seed(11)
nested_peak_ids <- lapply(seq_len(NNEST), function(i) {
    sample_by_group(rest_peaks, want_strata, "stratum")
})

# -- (c) nested null for the motifs: rest of the motifs, matched on word -----------
rest_motifs <- motifs[!id %in% k27_motif]

set.seed(2026)
nested_motif_ids <- lapply(seq_len(NNEST), function(i) {
    sample_by_group(rest_motifs, need_words, "word")
})


# =============================================================================
# 4. Statistics for one track
# =============================================================================

panel <- function(track, reader, stat) {

    # -- scores ---------------------------------------------------------------
    pk <- reader("peaks")
    ba <- reader("bg_atac")
    mo <- reader("motifs")
    bk <- reader("bg_kmer")

    ba[, draw := draw_id(name)]
    ba[, idx  := pair_id(name)]
    bk[, draw := draw_id(name)]

    # -- bar 1: all peaks vs the open-chromatin draws ----------------------------
    null_peaks <- ba[, stat(.SD), by = draw]$V1
    r_peaks    <- ztest(stat(pk), null_peaks)

    # -- bar 2: K27ac peaks vs the paired regions of the same draws --------------
    k27_peak_names <- pidx[idx %in% k27_idx, id]
    null_k27p      <- ba[idx %in% k27_idx, stat(.SD), by = draw]$V1
    r_k27p         <- ztest(stat(pk[name %in% k27_peak_names]), null_k27p)

    # -- bar 3: all motifs vs the word-matched draws -----------------------------
    null_motifs <- bk[, stat(.SD), by = draw]$V1
    r_mot       <- ztest(stat(mo), null_motifs)

    # -- bar 4: K27ac motifs vs the word-matched null at the subset's composition -
    null_k27m <- bk[name %in% k27_null_names, stat(.SD), by = draw]$V1
    r_k27m    <- ztest(stat(mo[name %in% k27_motif]), null_k27m)

    # -- assemble the four bars ---------------------------------------------------
    groups  <- c("peaks", "k27_peaks", "motifs", "k27_motifs")
    n_fg    <- c(nrow(pidx), length(k27_idx), nrow(motifs), length(k27_motif))
    results <- list(r_peaks, r_k27p, r_mot, r_k27m)

    V <- rbindlist(lapply(seq_along(groups), function(i) {
        data.table(panel = track, group = groups[i], n_fg = n_fg[i], as.data.table(results[[i]]))
    }))

    # -- nested brackets: subset vs composition-matched rest of the parent set ----
    nest_p <- sapply(nested_peak_ids,  function(ids) stat(pk[name %in% ids]))
    nest_m <- sapply(nested_motif_ids, function(ids) stat(mo[name %in% ids]))

    label_p <- sprintf("%d K27ac peaks vs %d other peaks",   length(k27_idx),   nrow(rest_peaks))
    label_m <- sprintf("%d K27ac motifs vs %d other motifs", length(k27_motif), nrow(rest_motifs))

    NEST <- rbind(
        data.table(panel = track, comparison = label_p, as.data.table(ztest(r_k27p$obs, nest_p))),
        data.table(panel = track, comparison = label_m, as.data.table(ztest(r_k27m$obs, nest_m))))

    list(V = V, NEST = NEST)
}


# =============================================================================
# 5. Run both tracks and write the tables
# =============================================================================

A <- panel("phastCons470way", read_cons, frac_hi)
B <- panel("LINSIGHT",        read_lin,  mean_sc)

V    <- rbind(A$V,    B$V)
NEST <- rbind(A$NEST, B$NEST)

fwrite(V,    file.path(RES, "final_comparison_values.tsv"), sep = "\t")
fwrite(NEST, file.path(RES, "final_comparison_nested.tsv"), sep = "\t")


# =============================================================================
# 6. Covariate balance of the peaks against their background
# =============================================================================

cov_pool <- fread(file.path(INT, "cov_pool.tsv"))
pool_bed <- fread(file.path(INT, "atac_pool.bed"), col.names = c("chrom", "start", "end", "id"))
bg_atac  <- fread(file.path(INT, "bg_atac.bed"),   col.names = c("chrom", "start", "end", "name"))

# background region -> pool id -> covariates
bg_cov <- merge(bg_atac, pool_bed, by = c("chrom", "start", "end"))
bg_cov <- merge(bg_cov, cov_pool, by = "id")
bg_cov[, idx := pair_id(name)]

# standardised mean difference of each covariate, focal minus background
balance <- function(focal, bg, label) {
    rbindlist(lapply(c("width", "gc", "cpg", "rep", "tssdist"), function(v) {
        data.table(comparison = label,
                   covariate  = v,
                   focal_mean = mean(focal[[v]]),
                   bg_mean    = mean(bg[[v]]),
                   smd        = (mean(focal[[v]]) - mean(bg[[v]])) / sd(focal[[v]]))
    }))
}

BAL <- rbind(
    balance(pidx,                   bg_cov,                   "all peaks vs full background"),
    balance(pidx[idx %in% k27_idx], bg_cov[idx %in% k27_idx], "K27ac peaks vs paired background"))

fwrite(BAL, file.path(RES, "covariate_balance.tsv"), sep = "\t")


# =============================================================================
# 7. Console summary
# =============================================================================

print(V[,    .(panel, group, n_fg, obs = round(obs, 4), bg_mean = round(bg_mean, 4), p = signif(p, 2))])
print(NEST[, .(panel, comparison, obs = round(obs, 4), rest = round(bg_mean, 4), p = signif(p, 2))])
