#!/usr/bin/env Rscript

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 3: background for the peaks.
#
# Other open CNCC regions of the same width, GC, CpG, repeat content and TSS
# distance as the SNAI2 peaks. Peaks are split into width-decile bands. Within each
# band nullranges::matchRanges draws one pool region per peak by propensity-score
# stratification on GC, CpG, repeat fraction and log10 TSS distance. Width is thus
# matched by construction and sequence composition by matching. This is repeated
# N_DRAWS_ATAC times with fixed seeds (1000 + draw number).
#
# Outputs, written to $INT:
#   peak_index.tsv   idx id band width gc cpg rep tssdist
#                    peaks ordered by width band, then by name; idx = row number
#   bg_atac.bed      chrom start end name, with name = atm{draw}_{idx}
#
# Every background region carries the idx of the peak it was drawn for. Restricting
# a draw to the idx values of the 834 K27ac peaks therefore gives the background
# already drawn for that subset, without a second matching run.

suppressMessages({
    library(data.table)
    library(GenomicRanges)
    library(nullranges)
})
setDTthreads(1)

INT   <- Sys.getenv("INT")
NDRAW <- as.integer(Sys.getenv("N_DRAWS_ATAC", "200"))


# =============================================================================
# 1. Load peaks and pool with their covariates
# =============================================================================

load_set <- function(bed_file, cov_file) {
    bed <- fread(file.path(INT, bed_file), header = FALSE,
                 col.names = c("chrom", "start", "end", "id"))
    cov <- fread(file.path(INT, cov_file))

    d <- merge(bed, cov, by = "id")

    # elements with 10% or more N bases are not usable for matching
    d <- d[nfrac < 0.1]

    d[, ltss := log10(tssdist + 1)]
    setorder(d, id)
    d
}

focal <- load_set("peaks.bed",     "cov_peaks.tsv")
pool  <- load_set("atac_pool.bed", "cov_pool.tsv")

# The pool is restricted to the width range of the peaks, so that the open-ended
# top and bottom bands do not contain regions far wider or narrower than any peak.
pool <- pool[width >= min(focal$width) & width <= max(focal$width)]


# =============================================================================
# 2. Width-decile bands and the band-major peak index
# =============================================================================

cuts <- quantile(focal$width, seq(0.1, 0.9, 0.1))

focal[, band := findInterval(width, cuts) + 1L]
pool[,  band := findInterval(width, cuts) + 1L]

# Order peaks by band, then by name; idx is the row number in that order.
setorder(focal, band, id)
focal[, idx := .I]

fwrite(focal[, .(idx, id, band, width, gc, cpg, rep, tssdist)],
       file.path(INT, "peak_index.tsv"), sep = "\t")


# =============================================================================
# 3. One matched draw
# =============================================================================

# data.table -> GRanges carrying the matching covariates. Only the covariates are
# needed here; the peak index is re-attached below, after matching.
as_gr <- function(d) {
    GRanges(d$chrom, IRanges(d$start + 1L, d$end),
            gc = d$gc, cpg = d$cpg, rep = d$rep, ltss = d$ltss)
}

# Match the peaks of one width band against the pool regions of that band.
match_band <- function(b) {
    f <- focal[band == b]
    p <- pool[band == b]

    # Sample with replacement only when the pool is small relative to the peaks.
    with_replacement <- nrow(p) < 3 * nrow(f)

    m <- matchRanges(focal = as_gr(f), pool = as_gr(p),
                     covar = ~ gc + cpg + rep + ltss,
                     method = "stratified", replace = with_replacement)
    m <- matched(m)

    data.table(chrom = as.character(seqnames(m)),
               start = start(m) - 1L,
               end   = end(m),
               idx   = f$idx)
}

draw_one <- function(k) {
    set.seed(1000L + k)

    bands <- sort(unique(focal$band))
    d <- rbindlist(lapply(bands, match_band))

    d[, name := sprintf("atm%d_%d", k, idx)]
    d[, .(chrom, start, end, name)]
}


# =============================================================================
# 4. All draws
# =============================================================================

draws <- vector("list", NDRAW)
for (k in seq_len(NDRAW)) {
    if (k %% 20 == 0) cat("draw", k, "\n")
    draws[[k]] <- draw_one(k)
}

bg <- rbindlist(draws)
setorder(bg, chrom, start, end)

fwrite(bg, file.path(INT, "bg_atac.bed"), sep = "\t", col.names = FALSE)

cat(sprintf("open-chromatin background: %d draws x %d regions\n", NDRAW, nrow(focal)))
