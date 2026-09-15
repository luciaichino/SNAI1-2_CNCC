#!/usr/bin/env bash

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 2: matching covariates for the peaks and for the open-chromatin pool.
#
# Outputs, written to $INT:
#   cov_peaks.tsv, cov_pool.tsv   columns: id width gc cpg rep nfrac tssdist
#
#   width    element length in bp
#   gc       fraction of G or C bases
#   cpg      CpG dinucleotides per bp
#   rep      soft-masked (lower-case, i.e. repeat) fraction
#   nfrac    fraction of N bases
#   tssdist  distance in bp to the nearest transcript start (tss.bed from step 01)

set -euo pipefail
: "${INT:?source config.sh first}"

cd "$INT"
mkdir -p tmp


# =============================================================================
# Function: covariates <BED4 in> <table out>
# =============================================================================

covariates () {
    local bed=$1
    local out=$2

    # -- sequence composition ------------------------------------------------
    # getfasta -tab -name gives two columns: "name::chrom:start-end" and sequence.
    bedtools getfasta -fi "$GENOME_FA" -bed "$bed" -tab -name > tmp/seq.tab

    # Count lower-case bases (repeats) on the raw sequence, then G/C, N and CpG on
    # the upper-cased sequence. gsub() returns the number of substitutions, which
    # is used here purely as a counter.
    awk -F'\t' -v OFS='\t' '{
            seq = $2
            L   = length(seq)
            if (L == 0) next

            rep = gsub(/[acgtn]/, "&", seq)

            u   = toupper($2)
            gc  = gsub(/[GC]/, "&", u)
            nn  = gsub(/N/,    "&", u)
            cg  = gsub(/CG/,   "&", u)

            split($1, a, "::")
            print a[1], L, gc / L, cg / L, rep / L, nn / L
        }' tmp/seq.tab > tmp/seq.unsorted

    sort -k1,1 tmp/seq.unsorted > tmp/seq.sorted

    # -- distance to the nearest TSS -----------------------------------------
    bedtools closest -a "$bed" -b tss.bed -d -t first > tmp/closest.tab
    awk -v OFS='\t' '{ print $4, $NF }' tmp/closest.tab > tmp/tss.unsorted
    sort -k1,1 tmp/tss.unsorted > tmp/tss.sorted

    # -- join the two tables on the element id -------------------------------
    printf 'id\twidth\tgc\tcpg\trep\tnfrac\ttssdist\n' > "$out"
    join -t $'\t' tmp/seq.sorted tmp/tss.sorted >> "$out"
}


# =============================================================================
# Run for the peaks and for the pool
# =============================================================================

covariates peaks.bed     cov_peaks.tsv
covariates atac_pool.bed cov_pool.tsv

n_peaks=$(( $(wc -l < cov_peaks.tsv) - 1 ))
n_pool=$((  $(wc -l < cov_pool.tsv)  - 1 ))
echo "covariates: $n_peaks peaks, $n_pool pool regions"

rm -rf tmp
