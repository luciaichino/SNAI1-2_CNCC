#!/usr/bin/env bash

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 5: hg38 -> hg19 coordinates for LINSIGHT, with a reciprocity check.
#
# An element is kept only if
#   - the forward lift (hg38 -> hg19) maps it exactly once,
#   - the reverse lift (hg19 -> hg38) of that result maps exactly once and lands on
#     the same chromosome within 5 bp of the original start and end,
#   - the hg19 position is on a primary chromosome (chr1-22, X, Y),
#   - the width is preserved: exactly for 9-bp motifs, within 5% for peaks and regions.
#
# Outputs, written to $INT:
#   hg19/<set>.bed          chrom start end name   (name = the hg38 element name)
#   hg19/liftover_stats.tsv set n_hg38 n_hg19

set -euo pipefail
: "${INT:?source config.sh first}"

# byte order for sort and join, so they agree with each other
export LC_ALL=C

cd "$INT"
mkdir -p hg19 tmp


# =============================================================================
# Function: lift <set name> <BED4 in> <width rule: exact | tol5>
# =============================================================================

lift () {
    local nm=$1
    local src=$2
    local rule=$3

    # -- tag every element with a unique key "u<row>|<name>" ------------------
    # liftOver keeps column 4, so the key survives both lifts and lets the
    # hg38, hg19 and round-trip coordinates be joined back together.
    awk -v OFS='\t' '{ print $1, $2, $3, "u" NR "|" $4 }' "$src" > tmp/$nm.in.bed

    # -- forward and reverse lifts ---------------------------------------------
    liftOver -minMatch=0.95 tmp/$nm.in.bed  "$CHAIN_38TO19" tmp/$nm.fwd.bed  tmp/$nm.unmapped_fwd  2>/dev/null
    liftOver -minMatch=0.95 tmp/$nm.fwd.bed "$CHAIN_19TO38" tmp/$nm.back.bed tmp/$nm.unmapped_back 2>/dev/null

    # -- one table per stage, keyed by the tag -------------------------------
    awk -v OFS='\t' '{ print $4, $1, $2, $3 }' tmp/$nm.in.bed   > tmp/$nm.orig.unsorted
    awk -v OFS='\t' '{ print $4, $1, $2, $3 }' tmp/$nm.fwd.bed  > tmp/$nm.fwd.unsorted
    awk -v OFS='\t' '{ print $4, $1, $2, $3 }' tmp/$nm.back.bed > tmp/$nm.back.unsorted

    sort -k1,1 tmp/$nm.orig.unsorted > tmp/$nm.orig.tab
    sort -k1,1 tmp/$nm.fwd.unsorted  > tmp/$nm.fwd.tab
    sort -k1,1 tmp/$nm.back.unsorted > tmp/$nm.back.tab

    # joined columns:
    #   1 tag  2-4 hg38 chrom start end  5-7 hg19 chrom start end  8-10 round-trip chrom start end
    join -t $'\t' tmp/$nm.orig.tab   tmp/$nm.fwd.tab  > tmp/$nm.orig_fwd.tab
    join -t $'\t' tmp/$nm.orig_fwd.tab tmp/$nm.back.tab > tmp/$nm.joined.tab

    # -- apply the filters -------------------------------------------------------
    # The file is read twice: the first pass counts how often each tag occurs
    # (more than once = ambiguous mapping), the second pass filters.
    awk -F'\t' -v OFS='\t' -v rule="$rule" -v tol=5 '
        NR == FNR { seen[$1]++; next }

        {
            tag = $1
            if (seen[tag] != 1) next

            c38 = $2; s38 = $3; e38 = $4
            c19 = $5; s19 = $6; e19 = $7
            cb  = $8; sb  = $9; eb  = $10

            # round trip must return to the same chromosome within tol bp
            if (cb != c38) next
            if (sb - s38 > tol || s38 - sb > tol) next
            if (eb - e38 > tol || e38 - eb > tol) next

            # primary hg19 chromosomes only
            if (c19 !~ /^chr([0-9]+|X|Y)$/) next

            # width rule
            w19 = e19 - s19
            w38 = e38 - s38
            if (rule == "exact" && w19 != w38) next
            if (rule != "exact" && (w19 < 0.95 * w38 || w19 > 1.05 * w38)) next

            # strip the "u<row>|" prefix to recover the original name
            split(tag, a, "|")
            print c19, s19, e19, a[2]
        }' tmp/$nm.joined.tab tmp/$nm.joined.tab > tmp/$nm.hg19.unsorted

    sort -k1,1 -k2,2n tmp/$nm.hg19.unsorted > hg19/$nm.bed

    # -- bookkeeping -------------------------------------------------------------
    printf '%s\t%d\t%d\n' "$nm" "$(wc -l < "$src")" "$(wc -l < hg19/$nm.bed)" >> hg19/liftover_stats.tsv

    rm -f tmp/$nm.*
}


# =============================================================================
# Lift the four element sets
# =============================================================================

printf 'set\tn_hg38\tn_hg19\n' > hg19/liftover_stats.tsv

# motifs_bound.bed has eight columns; liftOver needs BED4 here
cut -f1-4 motifs_bound.bed > tmp/motifs4.bed

lift peaks    peaks.bed        tol5
lift bg_atac  bg_atac.bed      tol5
lift motifs   tmp/motifs4.bed  exact
lift bg_kmer  bg_kmer.bed      exact

rm -rf tmp

column -t hg19/liftover_stats.tsv
