#!/usr/bin/env bash

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 6: per-element constraint scores for the four element sets.
#
# Outputs, written to $INT/scores:
#   <set>.phastCons470way.tsv   name size covered mean bp_ge0.8      (hg38)
#       covered  = bases with a phastCons value, i.e. aligned bases
#       mean     = mean phastCons470way over the covered bases
#       bp_ge0.8 = bases with phastCons470way >= 0.8
#   <set>.linsight.tsv          name size covered mean               (hg19)

set -euo pipefail
: "${INT:?source config.sh first}"

cd "$INT"
mkdir -p scores tmp


# =============================================================================
# 1. Genome-wide intervals with phastCons470way >= 0.8, computed once
# =============================================================================

# The bedGraph of the whole track is far too large to keep on disk, so this one
# command chain is deliberately streamed: bigWig -> bedGraph -> threshold ->
# merged intervals.
if [ ! -s phastCons470way_ge0.8.bed ]; then
    bigWigToBedGraph "$PHASTCONS_BW" stdout \
        | awk -v OFS='\t' '$4 >= 0.8 { print $1, $2, $3 }' \
        | bedtools merge -i - \
        | sort -k1,1 -k2,2n > phastCons470way_ge0.8.bed
fi


# =============================================================================
# Function: score_cons <set name> <BED4, hg38>
# =============================================================================

score_cons () {
    local nm=$1
    local bed=$2

    # bigWigAverageOverBed columns: name size covered sum mean0 mean
    bigWigAverageOverBed "$PHASTCONS_BW" "$bed" tmp/$nm.tab

    # bases per element inside the >= 0.8 intervals: -wo appends the overlap
    # length as the last column, summed per element name
    bedtools intersect -a "$bed" -b phastCons470way_ge0.8.bed -wo > tmp/$nm.overlaps
    awk -v OFS='\t' '{ bp[$4] += $NF } END { for (k in bp) print k, bp[k] }' tmp/$nm.overlaps > tmp/$nm.hi

    # combine: name size covered mean bp_ge0.8 (0 when the element has no such bases)
    awk -v OFS='\t' '
        NR == FNR { hi[$1] = $2; next }
        { print $1, $2, $3, $6, ($1 in hi ? hi[$1] : 0) }
        ' tmp/$nm.hi tmp/$nm.tab > scores/$nm.phastCons470way.tsv

    rm -f tmp/$nm.tab tmp/$nm.overlaps tmp/$nm.hi
}


# =============================================================================
# Function: score_lin <set name> <BED4, hg19>
# =============================================================================

score_lin () {
    local nm=$1
    local bed=$2

    bigWigAverageOverBed "$LINSIGHT_BW" "$bed" tmp/$nm.lin.tab
    awk -v OFS='\t' '{ print $1, $2, $3, $6 }' tmp/$nm.lin.tab > scores/$nm.linsight.tsv

    rm -f tmp/$nm.lin.tab
}


# =============================================================================
# 2. Score the four sets
# =============================================================================

# motifs_bound.bed has eight columns; the scorers need BED4
cut -f1-4 motifs_bound.bed > tmp/motifs.bed

score_cons peaks    peaks.bed
score_lin  peaks    hg19/peaks.bed

score_cons bg_atac  bg_atac.bed
score_lin  bg_atac  hg19/bg_atac.bed

score_cons motifs   tmp/motifs.bed
score_lin  motifs   hg19/motifs.bed

score_cons bg_kmer  bg_kmer.bed
score_lin  bg_kmer  hg19/bg_kmer.bed

for nm in peaks bg_atac motifs bg_kmer; do
    echo "scored $nm: $(wc -l < scores/$nm.phastCons470way.tsv) hg38, $(wc -l < scores/$nm.linsight.tsv) hg19"
done

rm -rf tmp
