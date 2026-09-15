#!/usr/bin/env bash

# Date: 9/9/2026
# Author: Jaaved Mohammed

# Step 1: define the element sets on hg38 primary chromosomes (chr1-22, X, Y).
#
# Outputs, all written to $INT (see config.sh):
#   chrom.sizes               primary chromosomes only
#   tss.bed                   one TSS per transcript (1 bp)
#   promoters.bed             TSS +/- PROMOTER_HALF bp, merged
#   exclude.bed               blacklist + assembly gaps, merged
#   peaks.bed                 SNAI2 distal peaks (BED4, original names)
#   peaks_k27ac.ids           names of the peaks overlapping H3K27ac-gain enhancers
#   motifs_all.bed            every FIMO occurrence outside excluded regions
#                             (chrom start end id score strand pval word)
#   motifs_bound.bed          FIMO occurrences inside peaks
#   motifs_unbound_distal.bed FIMO occurrences outside peaks and promoter windows
#   atac_pool.bed             CNCC accessible regions outside peaks and promoters
#                             (the pool the peak background is drawn from)

set -euo pipefail
: "${INT:?source config.sh first}"

cd "$INT"
mkdir -p tmp

PRIMARY='^chr([0-9]+|X|Y)$'


# =============================================================================
# 1. Chromosome sizes, primary chromosomes only
# =============================================================================

awk -v re="$PRIMARY" '$1 ~ re' "$CHROM_SIZES" > tmp/chrom.unsorted
sort -k1,1 tmp/chrom.unsorted > chrom.sizes


# =============================================================================
# 2. Transcription start sites and promoter windows
# =============================================================================

# One TSS per transcript from the GENCODE GTF. For + strand transcripts the TSS
# is the feature start, for - strand transcripts it is the feature end.
# GTF is 1-based, BED is 0-based, hence the "- 1".
awk -F'\t' -v re="$PRIMARY" '
    $1 ~ re && $3 == "transcript" {
        tss = ($7 == "+") ? $4 - 1 : $5 - 1
        print $1 "\t" tss "\t" tss + 1
    }' "$GENCODE_GTF" > tmp/tss.unsorted.bed

sort -k1,1 -k2,2n -u tmp/tss.unsorted.bed > tss.bed

# Promoter window = TSS extended by PROMOTER_HALF bp on each side, overlapping
# windows merged into one.
bedtools slop -i tss.bed -g chrom.sizes -b "$PROMOTER_HALF" > tmp/promoters.unmerged.bed
bedtools merge -i tmp/promoters.unmerged.bed > promoters.bed


# =============================================================================
# 3. Regions never used as elements or as background
# =============================================================================

# ENCODE blacklist plus assembly gaps, reduced to three columns and merged.
cat "$BLACKLIST_BED" "$GAPS_BED" > tmp/exclude.raw
cut -f1-3 tmp/exclude.raw > tmp/exclude.3col
awk -v re="$PRIMARY" '$1 ~ re' tmp/exclude.3col > tmp/exclude.primary
sort -k1,1 -k2,2n tmp/exclude.primary > tmp/exclude.sorted
bedtools merge -i tmp/exclude.sorted > exclude.bed


# =============================================================================
# 4. SNAI2 distal peaks and the H3K27ac-gain subset
# =============================================================================

# Peaks as BED4 with their original names.
cut -f1-4 "$PEAKS_BED" > tmp/peaks.4col
awk -v re="$PRIMARY" '$1 ~ re' tmp/peaks.4col > tmp/peaks.primary
sort -k1,1 -k2,2n tmp/peaks.primary > peaks.bed

# The K27ac subset is identified by peak name only; coordinates come from peaks.bed.
cut -f4 "$PEAKS_K27AC_BED" > tmp/k27ac.names
sort -u tmp/k27ac.names > peaks_k27ac.ids

echo "peaks: $(wc -l < peaks.bed)   K27ac subset: $(wc -l < peaks_k27ac.ids)"


# =============================================================================
# 5. Motif occurrences
# =============================================================================

# FIMO tsv columns:
#   1 motif_id  2 alt_id  3 sequence(chrom)  4 start(1-based)  5 stop  6 strand
#   7 score  8 p-value  9 q-value  10 matched_sequence
# Convert to BED (0-based start) and keep the 9-mer word in upper case.
awk -F'\t' -v OFS='\t' -v re="$PRIMARY" '
    NR > 1 && $3 ~ re {
        n++
        print $3, $4 - 1, $5, "fimo_" n, $7, $6, $8, toupper($10)
    }' "$FIMO_TSV" > tmp/motifs.unsorted.bed

sort -k1,1 -k2,2n tmp/motifs.unsorted.bed > tmp/motifs.sorted.bed

# Drop occurrences that touch the blacklist or a gap.
bedtools intersect -v -a tmp/motifs.sorted.bed -b exclude.bed > motifs_all.bed

# Bound = inside a peak.
bedtools intersect -u -a motifs_all.bed -b peaks.bed > motifs_bound.bed

# Unbound distal = outside every peak and outside every promoter window.
bedtools intersect -v -a motifs_all.bed -b peaks.bed > tmp/motifs.unbound.bed
bedtools intersect -v -a tmp/motifs.unbound.bed -b promoters.bed > motifs_unbound_distal.bed

echo "motifs: $(wc -l < motifs_all.bed) total, $(wc -l < motifs_bound.bed) in peaks, $(wc -l < motifs_unbound_distal.bed) unbound distal"


# =============================================================================
# 6. Open-chromatin pool for the peak background
# =============================================================================

# CNCC accessible regions, renamed atac_1 ... atac_n, then stripped of anything
# that overlaps a promoter window, an excluded region or a SNAI2 peak.
awk -v OFS='\t' -v re="$PRIMARY" '
    $1 ~ re {
        n++
        print $1, $2, $3, "atac_" n
    }' "$ATAC_PEAKS" > tmp/atac.unsorted.bed

sort -k1,1 -k2,2n tmp/atac.unsorted.bed > tmp/atac.sorted.bed
bedtools intersect -v -a tmp/atac.sorted.bed -b promoters.bed > tmp/atac.nopromoter.bed
bedtools intersect -v -a tmp/atac.nopromoter.bed -b exclude.bed > tmp/atac.noexclude.bed
bedtools intersect -v -a tmp/atac.noexclude.bed -b peaks.bed > atac_pool.bed

echo "open-chromatin pool: $(wc -l < atac_pool.bed)"


# =============================================================================
# 7. Clean up
# =============================================================================

rm -rf tmp
