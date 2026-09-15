#!/usr/bin/env bash
# Paths for the selection analysis pipeline. Edit this file only.

#=================================
# project inputs
#=================================
PROJ=./   # This is the base directory that contains all input files and output results
PEAKS_BED=$PROJ/peaks/s8_lib1_CNCCp4_SNAI2a3029_DMSO_r2_aSNAI2_reltodTag_peaks_minusINPUT_filtered_DISTALcomb_v7.bed
PEAKS_K27AC_BED=$PROJ/peaks/s8_lib1_CNCCp4_SNAI2a3029_DMSO_r2_aSNAI2_reltodTag_peaks_minusINPUT_filtered_DISTALcomb_v7_overlappingK27acUPregions3hor1d.bed
ATAC_PEAKS=$PROJ/peaks/human_cncc_hmmratac_accessible_regions.narrowPeak
FIMO_TSV=$PROJ/motifs/SNAI2_DBD_jolma2013_fimo/fimo.tsv

#=================================
# Public resources 
#=================================
# Downloadable from the web. 
# Please check the README.md for instructions on how to download these files.
GENOME_FA=$PROJ/data/genomes/hg38.fa
CHROM_SIZES=$PROJ/data/genomes/hg38.chrom.sizes
GAPS_BED=$PROJ/data/annotation/hg38.gap.bed
BLACKLIST_BED=$PROJ/data/annotation/hg38-blacklist.v2.bed
GENCODE_GTF=$PROJ/data/annotation/gencode.v49.basic.annotation.gtf
PHASTCONS_BW=$PROJ/data/tracks/phastCons470way/hg38.phastCons470way.bw
LINSIGHT_BW=$PROJ/data/tracks/LINSIGHT.bw
CHAIN_38TO19=$PROJ/data/chains/hg38ToHg19.over.chain.gz
CHAIN_19TO38=$PROJ/data/chains/hg19ToHg38.over.chain.gz

#=================================
# parameters
#=================================
N_DRAWS_ATAC=200      # covariate-matched open-chromatin draws
N_DRAWS_KMER=100      # word-matched unbound-motif draws
N_DRAWS_NESTED=1000   # subset-versus-rest resampling
PROMOTER_HALF=500     # promoter window = TSS +/- this many bp

#=================================
# project outputs
#=================================
ROOT=$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)
INT=$ROOT/intermediate
RES=$ROOT/results
FIG=$ROOT/figures


mkdir -p "$INT" "$RES" "$FIG"
export PROJ PEAKS_BED PEAKS_K27AC_BED FIMO_TSV ATAC_PEAKS GENOME_FA CHROM_SIZES GAPS_BED \
       BLACKLIST_BED GENCODE_GTF PHASTCONS_BW LINSIGHT_BW CHAIN_38TO19 CHAIN_19TO38 \
       N_DRAWS_ATAC N_DRAWS_KMER N_DRAWS_NESTED PROMOTER_HALF ROOT INT RES FIG
