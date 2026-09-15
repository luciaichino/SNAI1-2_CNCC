# Evolutionary constraint of SNAI2 distal peaks and motif occurrences

Self-contained pipeline that produces `figures/fig_final_comparison.pdf` (and `.png`)
from raw inputs: peak calls, a genome-wide FIMO motif scan, CNCC open-chromatin
regions, and two public constraint tracks (phastCons 470-way, LINSIGHT).

The FIMO scans need to be downloaded from here (fimo.tsv):
https://github.com/luciaichino/SNAI1-2_CNCC/releases/tag/SNAI2_DBD_jolma2013_fimo


```
config.sh              paths to inputs and tools; edit this file only
run_all.sh             runs the eight steps in order
scripts/1 ... 8        one step each (bash + R)
intermediate/          element sets, backgrounds, hg19 lifts, per-element scores
results/               statistics tables
figures/               fig_final_comparison.pdf / .png
```

## What the figure shows

Two element types, each as the full set and as the subset overlapping enhancers that
gain H3K27ac after SNAI2 induction, each against its own matched background:

| bar | set | matched background |
|---|---|---|
| red | all SNAI2 distal ChIP-seq peaks (3,877) | other open CNCC regions of the same width, GC, CpG, repeat content and TSS distance; 200 draws |
| green | the 834 peaks overlapping H3K27ac-gain enhancers | the same draws, restricted to the regions paired to those 834 peaks |
| purple | SNAI2 motif occurrences inside peaks, FIMO p < 1e-4 (4,493) | unbound genomic copies of the identical 9-mer word; 100 draws |
| yellow | the 982 occurrences inside the 834 peaks | the same rule applied at the subset's word composition |

Panel a: fraction of aligned bases with phastCons470way >= 0.8. Panel b: mean LINSIGHT.

## Inputs (glossary)

Set every path in `config.sh`. Project-specific files are listed first, public downloads after.

`PROJ` in `config.sh` is the directory holding the inputs, and defaults to the pipeline
directory itself. The other paths are written relative to it, so placing the files as
below means only `PROJ` has to be set:

```
$PROJ/peaks/       PEAKS_BED, PEAKS_K27AC_BED, ATAC_PEAKS
$PROJ/motifs/      FIMO_TSV
$PROJ/data/genomes/     GENOME_FA, CHROM_SIZES
$PROJ/data/annotation/  GAPS_BED, BLACKLIST_BED, GENCODE_GTF
$PROJ/data/tracks/      PHASTCONS_BW (in phastCons470way/), LINSIGHT_BW
$PROJ/data/chains/      CHAIN_38TO19, CHAIN_19TO38
```

| variable | file | what it provides | source |
|---|---|---|---|
| `PEAKS_BED` | `peaks/..._DISTALcomb_v7.bed` | 3,877 SNAI2 distal ChIP-seq peaks, hg38, BED6 with unique names | this study |
| `PEAKS_K27AC_BED` | `peaks/..._overlappingK27acUPregions3hor1d.bed` | the 834 of those peaks overlapping regions gaining H3K27ac at 3 h or 1 d of dTAG treatment; same names as `PEAKS_BED` | this study |
| `FIMO_TSV` | `motifs/SNAI2_DBD_jolma2013_fimo/fimo.tsv` | genome-wide FIMO scan of the Jolma 2013 SNAI2 DBD matrix against hg38 at the default p < 1e-4 | this study; regenerate with `fimo --o OUT --max-stored-scores 10000000 SNAI2_DBD_jolma2013.meme hg38.fa` (MEME suite 5.1.1) |
| `ATAC_PEAKS` | `human_cncc_hmmratac_accessible_regions.narrowPeak` | all accessible regions in the same CNCC cells (HMMRATAC calls); the pool for the peak background | this study |
| `GENOME_FA` | `hg38.fa` | soft-masked hg38 primary assembly (lower case = RepeatMasker/TRF) for GC, CpG and repeat content | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz |
| `CHROM_SIZES` | `hg38.chrom.sizes` | chromosome lengths | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes |
| `GAPS_BED` | `hg38.gap.bed` | assembly gaps (excluded from every element set) | UCSC table `gap`, e.g. `wget -qO- https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz \| zcat \| cut -f2-4` |
| `BLACKLIST_BED` | `hg38-blacklist.v2.bed` | ENCODE blacklist regions (excluded) | https://github.com/Boyle-Lab/Blacklist/raw/master/lists/hg38-blacklist.v2.bed.gz |
| `GENCODE_GTF` | `gencode.v49.basic.annotation.gtf` | transcript start sites for promoter windows (TSS +/- 500 bp) and TSS distance | https://www.gencodegenes.org/human/ (basic annotation, chr names as `chr1`) |
| `PHASTCONS_BW` | `hg38.phastCons470way.bw` | per-base phastCons over the 470-mammal Zoonomia alignment, hg38 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/phastCons470way/hg38.phastCons470way.bw |
| `LINSIGHT_BW` | `LINSIGHT.bw` | per-base LINSIGHT score, hg19 (non-commercial licence) | https://hgdownload.soe.ucsc.edu/gbdb/hg19/linsight/LINSIGHT.bw |
| `CHAIN_38TO19` | `hg38ToHg19.over.chain.gz` | liftOver chain hg38 -> hg19 | https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToHg19.over.chain.gz |
| `CHAIN_19TO38` | `hg19ToHg38.over.chain.gz` | liftOver chain hg19 -> hg38 (for the reciprocity check) | https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/hg19ToHg38.over.chain.gz |

## Software

- bedtools >= 2.30
- UCSC utilities: `bigWigAverageOverBed`, `bigWigToBedGraph`, `liftOver`
  (https://hgdownload.soe.ucsc.edu/admin/exe/)
- R >= 4.2 with `data.table`, `ggplot2`, `patchwork`, `GenomicRanges`, `nullranges`
  (the last two from Bioconductor: `BiocManager::install(c("GenomicRanges","nullranges"))`)

## Running

1. Edit `config.sh`.
2. Run `source config.sh`.
3. Run `bash run_all.sh`.

Steps 1-6 need the genome, tracks and chains and takes roughly 20 minutes, because these steps are dominated by scoring the background elements. 
Steps 7-8 run in seconds from the `intermediate/` files and can be rerun on their own, after `source config.sh`.

The random seeds are fixed in the scripts, so a rerun reproduces the tables and the
figure exactly on the same software versions.

## Step summary

| step | script | output |
|---|---|---|
| 1 | `1_define_elements.sh` | peaks, K27ac subset ids, motif occurrences (bound / unbound distal), CNCC open-chromatin pool |
| 2 | `2_covariates.sh` | width, GC, CpG, repeat fraction, N fraction, TSS distance for peaks and pool |
| 3 | `3_match_atac_background.R` | 200 width-stratified, covariate-matched draws from the pool; each element carries the index of the peak it is paired to |
| 4 | `4_match_kmer_background.R` | 100 draws of unbound occurrences matched word-for-word to the bound occurrences |
| 5 | `5_liftover_hg19.sh` | reciprocal hg38 -> hg19 -> hg38 liftOver of every set for LINSIGHT |
| 6 | `6_score.sh` | per-element phastCons470way mean, aligned bases, bases >= 0.8; per-element LINSIGHT mean |
| 7 | `7_statistics.R` | draw z-tests for the eight bars, nested subset-versus-rest tests, covariate balance |
| 8 | `8_figure.R` | `fig_final_comparison.pdf` / `.png` |

All statistics keep only elements with at least half of their bases aligned in the
track (`covered / size >= 0.5`).
