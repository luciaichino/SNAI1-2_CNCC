# Enhancer Tuning by Sequence-Dependent Repressors SNAI1 and SNAI2

Code accompanying the analyses in paper **link TBD

All new sequencing data have been deposited to GEO ** and SRA **

## Figure 1: human / chimp tetraployd hybrid analysis

--- 01_snpsplit_genome_prep.sh
Builds the SNP VCFs and N-masked,bowtie2-indexed dual-hybrid genomes needed for read splitting.

After this, ChIP-seq reads were mapped to the two bowtie2 indices H20961.hap1.C3649.hap1 (masked for human and chimp hap1 SNPs) and H20961.hap2.C3649.hap2 (masked for human and chimp hap2 SNPs)

--- 02_snpsplit_run_and_merge.sh
This script splits each sample's BAM into human and chimp reads, merges them by species, removes duplicates, converts to fastq.

After this, reads were mapped to the two bowtie2 indices H20961.hap1.H20961.hap2 (masked for all human SNPs) and C3649.hap1.C3649.hap2 (masked for all chimp SNPs)
Counts at CNCC regulatory regions were calculated with bedtools coverage

--- 03_human_chimp_ChIP_motifs_correlation.Rmd
R markdown of correlation analysis between human and chimpanzee ChIP changes and SNAI2 motif changes

## Figure 2/3: SNAI1/2, H3K27ac, ATAC, analyses

--- 04_Define_proximal_distal_v5.Rmd
R markdown including processing of data to classify regions as promoter proximal or distal and make panels shown in Figures 2 and 3 and related supplemental figures.
The data file containing all the processed informaion has been uploaded to this repo (data_v4.rds).

## Figure 4: TWIST1, NR2F1, TFAP2A ChIP analysis (motif distances)

--- 05_Fig4_analysis_TFs_motifsdist_v5_allpeaks.Rmd
R markdown including analysis to generate the kernel smoothened plots of TF activators ChIP-seq foldchanges relative to distance between motifs (Figure 4).



