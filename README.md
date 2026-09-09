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

--- 04_SNAI_RNAseq_DESeq2_and_downstream.Rmd

R markdown with processing of RNA-seq data 

--- 05_DESeq2_ChIP_ATAC.Rmd

R markdown with DEseq of H3K27ac and ATAC data, and MA plots

--- 06_Define_proximal_distal_v6.Rmd

R markdown including processing of data to classify regions as promoter proximal or distal and make panels shown in Figures 2 and 3 and related supplemental figures.
The data file containing all the processed information has been uploaded to this repo (data_v4.rds).

--- 07_Mixed_effect_model.Rmd

ATAC vs H3K27ac log2FC across dTAG time points — linear mixed-effects model

## Figure 4: TWIST1, NR2F1, TFAP2A ChIP analysis

--- 08_Fig4_analysis_TFs_motifsdist_v5_allpeaks.Rmd

R markdown including analysis to generate the kernel smoothened plots of TF activators ChIP-seq foldchanges relative to distance between activator motif and SNAI1/2 motif (Figure 4).

--- 09_Fig4_linear_model_v2.Rmd

R markdown for linear model predicting H3K27ac change based on TF binding changes (DMSO vs dTag)

## Figure 5 and 7: Fiber-seq and DAF-seq analyses
--- 10_Fiber_seq_plots.Rmd

Code for plotting Fiber-seq data (metaplot and heatmap)

--- 11_DAFseq_plots.Rmd

Nucleosome clustering analysis and TF footprint analysis DAFseq data

## Figure 6: DAF-seq codependency analysis

--- 12_DAFseq_codependency.Rmd

codependency analysis enhancers 2 and 3 DAF-seq data


## ChromBPnet analysis (Figures 1, S1, 4, S3)

20_curate_cncc_regions.Rmd

21_generate_nonpeaks.sh

22_chrombpnet_train_snail.sh

23_chrombpnet_train_fold_validation.sh

24_convert_chrombpnet_h5_to_npz.py

25_run_conversion.sh

26_run_modisco.py

27_modisco.sh

28_finemo.sh

29_extract_cwms.py

30_extract_cwms.sh

31_curate_motifs.Rmd

32_custom_sequence_generation.Rmd

33_predict_motif_effects.ipynb


## Evolutionary constraint (Figure S1)


