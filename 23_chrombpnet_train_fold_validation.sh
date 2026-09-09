#!/bin/bash

# Source conda environment

source ~/anaconda3/etc/profile.d/conda.sh

# activate chrombpnet

conda activate chrombpnet

# Specify which GPU to use (e.g., GPU 0, or GPU 1, etc.)
export CUDA_VISIBLE_DEVICES=0  # Use GPU 0

# Train bias models

chrombpnet bias pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_1.json \
    -b 0.5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_1/bias/ \
    -fp hg38_fold1

chrombpnet bias pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_2.json \
    -b 0.5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_2/bias/ \
    -fp hg38_fold2

chrombpnet bias pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_3.json \
    -b 0.5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_3/bias/ \
    -fp hg38_fold3

chrombpnet bias pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_4.json \
    -b 0.5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_4/bias/ \
    -fp hg38_fold4

chrombpnet bias pipeline \
    -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
    -d "ATAC" \
    -g ~/analysis/genomes/hg38.fa \
    -c ~/analysis/genomes/hg38.chrom.sizes \
    -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
    -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
    -fl ~/analysis/genomes/fold_5.json \
    -b 0.5 \
    -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_5/bias/ \
    -fp hg38_fold5

# Train single task ChromBPNet models

chrombpnet pipeline \
   -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
   -d "ATAC" \
   -g ~/analysis/genomes/hg38.fa \
   -c ~/analysis/genomes/hg38.chrom.sizes \
   -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
   -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
   -fl ~/analysis/genomes/fold_1.json \
   -b ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_1/bias/models/hg38_fold1_bias.h5 \
   -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_1/all_peaks_fold_1

chrombpnet pipeline \
   -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
   -d "ATAC" \
   -g ~/analysis/genomes/hg38.fa \
   -c ~/analysis/genomes/hg38.chrom.sizes \
   -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
   -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
   -fl ~/analysis/genomes/fold_2.json \
   -b ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_2/bias/models/hg38_fold2_bias.h5 \
   -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_2/all_peaks_fold_2

chrombpnet pipeline \
   -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
   -d "ATAC" \
   -g ~/analysis/genomes/hg38.fa \
   -c ~/analysis/genomes/hg38.chrom.sizes \
   -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
   -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
   -fl ~/analysis/genomes/fold_3.json \
   -b ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_3/bias/models/hg38_fold3_bias.h5 \
   -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_3/all_peaks_fold_3

chrombpnet pipeline \
   -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
   -d "ATAC" \
   -g ~/analysis/genomes/hg38.fa \
   -c ~/analysis/genomes/hg38.chrom.sizes \
   -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
   -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
   -fl ~/analysis/genomes/fold_4.json \
   -b ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_4/bias/models/hg38_fold4_bias.h5 \
   -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_4/all_peaks_fold_4

  chrombpnet pipeline \
     -ibam ~/analysis/lucia/bam/ATAC_SNAI12_combined_DMSO_merged.bam \
     -d "ATAC" \
     -g ~/analysis/genomes/hg38.fa \
     -c ~/analysis/genomes/hg38.chrom.sizes \
     -p ~/analysis/lucia/narrowpeak/20260607_cncc_training_atac_peaks.narrowpeak \
     -n ~/analysis/lucia/narrowpeak/NONPEAK_20260607_cncc_training_atac_peaks_negatives.bed \
     -fl ~/analysis/genomes/fold_5.json \
     -b ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_5/bias/models/hg38_fold5_bias.h5 \
     -o ~/analysis/lucia/chrombpnet_model/for_paper/snail_modisco_test/fold_validation/fold_5/all_peaks_fold_5
