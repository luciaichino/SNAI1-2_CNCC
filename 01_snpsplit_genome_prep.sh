#!/bin/bash
set -euo pipefail
source /etc/profile.d/modules.sh

# =============================================================================
# Human vs. chimp ChIP-seq read splitting -- STEP 1: GENOME PREPARATION
# =============================================================================
# Hybrid human x chimp cells (H20961 x C3649) were ChIP-seq'd. Each individual
# is diploid (hap1 + hap2), and the SNPs distinguishing every haplotype from
# every other haplotype are known.
#
# SNPsplit only compares 2 haplotypes/strains at a time, so human-vs-chimp
# read splitting is done in 2 independent rounds and combined afterwards
# (the combining happens in the second script):
#   Round 1: H20961 hap1 (human) vs C3649 hap1 (chimp)
#   Round 2: H20961 hap2 (human) vs C3649 hap2 (chimp)
# In both rounds SNPsplit's "genome1" = human, "genome2" = chimp.
#
# This script builds the 2 SNP VCFs and the 2 matching bowtie2-indexed,
# N-masked dual-hybrid genomes needed for those 2 rounds.
#

# =============================================================================

module load bowtie2
module load bcftools
module load bedtools
module load samtools

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
RAW_VCF_DIR=/.../... # path to VCF files H20961.hg38.hap1.vcf.gz and  C3649.hg38.hap1.vcf.gz
WORK_DIR=/.../... 
REF_GENOME=/.../... # path to hg38 reference fasta file for SNPsplit_genome_preparation
SNPSPLIT_DIR=${WORK_DIR}/...

mkdir -p "$WORK_DIR"
cd "$WORK_DIR"

# ---------------------------------------------------------------------------
# FUNCTION: prepare_diploid_snp_vcf <vcf_A> <vcf_B> <output_prefix>
#
# Builds the SNPsplit-ready VCF for a pair of haploid, single-sample VCFs.
#   1. bcftools merge   -> combine the 2 haploid VCFs into one 2-sample VCF
#   2. keep GT only     -> drop FORMAT/PL and FORMAT/AD, SNPsplit doesn't need them
#   3. diff sites only  -> keep positions where the 2 samples' genotypes differ
#                          (the only informative SNPs for splitting reads)
#   4. haploid->diploid -> recode "1" as 1/1 and anything else as 0/0, so each
#                          haplotype looks like a homozygous diploid sample
#                          (the format SNPsplit_genome_preparation expects)
#   5. add dummy FI     -> SNPsplit_genome_preparation requires an FI FORMAT
#                          field per sample; hardcoded to "1" (pass) throughout
# ---------------------------------------------------------------------------
prepare_diploid_snp_vcf() {
    local vcf_a="$1"
    local vcf_b="$2"
    local prefix="$3"

    # 1 -- merge
    bcftools merge \
      "$vcf_a" \
      "$vcf_b" \
      -Oz -o "${prefix}_merged.hg38.vcf.gz"
    bcftools index "${prefix}_merged.hg38.vcf.gz"

    # 2 -- keep GT only
    bcftools view -Ou "${prefix}_merged.hg38.vcf.gz" \
      | bcftools annotate -x FORMAT/PL,FORMAT/AD \
      | bcftools view -Oz -o "tmp_GTonly_${prefix}.vcf.gz"
    bcftools index "tmp_GTonly_${prefix}.vcf.gz"

    # 3 -- keep only sites where the 2 haplotypes differ
    bcftools view "tmp_GTonly_${prefix}.vcf.gz" \
      | awk 'BEGIN{OFS="\t"}
        /^#/ {print; next}
        {
          split($10,a,":");
          split($11,b,":");
          if(a[1] != b[1]) print
        }' \
      | bgzip > "${prefix}_diff.hg38.vcf.gz"
    bcftools index "${prefix}_diff.hg38.vcf.gz"

bcftools norm -m -any "${prefix}_diff.hg38.vcf.gz" \
  -Oz -o "${prefix}_biallelic.hg38.vcf.gz"
bcftools index "${prefix}_biallelic.hg38.vcf.gz"


    # 4 -- haploid GT -> diploid homozygous GT
    bcftools view "${prefix}_biallelic.hg38.vcf.gz" \
      | awk 'BEGIN{OFS="\t"}
        /^#/ {print; next}
        {
          split($10,a,":");
          split($11,b,":");

          if(a[1]=="1") $10="1/1";
          else $10="0/0";

          if(b[1]=="1") $11="1/1";
          else $11="0/0";

          $9="GT";

          print
        }' \
      | bgzip > "${prefix}_SNPsplit_ready.hg38.vcf.gz"
    bcftools index "${prefix}_SNPsplit_ready.hg38.vcf.gz"



    # 5 -- add dummy FI field
    zcat "${prefix}_SNPsplit_ready.hg38.vcf.gz" \
      | awk 'BEGIN{OFS="\t"}
        /^##/ {print; next}
        /^#CHROM/ {
            print "##FORMAT=<ID=FI,Number=1,Type=Integer,Description=\"Whether a sample was a Pass(1) or fail (0) based on FILTER values\">"
            print
            next
        }
        {
            # Append :FI to FORMAT column
            $9 = $9 ":FI"

            # Append :1 to all sample columns
            for(i=10; i<=NF; i++){
                $i = $i ":1"
            }
            print
        }' \
      | bgzip > "${prefix}_SNPsplit_ready_FIadded.hg38.vcf.gz"
    tabix -p vcf "${prefix}_SNPsplit_ready_FIadded.hg38.vcf.gz"
}

# ---------------------------------------------------------------------------
# Build the 2 SNP VCFs actually needed downstream
# ---------------------------------------------------------------------------
prepare_diploid_snp_vcf \
  "${RAW_VCF_DIR}/H20961.hg38.hap1.vcf.gz" \
  "${RAW_VCF_DIR}/C3649.hg38.hap1.vcf.gz" \
  "H20961_hap1_vs_C3649_hap1"

prepare_diploid_snp_vcf \
  "${RAW_VCF_DIR}/H20961.hg38.hap2.vcf.gz" \
  "${RAW_VCF_DIR}/C3649.hg38.hap2.vcf.gz" \
  "H20961_hap2_vs_C3649_hap2"

# ---------------------------------------------------------------------------
# SNPsplit genome preparation (dual hybrid mode) -- run on raccoon
# Requires SNPsplit (conda install -c bioconda snpsplit)
# ---------------------------------------------------------------------------
mkdir -p "${SNPSPLIT_DIR}/H20961_C3649_hap1_genomeprep"
cd "${SNPSPLIT_DIR}/H20961_C3649_hap1_genomeprep"
SNPsplit_genome_preparation \
  --genome_build hg38 \
  --dual_hybrid \
  --vcf_file "${WORK_DIR}/H20961_hap1_vs_C3649_hap1_SNPsplit_ready_FIadded.hg38.vcf.gz" \
  --reference_genome "$REF_GENOME" \
  --strain H20961.hg38.hap1.bam \
  --strain2 C3649.hg38.hap1.bam

mkdir -p "${SNPSPLIT_DIR}/H20961_C3649_hap2_genomeprep"
cd "${SNPSPLIT_DIR}/H20961_C3649_hap2_genomeprep"
SNPsplit_genome_preparation \
  --genome_build hg38 \
  --dual_hybrid \
  --vcf_file "${WORK_DIR}/H20961_hap2_vs_C3649_hap2_SNPsplit_ready_FIadded.hg38.vcf.gz" \
  --reference_genome "$REF_GENOME" \
  --strain H20961.hg38.hap2.bam \
  --strain2 C3649.hg38.hap2.bam

# ---------------------------------------------------------------------------
# Build bowtie2 indices from the N-masked genomes
# ---------------------------------------------------------------------------


cd "${SNPSPLIT_DIR}/H20961_C3649_hap1_genomeprep/H20961.hg38.hap1.bam_C3649.hg38.hap1.bam_dual_hybrid.based_on_hg38_N-masked"
bowtie2-build *N-masked.fa hg38_dual_masked_H20961_C3649_hap1

cd "${SNPSPLIT_DIR}/H20961_C3649_hap2_genomeprep/H20961.hg38.hap2.bam_C3649.hg38.hap2.bam_dual_hybrid.based_on_hg38_N-masked"
bowtie2-build *N-masked.fa hg38_dual_masked_H20961_C3649_hap2

echo "Genome prep done. Next: map ChIP-seq reads to the 2 bowtie2 indices above, then run the SNPsplit script."
