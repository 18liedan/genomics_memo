#!/usr/bin/env bash
set -euo pipefail

#This script automates variant calling from high quality samples, to gather reliable sites to be later used in BQSR.

###############################################################################
# RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="yourspecies"

# RESOURCE ALLOCATION (Optimized for NIG supercomputer; 1.5TB RAM / 192 Cores)
MAX_JOBS=12              # 12 samples in parallel
THREADS_PER_JOB=16       # 12 * 16 = 192 cores
JAVA_OPTS="-Xmx64g"      # 12 * 64GB = 768GB
COHORT_JAVA_OPTS="-Xmx256g" 

# PATHS
REF="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa" #your reference genome
SAMPLE_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples_hq.txt" #your list of high coverage, reliable samples
MASK_BED="${SPECIES_ID}_ref/${SPECIES_ID}_ref_masked_regions.bed" #your repeat masked sites
# Coverage TSV (Assumes format: SampleID [tab] CoverageValue)
COVERAGE_STATS="${SPECIES_ID}_stats/${SPECIES_ID}_coverage_summary.tsv"

BAM_DIR="${SPECIES_ID}_mapped" #directory with mapped and processed genomes
GVCF_DIR="${SPECIES_ID}_gvcf" #directory to save per-sample vcf.gz outputs from haplotypecaller
VCF_DIR="${SPECIES_ID}_vcf" #directory to save filtered vcf.gz
SAMPLE_VCF_DIR="${VCF_DIR}/individual_samples"

mkdir -p "$GVCF_DIR" "$VCF_DIR" "$SAMPLE_VCF_DIR"

# Print date on terminal log outputs
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }

###############################################################################
# 1. PARALLEL PER-SAMPLE VARIANT CALLING & DEPTH FILTERING
###############################################################################
# Export variables for GNU Parallel
export REF GVCF_DIR BAM_DIR SAMPLE_VCF_DIR THREADS_PER_JOB JAVA_OPTS COVERAGE_STATS

process_sample_to_filtered_vcf() {
    local sample="$1"
    local INPUT_BAM="${BAM_DIR}/${sample}_RGdup.bam"
    local OUTPUT_GVCF="${GVCF_DIR}/${sample}.g.vcf.gz"
    local DP_FILTERED_VCF="${SAMPLE_VCF_DIR}/${sample}.dp_filtered.g.vcf.gz"

    # Step A: HaplotypeCaller (Skip if DP filtered VCF already exists)
    if [[ ! -f "$DP_FILTERED_VCF" ]]; then
        
        # 1. HaplotypeCaller
        if [[ ! -f "$OUTPUT_GVCF" ]]; then
            gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
                -R "$REF" -I "$INPUT_BAM" -O "$OUTPUT_GVCF" -ERC GVCF \
                --native-pair-hmm-threads "$THREADS_PER_JOB"
        fi

        # 2. Depth Filtering (1/3 depth and 2*depth)
        # Parses the TSV produced after mapping
        DEPTH=$(grep -w "$sample" "$COVERAGE_STATS" | awk '{print $4}')
        if [[ -z "$DEPTH" ]]; then
            echo "Error: Depth for $sample not found in $COVERAGE_STATS" >&2
            return 1
        fi
        
        MIN_DP=$(echo "scale=2; $DEPTH / 3" | bc)
        MAX_DP=$(echo "scale=2; $DEPTH * 2" | bc)
        
        echo "Filtering $sample: Depth=$DEPTH (Min=$MIN_DP, Max=$MAX_DP)"
        
        vcftools --gzvcf "$OUTPUT_GVCF" --min-meanDP "$MIN_DP" --max-meanDP "$MAX_DP" --recode --stdout | \
            bgzip -c > "$DP_FILTERED_VCF"
        gatk IndexFeatureFile -I "$DP_FILTERED_VCF"
        
    else
        echo "Skipping $sample: Depth-filtered VCF already exists."
    fi
}

export -f process_sample_to_filtered_vcf

log "Step 1: Running Parallel HaplotypeCalling and Per-Sample Depth Filtering..."
tr -d '\r' < "$SAMPLE_LIST" | parallel --jobs "$MAX_JOBS" --progress process_sample_to_filtered_vcf {}

###############################################################################
# 2. COHORT MERGING
###############################################################################
MERGED_VCF="${VCF_DIR}/${SPECIES_ID}_merged.vcf.gz"
GENOTYPED_VCF="${VCF_DIR}/${SPECIES_ID}_genotyped.vcf.gz"

# 3a. Merging
if [[ ! -f "$MERGED_VCF" ]]; then
    log "Step 2: Merging per-sample depth-filtered VCFs into cohort..."
    VCF_LIST=()
    while read -r sample; do
        sample=$(echo "$sample" | tr -d '\r')
        [[ -f "${SAMPLE_VCF_DIR}/${sample}.dp_filtered.g.vcf.gz" ]] && VCF_LIST+=("${SAMPLE_VCF_DIR}/${sample}.dp_filtered.g.vcf.gz")
    done < "$SAMPLE_LIST"

    bcftools merge --threads 32 -O z -o "$MERGED_VCF" "${VCF_LIST[@]}"
    gatk IndexFeatureFile -I "$MERGED_VCF"
fi

# 3b. Genotyping
if [[ ! -f "$GENOTYPED_VCF" ]]; then
    log "Step 2b: Running GATK GenotypeGVCFs on merged cohort..."
    gatk --java-options "${COHORT_JAVA_OPTS}" GenotypeGVCFs \
        -R "$REF" \
        -V "$MERGED_VCF" \
        -O "$GENOTYPED_VCF"
else
    log "Step 2b: Genotyped cohort VCF already exists. Skipping."
fi

###############################################################################
# 3. FILTERING & EXTRACTION (SNPs and Indels)
###############################################################################
log "Step 3: Starting Filtering Pipeline..."

# 4a: SNP Hard Filtering
SNPS_HARD="${VCF_DIR}/${SPECIES_ID}_hardfilteredpass.snps.vcf.gz"
if [[ ! -f "$SNPS_HARD" ]]; then
    log "  - Extracting and Filtering SNPs"
    gatk SelectVariants -V "$GENOTYPED_VCF" -select-type SNP -O "${VCF_DIR}/temp_snps.vcf.gz"
    gatk VariantFiltration \
        -V "${VCF_DIR}/temp_snps.vcf.gz" \
        -filter "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
        --filter-name "snp_hard_filter" -O "${VCF_DIR}/temp_snps_filt.vcf.gz"
    bcftools view -f PASS "${VCF_DIR}/temp_snps_filt.vcf.gz" -Oz -o "$SNPS_HARD"
    gatk IndexFeatureFile -I "$SNPS_HARD"
fi

# 4b: INDEL Extraction and Filtering
INDELS_HARD="${VCF_DIR}/${SPECIES_ID}_hardfilteredpass.indels.vcf.gz"
if [[ ! -f "$INDELS_HARD" ]]; then
    log "  - 4b: Extracting and Filtering Indels"
    gatk SelectVariants -R "$REF" -V "$GENOTYPED_VCF" -select-type INDEL -O "${VCF_DIR}/temp_indels.vcf.gz"
    gatk VariantFiltration -R "$REF" -V "${VCF_DIR}/temp_indels.vcf.gz" -O "${VCF_DIR}/temp_indels_filt.vcf.gz" \
        -filter 'QD < 2.0' --filter-name 'QD2' \
        -filter 'QUAL < 30.0' --filter-name 'QUAL30' \
        -filter 'SOR > 4.0' --filter-name 'SOR4' \
        -filter 'FS > 60.0' --filter-name 'FS60'
    bcftools view -f PASS "${VCF_DIR}/temp_indels_filt.vcf.gz" -Oz -o "$INDELS_HARD"
    gatk IndexFeatureFile -I "$INDELS_HARD"
fi

# 4c: Softfiltering (Missingness/MAF)
SNPS_SOFT="${VCF_DIR}/${SPECIES_ID}_filtered.snps.vcf.gz"
if [[ ! -f "$SNPS_SOFT" ]]; then
    log "  - 4c: Softfiltering (Max-missing 0.85, MAF 0.05)"
    vcftools --gzvcf "$SNPS_HARD" --max-missing 0.85 --maf 0.05 --recode --stdout | bgzip -@ 16 -c > "$SNPS_SOFT"
    gatk IndexFeatureFile -I "$SNPS_SOFT"
fi

# 4d: Filter out repeat masker sites
NOREPEAT_VCF="${VCF_DIR}/${SPECIES_ID}_norepeat_filtered.snps.vcf.gz"
if [[ ! -f "$NOREPEAT_VCF" ]]; then
    if [[ -f "$MASK_BED" ]]; then
        log "  - 4d: Excluding Masked Regions"
        vcftools --gzvcf "$SNPS_SOFT" --exclude-bed "$MASK_BED" --recode --stdout | bgzip -@ 16 -c > "$NOREPEAT_VCF"
        gatk IndexFeatureFile -I "$NOREPEAT_VCF"
    else
        log "  - Warning: Mask BED not found. Skipping Step 4d."
    fi
fi

###############################################################################
# 4. VARIANT SUMMARY TABLE
###############################################################################
log "Final Step: Generating Variant Summary..."
SUMMARY_OUT="${SPECIES_ID}_stats/${SPECIES_ID}_variant_summary.txt"

{
    echo -e "STAGE\tCOUNT"
    echo -ne "Raw_Genotyped_VCF (DP_Filtered)\t"; bcftools view -H "$GENOTYPED_VCF" | wc -l
    echo -ne "HardFiltered_SNPs\t"; bcftools view -H "$SNPS_HARD" | wc -l
    echo -ne "HardFiltered_Indels\t"; bcftools view -H "$INDELS_HARD" | wc -l
    echo -ne "SoftFiltered_SNPs\t"; bcftools view -H "$SNPS_SOFT" | wc -l
    if [[ -f "$NOREPEAT_VCF" ]]; then
        echo -ne "Final_NoRepeat_SNPs\t"; bcftools view -H "$NOREPEAT_VCF" | wc -l
    fi
} > "$SUMMARY_OUT"

cat "$SUMMARY_OUT"

echo "Pre-BQSR variant calling and filtering for $SPECIES_ID completed. Variant stats exported to $SUMMARY_OUT"
