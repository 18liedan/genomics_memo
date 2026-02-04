#!/usr/bin/env bash
set -euo pipefail

# This script automates calling of reliable variants from high quality samples to be used in BQSR.

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
REF="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa"
SAMPLE_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples_hq.txt"
MASK_BED="${SPECIES_ID}_ref/${SPECIES_ID}_ref_masked_regions.bed"
# Coverage TSV (Assumes format: SampleID [tab] CoverageValue)
COVERAGE_STATS="${SPECIES_ID}_stats/${SPECIES_ID}_coverage_summary.tsv"

BAM_DIR="${SPECIES_ID}_mapped"
GVCF_DIR="${SPECIES_ID}_gvcf"
VCF_DIR="${SPECIES_ID}_vcf"
SAMPLE_VCF_DIR="${VCF_DIR}/individual_samples"

mkdir -p "$GVCF_DIR" "$VCF_DIR" "$SAMPLE_VCF_DIR"

# Print date on terminal log outputs
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }

# Export variables for GNU Parallel
export REF GVCF_DIR BAM_DIR SAMPLE_VCF_DIR THREADS_PER_JOB JAVA_OPTS COVERAGE_STATS

###############################################################################
# 1. PER-SAMPLE VARIANT CALLING (HaplotypeCaller)
###############################################################################
run_haplotype_caller() {
    local sample="$1"
    local INPUT_BAM="${BAM_DIR}/${sample}_RGdup.bam"
    local OUTPUT_GVCF="${GVCF_DIR}/${sample}.g.vcf.gz"

    if [[ ! -f "$OUTPUT_GVCF" ]]; then
        log "Starting HaplotypeCaller for $sample..."
        gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
            -R "$REF" -I "$INPUT_BAM" -O "$OUTPUT_GVCF" -ERC GVCF \
            --native-pair-hmm-threads "$THREADS_PER_JOB"
    else
        echo "Skipping HC: $sample gVCF already exists."
    fi
}
export -f run_haplotype_caller

log "Step 1: Running Parallel HaplotypeCalling..."
tr -d '\r' < "$SAMPLE_LIST" | parallel --jobs "$MAX_JOBS" --progress run_haplotype_caller {}

###############################################################################
# 2. PER-SAMPLE DEPTH FILTERING
###############################################################################

run_depth_filtering() {
    local sample="$1"
    local INPUT_GVCF="${GVCF_DIR}/${sample}.g.vcf.gz"
    local DP_FILTERED_VCF="${SAMPLE_VCF_DIR}/${sample}.dp_filtered.g.vcf.gz"

    if [[ ! -f "$DP_FILTERED_VCF" ]]; then
        # Ensure input exists
        if [[ ! -f "$INPUT_GVCF" ]]; then
            echo "Error: Cannot filter $sample, input $INPUT_GVCF not found." >&2
            return 1
        fi

        # Parses the TSV produced after mapping
        DEPTH=$(grep -w "$sample" "$COVERAGE_STATS" | awk '{print $4}')
        if [[ -z "$DEPTH" ]]; then
            echo "Error: Depth for $sample not found in $COVERAGE_STATS" >&2
            return 1
        fi
        
        # Calculate thresholds
        MIN_DP=$(echo "scale=2; $DEPTH / 3" | bc)
        MAX_DP=$(echo "scale=2; $DEPTH * 2" | bc)
        
        echo "Filtering $sample: Depth=$DEPTH (Min=$MIN_DP, Max=$MAX_DP)"
        
        vcftools --gzvcf "$INPUT_GVCF" --min-meanDP "$MIN_DP" --max-meanDP "$MAX_DP" --recode --stdout | \
            bgzip -c > "$DP_FILTERED_VCF"
        gatk IndexFeatureFile -I "$DP_FILTERED_VCF"
    else
        echo "Skipping Filter: $sample depth-filtered VCF already exists."
    fi
}
export -f run_depth_filtering

log "Step 2: Running Parallel Per-Sample Depth Filtering..."
tr -d '\r' < "$SAMPLE_LIST" | parallel --jobs "$MAX_JOBS" --progress run_depth_filtering {}

###############################################################################
# 3. COHORT MERGING & GENOTYPING
###############################################################################
MERGED_GVCF="${VCF_DIR}/${SPECIES_ID}_merged.vcf.gz"
GENOTYPED_VCF="${VCF_DIR}/${SPECIES_ID}_genotyped.vcf.gz"

# 2a. GATK CombineGVCFs
if [[ ! -f "$MERGED_GVCF" ]]; then
    log "Step 3a: Combining per-sample VCFs with GATK CombineGVCFs..."
    
    # Build the list of -V arguments
    VCF_INPUTS=()
    while IFS= read -r sample || [[ -n "$sample" ]]; do
        sample=${sample//$'\r'/} # strip CR if any
        [[ -z "$sample" ]] && continue
		VCF_FILE="${SAMPLE_VCF_DIR}/${sample}.dp_filtered.g.vcf.gz"
        if [[ -f "$VCF_FILE" ]]; then
            VCF_INPUTS+=("-V" "$VCF_FILE")
        else
			log "ERROR: missing filtered gVCF for $sample: $VCF_FILE"
			missing=1
		fi
    done < "$SAMPLE_LIST"

    gatk --java-options "${COHORT_JAVA_OPTS}" CombineGVCFs \
        -R "$REF" \
        "${VCF_INPUTS[@]}" \
        -O "$MERGED_GVCF"
else
    log "Step 3a: Combined gVCF already exists ($MERGED_GVCF). Skipping."
fi

# 3b. Genotyping
if [[ ! -f "$GENOTYPED_VCF" ]]; then
    log "Step 3b: Running GATK GenotypeGVCFs on merged cohort..."
    gatk --java-options "${COHORT_JAVA_OPTS}" GenotypeGVCFs \
        -R "$REF" \
        -V "$MERGED_GVCF" \
        -O "$GENOTYPED_VCF"
else
    log "Step 3b: Genotyped cohort VCF already exists. Skipping."
fi

###############################################################################
# 4. FILTERING & EXTRACTION (SNPs and Indels)
###############################################################################
log "Step 4: Starting Filtering Pipeline..."

# 4a: SNP Extraction and Hard Filtering
SNPS_HARD="${VCF_DIR}/${SPECIES_ID}_hardfilteredpass.snps.vcf.gz"
if [[ ! -f "$SNPS_HARD" ]]; then
    log "  - Extracting and Filtering SNPs"
    gatk SelectVariants -V "$GENOTYPED_VCF" -select-type SNP -O "${VCF_DIR}/temp_snps.vcf.gz"
    gatk VariantFiltration \
        -V "${VCF_DIR}/temp_snps.vcf.gz" \
        -filter "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || ReadPosRankSum < -8.0" \
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

# 4c: Filter out repeat masker sites
NOREPEAT_VCF="${VCF_DIR}/${SPECIES_ID}_norepeat.snps.vcf.gz"
if [[ ! -f "$NOREPEAT_VCF" ]]; then
    if [[ -f "$MASK_BED" ]]; then
        log "  - 4d: Excluding Masked Regions"
        vcftools --gzvcf "$SNPS_HARD" \
		--exclude-bed "$MASK_BED" \
		--recode --stdout | bgzip -@ 16 -c > "$NOREPEAT_VCF"
        gatk IndexFeatureFile -I "$NOREPEAT_VCF"
    else
        log "  - Warning: Mask BED not found. Skipping Step 4d."
    fi
fi

# 4d: SNP Softfiltering (Missingness/MAF) # adjust parameters as you wish
SNPS_SOFT="${VCF_DIR}/${SPECIES_ID}_softfiltered_norepeat_biallelic.snps.vcf.gz"
if [[ ! -f "$SNPS_SOFT" ]]; then
    log "  - 4c: Softfiltering (Max-missing 0.85, MAF 0.05)"
    vcftools --gzvcf "$NOREPEAT_VCF" \
	--max-missing 0.85 --maf 0.05 \
	--max-alleles 2 --min-alleles 2 \
	--recode --stdout | bgzip -@ 16 -c > "$SNPS_SOFT"
    gatk IndexFeatureFile -I "$SNPS_SOFT"
fi

###############################################################################
# 5. VARIANT SUMMARY TABLE
###############################################################################
log "Final Step: Generating Variant Summary..."
SUMMARY_OUT="${SPECIES_ID}_stats/${SPECIES_ID}_variant_summary.txt"

{
    echo -e "STAGE\tCOUNT"
    echo -ne "Raw_Genotyped_VCF (DP_Filtered)\t"; bcftools view -H "$GENOTYPED_VCF" | wc -l
    echo -ne "HardFiltered_SNPs\t"; bcftools view -H "$SNPS_HARD" | wc -l
    echo -ne "HardFiltered_Indels\t"; bcftools view -H "$INDELS_HARD" | wc -l
    echo -ne "Final_NoRepeat_SNPs\t"; bcftools view -H "$NOREPEAT_VCF" | wc -l
	echo -ne "SoftFiltered_SNPs\t"; bcftools view -H "$SNPS_SOFT" | wc -l
    fi
} > "$SUMMARY_OUT"

cat "$SUMMARY_OUT"

echo "Pre-BQSR variant calling and filtering for $SPECIES_ID completed. Variant stats exported to $SUMMARY_OUT"
