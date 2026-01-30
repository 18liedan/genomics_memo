#!/usr/bin/env bash
set -euo pipefail

# This script performs Final Haplotype Calling on BQSR-calibrated BAMs and 
# processes four different sample subsets for downstream analysis.

###############################################################################
# RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="yourspecies"
SUBSETS=("subset1" "subset2" "subset3") 
# Each subset file contains names of individuals per subpopulation.
# Add as many subpopulations as you need to analyze.
# Make sure to include sample list in your ref directory.

# RESOURCE ALLOCATION (Optimized for NIG supercomputer; 1.5TB RAM / 192 Cores)
MAX_JOBS=12              # 12 samples in parallel
THREADS_PER_JOB=16       # 12 * 16 = 192 cores
JAVA_OPTS="-Xmx64g"      # 12 * 64GB = 768GB
COHORT_JAVA_OPTS="-Xmx256g" 

# PATHS
REF="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa"
# Corresponds to BQSR script's sample list
SAMPLE_LIST_FULL="${SPECIES_ID}_ref/${SPECIES_ID}_samples.txt" 
MASK_BED="${SPECIES_ID}_ref/${SPECIES_ID}_ref_masked_regions.bed"
COVERAGE_STATS="${SPECIES_ID}_stats/${SPECIES_ID}_coverage_summary.tsv"
	# this is produced at the end of "01_trimming-mapping.sh"

# Input BAMs are the output of BQSR script
BAM_DIR="${SPECIES_ID}_bqsr" 
GVCF_DIR="${SPECIES_ID}_gvcf_bqsr"
VCF_DIR="${SPECIES_ID}_vcf_bqsr"
SAMPLE_VCF_DIR="${VCF_DIR}/individual_samples"

mkdir -p "$GVCF_DIR" "$VCF_DIR" "$SAMPLE_VCF_DIR"

# Export variables for GNU Parallel
export REF GVCF_DIR BAM_DIR SAMPLE_VCF_DIR THREADS_PER_JOB JAVA_OPTS COVERAGE_STATS

# Print date on terminal log outputs
log() { printf '[%s] %s\\n' "$(date +'%F %T')" "$*" >&2; }

###############################################################################
# 1. PER-SAMPLE FINAL VARIANT CALLING (HaplotypeCaller)
###############################################################################

run_final_hc() {
    local sample="$1"
    # Input BAMs are from BQSR, named _bqsr.bam
    local INPUT_BAM="${BAM_DIR}/${sample}_bqsr.bam"
    local OUTPUT_GVCF="${GVCF_DIR}/${sample}.bqsr.g.vcf.gz"

    if [[ ! -f "$OUTPUT_GVCF" ]]; then
        log "Starting Final HaplotypeCaller for $sample..."
        gatk --java-options "${JAVA_OPTS}" HaplotypeCaller \
            -R "$REF" -I "$INPUT_BAM" -O "$OUTPUT_GVCF" -ERC GVCF \
            --native-pair-hmm-threads "$THREADS_PER_JOB"
    else
        echo "Skipping Final HC: $sample gVCF already exists."
    fi
}
export -f run_final_hc

log "Step 1a: Running Parallel Final HaplotypeCalling..."
tr -d '\r' < "$SAMPLE_LIST_FULL" | parallel --jobs "$MAX_JOBS" --progress run_final_hc {}

###############################################################################
# 2. PER-SAMPLE FINAL DEPTH FILTERING
###############################################################################

run_final_filtering() {
    local sample="$1"
    local INPUT_GVCF="${GVCF_DIR}/${sample}.bqsr.g.vcf.gz"
    local DP_FILTERED_VCF="${SAMPLE_VCF_DIR}/${sample}.bqsr.dp_filtered.g.vcf.gz"

    if [[ ! -f "$DP_FILTERED_VCF" ]]; then
        # Ensure input exists
        if [[ ! -f "$INPUT_GVCF" ]]; then
            echo "Error: Cannot filter $sample, input $INPUT_GVCF not found." >&2
            return 1
        fi

        DEPTH=$(grep -w "$sample" "$COVERAGE_STATS" | awk '{print $4}')
        if [[ -z "$DEPTH" ]]; then
            echo "Error: Depth for $sample not found in $COVERAGE_STATS" >&2
            return 1
        fi
        
        # Calculate thresholds: Min = $$DEPTH / 3$$, Max = $$DEPTH \times 2$$
        MIN_DP=$(echo "scale=2; $DEPTH / 3" | bc)
        MAX_DP=$(echo "scale=2; $DEPTH * 2" | bc)
        
        echo "Filtering $sample: Depth=$DEPTH (Min=$MIN_DP, Max=$MAX_DP)"
        
        vcftools --gzvcf "$INPUT_GVCF" --min-meanDP "$MIN_DP" --max-meanDP "$MAX_DP" --recode --stdout | \
            bgzip -c > "$DP_FILTERED_VCF"
        gatk IndexFeatureFile -I "$DP_FILTERED_VCF"
    else
        echo "Skipping Filter: $sample final filtered VCF exists."
    fi
}
export -f run_final_filtering

log "Step 1b: Running Parallel Per-Sample Final Depth Filtering..."
tr -d '\r' < "$SAMPLE_LIST_FULL" | parallel --jobs "$MAX_JOBS" --progress run_final_filtering {}

###############################################################################
# 3. SUBSET PROCESSING (Merging, Genotyping, & Filtering)
###############################################################################

for SUBSET in "${SUBSETS[@]}"; do
    SUBSET_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples_${SUBSET}.txt"
    SUBSET_DIR="${VCF_DIR}/${SUBSET}"
    mkdir -p "$SUBSET_DIR"
    
    log "Processing SUBSET: $SUBSET"

    # --- 3a. CombineGVCFs for Subset ---
    MERGED_GVCF="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_merged.g.vcf.gz"
    if [[ ! -f "$MERGED_GVCF" ]]; then
        log "  - Combining gVCFs for $SUBSET"
        VCF_INPUTS=()
        while read -r sample; do
            sample=$(echo "$sample" | tr -d '\r')
            VCF_FILE="${SAMPLE_VCF_DIR}/${sample}.bqsr.dp_filtered.g.vcf.gz"
            [[ -f "$VCF_FILE" ]] && VCF_INPUTS+=("-V" "$VCF_FILE")
        done < "$SUBSET_LIST"

        gatk --java-options "${COHORT_JAVA_OPTS}" CombineGVCFs \
            -R "$REF" "${VCF_INPUTS[@]}" -O "$MERGED_GVCF"
    else
        log "  - Combined gVCF for $SUBSET already exists. Skipping."
    fi

    # --- 3b. GenotypeGVCFs for Subset ---
    GENOTYPED_VCF="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_genotyped.vcf.gz"
    if [[ ! -f "$GENOTYPED_VCF" ]]; then
        log "  - Genotyping $SUBSET"
        gatk --java-options "${COHORT_JAVA_OPTS}" GenotypeGVCFs \
            -R "$REF" -V "$MERGED_GVCF" -O "$GENOTYPED_VCF"
    else
        log "  - Genotyped VCF for $SUBSET already exists. Skipping."
    fi

    # --- 3c. SNP Extraction and Hard Filtering ---
    SNPS_HARD="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_hardfiltered.snps.vcf.gz"
    if [[ ! -f "$SNPS_HARD" ]]; then
        log "  - SNP Hard-Filtering ($SUBSET)"
        gatk SelectVariants -V "$GENOTYPED_VCF" -select-type SNP -O "${SUBSET_DIR}/temp_snps.vcf.gz"
        gatk VariantFiltration -V "${SUBSET_DIR}/temp_snps.vcf.gz" \
            -filter "QD < 2.0 || QUAL < 30.0 || SOR > 3.0 || FS > 60.0 || MQ < 40.0 || MQRankSum < -12.5 || ReadPosRankSum < -8.0" \
            --filter-name "snp_hard_filter" -O "${SUBSET_DIR}/temp_snps_filt.vcf.gz"
        bcftools view -f PASS "${SUBSET_DIR}/temp_snps_filt.vcf.gz" -Oz -o "$SNPS_HARD"
        gatk IndexFeatureFile -I "$SNPS_HARD"
        rm "${SUBSET_DIR}/temp_snps.vcf.gz" "${SUBSET_DIR}/temp_snps_filt.vcf.gz"* # Removed index file cleanup
    else
        log "  - Hard-filtered SNPs for $SUBSET already exists. Skipping."
    fi

    # --- 3d. Softmasking: Missingness & Repeat Masking ---
    SNPS_SOFT="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered.snps.vcf.gz"
    if [[ ! -f "$SNPS_SOFT" ]]; then
        log "  - First soft-filtering ($SUBSET: Missingness 0.8)"
		vcftools --gzvcf "$SNPS_HARD" --max-missing 0.8 --recode --stdout | \
			bgzip -c > "$SNPS_SOFT"
    else
        log "  - Soft-filtered VCF for $SUBSET already exists. Skipping."
    fi

	SNPS_BIALLELIC="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered_biallelic.snps.vcf.gz"
    if [[ ! -f "$SNPS_BIALLELIC" ]]; then
        log "  - Second soft-filtering ($SUBSET: Missingness 0.8, biallelic only)"
		vcftools --gzvcf "$SNPS_SOFT" --max-alleles 2 --min-alleles 2 --recode --stdout | \
			bgzip -c > "$SNPS_BIALLELIC"
    else
        log "  - Soft-filtered biallelic VCF for $SUBSET already exists. Skipping."
    fi

	SNPS_BIALLELIC_NOREPEAT ="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered_norepeat_biallelic.snps.vcf.gz"
    if [[ ! -f "$SNPS_NOREPEAT" ]]; then
        log "  - Third soft-filtering ($SUBSET: Missingness 0.8, biallelic, no repeats)"
        if [[ -f "$MASK_BED" ]]; then
            vcftools --gzvcf "$SNPS_BIALLELIC" --exclude-bed "$MASK_BED" --recode --stdout | \
                bgzip -c > "$SNPS_NOREPEAT"
			gatk IndexFeatureFile -I "$SNPS_BIALLELIC_NOREPEAT"
        else
			log "  - Hard masked sites not found in ref directory. Skipping."
        fi
    else
        log "  - No repeat VCF for $SUBSET already exists. Skipping."
    fi
done

###############################################################################
# 4. FINAL SUMMARY TABLE
###############################################################################
log "Final Step: Generating Subset Summary..."
SUMMARY_OUT="${SPECIES_ID}_stats/${SPECIES_ID}_final_variant_summary.txt"

{
    echo -e "SUBSET\tSTAGE\tCOUNT"
    for SUBSET in "${SUBSETS[@]}"; do
        SUBSET_DIR="${VCF_DIR}/${SUBSET}"
        GENO="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_genotyped.vcf.gz"
        HARD="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_hardfiltered.snps.vcf.gz"
        SOFT="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered.snps.vcf.gz"
        SOFT_BIALLELIC="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered_biallelic.snps.vcf.gz"
        SNPS_BIALLELIC_NOREPEAT ="${SUBSET_DIR}/${SPECIES_ID}_${SUBSET}_softfiltered_norepeat_biallelic.snps.vcf.gz"

        # Helper function to count variants safely
        count_vars() {
            local file=$1
            if [[ -f "$file" ]]; then
                # -H ignores the header for faster counting
                bcftools view -H "$file" 2>/dev/null | wc -l
            else
                echo "0"
            fi
        }

        log "  - Counting variants for subset: $SUBSET"
        
        COUNT_GENO=$(count_vars "$GENO")
        COUNT_HARD=$(count_vars "$HARD")
        COUNT_SOFT=$(count_vars "$SOFT")
        COUNT_SOFT_BIALLELIC=$(count_vars "$SOFT_BIALLELIC")
        COUNT_SOFT_NOREPEAT_BIALLELIC=$(count_vars "$SOFT_NOREPEAT_BIALLELIC")

        # Output to the summary file with correct STAGE labels
        echo -e "${SUBSET}\tRaw_Genotyped\t${COUNT_GENO}"
        echo -e "${SUBSET}\tHardFiltered_SNPs\t${COUNT_HARD}"
        echo -e "${SUBSET}\tSoftFiltered_SNPs\t${COUNT_SOFT}"
        echo -e "${SUBSET}\tBiallelic_SNPs\t${COUNT_SOFT_BIALLELIC}"
        echo -e "${SUBSET}\tNoRepeatBiallelic_SNPs\t${COUNT_SOFT_NOREPEAT_BIALLELIC}"
    done
} > "$SUMMARY_OUT"

cat "$SUMMARY_OUT"
log "Final variant calling completed for all subsets. Results are in $VCF_DIR"
