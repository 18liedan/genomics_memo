#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="yourspecies"

# RESOURCE ALLOCATION (Optimized for NIG supercomputer; 1.5TB RAM / 192 Cores)
MAX_JOBS=12              # 12 samples in parallel
THREADS_PER_JOB=16       # 12 * 16 = 192 cores
JAVA_OPTS="-Xmx64g"      # 12 * 64GB = 768GB

# PATHS
REF="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa"
SAMPLE_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples.txt"
BAM_DIR="${SPECIES_ID}_mapped"
VCF_DIR="${SPECIES_ID}_vcf"
BQSR_DIR="${SPECIES_ID}_bqsr"

# Known sites from previous script
KNOWN_SNPS="${VCF_DIR}/${SPECIES_ID}_filtered.snps.vcf.gz"
KNOWN_INDELS="${VCF_DIR}/${SPECIES_ID}_hardfilteredpass.indels.vcf.gz"

mkdir -p "$BQSR_DIR"

# Print date on terminal log outputs
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }

# Verify known sites exist before starting
if [[ ! -f "$KNOWN_SNPS" || ! -f "$KNOWN_INDELS" ]]; then
    log "ERROR: Known sites VCFs not found. Please ensure the previous script finished."
    exit 1
fi

###############################################################################
# 1. PARALLEL BQSR FUNCTION
###############################################################################
# Export variables for GNU Parallel
export REF BAM_DIR BQSR_DIR KNOWN_SNPS KNOWN_INDELS THREADS_PER_JOB JAVA_OPTS

process_bqsr_sample() {
    local sample="$1"
    local INPUT_BAM="${BAM_DIR}/${sample}_RGdup.bam"
    local RECAL_TABLE="${BQSR_DIR}/${sample}_recal_data.table"
    local OUTPUT_BAM="${BQSR_DIR}/${sample}_bqsr.bam"

    # Step A: BaseRecalibrator (Generates the recalibration table)
    if [[ ! -f "$RECAL_TABLE" ]]; then
        log "Starting BaseRecalibrator for $sample..."
        gatk --java-options "${JAVA_OPTS}" BaseRecalibrator \
            -R "$REF" \
            -I "$INPUT_BAM" \
            --known-sites "$KNOWN_SNPS" \
            --known-sites "$KNOWN_INDELS" \
            -O "$RECAL_TABLE"
    else
        log "Recalibration table for $sample already exists. Skipping."
    fi

    # Step B: ApplyBQSR (Applies the table to the BAM)
    if [[ ! -f "$OUTPUT_BAM" ]]; then
        log "Applying BQSR to $sample..."
        gatk --java-options "${JAVA_OPTS}" ApplyBQSR \
            -R "$REF" \
            -I "$INPUT_BAM" \
            --bqsr-recal-file "$RECAL_TABLE" \
            -O "$OUTPUT_BAM"
    else
        log "BQSR BAM for $sample already exists. Skipping."
    fi
}

export -f process_bqsr_sample

###############################################################################
# 2. EXECUTION
###############################################################################
log "Step 4: Running BQSR Recalibration in Parallel..."

# Clean the sample list of any Windows line endings and run parallel
tr -d '\r' < "$SAMPLE_LIST" | parallel --jobs "$MAX_JOBS" --progress process_bqsr_sample {}

log "BQSR process completed for all samples."

###############################################################################
# 3. OPTIONAL: SUMMARY
###############################################################################
echo "BQSR Summary for $SPECIES_ID:"
ls -lh "$BQSR_DIR"/*_bqsr.bam | awk '{print $9, $5}'

echo "BQSR for $SPECIES_ID completed. Move onto second haplotype calling for all samples."
