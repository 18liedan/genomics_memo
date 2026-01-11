#!/usr/bin/env bash
set -euo pipefail

#This script automates trimming and mapping of genomes, followed by processing in GATK in preparation for subsequent haplotype calling.

###############################################################################
# RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="yourspecies"

# RESOURCE ALLOCATION (Optimized for NIG supercomputer; 1.5TB RAM / 192 Cores)
MAX_JOBS=12              # 12 samples at a time
THREADS_PER_JOB=16       # 12 * 16 = 192 cores (Full utilization)
JAVA_OPTS="-Xmx64g"      # 12 * 64GB = 768GB (Uses half your RAM, very safe)
SAMTOOLS_MEM="4G"        # Memory per thread for sorting

# PATHS
REF="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa" #your reference genome
SAMPLE_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples.txt" #your list of samples
READS_DIR="${SPECIES_ID}_genomes" #directory with raw sample genomes
TRIM_DIR="${SPECIES_ID}_trimmed" #directory to save trimmed genome reads
BAM_DIR="${SPECIES_ID}_mapped" #directory to save mapped and processed genomes
STATS_DIR="${SPECIES_ID}_stats" #directory to save coverage and depth stats, etc.

mkdir -p "$TRIM_DIR" "$BAM_DIR" "$STATS_DIR"

# Print date on terminal log outputs
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }

# Header for summary (Only create if it doesn't exist)
COVERAGE_SUMMARY="${STATS_DIR}/${SPECIES_ID}_coverage_summary.tsv"
if [[ ! -f "$COVERAGE_SUMMARY" ]]; then
    printf "Sample\tNum_Reads\tCov_Breadth\tMean_Depth\n" > "$COVERAGE_SUMMARY"
fi

# Export variables for GNU Parallel
export SPECIES_ID REF READS_DIR TRIM_DIR BAM_DIR STATS_DIR THREADS_PER_JOB JAVA_OPTS SAMTOOLS_MEM

###############################################################################
# Define the Processing Function
###############################################################################
process_sample() {
    sample=$(echo "$1" | tr -d '\r' | xargs)
    [[ -z "$sample" || "$sample" =~ ^# ]] && return 0

    # Paths for each sample
    R1_IN="${READS_DIR}/${sample}_1.fastq.gz"
    R2_IN="${READS_DIR}/${sample}_2.fastq.gz"
    R1_TRIM="${TRIM_DIR}/trimmed.${sample}_1.fq.gz"
    R2_TRIM="${TRIM_DIR}/trimmed.${sample}_2.fq.gz"
    RAW_BAM="${BAM_DIR}/${sample}_sorted.bam"
    RG_BAM="${BAM_DIR}/${sample}_RG.bam"
    FINAL_BAM="${BAM_DIR}/${sample}_RGdup.bam"
    COV_FILE="${STATS_DIR}/${sample}.coverage.txt"

    echo "[$(date +'%T')] --- Starting: $sample ---"

    # 1. TRIMMING
    if [[ ! -f "$R1_TRIM" ]]; then
        fastp --thread "$THREADS_PER_JOB" -i "$R1_IN" -I "$R2_IN" -o "$R1_TRIM" -O "$R2_TRIM" \
              --detect_adapter_for_pe --html "${TRIM_DIR}/${sample}.html" || return 1
    fi

    # 2. MAPPING & SORTING (Added -m for Samtools)
    if [[ ! -f "$RAW_BAM" ]]; then
        bwa-mem2 mem -t "$THREADS_PER_JOB" "$REF" "$R1_TRIM" "$R2_TRIM" | \
        samtools view -@ "$THREADS_PER_JOB" -Sb - | \
        samtools sort -@ "$THREADS_PER_JOB" -m "$SAMTOOLS_MEM" -o "$RAW_BAM" || return 1
    fi

    # 3. READ GROUPS
    if [[ ! -f "$RG_BAM" ]]; then
        gatk --java-options "$JAVA_OPTS" AddOrReplaceReadGroups \
            -I "$RAW_BAM" -O "$RG_BAM" \
            -RGID "$sample" -RGPL illumina -RGLB "$sample" -RGPU "$sample" -RGSM "$sample" \
            -VALIDATION_STRINGENCY SILENT -SORT_ORDER coordinate || return 1
    fi

    # 4. MARK DUPLICATES
    if [[ ! -f "$FINAL_BAM" ]]; then
        gatk --java-options "$JAVA_OPTS" MarkDuplicates \
            -I "$RG_BAM" -O "$FINAL_BAM" \
            -METRICS_FILE "${STATS_DIR}/${sample}_dup_metrics.txt" \
            --CREATE_INDEX true \
            --VALIDATION_STRINGENCY SILENT || return 1
    fi

    # 5. COVERAGE
    if [[ ! -f "$COV_FILE" ]]; then
        samtools coverage "$FINAL_BAM" > "$COV_FILE"
        MEANDEPTHTOT=$(awk 'NR>1 {sum+=$7; count++} END {if (count > 0) print sum/count; else print 0}' "$COV_FILE")
        BREADTHTOT=$(awk 'NR>1 {sum+=$6; count++} END {if (count > 0) print sum/count; else print 0}' "$COV_FILE")
        READSTOT=$(samtools view -c "$FINAL_BAM")
        # Appending to a shared file is safe with small strings, but for perfect safety in parallel:
        printf "%s\t%s\t%s\t%s\n" "$sample" "$READSTOT" "$BREADTHTOT" "$MEANDEPTHTOT" >> "${STATS_DIR}/${SPECIES_ID}_coverage_summary.tsv"
    fi

    echo "[$(date +'%T')] --- Finished: $sample ---"
}

export -f process_sample

###############################################################################
# Run with GNU Parallel
###############################################################################
# --bar shows a progress bar
# --joblog keeps a record of which samples finished or failed
parallel --jobs "$MAX_JOBS" --joblog "${SPECIES_ID}_pipeline.log" --bar process_sample :::: "$SAMPLE_LIST"

echo "Trimming, mapping, and bam processing complete. Check ${STATS_DIR}/${SPECIES_ID}_coverage_summary.tsv"
