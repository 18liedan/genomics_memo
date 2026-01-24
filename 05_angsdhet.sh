#!/usr/bin/env bash
set -euo pipefail

# This script performs heterozygosity calculations in ANGSD using each sample's BQSR-calibrated BAM file.

# -------------------------------------------------------------------------
# SCRIPT CONFIGURATION
# -------------------------------------------------------------------------
SPECIES_ID="yourspecies"

# PATHS TO EXECUTABLES
ANGSD_PATH="angsd" # Directory containing 'misc/realSFS'

# RESOURCE ALLOCATION
MAX_JOBS=12              
THREADS_PER_JOB=8          

# PATH DEFINITIONS
REF_GENOME="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa"
SAMPLE_LIST="${SPECIES_ID}_ref/${SPECIES_ID}_samples.txt"

INPUT_BAM_DIR="${SPECIES_ID}_bqsr"
OUTPUT_SFS_DIR="${SPECIES_ID}_sfs_folded"

# Create output directory
mkdir -p "$OUTPUT_SFS_DIR"

# Print date on terminal log outputs
log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }

# Export variables for GNU Parallel
export ANGSD_PATH THREADS_PER_JOB REF_GENOME INPUT_BAM_DIR OUTPUT_SFS_DIR SPECIES_ID

# -------------------------------------------------------------------------
# VALIDATION
# -------------------------------------------------------------------------
if [ ! -f "$REF_GENOME" ]; then
    log "ERROR: Reference genome not found at '$REF_GENOME'"
    exit 1
fi

if [ ! -f "$SAMPLE_LIST" ]; then
    log "ERROR: Sample list not found at '$SAMPLE_LIST'"
    exit 1
fi

if [ ! -x "$ANGSD_PATH/misc/realSFS" ]; then
    log "ERROR: ANGSD executable not found at '$ANGSD_PATH'"
    exit 1
fi

if ! command -v samtools &> /dev/null; then
    log "ERROR: samtools could not be found. Please load the module or add to PATH."
    exit 1
fi

# -------------------------------------------------------------------------
# STEP 1: CALCULATE SAF AND SFS (Per Sample Function)
# -------------------------------------------------------------------------
run_angsd_het() {
    local sample="$1"
    local INPUT_BAM="${INPUT_BAM_DIR}/${sample}_bqsr.bam"
    local OUTPUT_PREFIX="${OUTPUT_SFS_DIR}/${sample}"
    local SFS_FILE="${OUTPUT_PREFIX}.sfs"
    local COV_FILE="${OUTPUT_SFS_DIR}/${sample}.samtools_cov.txt"

    # 1. CALCULATE MEAN DEPTH USING SAMTOOLS
    # Check if we already have the coverage file to save time
    if [[ ! -f "$COV_FILE" ]]; then
        echo "Calculating coverage for $sample using samtools..."
        samtools coverage "$INPUT_BAM" > "$COV_FILE"
    else
        echo "-> Coverage file exists for $sample. Reading depth..."
    fi

    # Parse Mean Depth:
    # samtools coverage outputs a row per scaffold. 
    # We calculate the weighted average depth across all scaffolds.
    local MEAN_DEPTH=$(awk 'NR>1 {sum+=$7*($3-$2); len+=$3-$2} END {if(len>0) print sum/len; else print 0}' "$COV_FILE")
    
    # Calculate thresholds using bc (1/3 and 2x)
    # printf "%.0f" rounds the float to the nearest integer for ANGSD
    local MIN_DP=$(printf "%.0f" $(echo "$MEAN_DEPTH / 3" | bc -l))
    local MAX_DP=$(printf "%.0f" $(echo "$MEAN_DEPTH * 2" | bc -l))

    echo "Sample: $sample | Mean Depth: $MEAN_DEPTH | MinDP: $MIN_DP | MaxDP: $MAX_DP"

    # 2. RUN ANGSD SAF
    if [[ ! -f "${OUTPUT_PREFIX}.saf.idx" ]]; then
        "${ANGSD_PATH}" \ # modify parameters as you wish
            -i "$INPUT_BAM" \
            -ref "$REF_GENOME" \
            -anc "$REF_GENOME" \
            -out "$OUTPUT_PREFIX" \
            -nThreads "$THREADS_PER_JOB" \
            -doSaf 1 \
            -GL 2 \
            -doCounts 1 \
            -minMapQ 20 \
            -minQ 30 \
            -remove_bads 1 \
            -uniqueOnly 1 \
            -only_proper_pairs 1 \
            -baq 1 \
            -C 50 \
            -setMinDepth "$MIN_DP" \
            -setMaxDepth "$MAX_DP"
    else
        echo "-> SAF index already exists for $sample. Skipping to realSFS."
    fi

    # 3. RUN REALSFS
    if [[ ! -f "$SFS_FILE" ]]; then
        "$ANGSD_PATH/misc/realSFS" \
            "$OUTPUT_PREFIX.saf.idx" \
            -P "$THREADS_PER_JOB" \
			-fold 1 \
            > "$SFS_FILE"
    else
        echo "-> SFS file already exists for $sample."
    fi

    echo "Finished $sample."
}
export -f run_angsd_het

log "Step 1: Running Parallel ANGSD Heterozygosity Calculation..."
tr -d '\r' < "$SAMPLE_LIST" | parallel --jobs "$MAX_JOBS" --progress run_angsd_het {}

# -------------------------------------------------------------------------
# STEP 2: SUMMARIZE RESULTS
# -------------------------------------------------------------------------
log "Step 2: Consolidating heterozygosity results..."

OUTPUT_SUMMARY_FILE="${OUTPUT_SFS_DIR}/${SPECIES_ID}_angsdhet_summary.txt"

# Header
echo -e "SampleID\tHomRef\tHeterozygous\tHomAlt\tTotalSites\tHeterozygosity" > "$OUTPUT_SUMMARY_FILE"

# Iterate through sample list to maintain order, rather than globbing
while read -r sample; do
    # Remove carriage returns if present
    sample=$(echo "$sample" | tr -d '\r')
    sfs_file="${OUTPUT_SFS_DIR}/${sample}.sfs"

    if [ -f "$sfs_file" ]; then
        # Read SFS values into array
        sfs_values=($(cat "$sfs_file"))

        # realSFS output for diploid: [0]=AA (HomRef), [1]=Aa (Het), [2]=aa (HomAlt)
        hom_ref=${sfs_values[0]:-0}
        het_sites=${sfs_values[1]:-0}
        hom_alt=${sfs_values[2]:-0}

        # Calculate using awk
        summary=$(awk -v h_ref="$hom_ref" -v het="$het_sites" -v h_alt="$hom_alt" '
            BEGIN {
                total = h_ref + het + h_alt;
                if (total > 0) {
                    het_rate = het / total;
                } else {
                    het_rate = 0;
                }
                # Output format: HomRef Het HomAlt Total HetRate
                printf "%.0f\t%.0f\t%.0f\t%.0f\t%0.8f", h_ref, het, h_alt, total, het_rate;
            }
        ')
        
        echo -e "${sample}\t${summary}" >> "$OUTPUT_SUMMARY_FILE"
    else
        echo "WARNING: No SFS file found for $sample during summary generation." >&2
    fi

done < "$SAMPLE_LIST"

log "Done. Summary saved to: $OUTPUT_SUMMARY_FILE"
cat "$OUTPUT_SUMMARY_FILE"
