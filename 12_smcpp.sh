#!/usr/bin/env bash
set -euo pipefail

###############################################################################
# 1. RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="mhe"       # Change this for different species
SUBSET_ID="subset1"      # Change this for different subpopulations
POP_NAME="JP"           # Short name for the subpopulation

# Paths
VCF_DIR="${SPECIES_ID}_vcf_bqsr/${SUBSET_ID}"
VCF="${VCF_DIR}/${SPECIES_ID}_${SUBSET_ID}_clean.vcf.gz"
REF_FAI="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa.fai"
SAMPLE_FILE="${SPECIES_ID}_ref/${SPECIES_ID}_samples_${SUBSET_ID}.txt" # One sample name per line
MASK_BED="${SPECIES_ID}_ref/${SPECIES_ID}_ref_masked_regions.bed.gz"

# Clean the sample file and the contig list of any Windows line endings
[[ -f "$SAMPLE_FILE" ]] && sed -i 's/\r//' "$SAMPLE_FILE"

# Parameters
MU="1.4e-8"
GEN_TIME="9.93"
KNOTS="30"
N_REPS="20"

# Parallelism settings
MAX_JOBS=12              # How many contigs/replicates to run at once
THREADS_PER_JOB=4        # Cores given to EACH smc++ process

# Singularity & Scripts
SIF="smcpp.sif"
SINGULARITY_BIN="${SINGULARITY_BIN:-singularity}"
BOOTSTRAP_SCRIPT="smcpp_bootstrap.py"
CSV_SCRIPT="smcpp_makecsv.py"

# Bootstrap parameters (Used in Step 3 & 4)
CHUNK_SIZE="5000000"    # 5Mb
CHUNKS_PER_CHR="10"     # How many chunks to sample per synthetic chromosome
NR_CHROMOSOMES="30"     # Number of synthetic chromosomes per replicate
N_REPS="20" 

###############################################################################
# 2. DYNAMIC SETUP & DIRECTORIES
###############################################################################
OUTDIR="${SPECIES_ID}_smcpp_${SUBSET_ID}"
SMC_DIR="${OUTDIR}/smc_per_contig"
FULL_EST_DIR="${OUTDIR}/estimate_full"
BOOT_DIR="${OUTDIR}/bootstrap_results"
FINAL_JSON_DIR="${OUTDIR}/final_jsons"
OUTPUT_CSV="${FINAL_JSON_DIR}/${OUTDIR}.csv"
CONTIG_LIST_FILE="${OUTDIR}/contig_list.txt" 
BOOT_PREFIX="${BOOT_DIR}/rep"

mkdir -p "$SMC_DIR" "$FULL_EST_DIR" "$BOOT_DIR" "$FINAL_JSON_DIR"

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

# Helper function to run SMC++ inside Singularity
smcpp_run() { "${SINGULARITY_BIN}" exec "$SIF" smc++ "$@"; }
export -f smcpp_run

###############################################################################
# 3. CONSTRUCT POPULATION SPECIFICATION
###############################################################################
if [[ ! -f "$SAMPLE_FILE" ]]; then die "Sample list $SAMPLE_FILE not found"; fi

# Clean Windows line endings if present
sed -i 's/\r//' "$SAMPLE_FILE"

SAMPLES_COMMA=$(paste -sd "," "$SAMPLE_FILE")
POP_SPEC="${POP_NAME}:${SAMPLES_COMMA}"
log "Population specification: $POP_SPEC"

# NOW export everything needed for 'parallel'
export SMC_DIR VCF POP_SPEC THREADS_PER_JOB SIF SINGULARITY_BIN MASK_BED

###############################################################################
# 4. PREPARE CONTIGS
###############################################################################
log "Preparing contig list (filtering for contigs > 1Mb)..."
CONTIG_REGEX="${CONTIG_REGEX:-}"

# Change {print $1 "\t" $2} to just {print $1}
if [[ -n "${CONTIG_REGEX:-}" ]]; then
    # We check if $2 > 1000000, but we ONLY print $1 (the name)
    awk '$2 > 1000000 {print $1}' "$REF_FAI" | grep -E "${CONTIG_REGEX}" > "$CONTIG_LIST_FILE" || die "No contigs match regex or meet size criteria."
else
    # We check if $2 > 1000000, but we ONLY print $1
    awk '$2 > 1000000 {print $1}' "$REF_FAI" > "$CONTIG_LIST_FILE"
fi

# Check if we actually found any contigs
if [[ ! -s "$CONTIG_LIST_FILE" ]]; then
    die "No contigs found after filtering by size (>1Mb). Check your .fai file paths."
fi

log "Found $(wc -l < "$CONTIG_LIST_FILE") contigs for analysis."

###############################################################################
# 5. STEP 1: vcf2smc (PARALLEL BY CONTIG) - WITH SKIP LOGIC
###############################################################################
log "Step 1: Converting VCF to SMC++ format (Parallel)..."

parallel --jobs "$MAX_JOBS" --progress \
    "if [ ! -f $SMC_DIR/{}.smc.gz ]; then \
        smcpp_run vcf2smc --mask $MASK_BED --cores $THREADS_PER_JOB $VCF $SMC_DIR/{}.smc.gz {} $POP_SPEC; \
     else \
        echo 'Skipping {}: file exists'; \
     fi" \
    :::: "$CONTIG_LIST_FILE"

###############################################################################
# 6. STEP 2: FULL DATA ESTIMATE
###############################################################################
log "Step 2: Running full-data estimate..."
FULL_POP_DIR="${FULL_EST_DIR}/${POP_NAME}"
mkdir -p "$FULL_POP_DIR"

if [[ ! -f "${FULL_POP_DIR}/model.final.json" ]]; then
	smcpp_run estimate \
        --cores $((MAX_JOBS * THREADS_PER_JOB)) \
        --knots "$KNOTS" \
        -o "$FULL_POP_DIR" \
        "$MU" \
        "$SMC_DIR"/*.smc.gz
fi

###############################################################################
# 7. STEP 3: Generate Bootstrap Replicates
###############################################################################
log "Step 3: Generating $N_REPS bootstrap replicates..."

# 1. Find the files and verify they exist
ALL_SMC_FILES=("${SMC_DIR}"/*.smc.gz)

# 2. Safety Check: If the first element doesn't exist, the glob failed
if [[ ! -f "${ALL_SMC_FILES[0]}" ]]; then
    die "No SMC files found in ${SMC_DIR}. Check Step 1 logs."
fi

log "Found ${#ALL_SMC_FILES[@]} SMC files. Starting bootstrap resampling..."

# 3. Run the bootstrap replicates
for i in $(seq 1 "$N_REPS"); do
    REP_DIR="${BOOT_DIR}/rep_${i}"
    mkdir -p "$REP_DIR"
    
    # Use the array we just created
    python3 "$BOOTSTRAP_SCRIPT" \
        --nr_bootstraps "$N_REPS" \
        --nr_chromosomes "$NR_CHROMOSOMES" \
        --chunks_per_chromosome "$CHUNKS_PER_CHR" \
        --chunk_size "$CHUNK_SIZE" \
        "${REP_DIR}/bootstrap" \
        "${ALL_SMC_FILES[@]}" 
done

###############################################################################
# 8. STEP 4: BOOTSTRAP ESTIMATION
###############################################################################
run_rep_est() {
    local i=$1
    local BOOT_PREFIX=$2
    local POP_NAME=$3
    local KNOTS=$4
    local MU=$5
    local THREADS=$6
    
    local REP_DIR="${BOOT_PREFIX}_${i}"
    local REP_POP_DIR="${REP_DIR}/${POP_NAME}"
    
    if [[ ! -f "${REP_POP_DIR}/model.final.json" ]]; then
        mkdir -p "$REP_POP_DIR"
        # Again, pass the directory ($REP_DIR) glob
        smcpp_run estimate \
            --cores "$THREADS" \
            --knots "$KNOTS" \
            -o "$REP_POP_DIR" \
            "$MU" \
            "$REP_DIR"/*.smc.gz
    fi
}
export -f run_rep_est

# Run bootstrap replicates in parallel
seq 0 $((N_REPS - 1)) | parallel --jobs "$MAX_JOBS" --progress \
    run_rep_est {} "$BOOT_PREFIX" "$POP_NAME" "$KNOTS" "$MU" "$THREADS_PER_JOB"

###############################################################################
# 9. STEP 5 & 6: AGGREGATE RESULTS
###############################################################################
log "Step 5: Collecting and Converting Results..."

# Main result
cp "${FULL_POP_DIR}/model.final.json" "${FINAL_JSON_DIR}/${POP_NAME}_main.json"

# Replicate results
for f in "${BOOT_DIR}"/rep_*/"${POP_NAME}"/model.final.json; do
    rep_num=$(echo "$f" | grep -oE 'rep_[0-9]+' | cut -d'_' -f2)
    cp "$f" "${FINAL_JSON_DIR}/${POP_NAME}_rep${rep_num}.json"
done

export GEN_TIME="$GEN_TIME"
python3 "$CSV_SCRIPT" "$OUTPUT_CSV" "$POP_NAME" "${FINAL_JSON_DIR}"/*.json

log "Workflow complete. Final CSV: $OUTPUT_CSV"
