#!/usr/bin/env bash
set -euo pipefail
###############################################################################
# 1. RESOURCE ALLOCATION & CONFIGURATION
###############################################################################
SPECIES_ID="yourspecies"       # Change this for different species
SUBSET_ID="subset1"      # Change this for different subpopulations
POP_NAME="popname"           # Short name for the subpopulation

# Paths
VCF_DIR="${SPECIES_ID}_vcf_bqsr/${SUBSET_ID}" # directory where post-BQSR VCF files are stored
VCF="${VCF_DIR}/${SPECIES_ID}_forsmcpp.vcf.gz" # your VCF file for SMC++
REF_FAI="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa.fai"
SAMPLE_FILE="${SPECIES_ID}_ref/${SPECIES_ID}_samples_${SUBSET_ID}.txt"

# Parameters
MU="1.4e-8" # Change this for your species
GEN_TIME="7.64" # Change this for your species
KNOTS="16"
T_START="1"
T_END="100000"
SPLINE="pchip"

# Bootstrap parameters
CHUNK_SIZE="5000000"
CHUNKS_PER_CHR="20"
NR_CHROMOSOMES="10"
N_BOOTSTRAPS="20" 

# Parallelism settings (Optimized for 192 cores)
MAX_JOBS=12
THREADS_PER_JOB=8

# Singularity & Scripts
SIF="smcpp.sif"
SINGULARITY_BIN="${SINGULARITY_BIN:-singularity}"
BOOTSTRAP_SCRIPT="smcpp_bootstrap.py"
CSV_SCRIPT="smcpp_makecsv.py"

###############################################################################
# 2. SETUP & DIRECTORIES
###############################################################################
OUTDIR="${SPECIES_ID}_smcpp_${SUBSET_ID}"
SMC_DIR="${OUTDIR}/smc_per_contig"
FULL_EST_DIR="${OUTDIR}/estimate_full"
BOOT_DIR="${OUTDIR}/bootstrap_results"
FINAL_JSON_DIR="${OUTDIR}/final_jsons"
OUTPUT_CSV="${FINAL_JSON_DIR}/${OUTDIR}.csv"
CONTIG_LIST_FILE="${OUTDIR}/contig_list.txt" 

mkdir -p "$SMC_DIR" "$FULL_EST_DIR" "$BOOT_DIR" "$FINAL_JSON_DIR"

log() { printf '[%s] %s\n' "$(date +'%F %T')" "$*" >&2; }
die() { log "ERROR: $*"; exit 1; }

smcpp_run() { "${SINGULARITY_BIN}" exec "$SIF" smc++ "$@"; }
export -f smcpp_run

###############################################################################
# 3. CONSTRUCT POPULATION SPECIFICATION
###############################################################################
if [[ ! -f "$SAMPLE_FILE" ]]; then die "Sample list $SAMPLE_FILE not found"; fi
sed -i 's/\r//' "$SAMPLE_FILE"
SAMPLES_COMMA=$(paste -sd "," "$SAMPLE_FILE")
POP_SPEC="${POP_NAME}:${SAMPLES_COMMA}"
export SMC_DIR VCF POP_SPEC THREADS_PER_JOB SIF SINGULARITY_BIN MASK_BED

###############################################################################
# 4. PREPARE CONTIGS
###############################################################################
log "Preparing contig list..."
awk '$2 > 1000000 {print $1}' "$REF_FAI" > "$CONTIG_LIST_FILE"
if [[ ! -s "$CONTIG_LIST_FILE" ]]; then die "No contigs >1Mb found."; fi

###############################################################################
# 5. STEP 1: vcf2smc
###############################################################################
log "Step 1: Converting VCF to SMC++ format..."
parallel --jobs "$MAX_JOBS" --progress \
    "if [ ! -f $SMC_DIR/{}.smc.gz ]; then \
        smcpp_run vcf2smc -c 50000 --cores $THREADS_PER_JOB $VCF $SMC_DIR/{}.smc.gz {} $POP_SPEC; \
     fi" \
    :::: "$CONTIG_LIST_FILE"

###############################################################################
# 6. STEP 2: FULL DATA ESTIMATE
###############################################################################
log "Step 2: Running full-data estimate..."
FULL_POP_DIR="${FULL_EST_DIR}/${POP_NAME}"
mkdir -p "$FULL_POP_DIR"

if [[ ! -f "${FULL_POP_DIR}/model.final.json" ]]; then
    # Expand wildcards here to ensure Singularity gets a clean list of files
    SMC_FILES=("$SMC_DIR"/*.smc.gz)
	smcpp_run estimate \
        --cores $((MAX_JOBS * THREADS_PER_JOB)) \
        --knots "$KNOTS" \
        --spline "$SPLINE" \
        --timepoints "$T_START" "$T_END" \
        -o "$FULL_POP_DIR" \
        "$MU" \
        "${SMC_FILES[@]}"
fi

###############################################################################
# 7. STEP 3: Generate Bootstrap Replicates
###############################################################################
log "Step 3: Checking bootstrap replicates..."
ALL_SMC_FILES=("$SMC_DIR"/*.smc.gz)
# Check for the last expected file to see if we need to regenerate
LAST_BOOT_FILE="${BOOT_DIR}/bootstrap_${N_BOOTSTRAPS}/bootstrap_chr${NR_CHROMOSOMES}.gz"

if [[ -f "$LAST_BOOT_FILE" ]]; then
    log "Skip Step 3: Bootstrap files already exist."
else
    log "Generating $N_BOOTSTRAPS bootstrap replicates..."
    python3 "$BOOTSTRAP_SCRIPT" \
        --nr_bootstraps "$N_BOOTSTRAPS" \
        --nr_chromosomes "$NR_CHROMOSOMES" \
        --chunks_per_chromosome "$CHUNKS_PER_CHR" \
        --chunk_size "$CHUNK_SIZE" \
        "${BOOT_DIR}/bootstrap" \
        "${ALL_SMC_FILES[@]}"
fi

###############################################################################
# 8. STEP 4: BOOTSTRAP ESTIMATION
###############################################################################
log "Step 4: Running bootstrap estimates..."

run_single_boot() {
    local BOOT_SUBDIR=$1
    local POP_NAME=$2
    local KNOTS=$3
    local MU=$4
    local THREADS=$5
    local SPLINE=$6
    local T_START=$7
    local T_END=$8
    
    local OUT_DIR="${BOOT_SUBDIR}/${POP_NAME}_estimate"
    
    if [[ -f "${OUT_DIR}/model.final.json" ]]; then
        return 0
    fi

    mkdir -p "$OUT_DIR"
    
    # Use nullglob so missing files don't result in a literal "*" string
    shopt -s nullglob
    local FILES=("$BOOT_SUBDIR"/*.gz)
    
    if [ ${#FILES[@]} -eq 0 ]; then
        echo "Warning: No .gz files found in $BOOT_SUBDIR"
        return 1
    fi

    smcpp_run estimate \
        --cores "$THREADS" \
        --knots "$KNOTS" \
        --spline "$SPLINE" \
        --timepoints "$T_START" "$T_END" \
        -o "$OUT_DIR" \
        "$MU" \
        "${FILES[@]}"
}

export -f run_single_boot

# FIX: Added -mindepth 1 so it doesn't try to run 'estimate' on the parent BOOT_DIR
find "${BOOT_DIR}" -mindepth 1 -maxdepth 1 -type d -name "bootstrap_*" | \
parallel --jobs "$MAX_JOBS" --progress \
    run_single_boot {} "$POP_NAME" "$KNOTS" "$MU" "$THREADS_PER_JOB" "$SPLINE" "$T_START" "$T_END"

###############################################################################
# 9. STEP 5: AGGREGATE RESULTS
###############################################################################
log "Step 5: Collecting and Converting Results..."

cp "${FULL_POP_DIR}/model.final.json" "${FINAL_JSON_DIR}/${POP_NAME}_main.json"

for f in "${BOOT_DIR}"/bootstrap_*/"${POP_NAME}_estimate"/model.final.json; do
    if [[ -f "$f" ]]; then
        boot_num=$(echo "$f" | grep -oE 'bootstrap_[0-9]+' | cut -d'_' -f2)
        cp "$f" "${FINAL_JSON_DIR}/${POP_NAME}_boot${boot_num}.json"
    fi
done

export GEN_TIME="$GEN_TIME"
python3 "$CSV_SCRIPT" "$OUTPUT_CSV" "$POP_NAME" "${FINAL_JSON_DIR}"/*.json

log "Workflow complete. Final CSV: $OUTPUT_CSV"

