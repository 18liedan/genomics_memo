#!/bin/bash

# This script performs Fst calculations in ANGSD using designated populations (subsets).

# -------------------------------------------------------------------------
# SCRIPT CONFIGURATION
# -------------------------------------------------------------------------
SPECIES_ID="yourspecies"
SUBSETS=("subset1" "subset2")
# Add as many subpopulations as you need to analyze.
# Make sure to include .filelists in your species' vcf_bqsr directory.
# Each file list should contain paths of post-bqsr bam files of the samples you want to analyze in each subset.

# PATHS TO EXECUTABLES
ANGSD_PATH="angsd" # Path to the angsd binary
REALSFS_PATH="angsd/misc/realSFS" # Path to realSFS

# RESOURCE ALLOCATION
MAX_JOBS=12              
THREADS_PER_JOB=8          

# PATH DEFINITIONS
REF_GENOME="${SPECIES_ID}_ref/${SPECIES_ID}_ref_softmasked_auto.fa"
OUTPUT_SFS_DIR="${SPECIES_ID}_sfspop_folded"

mkdir -p "$OUTPUT_SFS_DIR"

# Export variables for GNU Parallel subshells
export ANGSD_PATH REALSFS_PATH THREADS_PER_JOB REF_GENOME OUTPUT_SFS_DIR SPECIES_ID

# ------------------------------------------------------------------
# STEP 1: SAF Function
# ------------------------------------------------------------------
run_angsd_saf() {
    local SUBSET=$1
    local FILELIST="${SPECIES_ID}_bqsr/${SPECIES_ID}_${SUBSET}.filelist"
    local OUT_PREFIX="$OUTPUT_SFS_DIR/${SPECIES_ID}_$SUBSET"
    
    # 1. SAF Generation (modify parameters as you wish)
    if [[ ! -f "${OUT_PREFIX}.saf.idx" ]]; then
        echo "Processing SAF for $SUBSET..."
        # Ensure unix format
        dos2unix "$FILELIST" 2>/dev/null

        $ANGSD_PATH \
            -b "$FILELIST" \
            -ref "$REF_GENOME" \
            -anc "$REF_GENOME" \
            -out "$OUT_PREFIX" \
            -nThreads "$THREADS_PER_JOB" \
            -doSaf 1 \
			-GL 2 \
			-minMapQ 20 \
			-minQ 30 \
			-remove_bads 1 \
			-uniqueOnly 1 \
			-baq 1 \
			-C 50
    else
        echo "-> SAF index already exists for $SUBSET. Skipping."
    fi

    # 2. 1D SFS Generation
    if [[ ! -f "${OUT_PREFIX}.sfs" ]]; then
        echo "Processing 1D SFS for $SUBSET..."
        $REALSFS_PATH  \
			"${OUT_PREFIX}.saf.idx"  \
			-P "$THREADS_PER_JOB"  \
			-fold 1 \
			> "${OUT_PREFIX}.sfs"
    else
        echo "-> SFS index already exists for $SUBSET. Skipping."
	fi
}
export -f run_angsd_saf

echo "--- Step 1: Running SAF calculations in Parallel ---"
# Passing the array elements directly to parallel
parallel --jobs "$MAX_JOBS" run_angsd_saf ::: "${SUBSETS[@]}"

# ------------------------------------------------------------------
# STEP 2: Fst Function (Processes ONE PAIR at a time)
# ------------------------------------------------------------------
run_angsd_fst() {
    local SUBSET1=$1
    local SUBSET2=$2
    local PAIR="${SUBSET1}-${SUBSET2}"
    local OUT_P1="$OUTPUT_SFS_DIR/${SPECIES_ID}_$SUBSET1"
    local OUT_P2="$OUTPUT_SFS_DIR/${SPECIES_ID}_$SUBSET2"
    local OUT_PAIR="$OUTPUT_SFS_DIR/${SPECIES_ID}_$PAIR"

    if [[ ! -f "${OUT_PAIR}.fst_stats.txt" ]]; then
        echo "Calculating Fst for pair: $SUBSET1 vs $SUBSET2"

        # 1. Calculate 2D SFS (joint SFS)
        $REALSFS_PATH "${OUT_P1}.saf.idx" "${OUT_P2}.saf.idx" -P "$THREADS_PER_JOB" > "${OUT_PAIR}.ml"

        # 2. Index the Fst
        $REALSFS_PATH fst index "${OUT_P1}.saf.idx" "${OUT_P2}.saf.idx" \
            -sfs "${OUT_PAIR}.ml" -fstout "$OUT_PAIR" -P "$THREADS_PER_JOB"

        # 3. Get the global Fst stats
        $REALSFS_PATH fst stats "${OUT_PAIR}.fst.idx" -P "$THREADS_PER_JOB" > "${OUT_PAIR}.fst_stats.txt"
    else
        echo "-> Fst results already exist for $PAIR. Skipping."
    fi
}
export -f run_angsd_fst

echo "--- Step 2: Running Pairwise Fst in Parallel ---"

# Generate the unique pairs list and pipe into parallel
# This replaces the nested loop logic and ensures no duplicates (i < j)
LEN=${#SUBSETS[@]}
for (( i=0; i<LEN; i++ )); do
    for (( j=i+1; j<LEN; j++ )); do
        echo "${SUBSETS[i]} ${SUBSETS[j]}"
    done
done | parallel --colsep ' ' --jobs "$MAX_JOBS" run_angsd_fst {1} {2}

echo "All Fst calculations finished."
