#!/bin/bash
set -e  # Exit immediately if a command exits with a non-zero status.
set -u  # Treat unset variables as an error.
set -o pipefail # Causes a pipeline to return the exit status of the last command in the pipe.

# =========================================================================
#
# PSMC Pipeline
# DESCRIPTION: Runs a complete PSMC (Pairwise Sequentially Markovian Coalescent)
#              analysis for each sample to infer demographic history.
#
# =========================================================================

echo "================================================================"
echo "STARTING: PSMC Demographic History Pipeline"
echo "DATE: $(date)"
echo "================================================================"

# -------------------------------------------------------------------------
# Step 1: Define Variables and Set Up Environment
# -------------------------------------------------------------------------
echo "Step 1: Setting up variables and environment..."

# --- User-Defined Variables ---
# define your samples and species here
declare -a SAMPLES=(
	"sample1"
	"sample2"
	"sample3"
)

SPECIES_ID="yourspecies"

# Name of the Conda environment containing bcftools, samtools, and psmc
CONDA_ENV_NAME="yourenv"

# --- Project Structure Variables ---

# Base directory for the project
BASE_DIR=$(pwd)

# Input directory containing the final, recalibrated BAM files (see bqsr script in same repository for detials)
INPUT_BAM_DIR="$BASE_DIR/${SPECIES_ID}_bqsr"

# Reference genome file
REF_GENOME="$BASE_DIR/${SPECIES_ID}_ref/${SPECIES_ID}_ref.fa"

# Output directory for all PSMC results
PSMC_OUT_DIR="$BASE_DIR/${SPECIES_ID}_psmc"

# Directory to store calculated mean coverage values
COVERAGE_STATS_DIR="$BASE_DIR/${SPECIES_ID}_covstats"

# --- PSMC Analysis Parameters ---

# bcftools mpileup: coefficient for downgrading mapping quality
MPILEUP_C=50

# Minimum mapping quality
MIN_MAPPING_QUALITY=20
# fq2psmcfa: minimum base quality
PSMC_Q=20

# --- Depth Filtering Factors ---
MIN_DEPTH_FACTOR=0.33333
MAX_DEPTH_FACTOR=2.0

# psmc parameters
PSMC_N=25 # Number of iterations
PSMC_T=15 # Max number of hidden states in the HMM (t_max)
PSMC_R=5  # Ratio of recombination to mutation rate

# psmc: different time interval patterns to test
# This one remains an associative array as it uses pattern names as keys
declare -A PSMC_PATTERNS
PSMC_PATTERNS["a"]="1+1+1+1+25*2+4+6"
# try other patterns if you wish, like: PSMC_PATTERNS["b"]="2+2+25*2+4+6"

# psmc_plot.pl parameters
  # replace "generation_time" with your species' actual generation time as a number
  # replace "mutation_rate" with your species' estimated per-generation mutation rate as a number
declare -a PLOT_PARAMS=(
    "generation_time mutation_rate"
)

# -------------------------------------------------------------------------
# Step 2: Activate Conda Environment and Validate Inputs
# -------------------------------------------------------------------------
echo "Step 2: Activating Conda environment and validating inputs..."

# --- Activate Conda ---
source "$(conda info --base)/etc/profile.d/conda.sh" || { echo "ERROR: Could not source conda.sh. Is Conda installed?" >&2; exit 1; }
conda activate "$CONDA_ENV_NAME" || { echo "ERROR: Could not activate Conda env '$CONDA_ENV_NAME'." >&2; exit 1; }
echo "Conda environment '$CONDA_ENV_NAME' activated."

# --- Create Output Directory ---
mkdir -p "$PSMC_OUT_DIR"
mkdir -p "$COVERAGE_STATS_DIR"

# --- Input Validation ---
if [ ! -f "$REF_GENOME" ]; then echo "ERROR: Ref genome not found at '$REF_GENOME'" >&2; exit 1; fi
if [ ! -d "$INPUT_BAM_DIR" ]; then echo "ERROR: Input BAM directory '$INPUT_BAM_DIR' not found." >&2; exit 1; fi

for tool in bcftools vcfutils.pl fq2psmcfa psmc psmc_plot.pl samtools awk bc; do
    if ! command -v "$tool" &> /dev/null; then
        echo "ERROR: Tool '$tool' not found in PATH. Is the Conda environment correct?" >&2; exit 1;
    fi
done
echo "All checks passed. Tools and inputs are available."

# -------------------------------------------------------------------------
# Step 3: Loop Through Samples and Run PSMC Pipeline
# -------------------------------------------------------------------------
echo "Step 3: Starting main PSMC analysis loop for each sample..."

# Loop through the values of the SAMPLES array (removed the !)
for SAMPLE_ID in "${SAMPLES[@]}"; do

    echo -e "\n----------------------------------------------------------------"
    echo "Processing Sample: $SAMPLE_ID"
    echo "----------------------------------------------------------------"

    # --- Define sample-specific files ---
    INPUT_BAM="${INPUT_BAM_DIR}/${SAMPLE_ID}_bqsr.bam"
    PSMC_FQ_GZ="${PSMC_OUT_DIR}/${SAMPLE_ID}.psmc.fq.gz"
    PSMC_FASTA="${PSMC_OUT_DIR}/${SAMPLE_ID}.psmcfa"
    MEAN_COVERAGE_FILE="${COVERAGE_STATS_DIR}/${SAMPLE_ID}.mean_coverage.txt"

    # Check if the input BAM exists
    if [ ! -f "$INPUT_BAM" ]; then
        echo "WARNING: Recalibrated BAM file not found for '$SAMPLE_ID'. Skipping." >&2
        continue
    fi

    # --- Part A: Calculate Mean Depth and Dynamic Filters ---
    echo "-> Calculating depth filters for $SAMPLE_ID..."
    
    local_min_depth=0
    local_max_depth=0

    # ADDED SKIP LOGIC: Check if the mean coverage file already exists and is not empty
    if [ ! -s "$MEAN_COVERAGE_FILE" ]; then
        echo "  - Calculating mean depth using samtools depth on '$INPUT_BAM' with minMapQ $MIN_MAPPING_QUALITY..."
        MEAN_DP=$(samtools depth -a -q "$MIN_MAPPING_QUALITY" -Q "$MIN_MAPPING_QUALITY" "$INPUT_BAM" | awk '{sum+=$3} END { if (NR>0) print sum/NR; else print 0}')
        
        # Save the mean depth to a file
        echo "$MEAN_DP" > "$MEAN_COVERAGE_FILE"
        echo "  - Calculated mean depth: $MEAN_DP (saved to $MEAN_COVERAGE_FILE)"
    else
        # Read the existing depth from the file
        MEAN_DP=$(cat "$MEAN_COVERAGE_FILE")
        echo "  - Part A: Mean coverage file $MEAN_COVERAGE_FILE already exists. Read mean depth: $MEAN_DP. Skipping calculation."
    fi

    # Calculate min and max depth using the defined factors and ensure they are integers
    local_min_depth=$(awk -v mean="$MEAN_DP" -v factor="$MIN_DEPTH_FACTOR" 'BEGIN {printf "%.0f", mean * factor}')
    local_max_depth=$(awk -v mean="$MEAN_DP" -v factor="$MAX_DEPTH_FACTOR" 'BEGIN {printf "%.0f", mean * factor}')

    # Fallback for very low coverage samples
    if [ "$local_min_depth" -eq 0 ]; then
        local_min_depth=1
        echo "  - WARNING: Calculated min depth was 0. Setting to a default of 1." >&2
    fi
    if [ "$local_max_depth" -eq 0 ]; then
        local_max_depth=100
        echo "  - WARNING: Calculated max depth was 0. Setting to a default of 100." >&2
    fi
    # Ensure min depth is not greater than max depth
    if (( $(echo "$local_min_depth > $local_max_depth" | bc -l) )); then
        echo "  - WARNING: Calculated min depth ($local_min_depth) was greater than max depth ($local_max_depth). Adjusting min depth to be equal to max depth." >&2
        local_min_depth="$local_max_depth"
    fi

    echo "  - Using dynamic Min Depth for vcf2fq: $local_min_depth | Max Depth: $local_max_depth"


    # --- Part B: Generate consensus sequence (fq.gz) ---
    if [ ! -s "$PSMC_FQ_GZ" ]; then
        echo "-> Part B: Generating PSMC consensus fastq for $SAMPLE_ID..."
        bcftools mpileup -Ou -C"$MPILEUP_C" -f "$REF_GENOME" "$INPUT_BAM" | \
        bcftools call -c -Ov | \
        vcfutils.pl vcf2fq -d "$local_min_depth" -D "$local_max_depth" | \
        gzip > "$PSMC_FQ_GZ"
        echo "-> Consensus sequence created: $PSMC_FQ_GZ"
    else
        echo "-> Part B: Consensus fastq $PSMC_FQ_GZ already exists. Skipping generation."
    fi

    # --- Part C: Convert fq.gz to psmcfa format ---
    if [ ! -s "$PSMC_FASTA" ]; then
        echo "-> Part C: Converting to .psmcfa format..."
        fq2psmcfa -q"$PSMC_Q" "$PSMC_FQ_GZ" > "$PSMC_FASTA"
        if [ ! -s "$PSMC_FASTA" ]; then
            echo "ERROR: PSMC FASTA file is empty. The depth filters might be too strict, or the input FQ file was empty. Skipping rest of analysis for $SAMPLE_ID." >&2
            continue
        fi
        echo "-> PSMC FASTA file created: $PSMC_FASTA"
    else
        echo "-> Part C: PSMC FASTA file $PSMC_FASTA already exists. Skipping conversion."
    fi

    # --- Part D: Run PSMC with different patterns (if you wish) ---
    echo "-> Part D: Running PSMC with different patterns..."
    for suffix in "${!PSMC_PATTERNS[@]}"; do
        pattern="${PSMC_PATTERNS[$suffix]}"
        PSMC_OUT_FILE="${PSMC_OUT_DIR}/${SAMPLE_ID}.${suffix}.psmc"

        if [ ! -s "$PSMC_OUT_FILE" ]; then
            echo "  - Running pattern '$suffix' ($pattern)..."
            psmc -N"$PSMC_N" -t"$PSMC_T" -r"$PSMC_R" -p "$pattern" -o "$PSMC_OUT_FILE" "$PSMC_FASTA"
        else
            echo "  - Part D: PSMC output $PSMC_OUT_FILE already exists. Skipping run."
        fi
    done

    # --- Part E: Generate plots for each PSMC output ---
    echo "-> Part E: Generating plots for each completed run..."
    for suffix in "${!PSMC_PATTERNS[@]}"; do
        PSMC_OUT_FILE="${PSMC_OUT_DIR}/${SAMPLE_ID}.${suffix}.psmc"
        if [ ! -s "$PSMC_OUT_FILE" ]; then
            echo "  - WARNING: Cannot plot for pattern '$suffix' as output file is missing. Skipping." >&2
            continue
        fi

        for param_pair in "${PLOT_PARAMS[@]}"; do
            read -r gen_time mut_rate <<< "$param_pair"
            PLOT_PREFIX="${PSMC_OUT_DIR}/${SAMPLE_ID}.${suffix}.g${gen_time}.mu${mut_rate}"
            
            if [ ! -f "${PLOT_PREFIX}.plot.eps" ]; then
                echo "  - Plotting pattern '$suffix' with g=$gen_time, u=$mut_rate"
                psmc_plot.pl -g "$gen_time" -u "$mut_rate" "$PLOT_PREFIX" "$PSMC_OUT_FILE"
            else
                echo "  - Part E: Plot for pattern '$suffix' with g=$gen_time, u=$mut_rate already exists. Skipping."
            fi
        done
    done

    echo "Finished processing sample: $SAMPLE_ID"

done

# -------------------------------------------------------------------------
# Step 4: Convert results to txt for visualization in R
# -------------------------------------------------------------------------
# --- Parameters (Must match your downstream R script) ---
G_TIME=5 #replace this with your actual generation time
MUT_RATE="1e-8" #replace this with your actual generation time
PSMC_DIR="./${SPECIES_ID}_psmc"

# Check if directory exists
if [ ! -d "$PSMC_DIR" ]; then
    echo "ERROR: Directory $PSMC_DIR not found."
    exit 1
fi

echo "----------------------------------------------------------------"
echo "Starting conversion for R plotting..."
echo "Parameters: g=$G_TIME, u=$MUT_RATE"
echo "----------------------------------------------------------------"

# Loop through all .psmc files in the directory
# This handles "n" samples automatically
for PSMC_FILE in "${PSMC_DIR}"/*.psmc; do
    
    # Get the base name of the file (e.g., bobu_01.a)
    BASE_NAME=$(basename "$PSMC_FILE" .psmc)
    
    # Define the output prefix and the expected data file name
    # psmc_plot.pl -R adds '.0.txt' to the prefix you provide
    OUT_PREFIX="${PSMC_DIR}/${BASE_NAME}.scaled"
    DATA_FILE="${OUT_PREFIX}.0.txt"

    # --- Skip Logic ---
    if [ -s "$DATA_FILE" ]; then
        echo "[SKIP] Data file already exists for: $BASE_NAME"
    else
        echo "[RUN] Generating R data for: $BASE_NAME"
        
        # -R: Generate R-readable text file instead of plotting with gnuplot
        # -g: Generation time
        # -u: Mutation rate
        psmc_plot.pl -R -g "$G_TIME" -u "$MUT_RATE" "$OUT_PREFIX" "$PSMC_FILE"
        
        # Check if the command succeeded
        if [ $? -eq 0 ]; then
            echo "      Successfully created $DATA_FILE"
        else
            echo "      ERROR: Failed to process $BASE_NAME"
        fi
    fi
done

echo "----------------------------------------------------------------"
echo "Conversion complete. You can now run your R script."
echo "----------------------------------------------------------------"

# --- Deactivate Conda Environment ---
conda deactivate
echo "Conda environment deactivated."

echo "================================================================"
echo "SCRIPT COMPLETE: PSMC Demographic History Pipeline"
echo "DATE: $(date)"
echo "All outputs are located in: $PSMC_OUT_DIR"
echo "================================================================"
