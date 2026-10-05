#!/bin/bash
#SBATCH -o slurm-%j-minimap2.out
#SBATCH --job-name=minimap2_ont
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --mem=90G

# ==============================
# User-defined paths
# ==============================

BASE_DIR="/path/to/project"
REFERENCE_DIR="/path/to/reference"

FASTQ_DIR="${BASE_DIR}/fastplong_output"
ALIGN_DIR="${BASE_DIR}/minimap2_alignment"
BAM_DIR="${BASE_DIR}/bam"
LOG_DIR="${BASE_DIR}/logs_minimap2"

PROCESSED_LOG="${LOG_DIR}/processed_minimap2.log"

REFERENCE="${REFERENCE_DIR}/GRCh38.p14.genome.fa"
ANNOTATION_BED="${REFERENCE_DIR}/gencode_v49.bed"

# ==============================
# Conda environment
# ==============================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate minimap2_env

# ==============================
# Parallelisation settings
# ==============================

MAX_PARALLEL_JOBS=3
CORES_PER_JOB=8

mkdir -p "$ALIGN_DIR" "$BAM_DIR" "$LOG_DIR"
touch "$PROCESSED_LOG"

cd "$BASE_DIR" || exit 1

# ==============================
# Function to process one sample
# ==============================

process_sample() {
    local FASTQ=$1
    local CORES=$2

    local BASENAME
    BASENAME=$(basename "$FASTQ" .filtered.fastq.gz)

    local LOG_FILE="${LOG_DIR}/${BASENAME}.log"
    local SAM_FILE="${ALIGN_DIR}/${BASENAME}.sam"
    local BAM_FILE="${BAM_DIR}/${BASENAME}.sorted.bam"

    {
        echo "--- Starting minimap2 for ${BASENAME} at $(date) ---"

        minimap2 \
            -a \
            -x splice \
            --junc-bed "$ANNOTATION_BED" \
            -k14 \
            -t "$CORES" \
            "$REFERENCE" \
            "$FASTQ" > "$SAM_FILE"

        if [ $? -ne 0 ]; then
            echo "ERROR: minimap2 failed for ${BASENAME}"
            return 1
        fi

        echo "--- Sorting BAM for ${BASENAME} ---"

        samtools sort \
            -@ "$CORES" \
            -o "$BAM_FILE" \
            "$SAM_FILE"

        if [ $? -ne 0 ]; then
            echo "ERROR: samtools sort failed for ${BASENAME}"
            return 1
        fi

        samtools index "$BAM_FILE"

        if [ $? -ne 0 ]; then
            echo "ERROR: samtools index failed for ${BASENAME}"
            return 1
        fi

        rm "$SAM_FILE"

        echo "$BASENAME" >> "$PROCESSED_LOG"
        echo "--- Finished ${BASENAME} at $(date) ---"

    } >> "$LOG_FILE" 2>&1
}

export -f process_sample
export FASTQ_DIR ALIGN_DIR BAM_DIR LOG_DIR PROCESSED_LOG REFERENCE ANNOTATION_BED

# ==============================
# Control number of parallel jobs
# ==============================

wait_for_slot() {
    while true; do
        local job_count
        job_count=$(jobs -p | wc -l)

        if [ "$job_count" -lt "$MAX_PARALLEL_JOBS" ]; then
            break
        fi

        sleep 2
    done
}

# ==============================
# Process all filtered FASTQ files
# ==============================

for FASTQ in "$FASTQ_DIR"/*.filtered.fastq.gz; do

    if [[ ! -f "$FASTQ" ]]; then
        echo "No filtered FASTQ files found"
        break
    fi

    BASENAME=$(basename "$FASTQ" .filtered.fastq.gz)

    if grep -Fxq "$BASENAME" "$PROCESSED_LOG"; then
        echo "Skipping ${BASENAME}, already processed."
        continue
    fi

    echo "Starting minimap2 for ${BASENAME}"

    wait_for_slot
    process_sample "$FASTQ" "$CORES_PER_JOB" &

done

wait

echo "All minimap2 alignment jobs completed."
