#!/bin/bash
#SBATCH -o slurm-%j-fastplong.out
#SBATCH --job-name=fastplong_ont
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --mem=90G

cd /data/scratch/ha25158/sra_work || exit 1

source ~/.bashrc
conda activate fastplong_env

BASE_DIR="/path/to/project"
FASTQ_DIR="${BASE_DIR}/fastq"
OUTPUT_DIR="${BASE_DIR}/fastplong_output"
REPORT_DIR="${BASE_DIR}/fastplong_reports"
LOG_DIR="${BASE_DIR}/logs_fastplong"
PROCESSED_LOG="${LOG_DIR}/processed_fastplong.log"

MAX_PARALLEL_JOBS=3
CORES_PER_JOB=8

mkdir -p "$OUTPUT_DIR" "$REPORT_DIR" "$LOG_DIR"
touch "$PROCESSED_LOG"

process_sample() {
    local INPUT_FASTQ=$1
    local CORES=$2

    local BASENAME
    BASENAME=$(basename "$INPUT_FASTQ" .fastq.gz)

    local LOG_FILE="${LOG_DIR}/${BASENAME}.log"

    {
        echo "--- Starting fastplong for ${BASENAME} at $(date) ---"

        fastplong \
            --disable_quality_filtering \
            -i "$INPUT_FASTQ" \
            -o "${OUTPUT_DIR}/${BASENAME}.filtered.fastq.gz" \
            -h "${REPORT_DIR}/${BASENAME}.html" \
            -j "${REPORT_DIR}/${BASENAME}.json" \
            -w "$CORES"

        if [ $? -ne 0 ]; then
            echo "ERROR: fastplong failed for ${BASENAME}"
            return 1
        fi

        echo "$BASENAME" >> "$PROCESSED_LOG"
        echo "--- Finished fastplong for ${BASENAME} at $(date) ---"

    } >> "$LOG_FILE" 2>&1
}

export -f process_sample
export FASTQ_DIR OUTPUT_DIR REPORT_DIR LOG_DIR PROCESSED_LOG

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

for FASTQ in "$FASTQ_DIR"/*.fastq.gz; do
    if [[ ! -f "$FASTQ" ]]; then
        echo "No FASTQ files found"
        break
    fi

    BASENAME=$(basename "$FASTQ" .fastq.gz)

    if grep -Fxq "$BASENAME" "$PROCESSED_LOG"; then
        echo "Skipping ${BASENAME}, already processed."
        continue
    fi

    echo "Starting fastplong for ${BASENAME}"
    wait_for_slot
    process_sample "$FASTQ" "$CORES_PER_JOB" &

done

wait
echo "All fastplong jobs completed."
