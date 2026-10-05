#!/bin/bash
#SBATCH -o slurm-%j-fastp.out
#SBATCH --job-name=fastp_illumina
#SBATCH --time=120:00:00
#SBATCH --cpus-per-task=24
#SBATCH --nodes=1
#SBATCH --mem=90G

cd /data/scratch/ha25158/sra_illumina || exit 1

source ~/.bashrc
conda activate fastp_env

BASE_DIR="/data/scratch/ha25158/sra_illumina"
FASTQ_DIR="${BASE_DIR}/fastq"
TRIMMED_DIR="${BASE_DIR}/fastp_trimmed"
REPORT_DIR="${BASE_DIR}/fastp_reports"
LOG_DIR="${BASE_DIR}/logs_fastp"
PROCESSED_LOG="${LOG_DIR}/processed_fastp.log"

MAX_PARALLEL_JOBS=3
CORES_PER_JOB=8

mkdir -p "$TRIMMED_DIR" "$REPORT_DIR" "$LOG_DIR"
touch "$PROCESSED_LOG"

process_sample() {
    local R1=$1
    local CORES=$2

    local BASENAME
    BASENAME=$(basename "$R1" _1.fastq.gz)

    local R2="${FASTQ_DIR}/${BASENAME}_2.fastq.gz"

    local LOG_FILE="${LOG_DIR}/${BASENAME}.log"

    {
        echo "--- Starting fastp for ${BASENAME} at $(date) ---"

        fastp \
            --length_required 20 \
            --detect_adapter_for_pe \
            -x \
            -z 6 \
            -i "$R1" \
            -I "$R2" \
            -o "${TRIMMED_DIR}/${BASENAME}_R1.trimmed.fastq.gz" \
            -O "${TRIMMED_DIR}/${BASENAME}_R2.trimmed.fastq.gz" \
            -h "${REPORT_DIR}/${BASENAME}.html" \
            -j "${REPORT_DIR}/${BASENAME}.json" \
            -w "$CORES"

        if [ $? -ne 0 ]; then
            echo "ERROR: fastp failed for ${BASENAME}"
            return 1
        fi

        echo "$BASENAME" >> "$PROCESSED_LOG"
        echo "--- Finished fastp for ${BASENAME} at $(date) ---"

    } >> "$LOG_FILE" 2>&1
}

export -f process_sample
export FASTQ_DIR TRIMMED_DIR REPORT_DIR LOG_DIR PROCESSED_LOG

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

for R1 in "$FASTQ_DIR"/*_1.fastq.gz; do
    if [[ ! -f "$R1" ]]; then
        echo "No FASTQ files found"
        break
    fi

    BASENAME=$(basename "$R1" _1.fastq.gz)

    if grep -Fxq "$BASENAME" "$PROCESSED_LOG"; then
        echo "Skipping ${BASENAME}, already processed."
        continue
    fi

    echo "Starting fastp for ${BASENAME}"
    wait_for_slot
    process_sample "$R1" "$CORES_PER_JOB" &

done

wait
echo "All fastp jobs completed."
