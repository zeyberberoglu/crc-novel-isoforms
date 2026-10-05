#!/bin/bash
#SBATCH --job-name=fastp_salmon
#SBATCH --output=logs/fastp_salmon_%j.out
#SBATCH --error=logs/fastp_salmon_%j.err
#SBATCH --time=240:00:00
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=90G
#SBATCH --partition=compute

set -euo pipefail

module load miniforge
conda activate chitebuster_py311
module load sra-tools

# ==============================
# User-defined paths
# ==============================

BASE_DIR="/path/to/PRJNA658001"
SAMPLE_INFO="${BASE_DIR}/run_to_new_name.tsv"

RAW_DIR="${BASE_DIR}/raw_reads"
TRIM_DIR="${BASE_DIR}/trimmed_reads"
FASTP_DIR="${BASE_DIR}/fastp_reports"
SALMON_DIR="${BASE_DIR}/salmon_quant"
LOG_DIR="${BASE_DIR}/logs"

SALMON_INDEX="/path/to/salmon_index"

mkdir -p \
    "$RAW_DIR" \
    "$TRIM_DIR" \
    "$FASTP_DIR" \
    "$SALMON_DIR" \
    "$LOG_DIR"

cd "$BASE_DIR" || exit 1

# ==============================
# Process samples
# ==============================

while IFS=$'\t' read -r RUN SAMPLE
do
    echo "======================================"
    echo "Processing $RUN as $SAMPLE"
    echo "======================================"

    if [ -z "$RUN" ] || [ -z "$SAMPLE" ]; then
        echo "Skipping empty line"
        continue
    fi

    if [ -f "${SALMON_DIR}/${SAMPLE}/quant.sf" ]; then
        echo "Salmon already completed for $SAMPLE, skipping."
        continue
    fi

    echo "Downloading $RUN..."
    prefetch "$RUN" \
        --output-directory "$RAW_DIR"

    echo "Converting $RUN to FASTQ..."
    fasterq-dump "$RAW_DIR/$RUN" \
        --split-files \
        --threads 8 \
        --outdir "$RAW_DIR"

    echo "Compressing FASTQs..."
    pigz -p 8 \
        "${RAW_DIR}/${RUN}_1.fastq" \
        "${RAW_DIR}/${RUN}_2.fastq"

    echo "Running fastp..."
    fastp \
        -i "${RAW_DIR}/${RUN}_1.fastq.gz" \
        -I "${RAW_DIR}/${RUN}_2.fastq.gz" \
        -o "${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz" \
        -O "${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz" \
        --html "${FASTP_DIR}/${SAMPLE}.fastp.html" \
        --json "${FASTP_DIR}/${SAMPLE}.fastp.json" \
        --thread 8

    echo "Running Salmon..."
    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz" \
        -2 "${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz" \
        -p 16 \
        --validateMappings \
        -o "${SALMON_DIR}/${SAMPLE}"

    echo "Checking Salmon output..."

    if [ -s "${SALMON_DIR}/${SAMPLE}/quant.sf" ]; then

        echo "Salmon completed for $SAMPLE. Deleting FASTQ files..."

        rm -f "${RAW_DIR}/${RUN}_1.fastq.gz"
        rm -f "${RAW_DIR}/${RUN}_2.fastq.gz"
        rm -rf "${RAW_DIR}/${RUN}"

        rm -f "${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
        rm -f "${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"

        echo "$SAMPLE completed and FASTQs deleted."

    else

        echo "ERROR: Salmon failed for $SAMPLE. FASTQs were NOT deleted."
        exit 1

    fi

done < "$SAMPLE_INFO"

echo "All samples finished."
