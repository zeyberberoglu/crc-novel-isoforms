#!/bin/bash
#SBATCH --job-name=cellline_fastp_salmon
#SBATCH --output=logs/cellline_fastp_salmon_%j.out
#SBATCH --error=logs/cellline_fastp_salmon_%j.err
#SBATCH --time=240:00:00
#SBATCH --ntasks=32
#SBATCH --nodes=1
#SBATCH --mem=120G
#SBATCH --partition=compute

set -euo pipefail

cd /path/to/cell_lines || exit 1

module load miniforge
conda activate chitebuster_py311
module load sra-tools

BASE_DIR="/path/to/cell_lines"
CSV_FILE="${BASE_DIR}/comprehensive_sample_table_cell_lines.csv"

RAW_DIR="${BASE_DIR}/raw_reads"
TRIM_DIR="${BASE_DIR}/trimmed_reads"
FASTP_DIR="${BASE_DIR}/fastp_reports"
SALMON_DIR="${BASE_DIR}/salmon_quant"
LOG_DIR="${BASE_DIR}/logs"

SALMON_INDEX="/path/to/salmon_index"

THREADS=32

mkdir -p "$RAW_DIR" "$TRIM_DIR" "$FASTP_DIR" "$SALMON_DIR" "$LOG_DIR"

MANIFEST="${BASE_DIR}/cell_line_manifest.tsv"

python - <<EOF
import csv

csv_file = "$CSV_FILE"
out_file = "$MANIFEST"

with open(csv_file, newline="") as f, open(out_file, "w") as out:
    reader = csv.DictReader(f)
    out.write("Run\tSample\tLayout\n")
    for row in reader:
        run = row["Run"].strip()
        sample = row["Unified_Sample_Name"].strip()
        layout = row["LibraryLayout"].strip()
        if run and sample:
            out.write(f"{run}\\t{sample}\\t{layout}\\n")
EOF

tail -n +2 "$MANIFEST" | while IFS=$'\t' read -r RUN SAMPLE LAYOUT
do
    echo "======================================"
    echo "Processing ${RUN} as ${SAMPLE}"
    echo "Layout: ${LAYOUT}"
    echo "======================================"

    if [ "$LAYOUT" != "PAIRED" ]; then
        echo "Skipping ${SAMPLE}: not paired-end"
        continue
    fi

    RAW_1="${RAW_DIR}/${SAMPLE}_1.fastq.gz"
    RAW_2="${RAW_DIR}/${SAMPLE}_2.fastq.gz"

    TRIM_1="${TRIM_DIR}/${SAMPLE}_1.trimmed.fastq.gz"
    TRIM_2="${TRIM_DIR}/${SAMPLE}_2.trimmed.fastq.gz"

    SALMON_OUT="${SALMON_DIR}/${SAMPLE}"

    if [ -f "${SALMON_OUT}/quant.sf" ]; then
        echo "Salmon already completed for ${SAMPLE}, skipping."
        continue
    fi

    echo "Downloading ${RUN}"

    fasterq-dump "$RUN" \
        --split-files \
        --threads "$THREADS" \
        --outdir "$RAW_DIR"

    gzip -f "${RAW_DIR}/${RUN}_1.fastq"
    gzip -f "${RAW_DIR}/${RUN}_2.fastq"

    mv "${RAW_DIR}/${RUN}_1.fastq.gz" "$RAW_1"
    mv "${RAW_DIR}/${RUN}_2.fastq.gz" "$RAW_2"

    echo "Running fastp for ${SAMPLE}"

    fastp \
        -i "$RAW_1" \
        -I "$RAW_2" \
        -o "$TRIM_1" \
        -O "$TRIM_2" \
        --thread "$THREADS" \
        --detect_adapter_for_pe \
        --html "${FASTP_DIR}/${SAMPLE}.fastp.html" \
        --json "${FASTP_DIR}/${SAMPLE}.fastp.json"

    echo "Running Salmon for ${SAMPLE}"

    salmon quant \
        -i "$SALMON_INDEX" \
        -l A \
        -1 "$TRIM_1" \
        -2 "$TRIM_2" \
        -p "$THREADS" \
        --seqBias \
        --gcBias \
        --posBias \
        -o "$SALMON_OUT"

    if [ -f "${SALMON_OUT}/quant.sf" ]; then
        echo "Salmon successful for ${SAMPLE}. Deleting raw FASTQs."
        rm -f "$RAW_1" "$RAW_2"
    else
        echo "ERROR: Salmon failed for ${SAMPLE}. Raw FASTQs kept."
        exit 1
    fi

    echo "Done with ${SAMPLE}"

done

echo "All cell line samples completed."
