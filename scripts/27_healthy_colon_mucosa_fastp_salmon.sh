#!/bin/bash
#SBATCH --job-name=colon_mucosa
#SBATCH --output=logs/colon_mucosa_%A_%a.out
#SBATCH --error=logs/colon_mucosa_%A_%a.err
#SBATCH --partition=compute
#SBATCH --time=18:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --array=1-8%4

set -euo pipefail

module load miniforge
conda activate chitebuster_py311

BASE="/path/to/healthy_tissues/PRJNA1085261"

MANIFEST="${BASE}/metadata/healthy_colon_mucosa_manifest.tsv"
RAW_DIR="${BASE}/raw_reads"
TRIM_DIR="${BASE}/trimmed_reads"
FASTP_DIR="${BASE}/fastp_reports"
SALMON_DIR="${BASE}/salmon_quant"
LOG_DIR="${BASE}/logs"

# Replace this with your actual existing Salmon index
SALMON_INDEX="/path/to/salmon_index"

mkdir -p \
    "${RAW_DIR}" \
    "${TRIM_DIR}" \
    "${FASTP_DIR}" \
    "${SALMON_DIR}" \
    "${LOG_DIR}"

command -v wget >/dev/null 2>&1 || {
    echo "ERROR: wget not found"
    exit 1
}

command -v fastp >/dev/null 2>&1 || {
    echo "ERROR: fastp not found in chitebuster_py311"
    exit 1
}

command -v salmon >/dev/null 2>&1 || {
    echo "ERROR: salmon not found in chitebuster_py311"
    exit 1
}

[[ -s "${MANIFEST}" ]] || {
    echo "ERROR: Manifest not found: ${MANIFEST}"
    exit 1
}

[[ -d "${SALMON_INDEX}" ]] || {
    echo "ERROR: Salmon index not found: ${SALMON_INDEX}"
    exit 1
}

LINE=$(tail -n +2 "${MANIFEST}" | sed -n "${SLURM_ARRAY_TASK_ID}p")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No manifest row for task ${SLURM_ARRAY_TASK_ID}"
    exit 1
fi

IFS=$'\t' read -r SAMPLE RUN LAYOUT FASTQ_FTP <<< "${LINE}"

if [[ "${LAYOUT}" != "PAIRED" ]]; then
    echo "ERROR: Expected PAIRED layout but found: ${LAYOUT}"
    exit 1
fi

IFS=';' read -r FTP_R1 FTP_R2 <<< "${FASTQ_FTP}"

if [[ -z "${FTP_R1}" || -z "${FTP_R2}" ]]; then
    echo "ERROR: Missing paired FASTQ URLs for ${RUN}"
    exit 1
fi

URL_R1="https://${FTP_R1}"
URL_R2="https://${FTP_R2}"

SAMPLE_RAW="${RAW_DIR}/${SAMPLE}"
SAMPLE_TRIM="${TRIM_DIR}/${SAMPLE}"
SAMPLE_QUANT="${SALMON_DIR}/${SAMPLE}"

mkdir -p "${SAMPLE_RAW}" "${SAMPLE_TRIM}"

RAW_R1="${SAMPLE_RAW}/${RUN}_1.fastq.gz"
RAW_R2="${SAMPLE_RAW}/${RUN}_2.fastq.gz"

TRIM_R1="${SAMPLE_TRIM}/${SAMPLE}_trimmed_R1.fastq.gz"
TRIM_R2="${SAMPLE_TRIM}/${SAMPLE}_trimmed_R2.fastq.gz"

FASTP_HTML="${FASTP_DIR}/${SAMPLE}_fastp.html"
FASTP_JSON="${FASTP_DIR}/${SAMPLE}_fastp.json"
FASTP_LOG="${FASTP_DIR}/${SAMPLE}_fastp.log"

echo "========================================"
echo "Sample: ${SAMPLE}"
echo "Run: ${RUN}"
echo "Layout: ${LAYOUT}"
echo "R1: ${URL_R1}"
echo "R2: ${URL_R2}"
echo "========================================"

# ---------------------------------------------------------
# Download raw paired FASTQs
# ---------------------------------------------------------

if [[ ! -s "${RAW_R1}" ]]; then
    wget \
        --continue \
        --tries=10 \
        --timeout=60 \
        --retry-connrefused \
        -O "${RAW_R1}" \
        "${URL_R1}"
else
    echo "R1 already exists; skipping download."
fi

if [[ ! -s "${RAW_R2}" ]]; then
    wget \
        --continue \
        --tries=10 \
        --timeout=60 \
        --retry-connrefused \
        -O "${RAW_R2}" \
        "${URL_R2}"
else
    echo "R2 already exists; skipping download."
fi

echo "Checking downloaded gzip files..."

gzip -t "${RAW_R1}"
gzip -t "${RAW_R2}"

# ---------------------------------------------------------
# Trim with fastp
# ---------------------------------------------------------

if [[ ! -s "${TRIM_R1}" || ! -s "${TRIM_R2}" ]]; then

    fastp \
        --in1 "${RAW_R1}" \
        --in2 "${RAW_R2}" \
        --out1 "${TRIM_R1}" \
        --out2 "${TRIM_R2}" \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 30 \
        --cut_front \
        --cut_tail \
        --cut_window_size 4 \
        --cut_mean_quality 20 \
        --thread "${SLURM_CPUS_PER_TASK}" \
        --html "${FASTP_HTML}" \
        --json "${FASTP_JSON}" \
        2> "${FASTP_LOG}"

else
    echo "Trimmed FASTQs already exist; skipping fastp."
fi

gzip -t "${TRIM_R1}"
gzip -t "${TRIM_R2}"

# ---------------------------------------------------------
# Salmon quantification
# ---------------------------------------------------------

if [[ ! -s "${SAMPLE_QUANT}/quant.sf" ]]; then

    salmon quant \
        --index "${SALMON_INDEX}" \
        --libType A \
        --mates1 "${TRIM_R1}" \
        --mates2 "${TRIM_R2}" \
        --threads "${SLURM_CPUS_PER_TASK}" \
        --validateMappings \
        --seqBias \
        --gcBias \
        --posBias \
        --numGibbsSamples 50 \
        --output "${SAMPLE_QUANT}"

else
    echo "quant.sf already exists; skipping Salmon."
fi
    
# ---------------------------------------------------------
# Validate Salmon output and delete FASTQ files
# ---------------------------------------------------------

if [[ ! -s "${SAMPLE_QUANT}/quant.sf" ]]; then
    echo "ERROR: Salmon quantification failed for ${SAMPLE}"
    echo "FASTQ files will be retained."
    exit 1
fi

echo "Salmon quantification completed successfully."

# Delete both raw and trimmed reads to save storage
echo "Deleting raw FASTQ files..."
rm -f "${RAW_R1}" "${RAW_R2}"

echo "Deleting trimmed FASTQ files..."
rm -f "${TRIM_R1}" "${TRIM_R2}"

# Remove empty per-sample directories
rmdir "${SAMPLE_RAW}" 2>/dev/null || true
rmdir "${SAMPLE_TRIM}" 2>/dev/null || true

echo -e "${SAMPLE}\t${RUN}" >> "${BASE}/completed_samples.tsv"

echo "Successfully completed ${SAMPLE}"
echo "Raw and trimmed FASTQ files were deleted."
