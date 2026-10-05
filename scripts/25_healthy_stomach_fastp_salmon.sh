#!/bin/bash
#SBATCH --job-name=gastric_salmon
#SBATCH --output=logs/gastric_%A_%a.out
#SBATCH --error=logs/gastric_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --ntasks=8
#SBATCH --mem=40G
#SBATCH --partition=compute

set -euo pipefail

###############################################################################
# PRJNA531171
#
# Included:
#   - Low-BMI, normal gastric motility
#   - High-BMI, normal gastric motility
#
# Excluded:
#   - Low-BMI gastroparesis
#   - High-BMI gastroparesis
###############################################################################

module load miniforge
conda activate chitebuster_py311

THREADS="${SLURM_NTASKS:-8}"

BASE_DIR="/path/to/healthy_tissues/PRJNA531171_gastric"

MANIFEST="${BASE_DIR}/PRJNA531171_gastric_manifest.tsv"

RAW_DIR="${BASE_DIR}/raw_reads"
TRIM_DIR="${BASE_DIR}/trimmed_reads"
FASTP_DIR="${BASE_DIR}/fastp_reports"
SALMON_DIR="${BASE_DIR}/salmon_quant"
LOG_DIR="${BASE_DIR}/logs"

###############################################################################
# EDIT THIS PATH
###############################################################################

SALMON_INDEX="/path/to/salmon_index"

###############################################################################
# Checks
###############################################################################

mkdir -p \
    "${RAW_DIR}" \
    "${TRIM_DIR}" \
    "${FASTP_DIR}" \
    "${SALMON_DIR}" \
    "${LOG_DIR}"

for PROGRAM in wget fastp salmon; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} is not available." >&2
        exit 1
    fi
done

if [[ ! -d "${SALMON_INDEX}" ]]; then
    echo "ERROR: Salmon index does not exist:"
    echo "${SALMON_INDEX}"
    exit 1
fi

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Manifest does not exist or is empty:"
    echo "${MANIFEST}"
    exit 1
fi

if [[ -z "${SLURM_ARRAY_TASK_ID:-}" ]]; then
    echo "ERROR: Submit this script as a SLURM array job."
    exit 1
fi

###############################################################################
# Read one manifest row
#
# Expected columns:
# run_accession  sample_name  group  library_layout  fastq_ftp
###############################################################################

LINE=$(awk -F '\t' \
    -v row="$((SLURM_ARRAY_TASK_ID + 1))" \
    'NR == row {print; exit}' \
    "${MANIFEST}")

if [[ -z "${LINE}" ]]; then
    echo "ERROR: No manifest row for array task ${SLURM_ARRAY_TASK_ID}."
    exit 1
fi

IFS=$'\t' read -r \
    RUN \
    SAMPLE_NAME \
    GROUP \
    LIBRARY_LAYOUT \
    FASTQ_FTP \
    <<< "${LINE}"

RUN=$(echo "${RUN}" | tr -d '\r')
SAMPLE_NAME=$(echo "${SAMPLE_NAME}" | tr -d '\r')
GROUP=$(echo "${GROUP}" | tr -d '\r')
LIBRARY_LAYOUT=$(echo "${LIBRARY_LAYOUT}" | tr -d '\r')
FASTQ_FTP=$(echo "${FASTQ_FTP}" | tr -d '\r')

if [[ -z "${RUN}" || -z "${SAMPLE_NAME}" || -z "${FASTQ_FTP}" ]]; then
    echo "ERROR: Missing required manifest information."
    echo "RUN=${RUN}"
    echo "SAMPLE_NAME=${SAMPLE_NAME}"
    echo "GROUP=${GROUP}"
    echo "LIBRARY_LAYOUT=${LIBRARY_LAYOUT}"
    echo "FASTQ_FTP=${FASTQ_FTP}"
    exit 1
fi

case "${GROUP}" in
    Gastric_LowBMI|Gastric_HighBMI)
        ;;
    *)
        echo "ERROR: Unexpected group: ${GROUP}"
        echo "Allowed groups: Gastric_LowBMI or Gastric_HighBMI"
        exit 1
        ;;
esac

if [[ "${LIBRARY_LAYOUT^^}" != "PAIRED" ]]; then
    echo "ERROR: Expected paired-end data but found ${LIBRARY_LAYOUT}"
    exit 1
fi

###############################################################################
# Extract paired ENA URLs
###############################################################################

IFS=';' read -ra URL_ARRAY <<< "${FASTQ_FTP}"

if [[ "${#URL_ARRAY[@]}" -ne 2 ]]; then
    echo "ERROR: Expected two FASTQ URLs but found ${#URL_ARRAY[@]}."
    echo "${FASTQ_FTP}"
    exit 1
fi

URL1="${URL_ARRAY[0]}"
URL2="${URL_ARRAY[1]}"

[[ "${URL1}" =~ ^https?://|^ftp:// ]] || URL1="https://${URL1}"
[[ "${URL2}" =~ ^https?://|^ftp:// ]] || URL2="https://${URL2}"

RAW_R1="${RAW_DIR}/${SAMPLE_NAME}_${RUN}_1.fastq.gz"
RAW_R2="${RAW_DIR}/${SAMPLE_NAME}_${RUN}_2.fastq.gz"

TRIM_R1="${TRIM_DIR}/${SAMPLE_NAME}_${RUN}_trimmed_1.fastq.gz"
TRIM_R2="${TRIM_DIR}/${SAMPLE_NAME}_${RUN}_trimmed_2.fastq.gz"

HTML_REPORT="${FASTP_DIR}/${SAMPLE_NAME}_${RUN}_fastp.html"
JSON_REPORT="${FASTP_DIR}/${SAMPLE_NAME}_${RUN}_fastp.json"

SAMPLE_QUANT_DIR="${SALMON_DIR}/${SAMPLE_NAME}_${RUN}"

###############################################################################
# Skip completed samples
###############################################################################

if [[ -s "${SAMPLE_QUANT_DIR}/quant.sf" ]]; then
    echo "Sample already completed: ${SAMPLE_NAME}_${RUN}"
    exit 0
fi

echo "============================================================"
echo "Dataset: PRJNA531171"
echo "Run: ${RUN}"
echo "Sample: ${SAMPLE_NAME}"
echo "Group: ${GROUP}"
echo "Array task: ${SLURM_ARRAY_TASK_ID}"
echo "============================================================"

###############################################################################
# Download
###############################################################################

echo "Downloading read 1:"
echo "${URL1}"

wget \
    --continue \
    --tries=5 \
    --timeout=60 \
    --waitretry=10 \
    --output-document="${RAW_R1}" \
    "${URL1}"

echo "Downloading read 2:"
echo "${URL2}"

wget \
    --continue \
    --tries=5 \
    --timeout=60 \
    --waitretry=10 \
    --output-document="${RAW_R2}" \
    "${URL2}"

if [[ ! -s "${RAW_R1}" || ! -s "${RAW_R2}" ]]; then
    echo "ERROR: One or both downloaded FASTQ files are empty."
    exit 1
fi

gzip -t "${RAW_R1}"
gzip -t "${RAW_R2}"

###############################################################################
# Trim
###############################################################################

echo "Running fastp..."

fastp \
    --in1 "${RAW_R1}" \
    --in2 "${RAW_R2}" \
    --out1 "${TRIM_R1}" \
    --out2 "${TRIM_R2}" \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 30 \
    --thread "${THREADS}" \
    --html "${HTML_REPORT}" \
    --json "${JSON_REPORT}"

if [[ ! -s "${TRIM_R1}" || ! -s "${TRIM_R2}" ]]; then
    echo "ERROR: fastp did not produce both trimmed FASTQs."
    exit 1
fi

###############################################################################
# Salmon
###############################################################################

echo "Running Salmon..."

mkdir -p "${SAMPLE_QUANT_DIR}"

salmon quant \
    --index "${SALMON_INDEX}" \
    --libType A \
    --mates1 "${TRIM_R1}" \
    --mates2 "${TRIM_R2}" \
    --threads "${THREADS}" \
    --validateMappings \
    --seqBias \
    --gcBias \
    --posBias \
    --numGibbsSamples 50 \
    --output "${SAMPLE_QUANT_DIR}"

###############################################################################
# Verify output and clean FASTQs
###############################################################################

if [[ -s "${SAMPLE_QUANT_DIR}/quant.sf" ]]; then
    echo "Salmon quantification completed successfully."

    rm -f \
        "${RAW_R1}" \
        "${RAW_R2}" \
        "${TRIM_R1}" \
        "${TRIM_R2}"

    touch "${SAMPLE_QUANT_DIR}/PROCESSING_COMPLETE.txt"

    echo "FASTQ files deleted successfully."
else
    echo "ERROR: quant.sf was not produced."
    echo "FASTQ files have been retained."
    exit 1
fi

echo "Completed: ${SAMPLE_NAME}_${RUN}"
