#!/bin/bash
#SBATCH --job-name=frac_salmon
#SBATCH --output=logs/fractionation_fastp_salmon_%j.out
#SBATCH --error=logs/fractionation_fastp_salmon_%j.err
#SBATCH --time=48:00:00
#SBATCH --ntasks=16
#SBATCH --nodes=1
#SBATCH --mem=64G
#SBATCH --partition=compute

set -euo pipefail


# ============================================================
# PATHS AND SETTINGS
# ============================================================

BASE_DIR="/path/to/PRJNA397552_fractionation"

MANIFEST="${BASE_DIR}/metadata/fractionation_manifest.tsv"

RAW_DIR="${BASE_DIR}/raw_reads"
TRIM_DIR="${BASE_DIR}/trimmed_reads"
FASTP_DIR="${BASE_DIR}/fastp_reports"
SALMON_DIR="${BASE_DIR}/salmon_quant"
LOG_DIR="${BASE_DIR}/logs"

SALMON_INDEX="/path/to/salmon_index"

THREADS=16


# ============================================================
# ENVIRONMENT
# ============================================================

module load miniforge
conda activate chitebuster_py311


# ============================================================
# CREATE OUTPUT DIRECTORIES
# ============================================================

mkdir -p \
    "${RAW_DIR}" \
    "${TRIM_DIR}" \
    "${FASTP_DIR}" \
    "${SALMON_DIR}" \
    "${LOG_DIR}"


# ============================================================
# CHECK INPUTS AND PROGRAMS
# ============================================================

if [[ ! -s "${MANIFEST}" ]]; then
    echo "ERROR: Manifest not found or empty:"
    echo "${MANIFEST}"
    exit 1
fi

if [[ ! -d "${SALMON_INDEX}" ]]; then
    echo "ERROR: Salmon index not found:"
    echo "${SALMON_INDEX}"
    exit 1
fi

for PROGRAM in wget fastp salmon gzip; do
    if ! command -v "${PROGRAM}" >/dev/null 2>&1; then
        echo "ERROR: ${PROGRAM} is not available in longread_de."
        exit 1
    fi
done


echo "============================================================"
echo "Fractionation fastp + Salmon pipeline"
echo "Start time:    $(date)"
echo "Manifest:      ${MANIFEST}"
echo "Salmon index:  ${SALMON_INDEX}"
echo "Environment:   ${CONDA_DEFAULT_ENV:-unknown}"
echo "Threads:       ${THREADS}"
echo "============================================================"


# ============================================================
# COUNT MANIFEST SAMPLES
# ============================================================

TOTAL_SAMPLES=$(
    awk 'NR > 1 && NF > 0 {count++} END {print count+0}' \
        "${MANIFEST}"
)

echo "Samples in manifest: ${TOTAL_SAMPLES}"


# ============================================================
# PROCESS EACH SAMPLE
# ============================================================
#
# Manifest columns:
#
# 1  Run
# 2  Sample
# 3  LibraryType
# 4  CellLine
# 5  Fraction
# 6  Replicate
# 7  Layout
# 8  FastqFTP
# 9  FastqMD5
# 10 SampleTitle
# 11 ExperimentTitle
#
# All selected samples are single-end.
# ============================================================

while IFS=$'\t' read -r \
    RUN \
    SAMPLE \
    LIBRARY_TYPE \
    CELL_LINE \
    FRACTION \
    REPLICATE \
    LAYOUT \
    FASTQ_FTP \
    FASTQ_MD5 \
    SAMPLE_TITLE \
    EXPERIMENT_TITLE
do

    # Remove possible carriage-return characters.
    RUN="${RUN//$'\r'/}"
    SAMPLE="${SAMPLE//$'\r'/}"
    LIBRARY_TYPE="${LIBRARY_TYPE//$'\r'/}"
    CELL_LINE="${CELL_LINE//$'\r'/}"
    FRACTION="${FRACTION//$'\r'/}"
    REPLICATE="${REPLICATE//$'\r'/}"
    LAYOUT="${LAYOUT//$'\r'/}"
    FASTQ_FTP="${FASTQ_FTP//$'\r'/}"
    FASTQ_MD5="${FASTQ_MD5//$'\r'/}"

    if [[ -z "${RUN}" || -z "${SAMPLE}" ]]; then
        echo "WARNING: Skipping incomplete manifest row."
        continue
    fi

    echo
    echo "============================================================"
    echo "Run:          ${RUN}"
    echo "Sample:       ${SAMPLE}"
    echo "Library type: ${LIBRARY_TYPE}"
    echo "Cell line:    ${CELL_LINE}"
    echo "Fraction:     ${FRACTION}"
    echo "Replicate:    ${REPLICATE}"
    echo "Layout:       ${LAYOUT}"
    echo "Start time:   $(date)"
    echo "============================================================"


    # --------------------------------------------------------
    # Confirm single-end layout
    # --------------------------------------------------------

    if [[ "${LAYOUT}" != "SINGLE" ]]; then
        echo "ERROR: Sample is not marked as SINGLE:"
        echo "${SAMPLE}: ${LAYOUT}"
        exit 1
    fi


    # --------------------------------------------------------
    # Check ENA FASTQ URL
    # --------------------------------------------------------

    if [[ -z "${FASTQ_FTP}" ]]; then
        echo "ERROR: FastqFTP is empty for ${SAMPLE}."
        exit 1
    fi

    # Single-end runs should have one URL.
    # Stop if multiple semicolon-separated URLs are present.
    if [[ "${FASTQ_FTP}" == *";"* ]]; then
        echo "ERROR: Multiple FASTQ URLs detected for single-end sample:"
        echo "${SAMPLE}: ${FASTQ_FTP}"
        exit 1
    fi

    # ENA often supplies:
    # ftp.sra.ebi.ac.uk/vol1/fastq/...
    #
    # Convert it to:
    # https://ftp.sra.ebi.ac.uk/vol1/fastq/...
    if [[ "${FASTQ_FTP}" =~ ^https?:// ]]; then
        DOWNLOAD_URL="${FASTQ_FTP}"
    elif [[ "${FASTQ_FTP}" =~ ^ftp:// ]]; then
        DOWNLOAD_URL="${FASTQ_FTP/ftp:\/\//https://}"
    else
        DOWNLOAD_URL="https://${FASTQ_FTP}"
    fi


    # --------------------------------------------------------
    # Output filenames
    # --------------------------------------------------------

    SOURCE_FILENAME=$(basename "${FASTQ_FTP}")

    RAW_FASTQ_GZ="${RAW_DIR}/${RUN}.fastq.gz"
    DOWNLOADED_FILE="${RAW_DIR}/${SOURCE_FILENAME}"

    TRIMMED_FASTQ="${TRIM_DIR}/${SAMPLE}.trimmed.fastq.gz"

    FASTP_HTML="${FASTP_DIR}/${SAMPLE}_fastp.html"
    FASTP_JSON="${FASTP_DIR}/${SAMPLE}_fastp.json"
    FASTP_LOG="${FASTP_DIR}/${SAMPLE}_fastp.log"

    SAMPLE_SALMON_DIR="${SALMON_DIR}/${SAMPLE}"
    QUANT_FILE="${SAMPLE_SALMON_DIR}/quant.sf"


    # --------------------------------------------------------
    # Skip fully completed samples
    # --------------------------------------------------------

    if [[ -s "${QUANT_FILE}" ]]; then
        echo "Valid quant.sf already exists. Skipping ${SAMPLE}."
        continue
    fi


    # --------------------------------------------------------
    # Download FASTQ with wget
    # --------------------------------------------------------

    if [[ ! -s "${RAW_FASTQ_GZ}" ]]; then

        echo "Downloading FASTQ with wget..."
        echo "${DOWNLOAD_URL}"

        wget \
            --continue \
            --tries=5 \
            --timeout=60 \
            --read-timeout=60 \
            --output-document="${DOWNLOADED_FILE}" \
            "${DOWNLOAD_URL}"

        if [[ ! -s "${DOWNLOADED_FILE}" ]]; then
            echo "ERROR: Download failed or created an empty file:"
            echo "${DOWNLOADED_FILE}"
            exit 1
        fi

        # Standardise the raw filename to RUN.fastq.gz.
        if [[ "${DOWNLOADED_FILE}" != "${RAW_FASTQ_GZ}" ]]; then
            mv -f "${DOWNLOADED_FILE}" "${RAW_FASTQ_GZ}"
        fi

    else
        echo "Raw FASTQ already exists:"
        echo "${RAW_FASTQ_GZ}"
    fi


    # --------------------------------------------------------
    # Validate gzip file
    # --------------------------------------------------------

    echo "Checking FASTQ gzip integrity..."

    if ! gzip -t "${RAW_FASTQ_GZ}"; then
        echo "ERROR: Corrupted gzip FASTQ:"
        echo "${RAW_FASTQ_GZ}"
        rm -f "${RAW_FASTQ_GZ}"
        exit 1
    fi


    # --------------------------------------------------------
    # Optional ENA MD5 validation
    # --------------------------------------------------------

    if [[ -n "${FASTQ_MD5}" ]]; then

        OBSERVED_MD5=$(md5sum "${RAW_FASTQ_GZ}" | awk '{print $1}')

        if [[ "${OBSERVED_MD5}" != "${FASTQ_MD5}" ]]; then
            echo "ERROR: MD5 mismatch for ${SAMPLE}"
            echo "Expected: ${FASTQ_MD5}"
            echo "Observed: ${OBSERVED_MD5}"
            rm -f "${RAW_FASTQ_GZ}"
            exit 1
        fi

        echo "MD5 check passed."
    else
        echo "No MD5 value available; skipping MD5 validation."
    fi


    # --------------------------------------------------------
    # Run fastp
    # --------------------------------------------------------

    if [[ ! -s "${TRIMMED_FASTQ}" ]]; then

        echo "Running fastp..."

        fastp \
            --in1 "${RAW_FASTQ_GZ}" \
            --out1 "${TRIMMED_FASTQ}" \
            --thread "${THREADS}" \
            --qualified_quality_phred 20 \
            --length_required 20 \
            --trim_poly_g \
            --trim_poly_x \
            --html "${FASTP_HTML}" \
            --json "${FASTP_JSON}" \
            > "${FASTP_LOG}" 2>&1

    else
        echo "Trimmed FASTQ already exists:"
        echo "${TRIMMED_FASTQ}"
    fi


    if [[ ! -s "${TRIMMED_FASTQ}" ]]; then
        echo "ERROR: fastp did not create a valid trimmed FASTQ:"
        echo "${TRIMMED_FASTQ}"
        exit 1
    fi


    # --------------------------------------------------------
    # Run Salmon
    # --------------------------------------------------------

    echo "Running Salmon quantification..."

    rm -rf "${SAMPLE_SALMON_DIR}"

    salmon quant \
        --index "${SALMON_INDEX}" \
        --libType A \
        --unmatedReads "${TRIMMED_FASTQ}" \
        --threads "${THREADS}" \
        --validateMappings \
        --seqBias \
        --gcBias \
        --numGibbsSamples 50 \
        --output "${SAMPLE_SALMON_DIR}"


    if [[ ! -s "${QUANT_FILE}" ]]; then
        echo "ERROR: Salmon failed for ${SAMPLE}."
        exit 1
    fi


    # --------------------------------------------------------
    # Save sample metadata
    # --------------------------------------------------------

    cat > "${SAMPLE_SALMON_DIR}/sample_metadata.tsv" <<EOF
Run	Sample	LibraryType	CellLine	Fraction	Replicate	Layout
${RUN}	${SAMPLE}	${LIBRARY_TYPE}	${CELL_LINE}	${FRACTION}	${REPLICATE}	${LAYOUT}
EOF


    # --------------------------------------------------------
    # Delete raw FASTQ after successful Salmon quantification
    # --------------------------------------------------------

    echo "Salmon completed successfully."

    rm -f "${RAW_FASTQ_GZ}"

    echo "Finished ${SAMPLE}: $(date)"

done < <(tail -n +2 "${MANIFEST}")


# ============================================================
# FINAL CHECK
# ============================================================

COMPLETED=$(
    find "${SALMON_DIR}" \
        -mindepth 2 \
        -maxdepth 2 \
        -type f \
        -name "quant.sf" \
        -size +0c |
    wc -l
)

echo
echo "============================================================"
echo "Pipeline finished:          $(date)"
echo "Expected quantifications:   ${TOTAL_SAMPLES}"
echo "Completed quantifications:  ${COMPLETED}"
echo "============================================================"

if [[ "${COMPLETED}" -ne "${TOTAL_SAMPLES}" ]]; then
    echo "ERROR: Not all samples have a valid quant.sf."
    exit 1
fi

echo "All 12 fractionation samples were quantified successfully."
