#!/bin/bash
#SBATCH --job-name=novel_transcripts
#SBATCH --output=slurm-%j-novel_transcripts.out
#SBATCH --error=slurm-%j-novel_transcripts.err
#SBATCH --time=240:00:00
#SBATCH --cpus-per-task=10
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4G
#SBATCH --partition=compute

set -euo pipefail

# ==============================
# Conda environment
# ==============================

source ~/miniconda3/etc/profile.d/conda.sh
conda activate novel_annot

# ==============================
# User-defined paths
# ==============================

WORKDIR="/path/to/sqanti3_output"
REFDIR="/path/to/reference"
CPATDIR="/path/to/cpat_hg38"

cd "$WORKDIR"

# ==============================
# Compare transcripts to reference
# ==============================

echo "Running gffcompare..."

gffcompare \
  -r "${REFDIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3" \
  -o "${WORKDIR}/gffcomp" \
  "${WORKDIR}/rescue_refined/curated_transcriptome_refined_rescued.gtf"

# ==============================
# Extract multi-exon novel transcripts
# ==============================

echo "Extracting multi-exon novel transcript IDs..."

awk '$4 !~ /^(=|c|e|p|s|k|m|n)$/ {
    split($5,a,":");
    split(a[2],b,"|");
    if (b[3] > 1) print b[2]
}' "${WORKDIR}/gffcomp.tracking" \
  > "${WORKDIR}/multi_exon_novel_ids.txt"

echo "Number of multi-exon novel transcripts:"
wc -l "${WORKDIR}/multi_exon_novel_ids.txt"

echo "Extracting novel transcript GTF..."

awk 'FNR==NR {ids[$1]; next}
{
    tid="";
    for(i=9;i<=NF;i++) {
        if($i=="transcript_id") {
            tid=$(i+1);
            gsub(/[";]/, "", tid);
            break
        }
    }
    if(tid in ids) print $0
}' \
  "${WORKDIR}/multi_exon_novel_ids.txt" \
  "${WORKDIR}/gffcomp.annotated.gtf" \
  > "${WORKDIR}/novel_transcripts.gtf"

# ==============================
# ORFanage
# ==============================

echo "Running ORFanage..."

orfanage \
  --query "${WORKDIR}/novel_transcripts.gtf" \
  --reference "${REFDIR}/GRCh38.p14.genome.fa" \
  --output "${WORKDIR}/novel_transcripts_orfanage.gtf" \
  --threads 10 \
  --cleanq \
  --cleant \
  --rescue \
  --overhang 150 \
  --spliced_overhang \
  --mode LONGEST_MATCH \
  "${REFDIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3"

echo "Extracting transcripts with ORFanage ORFs..."

awk '$3=="CDS"' "${WORKDIR}/novel_transcripts_orfanage.gtf" \
  | sed 's/.*transcript_id "\([^"]*\)".*/\1/' \
  | sort | uniq \
  > "${WORKDIR}/with_orf_ids.txt"

# ==============================
# Extract transcripts without ORFs
# ==============================

awk 'FNR==NR {ids[$1]; next}
{
    tid="";
    for(i=9;i<=NF;i++) {
        if($i=="transcript_id") {
            tid=$(i+1);
            gsub(/[";]/, "", tid);
            break
        }
    }
    if(!(tid in ids)) print $0
}' \
  "${WORKDIR}/with_orf_ids.txt" \
  "${WORKDIR}/novel_transcripts.gtf" \
  > "${WORKDIR}/no_orf_transcripts.gtf"

# ==============================
# CPAT
# ==============================

if [[ ! -s "${WORKDIR}/no_orf_transcripts.gtf" ]]; then

  echo "No transcripts without ORFanage ORFs. Skipping CPAT."
  touch "${WORKDIR}/cpat_annotated_no_orf.gtf"

else

  echo "Extracting FASTA..."

  gffread \
    -w "${WORKDIR}/no_orf_transcripts.fasta" \
    -g "${REFDIR}/GRCh38.p14.genome.fa" \
    "${WORKDIR}/no_orf_transcripts.gtf"

  echo "Running CPAT..."

  cpat \
    --gene "${WORKDIR}/no_orf_transcripts.fasta" \
    --outfile "${WORKDIR}/cpat_novel_results" \
    --logitModel "${CPATDIR}/CPAT_logitModel.RData" \
    --hex "${CPATDIR}/CPAT_Hexamer.tsv" \
    --ref "${REFDIR}/GRCh38.p14.genome.fa" \
    --best-orf=p \
    --min-orf=75

  python add_cpat_cds.py \
    "${WORKDIR}/cpat_novel_results.ORF_prob.best.tsv" \
    "${WORKDIR}/no_orf_transcripts.gtf" \
    "${WORKDIR}/cpat_annotated_no_orf.gtf"

fi

# ==============================
# Final formatting
# ==============================

awk 'FNR==NR {ids[$1]; next}
{
    tid="";
    for(i=9;i<=NF;i++) {
        if($i=="transcript_id") {
            tid=$(i+1);
            gsub(/[";]/, "", tid);
            break
        }
    }
    if(tid in ids) print $0
}' \
  "${WORKDIR}/with_orf_ids.txt" \
  "${WORKDIR}/novel_transcripts_orfanage.gtf" \
  > "${WORKDIR}/orfanage_filtered.gtf"

python format_annotation.py \
  "${WORKDIR}/merged_assembly_rescued_qc_classification.txt" \
  "${WORKDIR}/orfanage_filtered.gtf" \
  "${WORKDIR}/cpat_annotated_no_orf.gtf" \
  "${REFDIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3" \
  "${WORKDIR}/gffcomp.annotated.gtf" \
  "${WORKDIR}/formatted_novel_transcripts.gff"

cat \
  "${REFDIR}/gencode.v49.chr_patch_hapl_scaff.annotation.gff3" \
  "${WORKDIR}/formatted_novel_transcripts.gff" \
  > "${WORKDIR}/final_extended_annotation.gff"

echo "Done."
echo "Final output:"
ls -lh "${WORKDIR}/final_extended_annotation.gff"
