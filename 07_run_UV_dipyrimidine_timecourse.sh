#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH -t 12:00:00
#SBATCH --job-name=UV_dipyrimidine
#SBATCH --output=logs/UV_dipyrimidine_%j.out
#SBATCH --error=logs/UV_dipyrimidine_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

mkdir -p logs

module purge
module load python/3.9
module load samtools
module load bedtools

# ============================================================
# SETTINGS
# ============================================================

BASE="/work/users/a/r/arslank"

REFERENCE="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

EXTRACT_RESCUE="${BASE}/uvv/07_extract_rescued_full_read_bed.py"

COUNT_DIMER="${BASE}/uvv/08_count_7th_dipyrimidine.py"

ROOT_OUT="${BASE}/uvv/XPC_UV_dipyrimidine_7th_from_3prime"

mkdir -p "$ROOT_OUT"

LEN_MIN=26
LEN_MAX=30

# This mirrors the ATL UMI-based deduplication.
HAS_UMI="no"
UMI_LEN=10

# Dimer 1 = final two bases.
# Dimer 7 = bases 8 and 7 from the 3' end.
DIMER_FROM_3PRIME=7

COMBINED_SUMMARY="${ROOT_OUT}/UV_dipyrimidine_7th_from_3prime_summary.tsv"
COMBINED_QC="${ROOT_OUT}/UV_dipyrimidine_7th_from_3prime_QC.tsv"

rm -f \
  "$COMBINED_SUMMARY" \
  "$COMBINED_QC"

first_summary="yes"
first_qc="yes"

# ============================================================
# SAMPLE TABLE
# ============================================================

while IFS=$'\t' read -r SAMPLE TIMEPOINT TIME_H; do

  [[ -z "${SAMPLE:-}" ]] && continue
  [[ "$SAMPLE" =~ ^# ]] && continue

  echo
  echo "============================================================"
  echo "Sample: ${SAMPLE}"
  echo "Time:   ${TIMEPOINT}"
  echo "============================================================"

  SAMPLE_OUT="${ROOT_OUT}/${SAMPLE}"

  mkdir -p "$SAMPLE_OUT"

  REGULAR_BED="${BASE}/${SAMPLE}_atl_output/genome/${SAMPLE}.dedup.bed"

  RESCUE_DIR="${BASE}/${SAMPLE}_CC_to_TT_rescue"

  RESCUE_BAM="${RESCUE_DIR}/${SAMPLE}.single_TT_to_CC.v0.sorted.bam"

  BEST_EVENTS="${RESCUE_DIR}/full_reads_and_NCCN/${SAMPLE}.CCTT_best_per_original_read.tsv"

  RESCUED_BED="${SAMPLE_OUT}/${SAMPLE}.accepted_CCTT_rescued_full_reads.bed"

  RESCUE_SUMMARY="${SAMPLE_OUT}/${SAMPLE}.accepted_CCTT_rescued_full_reads.summary.tsv"

  CANDIDATE_FILE="${SAMPLE_OUT}/${SAMPLE}.regular_plus_rescued.candidates.tsv"

  COMBINED_BED="${SAMPLE_OUT}/${SAMPLE}.regular_plus_rescued.dedup.bed"

  LENGTH_BED="${SAMPLE_OUT}/${SAMPLE}.regular_plus_rescued.len${LEN_MIN}_${LEN_MAX}.bed"

  SEQUENCE_TSV="${SAMPLE_OUT}/${SAMPLE}.regular_plus_rescued.len${LEN_MIN}_${LEN_MAX}.sequence.tsv"

  EVENT_TSV="${SAMPLE_OUT}/${SAMPLE}.dipyrimidine7.events.tsv"

  SUMMARY_TSV="${SAMPLE_OUT}/${SAMPLE}.dipyrimidine7.summary.tsv"

  QC_TSV="${SAMPLE_OUT}/${SAMPLE}.dipyrimidine7.QC.tsv"

  # ==========================================================
  # CHECK INPUTS
  # ==========================================================

  for file in \
    "$REGULAR_BED" \
    "$RESCUE_BAM" \
    "$BEST_EVENTS" \
    "$REFERENCE"
  do
    if [[ ! -s "$file" ]]; then
      echo "ERROR: missing input:"
      echo "$file"
      exit 1
    fi
  done

  # ==========================================================
  # 1. EXTRACT ACCEPTED RESCUED FULL-READ BED
  # ==========================================================

  echo "[1/6] Extracting accepted CC>TT rescued alignments"

  python "$EXTRACT_RESCUE" \
    --best-events "$BEST_EVENTS" \
    --bam "$RESCUE_BAM" \
    --output-bed "$RESCUED_BED" \
    --summary "$RESCUE_SUMMARY"

  N_REGULAR=$(wc -l < "$REGULAR_BED")
  N_RESCUED=$(wc -l < "$RESCUED_BED")

  echo "Regular dedup reads: ${N_REGULAR}"
  echo "Accepted rescued reads: ${N_RESCUED}"

  # ==========================================================
  # 2. COMBINE CANDIDATES
  #
  # Columns:
  # chrom start end read_id score strand source source_priority
  # ==========================================================

  echo "[2/6] Combining regular and rescued BED records"

  {
    awk '
    BEGIN {
      FS = OFS = "\t"
    }

    NF >= 6 {
      print $1, $2, $3, $4, $5, $6, "regular", 0
    }
    ' "$REGULAR_BED"

    awk '
    BEGIN {
      FS = OFS = "\t"
    }

    NF >= 6 {
      print $1, $2, $3, $4, $5, $6, "rescued", 1
    }
    ' "$RESCUED_BED"

  } > "$CANDIDATE_FILE"

# ==========================================================
# 3. DEDUPLICATE ACROSS REGULAR AND RESCUED READS
#
# UMI mode:
# chrom + start + end + strand + UMI
#
# No-UMI mode:
# chrom + start + end + strand + read length
#
# Regular reads receive priority when the same molecular
# key exists in both sources.
# ==========================================================

echo "[3/6] Deduplicating across regular and rescued reads"

awk \
  -v has_umi="$HAS_UMI" \
  -v umi_len="$UMI_LEN" '
BEGIN {
  FS = OFS = "\t"
}

{
  read_id = $4
  read_length = $3 - $2

  if (has_umi == "yes") {
    if (length(read_id) < umi_len) {
      next
    }

    dedup_component = substr(read_id, length(read_id) - umi_len + 1, umi_len)
  } else {
    dedup_component = read_length
  }

  print $1, $2, $3, $4, $5, $6, $7, $8, dedup_component
}
' "$CANDIDATE_FILE" |
LC_ALL=C sort \
  -k1,1 \
  -k2,2n \
  -k3,3n \
  -k6,6 \
  -k9,9 \
  -k8,8n |
awk '
BEGIN {
  FS = OFS = "\t"
}

{
  key = $1 FS $2 FS $3 FS $6 FS $9

  if (!(key in seen)) {
    seen[key] = 1
    output_name = $7 "|" $4

    print $1, $2, $3, output_name, $5, $6
  }
}
' > "$COMBINED_BED"

N_COMBINED=$(wc -l < "$COMBINED_BED")

echo "Combined deduplicated reads: ${N_COMBINED}"

# ==========================================================
# 4. LENGTH FILTER
# ==========================================================

echo "[4/6] Applying ${LEN_MIN}-${LEN_MAX} nt filter"

awk \
  -v minimum="$LEN_MIN" \
  -v maximum="$LEN_MAX" '
BEGIN {
  FS = OFS = "\t"
}

{
  read_length = $3 - $2

  if (read_length >= minimum && read_length <= maximum) {
    print
  }
}
' "$COMBINED_BED" > "$LENGTH_BED"

N_LENGTH=$(wc -l < "$LENGTH_BED")

echo "Length-filtered reads: ${N_LENGTH}"

  # ==========================================================
  # 5. STRAND-ORIENTED REFERENCE SEQUENCES
  # ==========================================================

  echo "[5/6] Extracting reference sequences"

  bedtools getfasta \
    -fi "$REFERENCE" \
    -bed "$LENGTH_BED" \
    -s \
    -name \
    -tab \
    > "$SEQUENCE_TSV"

  N_SEQUENCES=$(wc -l < "$SEQUENCE_TSV")

  echo "Sequences extracted: ${N_SEQUENCES}"

  if [[ "$N_LENGTH" -ne "$N_SEQUENCES" ]]; then
    echo "ERROR: BED and sequence row counts differ."
    exit 1
  fi

  # ==========================================================
  # 6. COUNT THE 7TH DINUCLEOTIDE
  # ==========================================================

  echo "[6/6] Counting CC, CT, TC and TT"

  python "$COUNT_DIMER" \
    --input "$SEQUENCE_TSV" \
    --sample "$SAMPLE" \
    --timepoint "$TIMEPOINT" \
    --time-h "$TIME_H" \
    --dimer-from-3prime "$DIMER_FROM_3PRIME" \
    --events "$EVENT_TSV" \
    --summary "$SUMMARY_TSV" \
    --qc "$QC_TSV"

  # ==========================================================
  # COMBINE SAMPLE TABLES
  # ==========================================================

  if [[ "$first_summary" == "yes" ]]; then
    cat "$SUMMARY_TSV" \
      > "$COMBINED_SUMMARY"

    first_summary="no"
  else
    tail -n +2 "$SUMMARY_TSV" \
      >> "$COMBINED_SUMMARY"
  fi

  if [[ "$first_qc" == "yes" ]]; then
    cat "$QC_TSV" \
      > "$COMBINED_QC"

    first_qc="no"
  else
    tail -n +2 "$QC_TSV" \
      >> "$COMBINED_QC"
  fi

done <<'SAMPLES'
XPC-UVCPD-30m	0.5h	0.5
XPC-UVCPD-2h	2h	2
XPC-UVCPD-4h	4h	4
XPC-UVCPD-8h	8h	8
SAMPLES

echo
echo "============================================================"
echo "Done"
echo "============================================================"
echo
echo "Combined summary:"
echo "$COMBINED_SUMMARY"
echo
echo "Combined QC:"
echo "$COMBINED_QC"