#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=128g
#SBATCH -t 24:00:00
#SBATCH --job-name=CC_TT
#SBATCH --array=0-14%2
#SBATCH --output=logs/mismatch_CCtoTT_%A_%a.out
#SBATCH --error=logs/mismatch_CCtoTT_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

# ============================================================
# MODULES
# ============================================================

module purge
module load python/3.9
module load samtools
module load bowtie2
module load bowtie

# ============================================================
# SETTINGS
# ============================================================

THREADS="${SLURM_CPUS_PER_TASK:-4}"
MAPQ=10

# samtools sort memory is per thread
SORT_MEM="2G"

BOWTIE2_IND="/proj/seq/data/hg38_UCSC/Sequence/Bowtie2Index/genome"
BOWTIE_INDEX="/proj/seq/data/hg38_UCSC/Sequence/BowtieIndex/genome"
REF_FASTA="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

PYTHON_SCRIPT="./cc_tt_rescue_tools.py"

SAMPLES=(
  "XPC-UVCPD-30m"
  "XPC-UVCPD-2h"
  "XPC-UVCPD-4h"
  "XPC-UVCPD-8h"
  "XPC-UVCPD-30m-r2"
  "XPC-UVCPD-2h-r2"
  "XPC-UVCPD-4h-r2"
  "XPC-UVCPD-8h-r2"
)

# ============================================================
# SELECT ARRAY SAMPLE
# ============================================================

TASK_ID="${SLURM_ARRAY_TASK_ID:?SLURM_ARRAY_TASK_ID is not defined}"

if (( TASK_ID < 0 || TASK_ID >= ${#SAMPLES[@]} )); then
  echo "ERROR: Invalid array index: ${TASK_ID}" >&2
  exit 1
fi

SAMPLE="${SAMPLES[$TASK_ID]}"

# ============================================================
# PATHS
# ============================================================

OUT="${SAMPLE}_atl_output"
CLEAN_FASTQ="${OUT}/${SAMPLE}.clean.fastq"

OUTDIR="${SAMPLE}_CC_to_TT_rescue"
QCDIR="${OUTDIR}/qc"
GENOMEDIR="${OUTDIR}/genome"

mkdir -p "${OUTDIR}" "${QCDIR}" "${GENOMEDIR}"

INITIAL_UNMAPPED="${OUTDIR}/${SAMPLE}.initial.unmapped.fastq"
INITIAL_BAM="${GENOMEDIR}/${SAMPLE}.initial.sorted.bam"

VARIANT_FASTQ="${OUTDIR}/${SAMPLE}.single_TT_to_CC_variants.fastq"
VARIANT_MAP="${OUTDIR}/${SAMPLE}.single_TT_to_CC_variant_map.tsv"

RAW_SAM="${OUTDIR}/${SAMPLE}.single_TT_to_CC.v0.sam"
CLEAN_SAM="${OUTDIR}/${SAMPLE}.single_TT_to_CC.v0.clean.sam"
BAD_SAM="${QCDIR}/${SAMPLE}.single_TT_to_CC.v0.bad_sam_lines.tsv"

RESCUE_BAM="${OUTDIR}/${SAMPLE}.single_TT_to_CC.v0.sorted.bam"

EVENT_TSV="${OUTDIR}/${SAMPLE}.CC_to_TT_tandem_events.tsv"
EVENT_BED="${OUTDIR}/${SAMPLE}.CC_to_TT_tandem_events.bed"
SORTED_EVENT_BED="${OUTDIR}/${SAMPLE}.CC_to_TT_tandem_events.sorted.bed"

STEP5_LOG="${QCDIR}/${SAMPLE}.call_events.time_and_error.log"

# ============================================================
# FUNCTIONS
# ============================================================

bam_is_valid() {
  local bam="$1"

  [[ -s "${bam}" ]] || return 1
  samtools quickcheck -q "${bam}" || return 1

  return 0
}

remove_sort_temp() {
  rm -f \
    "${GENOMEDIR}/${SAMPLE}.initial.sorttmp."*.bam \
    "${OUTDIR}/${SAMPLE}.rescue.sorttmp."*.bam \
    2>/dev/null || true
}

print_file_size() {
  local file="$1"

  if [[ -e "${file}" ]]; then
    ls -lh "${file}"
  else
    echo "Missing: ${file}"
  fi
}

trap remove_sort_temp EXIT

# ============================================================
# INITIAL CHECKS
# ============================================================

echo "============================================================"
echo "Sample:          ${SAMPLE}"
echo "Array index:     ${TASK_ID}"
echo "SLURM job:       ${SLURM_JOB_ID:-NA}"
echo "Threads:         ${THREADS}"
echo "Requested RAM:   ${SLURM_MEM_PER_NODE:-unknown} MB"
echo "Started:         $(date)"
echo "============================================================"

if [[ ! -s "${PYTHON_SCRIPT}" ]]; then
  echo "ERROR: Missing Python script: ${PYTHON_SCRIPT}" >&2
  exit 1
fi

if [[ ! -s "${REF_FASTA}" ]]; then
  echo "ERROR: Missing reference FASTA: ${REF_FASTA}" >&2
  exit 1
fi

if [[ ! -s "${REF_FASTA}.fai" ]]; then
  echo "Reference FASTA index not found. Creating it."
  samtools faidx "${REF_FASTA}"
fi

# ============================================================
# STEP 1
# Initial Bowtie2 alignment
# ============================================================

echo
echo "Step 1: Initial Bowtie2 alignment and unmapped-read recovery"

if bam_is_valid "${INITIAL_BAM}" && \
   [[ -s "${INITIAL_BAM}.bai" ]] && \
   [[ -s "${INITIAL_UNMAPPED}" ]]; then

  echo "Step 1 already completed."
  echo "Reusing:"
  print_file_size "${INITIAL_BAM}"
  print_file_size "${INITIAL_UNMAPPED}"

else
  if [[ ! -s "${CLEAN_FASTQ}" ]]; then
    echo "ERROR: Missing or empty clean FASTQ:" >&2
    echo "       ${CLEAN_FASTQ}" >&2
    exit 1
  fi

  rm -f \
    "${INITIAL_BAM}" \
    "${INITIAL_BAM}.bai" \
    "${INITIAL_BAM}.csi"

  remove_sort_temp

  bowtie2 \
    -x "${BOWTIE2_IND}" \
    -U "${CLEAN_FASTQ}" \
    --very-sensitive \
    -p "${THREADS}" \
    --un "${INITIAL_UNMAPPED}" \
    2> "${QCDIR}/${SAMPLE}.initial.bt2.log" \
  | samtools view \
      -@ "${THREADS}" \
      -b \
      -F 0x904 \
      -q "${MAPQ}" \
      - \
  | samtools sort \
      -@ "${THREADS}" \
      -m "${SORT_MEM}" \
      -T "${GENOMEDIR}/${SAMPLE}.initial.sorttmp" \
      -o "${INITIAL_BAM}" \
      -

  samtools quickcheck -v "${INITIAL_BAM}"
  samtools index -@ "${THREADS}" "${INITIAL_BAM}"

  if [[ ! -s "${INITIAL_UNMAPPED}" ]]; then
    echo "WARNING: No unmapped reads were produced."
    exit 0
  fi
fi

echo "Initial alignment count retained after filtering:"
samtools view -@ "${THREADS}" -c "${INITIAL_BAM}"

echo "Initial unmapped-read summary:"
grep "aligned 0 times" \
  "${QCDIR}/${SAMPLE}.initial.bt2.log" 2>/dev/null || true

# ============================================================
# STEP 2
# Generate TT-to-CC variants
# ============================================================

echo
echo "Step 2: Generate one TT-to-CC variant per TT occurrence"

if [[ -s "${VARIANT_FASTQ}" && -s "${VARIANT_MAP}" ]]; then
  echo "Step 2 already completed."
  echo "Reusing:"
  print_file_size "${VARIANT_FASTQ}"
  print_file_size "${VARIANT_MAP}"

else
  rm -f "${VARIANT_FASTQ}" "${VARIANT_MAP}"

  python3 "${PYTHON_SCRIPT}" make-variants \
    --in-fastq "${INITIAL_UNMAPPED}" \
    --out-fastq "${VARIANT_FASTQ}" \
    --map-tsv "${VARIANT_MAP}"

  if [[ ! -s "${VARIANT_FASTQ}" ]]; then
    echo "WARNING: No TT-to-CC variants were generated."
    exit 0
  fi

  if [[ ! -s "${VARIANT_MAP}" ]]; then
    echo "ERROR: Variant map was not generated." >&2
    exit 1
  fi
fi

echo "Variant FASTQ reads:"
awk 'END {print int(NR / 4)}' "${VARIANT_FASTQ}"

echo "Variant-map lines:"
wc -l "${VARIANT_MAP}"

# ============================================================
# STEPS 3–4
# Strict Bowtie rescue and BAM conversion
# ============================================================

echo
echo "Steps 3-4: Strict Bowtie rescue and sorted BAM generation"

if bam_is_valid "${RESCUE_BAM}" && [[ -s "${RESCUE_BAM}.bai" ]]; then
  echo "Rescue BAM already completed."
  echo "Reusing:"
  print_file_size "${RESCUE_BAM}"

else
  # ----------------------------------------------------------
  # Step 3
  # ----------------------------------------------------------

  if [[ -s "${RAW_SAM}" ]]; then
    echo "Existing rescue SAM found. Reusing it."
    print_file_size "${RAW_SAM}"
  else
    echo "Step 3: Strict Bowtie alignment with zero mismatches"

    bowtie \
      -q \
      -v 0 \
      -m 1 \
      --best \
      --strata \
      --sam \
      -p "${THREADS}" \
      "${BOWTIE_INDEX}" \
      "${VARIANT_FASTQ}" \
      > "${RAW_SAM}" \
      2> "${QCDIR}/${SAMPLE}.single_TT_to_CC.v0.bowtie.log"
  fi

  if [[ ! -s "${RAW_SAM}" ]]; then
    echo "ERROR: Rescue SAM is missing or empty." >&2
    exit 1
  fi

  # ----------------------------------------------------------
  # Step 4
  # ----------------------------------------------------------

  echo "Step 4: Validate SAM records and create sorted BAM"

  rm -f \
    "${CLEAN_SAM}" \
    "${BAD_SAM}" \
    "${RESCUE_BAM}" \
    "${RESCUE_BAM}.bai" \
    "${RESCUE_BAM}.csi"

  awk \
    -v clean="${CLEAN_SAM}" \
    -v bad="${BAD_SAM}" '
    BEGIN {
      FS = OFS = "\t"
    }

    /^@/ {
      print > clean
      next
    }

    NF < 11 {
      print NR, "NF_lt_11", NF, $0 > bad
      next
    }

    $10 != "*" && $11 != "*" && length($10) != length($11) {
      print NR,
            "SEQ_QUAL_length_mismatch",
            length($10),
            length($11),
            $1 > bad
      next
    }

    {
      print > clean
    }
  ' "${RAW_SAM}"

  echo "Malformed SAM lines removed:"

  if [[ -s "${BAD_SAM}" ]]; then
    wc -l "${BAD_SAM}"
    echo "First malformed records:"
    head -n 10 "${BAD_SAM}" || true
  else
    echo "0"
  fi

  remove_sort_temp

  samtools view \
    -@ "${THREADS}" \
    -b \
    "${CLEAN_SAM}" \
  | samtools sort \
      -@ "${THREADS}" \
      -m "${SORT_MEM}" \
      -T "${OUTDIR}/${SAMPLE}.rescue.sorttmp" \
      -o "${RESCUE_BAM}" \
      -

  samtools quickcheck -v "${RESCUE_BAM}"
  samtools index -@ "${THREADS}" "${RESCUE_BAM}"
fi

echo "Rescued alignment count:"
samtools view -@ "${THREADS}" -c "${RESCUE_BAM}"

# ============================================================
# REMOVE LARGE INTERMEDIATE SAM FILES
# ============================================================

echo
echo "Removing large SAM intermediates before event calling"

if bam_is_valid "${RESCUE_BAM}"; then
  rm -f "${RAW_SAM}" "${CLEAN_SAM}"
else
  echo "ERROR: Rescue BAM validation failed; SAM files will not be removed." >&2
  exit 1
fi

remove_sort_temp

# ============================================================
# STEP 5
# Call tandem events
# ============================================================

echo
echo "Step 5: Call rescued CC>TT tandem events"
echo "Step 5 started: $(date)"

if [[ -s "${EVENT_TSV}" && -s "${EVENT_BED}" ]]; then
  echo "Step 5 already completed."
  echo "Reusing:"
  print_file_size "${EVENT_TSV}"
  print_file_size "${EVENT_BED}"

else
  # Remove potentially truncated files left by a previous OOM kill.
  rm -f "${EVENT_TSV}" "${EVENT_BED}" "${STEP5_LOG}"

  set +e

  /usr/bin/time -v \
    python3 "${PYTHON_SCRIPT}" call-events \
      --variant-map "${VARIANT_MAP}" \
      --bam "${RESCUE_BAM}" \
      --ref-fasta "${REF_FASTA}" \
      --out-tsv "${EVENT_TSV}" \
      --out-bed "${EVENT_BED}" \
      2> "${STEP5_LOG}"

  STEP5_STATUS=$?

  set -e

  echo
  echo "Step 5 resource summary:"

  grep -E \
    "Maximum resident set size|Elapsed|User time|System time|Exit status|Command terminated" \
    "${STEP5_LOG}" || true

  if [[ "${STEP5_STATUS}" -ne 0 ]]; then
    echo "ERROR: Step 5 failed with exit status ${STEP5_STATUS}." >&2
    echo "See: ${STEP5_LOG}" >&2

    rm -f "${EVENT_TSV}" "${EVENT_BED}"
    exit "${STEP5_STATUS}"
  fi

  if [[ ! -s "${EVENT_TSV}" ]]; then
    echo "ERROR: Step 5 completed but event TSV is missing or empty." >&2
    exit 1
  fi

  if [[ ! -e "${EVENT_BED}" ]]; then
    echo "ERROR: Step 5 completed but event BED was not created." >&2
    exit 1
  fi
fi

echo "Step 5 completed: $(date)"

# ============================================================
# STEP 6
# Sort BED
# ============================================================

echo
echo "Step 6: Sort event BED"

if [[ -s "${SORTED_EVENT_BED}" ]]; then
  echo "Sorted BED already exists."
else
  rm -f "${SORTED_EVENT_BED}"

  if [[ -s "${EVENT_BED}" ]]; then
    LC_ALL=C sort \
      -S 8G \
      -T "${OUTDIR}" \
      -k1,1 \
      -k2,2n \
      "${EVENT_BED}" \
      > "${SORTED_EVENT_BED}"
  else
    echo "WARNING: Event BED is empty."
    touch "${SORTED_EVENT_BED}"
  fi
fi

# ============================================================
# FINAL SUMMARY
# ============================================================

echo
echo "============================================================"
echo "Completed sample: ${SAMPLE}"
echo "Finished:         $(date)"
echo "============================================================"

echo "Main outputs:"
print_file_size "${EVENT_TSV}"
print_file_size "${SORTED_EVENT_BED}"

echo
echo "Event TSV lines:"
wc -l "${EVENT_TSV}"

echo
echo "Sorted BED lines:"
wc -l "${SORTED_EVENT_BED}"

echo
echo "Step 5 memory record:"
echo "${STEP5_LOG}"

echo "============================================================"