#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32g
#SBATCH -t 18:00:00
#SBATCH --job-name=UV_7th_dinuc
#SBATCH --output=logs/UV_7th_dinuc_%j.out
#SBATCH --error=logs/UV_7th_dinuc_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

mkdir -p logs

module purge
module load python/3.9
module load bedtools

# ============================================================
# SETTINGS
# ============================================================

BASE="/work/users/a/r/arslank"

REF_FASTA="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

PYTHON_SCRIPT="./26_build_combined_dedup_7th_dinuc.py"

OUTDIR="${BASE}/XPC_UV_TFBS_7th_dinucleotide_repair"

mkdir -p "$OUTDIR"

SAMPLES=(
  "XPC-UVCPD-30m|0.5h|0.5"
  "XPC-UVCPD-2h|2h|2"
  "XPC-UVCPD-4h|4h|4"
  "XPC-UVCPD-8h|8h|8"  
  "XPC-UVCPD-30m-r2|0.5h|0.5"
  "XPC-UVCPD-2h-r2|2h|2"
  "XPC-UVCPD-4h-r2|4h|4"
  "XPC-UVCPD-8h-r2|8h|8"
)

SOURCES=(
  "regular"
  "CC>TT_rescue"
)

DINUCLEOTIDES=(
  "CC"
  "CT"
  "TC"
  "TT"
)

# ============================================================
# CHECK STATIC INPUTS
# ============================================================

for path in \
  "$REF_FASTA" \
  "$PYTHON_SCRIPT"
do
  if [[ ! -s "$path" ]]; then
    echo "ERROR: missing or empty file:"
    echo "$path"
    exit 1
  fi
done

# ============================================================
# ARRAYS FOR COMBINED OUTPUTS
# ============================================================

ALL_DIPYRIMIDINE_FILES=()
ALL_COMBINED_DEDUP_FILES=()

# ============================================================
# PROCESS EACH SAMPLE
# ============================================================

for entry in "${SAMPLES[@]}"; do

  IFS='|' read -r SAMPLE TIMEPOINT TIME_H <<< "$entry"

  echo
  echo "============================================================"
  echo "Processing: $SAMPLE"
  echo "Time point: $TIMEPOINT"
  echo "============================================================"
  echo

  REGULAR_BED="${BASE}/${SAMPLE}_atl_output/genome/${SAMPLE}.dedup.bed"

  RESCUED_TABLE="${BASE}/${SAMPLE}_CC_to_TT_rescue/full_reads_and_NCCN/${SAMPLE}.CCTT_best_per_original_read.tsv"

  SAMPLE_OUTDIR="${OUTDIR}/${SAMPLE}"

  mkdir -p "$SAMPLE_OUTDIR"

  # ----------------------------------------------------------
  # Main output files
  # ----------------------------------------------------------

  COMBINED_BED="${SAMPLE_OUTDIR}/${SAMPLE}.regular_plus_CCTT.final_dedup.bed"

  DINUC_BED="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.bed"

  DINUC_SEQUENCE="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.sequence.tsv"

  DIPYRIMIDINE_BED="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.dipyrimidines.bed"

  COMBINE_QC="${SAMPLE_OUTDIR}/${SAMPLE}.combined_dedup.QC.tsv"

  SEQUENCE_QC="${SAMPLE_OUTDIR}/${SAMPLE}.7th_dinucleotide.QC.tsv"

  SAMPLE_SUMMARY="${SAMPLE_OUTDIR}/${SAMPLE}.processing_summary.tsv"

  # ----------------------------------------------------------
  # Temporary sequence-extraction files
  #
  # DINUC_BED has 12 custom columns. BEDTools would interpret
  # that directly as formal BED12, so a BED6 copy is used.
  # ----------------------------------------------------------

  DINUC_GETFASTA_BED6="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.getfasta.bed6"

  DINUC_GETFASTA_OUTPUT="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.getfasta.bed6.sequence.tsv"

  DINUC_SEQUENCE_ONLY="${SAMPLE_OUTDIR}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.sequence_only.txt"

  # ==========================================================
  # CHECK SAMPLE INPUTS
  # ==========================================================

  for path in \
    "$REGULAR_BED" \
    "$RESCUED_TABLE"
  do
    if [[ ! -s "$path" ]]; then
      echo "ERROR: missing or empty sample input:"
      echo "$path"
      exit 1
    fi
  done

  # ==========================================================
  # 1. ADD RESCUED CC>TT READS AND DEDUPLICATE AGAIN
  # ==========================================================

  echo "Combining regular and rescued reads..."

  python "$PYTHON_SCRIPT" combine \
    --sample "$SAMPLE" \
    --timepoint "$TIMEPOINT" \
    --time-h "$TIME_H" \
    --regular-bed "$REGULAR_BED" \
    --rescued-table "$RESCUED_TABLE" \
    --combined-bed "$COMBINED_BED" \
    --dinucleotide-bed "$DINUC_BED" \
    --qc "$COMBINE_QC"

  if [[ ! -s "$COMBINED_BED" ]]; then
    echo "ERROR: combined deduplicated BED was not created:"
    echo "$COMBINED_BED"
    exit 1
  fi

  if [[ ! -s "$DINUC_BED" ]]; then
    echo "ERROR: seventh-dinucleotide BED was not created:"
    echo "$DINUC_BED"
    exit 1
  fi

  # ==========================================================
  # 2. VALIDATE THE 12-COLUMN DINUCLEOTIDE BED
  # ==========================================================

  echo "Validating seventh-dinucleotide intervals..."

  awk '
  BEGIN {
    FS = OFS = "\t"
  }

  NF != 12 {
    print "ERROR: expected 12 columns at line " NR \
          ", observed " NF > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  $1 == "" {
    print "ERROR: empty chromosome at line " NR > "/dev/stderr"
    exit 1
  }

  $2 !~ /^[0-9]+$/ || $3 !~ /^[0-9]+$/ {
    print "ERROR: invalid coordinates at line " NR > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  ($3 - $2) != 2 {
    print "ERROR: dinucleotide interval is not 2 bp at line " NR \
          > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  $6 != "+" && $6 != "-" {
    print "ERROR: invalid strand at line " NR > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  $11 !~ /^[0-9]+$/ {
    print "ERROR: invalid read length at line " NR > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  $11 < 26 || $11 > 30 {
    print "ERROR: read length outside 26-30 nt at line " NR \
          > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }
  ' "$DINUC_BED"

  # ==========================================================
  # 3. CREATE BED6 FOR BEDTOOLS GETFASTA
  # ==========================================================

  echo "Creating temporary BED6..."

  awk '
  BEGIN {
    FS = OFS = "\t"
  }

  {
    print $1, $2, $3, $4, $5, $6
  }
  ' "$DINUC_BED" > "$DINUC_GETFASTA_BED6"

  if [[ ! -s "$DINUC_GETFASTA_BED6" ]]; then
    echo "ERROR: BED6 extraction file was not created:"
    echo "$DINUC_GETFASTA_BED6"
    exit 1
  fi

  # ==========================================================
  # 4. EXTRACT STRAND-ORIENTED DINUCLEOTIDE SEQUENCES
  # ==========================================================

  echo "Extracting strand-oriented reference dinucleotides..."

  bedtools getfasta \
    -fi "$REF_FASTA" \
    -bed "$DINUC_GETFASTA_BED6" \
    -s \
    -bedOut \
    > "$DINUC_GETFASTA_OUTPUT"

  if [[ ! -s "$DINUC_GETFASTA_OUTPUT" ]]; then
    echo "ERROR: bedtools getfasta produced no output:"
    echo "$DINUC_GETFASTA_OUTPUT"
    exit 1
  fi

  # BED6 plus the appended sequence gives seven columns.

  awk '
  BEGIN {
    FS = OFS = "\t"
  }

  NF != 7 {
    print "ERROR: expected 7 columns from getfasta at line " NR \
          ", observed " NF > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }

  {
    sequence = toupper($7)

    if (sequence !~ /^[ACGTN][ACGTN]$/) {
      print "ERROR: invalid extracted dinucleotide at line " NR \
            ": " sequence > "/dev/stderr"
      exit 1
    }

    print sequence
  }
  ' "$DINUC_GETFASTA_OUTPUT" > "$DINUC_SEQUENCE_ONLY"

  # ==========================================================
  # 5. VERIFY ROW COUNTS
  # ==========================================================

  N_METADATA=$(wc -l < "$DINUC_BED")
  N_BED6=$(wc -l < "$DINUC_GETFASTA_BED6")
  N_GETFASTA=$(wc -l < "$DINUC_GETFASTA_OUTPUT")
  N_SEQUENCE=$(wc -l < "$DINUC_SEQUENCE_ONLY")

  echo "Metadata rows: $N_METADATA"
  echo "BED6 rows:     $N_BED6"
  echo "Getfasta rows: $N_GETFASTA"
  echo "Sequence rows: $N_SEQUENCE"

  if [[ "$N_METADATA" -ne "$N_BED6" ]]; then
    echo "ERROR: metadata and BED6 row counts differ."
    exit 1
  fi

  if [[ "$N_METADATA" -ne "$N_GETFASTA" ]]; then
    echo "ERROR: metadata and getfasta row counts differ."
    exit 1
  fi

  if [[ "$N_METADATA" -ne "$N_SEQUENCE" ]]; then
    echo "ERROR: metadata and extracted sequence row counts differ."
    exit 1
  fi

  # ==========================================================
  # 6. APPEND SEQUENCE TO THE ORIGINAL 12-COLUMN TABLE
  # ==========================================================

  echo "Appending extracted sequence to metadata..."

  paste \
    "$DINUC_BED" \
    "$DINUC_SEQUENCE_ONLY" \
    > "$DINUC_SEQUENCE"

  if [[ ! -s "$DINUC_SEQUENCE" ]]; then
    echo "ERROR: sequence-annotated table was not created:"
    echo "$DINUC_SEQUENCE"
    exit 1
  fi

  awk '
  BEGIN {
    FS = OFS = "\t"
  }

  NF != 13 {
    print "ERROR: expected 13 columns at line " NR \
          ", observed " NF > "/dev/stderr"
    print $0 > "/dev/stderr"
    exit 1
  }
  ' "$DINUC_SEQUENCE"

  # ==========================================================
  # 7. RETAIN CC, CT, TC, AND TT
  # ==========================================================

  echo "Selecting dipyrimidines..."

  python "$PYTHON_SCRIPT" filter-sequences \
    --input "$DINUC_SEQUENCE" \
    --output "$DIPYRIMIDINE_BED" \
    --qc "$SEQUENCE_QC"

  if [[ ! -s "$DIPYRIMIDINE_BED" ]]; then
    echo "ERROR: no dipyrimidine output was created:"
    echo "$DIPYRIMIDINE_BED"
    exit 1
  fi

  # ==========================================================
  # 8. SAMPLE-LEVEL QC SUMMARY
  # ==========================================================

  TOTAL_FINAL_READS=$(wc -l < "$COMBINED_BED")
  TOTAL_26TO30_READS=$(wc -l < "$DINUC_BED")
  TOTAL_DIPYRIMIDINES=$(wc -l < "$DIPYRIMIDINE_BED")

  REGULAR_DIPYRIMIDINES=$(
    awk -F $'\t' '
    $10 == "regular" {
      n++
    }

    END {
      print n + 0
    }
    ' "$DIPYRIMIDINE_BED"
  )

  RESCUED_DIPYRIMIDINES=$(
    awk -F $'\t' '
    $10 == "CC>TT_rescue" {
      n++
    }

    END {
      print n + 0
    }
    ' "$DIPYRIMIDINE_BED"
  )

  {
    echo -e "metric\tvalue"
    echo -e "sample\t${SAMPLE}"
    echo -e "timepoint\t${TIMEPOINT}"
    echo -e "final_combined_deduplicated_reads\t${TOTAL_FINAL_READS}"
    echo -e "final_26to30nt_reads\t${TOTAL_26TO30_READS}"
    echo -e "selected_dipyrimidine_events\t${TOTAL_DIPYRIMIDINES}"
    echo -e "selected_regular_events\t${REGULAR_DIPYRIMIDINES}"
    echo -e "selected_CCTT_rescue_events\t${RESCUED_DIPYRIMIDINES}"
  } > "$SAMPLE_SUMMARY"

  # ==========================================================
  # 9. COLLECT SAMPLE OUTPUTS
  # ==========================================================

  ALL_DIPYRIMIDINE_FILES+=(
    "$DIPYRIMIDINE_BED"
  )

  ALL_COMBINED_DEDUP_FILES+=(
    "$COMBINED_BED"
  )

  # Temporary sequence-only file is no longer needed.

  rm -f "$DINUC_SEQUENCE_ONLY"

  echo
  echo "Completed: $SAMPLE"
  echo "Final combined reads: $TOTAL_FINAL_READS"
  echo "26-30 nt reads: $TOTAL_26TO30_READS"
  echo "Dipyrimidines: $TOTAL_DIPYRIMIDINES"
  echo

done

# ============================================================
# COMBINE ALL TIME POINTS
# ============================================================

ALL_DIPYRIMIDINES="${OUTDIR}/UV_all_timepoints.26to30nt.7th_dinucleotide_from_3prime.dipyrimidines.bed"

ALL_COMBINED_DEDUP="${OUTDIR}/UV_all_timepoints.regular_plus_CCTT.final_dedup.bed"

echo
echo "Combining all time points..."
echo

cat "${ALL_DIPYRIMIDINE_FILES[@]}" \
| LC_ALL=C sort \
  -k1,1 \
  -k2,2n \
  -k3,3n \
  -k7,7 \
  -k8,8 \
> "$ALL_DIPYRIMIDINES"

cat "${ALL_COMBINED_DEDUP_FILES[@]}" \
| LC_ALL=C sort \
  -k1,1 \
  -k2,2n \
  -k3,3n \
  -k7,7 \
  -k8,8 \
> "$ALL_COMBINED_DEDUP"

if [[ ! -s "$ALL_DIPYRIMIDINES" ]]; then
  echo "ERROR: combined dipyrimidine file is empty."
  exit 1
fi

if [[ ! -s "$ALL_COMBINED_DEDUP" ]]; then
  echo "ERROR: combined deduplicated-read file is empty."
  exit 1
fi

# ============================================================
# FINAL DINUCLEOTIDE SUMMARY
# ============================================================

SUMMARY="${OUTDIR}/UV_7th_dinucleotide_summary.tsv"

{
  echo -e "sample\ttimepoint\ttime_h\tdinucleotide\tcount\ttotal_dipyrimidines\tpercent_of_dipyrimidines"

  for entry in "${SAMPLES[@]}"; do

    IFS='|' read -r SAMPLE TIMEPOINT TIME_H <<< "$entry"

    FILE="${OUTDIR}/${SAMPLE}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.dipyrimidines.bed"

    TOTAL=$(wc -l < "$FILE")

    for DINUC in "${DINUCLEOTIDES[@]}"; do

      COUNT=$(
        awk \
          -F $'\t' \
          -v d="$DINUC" '
          $13 == d {
            n++
          }

          END {
            print n + 0
          }
          ' "$FILE"
      )

      PERCENT=$(
        awk \
          -v count="$COUNT" \
          -v total="$TOTAL" '
          BEGIN {
            if (total > 0) {
              printf "%.6f", 100 * count / total
            } else {
              printf "0"
            }
          }
          '
      )

      echo -e "${SAMPLE}\t${TIMEPOINT}\t${TIME_H}\t${DINUC}\t${COUNT}\t${TOTAL}\t${PERCENT}"

    done

  done

} > "$SUMMARY"

# ============================================================
# FINAL SOURCE SUMMARY
# ============================================================

SOURCE_SUMMARY="${OUTDIR}/UV_7th_dinucleotide_source_summary.tsv"

{
  echo -e "sample\ttimepoint\tsource\tdinucleotide\tcount"

  for entry in "${SAMPLES[@]}"; do

    IFS='|' read -r SAMPLE TIMEPOINT TIME_H <<< "$entry"

    FILE="${OUTDIR}/${SAMPLE}/${SAMPLE}.26to30nt.7th_dinucleotide_from_3prime.dipyrimidines.bed"

    for SOURCE in "${SOURCES[@]}"; do

      for DINUC in "${DINUCLEOTIDES[@]}"; do

        COUNT=$(
          awk \
            -F $'\t' \
            -v source="$SOURCE" \
            -v d="$DINUC" '
            $10 == source && $13 == d {
              n++
            }

            END {
              print n + 0
            }
            ' "$FILE"
        )

        echo -e "${SAMPLE}\t${TIMEPOINT}\t${SOURCE}\t${DINUC}\t${COUNT}"

      done

    done

  done

} > "$SOURCE_SUMMARY"

# ============================================================
# FINAL QC CHECK
#
# Sum the rows in all sample files and compare them with the
# combined all-time-point file.
# ============================================================

EXPECTED_ALL_DIPYRIMIDINES=0

for file in "${ALL_DIPYRIMIDINE_FILES[@]}"; do

  if [[ ! -f "$file" ]]; then
    echo "ERROR: missing dipyrimidine file:"
    echo "$file"
    exit 1
  fi

  N=$(wc -l < "$file")

  EXPECTED_ALL_DIPYRIMIDINES=$((EXPECTED_ALL_DIPYRIMIDINES + N))

done

OBSERVED_ALL_DIPYRIMIDINES=$(wc -l < "$ALL_DIPYRIMIDINES")

echo
echo "Expected combined dipyrimidine rows: $EXPECTED_ALL_DIPYRIMIDINES"
echo "Observed combined dipyrimidine rows: $OBSERVED_ALL_DIPYRIMIDINES"
echo

if [[ "$EXPECTED_ALL_DIPYRIMIDINES" -ne "$OBSERVED_ALL_DIPYRIMIDINES" ]]; then
  echo "ERROR: combined dipyrimidine row count is incorrect."
  echo "Expected: $EXPECTED_ALL_DIPYRIMIDINES"
  echo "Observed: $OBSERVED_ALL_DIPYRIMIDINES"
  exit 1
fi

# ============================================================
# FINAL COMBINED QC SUMMARY
# ============================================================

FINAL_QC="${OUTDIR}/UV_7th_dinucleotide_final_QC.tsv"

TOTAL_COMBINED_READS=$(wc -l < "$ALL_COMBINED_DEDUP")
TOTAL_COMBINED_DIPYRIMIDINES=$(wc -l < "$ALL_DIPYRIMIDINES")

{
  echo -e "metric\tvalue"
  echo -e "number_of_samples\t${#SAMPLES[@]}"
  echo -e "combined_final_deduplicated_reads\t${TOTAL_COMBINED_READS}"
  echo -e "expected_combined_dipyrimidines\t${EXPECTED_ALL_DIPYRIMIDINES}"
  echo -e "observed_combined_dipyrimidines\t${OBSERVED_ALL_DIPYRIMIDINES}"
  echo -e "combined_dipyrimidines\t${TOTAL_COMBINED_DIPYRIMIDINES}"
} > "$FINAL_QC"

# ============================================================
# FINISHED
# ============================================================

echo
echo "============================================================"
echo "Done"
echo "============================================================"
echo

echo "Final combined deduplicated reads:"
echo "$ALL_COMBINED_DEDUP"
echo

echo "Final seventh-dinucleotide dipyrimidine events:"
echo "$ALL_DIPYRIMIDINES"
echo

echo "Dinucleotide summary:"
echo "$SUMMARY"
echo

echo "Source summary:"
echo "$SOURCE_SUMMARY"
echo

echo "Final QC:"
echo "$FINAL_QC"
echo