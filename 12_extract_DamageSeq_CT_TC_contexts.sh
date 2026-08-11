#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=8g
#SBATCH -t 04:00:00
#SBATCH --job-name=damage_CT_TC_ctx
#SBATCH --output=logs/damage_CT_TC_ctx_%j.out
#SBATCH --error=logs/damage_CT_TC_ctx_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

mkdir -p logs

# ============================================================
# SETTINGS
# ============================================================

BASE="/work/users/a/r/arslank"

INPUT="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_CC_NCCN_from_read_starts/NHF1_CPD_0h.minus3_to_0.raw_sequence.tsv"

OUTDIR="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_CT_TC_contexts"

EVENTS="${OUTDIR}/DamageSeq_NCTN_NTCN_events.tsv"
NCTN_COUNTS="${OUTDIR}/DamageSeq_NCTN_counts.tsv"
NTCN_COUNTS="${OUTDIR}/DamageSeq_NTCN_counts.tsv"
TRINUC_COUNTS="${OUTDIR}/DamageSeq_NCT_TCN_counts.tsv"
QC="${OUTDIR}/DamageSeq_NCTN_NTCN_QC.tsv"

mkdir -p "$OUTDIR"

# ============================================================
# CHECK INPUT
# ============================================================

if [[ ! -s "$INPUT" ]]; then
  echo "ERROR: missing input:"
  echo "$INPUT"
  exit 1
fi

# ============================================================
# EXTRACT
#
# Input sequence order:
#
#   sequence base:       1   2   3   4
#   relative position:  -3  -2  -1   0
#
# Central dipyrimidine:
#
#   substr(sequence, 2, 2)
#
# CT:
#   full context = NCTN = bases 1-4
#   trinucleotide = NCT = bases 1-3
#
# TC:
#   full context = NTCN = bases 1-4
#   trinucleotide = TCN = bases 2-4
# ============================================================

echo "Extracting Damage-seq CT and TC sequence contexts"

awk \
  -v events="$EVENTS" \
  -v nctn_counts="$NCTN_COUNTS" \
  -v ntcn_counts="$NTCN_COUNTS" \
  -v trinuc_counts="$TRINUC_COUNTS" \
  -v qc="$QC" '
BEGIN {
  FS = OFS = "\t"

  base[1] = "A"
  base[2] = "C"
  base[3] = "G"
  base[4] = "T"

  print \
    "record_id", \
    "sequence_minus3_to_0", \
    "lesion_type", \
    "tetranucleotide_context", \
    "trinucleotide_context", \
    "base_minus3", \
    "base_minus2", \
    "base_minus1", \
    "base_position0" \
    > events
}

{
  input_rows++

  if (NF < 2) {
    malformed_rows++
    next
  }

  record_id = $1
  sequence = toupper($2)

  if (length(sequence) != 4) {
    non_4nt_rows++
    next
  }

  if (sequence !~ /^[ACGT][ACGT][ACGT][ACGT]$/) {
    non_acgt_rows++
    next
  }

  valid_rows++

  base_minus3 = substr(sequence, 1, 1)
  base_minus2 = substr(sequence, 2, 1)
  base_minus1 = substr(sequence, 3, 1)
  base_0 = substr(sequence, 4, 1)

  dimer = substr(sequence, 2, 2)

  if (dimer == "CT") {
    lesion_type = "CT"
    tetranucleotide = sequence
    trinucleotide = substr(sequence, 1, 3)

    nctn_count[tetranucleotide]++
    nct_count[trinucleotide]++

    total_ct++
    selected_rows++

    print \
      record_id, \
      sequence, \
      lesion_type, \
      tetranucleotide, \
      trinucleotide, \
      base_minus3, \
      base_minus2, \
      base_minus1, \
      base_0 \
      > events

  } else if (dimer == "TC") {
    lesion_type = "TC"
    tetranucleotide = sequence
    trinucleotide = substr(sequence, 2, 3)

    ntcn_count[tetranucleotide]++
    tcn_count[trinucleotide]++

    total_tc++
    selected_rows++

    print \
      record_id, \
      sequence, \
      lesion_type, \
      tetranucleotide, \
      trinucleotide, \
      base_minus3, \
      base_minus2, \
      base_minus1, \
      base_0 \
      > events

  } else {
    other_dimer_rows++
  }
}

END {
  # ==========================================================
  # NCTN
  # ==========================================================

  print \
    "NCTN", \
    "count", \
    "total_CT_damage", \
    "percent_within_CT" \
    > nctn_counts

  for (left = 1; left <= 4; left++) {
    for (right = 1; right <= 4; right++) {
      context = base[left] "CT" base[right]
      value = nctn_count[context] + 0

      if (total_ct > 0) {
        percent = 100 * value / total_ct
      } else {
        percent = 0
      }

      print \
        context, \
        value, \
        total_ct + 0, \
        percent \
        > nctn_counts
    }
  }

  # ==========================================================
  # NTCN
  # ==========================================================

  print \
    "NTCN", \
    "count", \
    "total_TC_damage", \
    "percent_within_TC" \
    > ntcn_counts

  for (left = 1; left <= 4; left++) {
    for (right = 1; right <= 4; right++) {
      context = base[left] "TC" base[right]
      value = ntcn_count[context] + 0

      if (total_tc > 0) {
        percent = 100 * value / total_tc
      } else {
        percent = 0
      }

      print \
        context, \
        value, \
        total_tc + 0, \
        percent \
        > ntcn_counts
    }
  }

  # ==========================================================
  # NCT AND TCN
  # ==========================================================

  print \
    "lesion_type", \
    "trinucleotide", \
    "count", \
    "lesion_total", \
    "total_CT_TC_damage", \
    "fraction_within_lesion_type", \
    "fraction_of_all_CT_TC_damage", \
    "percent_within_lesion_type", \
    "percent_of_all_CT_TC_damage" \
    > trinuc_counts

  total_ct_tc = total_ct + total_tc

  for (left = 1; left <= 4; left++) {
    context = base[left] "CT"
    value = nct_count[context] + 0

    if (total_ct > 0) {
      fraction_class = value / total_ct
    } else {
      fraction_class = 0
    }

    if (total_ct_tc > 0) {
      fraction_all = value / total_ct_tc
    } else {
      fraction_all = 0
    }

    print \
      "CT", \
      context, \
      value, \
      total_ct + 0, \
      total_ct_tc + 0, \
      fraction_class, \
      fraction_all, \
      100 * fraction_class, \
      100 * fraction_all \
      > trinuc_counts
  }

  for (right = 1; right <= 4; right++) {
    context = "TC" base[right]
    value = tcn_count[context] + 0

    if (total_tc > 0) {
      fraction_class = value / total_tc
    } else {
      fraction_class = 0
    }

    if (total_ct_tc > 0) {
      fraction_all = value / total_ct_tc
    } else {
      fraction_all = 0
    }

    print \
      "TC", \
      context, \
      value, \
      total_tc + 0, \
      total_ct_tc + 0, \
      fraction_class, \
      fraction_all, \
      100 * fraction_class, \
      100 * fraction_all \
      > trinuc_counts
  }

  # ==========================================================
  # QC
  # ==========================================================

  print \
    "input_rows", \
    "valid_4nt_rows", \
    "malformed_rows", \
    "non_4nt_rows", \
    "non_ACGT_rows", \
    "CT_rows", \
    "TC_rows", \
    "selected_CT_TC_rows", \
    "other_dimer_rows" \
    > qc

  print \
    input_rows + 0, \
    valid_rows + 0, \
    malformed_rows + 0, \
    non_4nt_rows + 0, \
    non_acgt_rows + 0, \
    total_ct + 0, \
    total_tc + 0, \
    selected_rows + 0, \
    other_dimer_rows + 0 \
    > qc
}
' "$INPUT"

echo
echo "Done."
echo
echo "Events:"
echo "$EVENTS"
echo
echo "NCTN counts:"
echo "$NCTN_COUNTS"
echo
echo "NTCN counts:"
echo "$NTCN_COUNTS"
echo
echo "NCT/TCN counts:"
echo "$TRINUC_COUNTS"
echo
echo "QC:"
echo "$QC"