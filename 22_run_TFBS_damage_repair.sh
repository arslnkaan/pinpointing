#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --mem=32g
#SBATCH -t 18:00:00
#SBATCH --job-name=TFBS_UV
#SBATCH --output=logs/TFBS_UV_%j.out
#SBATCH --error=logs/TFBS_UV_%j.err
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

TFBS="${BASE}/active_TFBS.bed"

DAMAGE="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_CC_NCCN_from_read_starts/NHF1_CPD_0h.minus3_to_0.raw_sequence.tsv"

PYTHON_SCRIPT="./22_prepare_TFBS_damage_repair.py"

OUTDIR="${BASE}/UV_TFBS_damage_repair"

mkdir -p "$OUTDIR"

REPAIR_MAP="${OUTDIR}/UV_repair_sample_map.tsv"

# ============================================================
# REPAIR SAMPLE MAP
# ============================================================

cat > "$REPAIR_MAP" <<EOF
sample	timepoint	time_h	replicate	input
NHF1-UVCPD-30m	0.5h	0.5	R1	${BASE}/NHF1-UVCPD-30m_mismatch_pipeline_CtoT/02_filtered_events/NHF1-UVCPD-30m_singleMismatch_C_to_T_6to13nt_from3prime_20to30mers.tsv
NHF1-UVCPD-2h	2h	2	R1	${BASE}/NHF1-UVCPD-2h_mismatch_pipeline_CtoT/02_filtered_events/NHF1-UVCPD-2h_singleMismatch_C_to_T_6to13nt_from3prime_20to30mers.tsv
NHF1-UVCPD-4h	4h	4	R1	${BASE}/NHF1-UVCPD-4h_mismatch_pipeline_CtoT/02_filtered_events/NHF1-UVCPD-4h_singleMismatch_C_to_T_6to13nt_from3prime_20to30mers.tsv
NHF1-UVCPD-8h	8h	8	R1	${BASE}/NHF1-UVCPD-8h_mismatch_pipeline_CtoT/02_filtered_events/NHF1-UVCPD-8h_singleMismatch_C_to_T_6to13nt_from3prime_20to30mers.tsv
EOF

# ============================================================
# CHECK INPUTS
# ============================================================

for path in \
  "$TFBS" \
  "$DAMAGE" \
  "$PYTHON_SCRIPT"
do
  if [[ ! -s "$path" ]]; then
    echo "ERROR: missing input:"
    echo "$path"
    exit 1
  fi
done

# ============================================================
# PREPARE EVENT FILES
# ============================================================

python "$PYTHON_SCRIPT" \
  --tfbs "$TFBS" \
  --damage "$DAMAGE" \
  --repair-map "$REPAIR_MAP" \
  --outdir "$OUTDIR"

TFBS_EXPLODED="${OUTDIR}/active_TFBS.exploded.bed"
TFBS_UNIQUE="${OUTDIR}/active_TFBS.unique_loci.bed"

DAMAGE_BED="${OUTDIR}/DamageSeq_minus2_minus1_dipyrimidines.bed"
REPAIR_BED="${OUTDIR}/UV_repair_pinpointed_1bp.bed"

# ============================================================
# CREATE MERGED ALL-TFBS BED
#
# This counts each active TFBS genomic region only once for
# the pooled active-TFBS analysis.
# ============================================================

TFBS_MERGED="${OUTDIR}/active_TFBS.all_loci_merged.bed"

LC_ALL=C sort \
  -k1,1 \
  -k2,2n \
  "$TFBS_UNIQUE" \
| bedtools merge \
  -i - \
> "$TFBS_MERGED"

# ============================================================
# TF-SPECIFIC OVERLAPS
#
# Repair A columns: 10
# TFBS B columns:    4
#
# Damage A columns: 10
# TFBS B columns:    4
# ============================================================

bedtools intersect \
  -a "$REPAIR_BED" \
  -b "$TFBS_EXPLODED" \
  -wa \
  -wb \
  > "${OUTDIR}/UV_repair_pinpointed_TFBS_overlaps.tsv"

bedtools intersect \
  -a "$DAMAGE_BED" \
  -b "$TFBS_EXPLODED" \
  -wa \
  -wb \
  > "${OUTDIR}/DamageSeq_dipyrimidine_TFBS_overlaps.tsv"

# ============================================================
# POOLED ACTIVE-TFBS OVERLAPS
#
# -u keeps each event only once, regardless of how many TFBS
# intervals it overlaps.
# ============================================================

bedtools intersect \
  -a "$REPAIR_BED" \
  -b "$TFBS_MERGED" \
  -u \
  > "${OUTDIR}/UV_repair_pinpointed_any_active_TFBS.bed"

bedtools intersect \
  -a "$DAMAGE_BED" \
  -b "$TFBS_MERGED" \
  -u \
  > "${OUTDIR}/DamageSeq_dipyrimidine_any_active_TFBS.bed"

# ============================================================
# BASIC COUNTS
# ============================================================

{
  echo -e "file\tcount"
  echo -e "all_repair_events\t$(wc -l < "$REPAIR_BED")"
  echo -e "repair_events_at_any_TFBS\t$(wc -l < "${OUTDIR}/UV_repair_pinpointed_any_active_TFBS.bed")"
  echo -e "all_damage_events\t$(wc -l < "$DAMAGE_BED")"
  echo -e "damage_events_at_any_TFBS\t$(wc -l < "${OUTDIR}/DamageSeq_dipyrimidine_any_active_TFBS.bed")"
  echo -e "exploded_TFBS_assignments\t$(wc -l < "$TFBS_EXPLODED")"
  echo -e "unique_original_TFBS_loci\t$(wc -l < "$TFBS_UNIQUE")"
  echo -e "merged_active_TFBS_intervals\t$(wc -l < "$TFBS_MERGED")"
} > "${OUTDIR}/TFBS_overlap_basic_QC.tsv"

echo
echo "============================================================"
echo "Done"
echo "============================================================"
echo
echo "Output directory:"
echo "$OUTDIR"