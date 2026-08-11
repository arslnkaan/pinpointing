#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH -t 12:00:00
#SBATCH --job-name=CCTT_full_reads
#SBATCH --output=logs/CCTT_full_reads_%j.out
#SBATCH --error=logs/CCTT_full_reads_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

mkdir -p logs

module purge
module load python/3.9
module load samtools

BASE="/work/users/a/r/arslank"

PYTHON_SCRIPT="${BASE}/uvv/01_reconstruct_CCTT_reads_and_NCCN.py"

MAPQ=10

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

for sample in "${SAMPLES[@]}"; do

  rescue_dir="${BASE}/${sample}_CC_to_TT_rescue"

  event_file="${rescue_dir}/${sample}.CC_to_TT_tandem_events.tsv"

  rescue_bam="${rescue_dir}/${sample}.single_TT_to_CC.v0.sorted.bam"

  sample_outdir="${rescue_dir}/full_reads_and_NCCN"

  echo
  echo "============================================================"
  echo "Sample: ${sample}"
  echo "============================================================"

  if [[ ! -s "$event_file" ]]; then
    echo "ERROR: Missing event table:"
    echo "$event_file"
    exit 1
  fi

  if [[ ! -s "$rescue_bam" ]]; then
    echo "ERROR: Missing rescued BAM:"
    echo "$rescue_bam"
    exit 1
  fi

  mkdir -p "$sample_outdir"

  python "$PYTHON_SCRIPT" \
    --sample "$sample" \
    --events "$event_file" \
    --rescue-bam "$rescue_bam" \
    --outdir "$sample_outdir" \
    --mapq "$MAPQ"

done

echo
echo "All samples completed."