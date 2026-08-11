#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=64g
#SBATCH -t 12:00:00
#SBATCH --job-name=mismatch_CtoT
#SBATCH --output=logs/mismatch_CtoT_%A_%a.out
#SBATCH --error=logs/mismatch_CtoT_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu
#SBATCH --array=0-12

set -euo pipefail

module load samtools
module load python/3.9

mkdir -p logs

samples=(
"XPC-UVCPD-30m"
"XPC-UVCPD-30m-r2"
"XPC-UVCPD-2h"
"XPC-UVCPD-2h-r2"
"XPC-UVCPD-4h"
"XPC-UVCPD-4h-r2"
"XPC-UVCPD-8h"
"XPC-UVCPD-8h-r2"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

ref_fasta="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"
scripts_dir="/work/users/a/r/arslank/mismatch_pipeline"

input_bam="/work/users/a/r/arslank/${sample}_atl_output/genome/${sample}.sorted.bam"

out="/work/users/a/r/arslank/${sample}_mismatch_pipeline_CtoT"
dedup_dir="${out}/00_dedup_bam"

mkdir -p "$out" "$dedup_dir"

mapq=20
min_len=20
max_len=30
from3_min=6
from3_max=13
changes="C>T"

dedup_bam="${dedup_dir}/${sample}.dedup.bam"

echo "============================================================"
echo "Sample: $sample"
echo "Input BAM: $input_bam"
echo "Dedup BAM: $dedup_bam"
echo "Output: $out"
echo "Change: $changes"
echo "============================================================"

if [[ ! -f "$input_bam" ]]; then
  echo "ERROR: input BAM does not exist:"
  echo "$input_bam"
  exit 1
fi

echo "[0/4] Deduplicating BAM with samtools markdup"

if [[ ! -f "$dedup_bam" ]]; then

  name_sorted_bam="${dedup_dir}/${sample}.namesort.bam"
  fixmate_bam="${dedup_dir}/${sample}.fixmate.bam"
  coord_sorted_bam="${dedup_dir}/${sample}.fixmate.coordsort.bam"

  samtools sort -n -@ 8 \
    -o "$name_sorted_bam" \
    "$input_bam"

  samtools fixmate -m \
    "$name_sorted_bam" \
    "$fixmate_bam"

  samtools sort -@ 8 \
    -o "$coord_sorted_bam" \
    "$fixmate_bam"

  samtools markdup -r -@ 8 \
    "$coord_sorted_bam" \
    "$dedup_bam"

  samtools index "$dedup_bam"

  rm -f "$name_sorted_bam" "$fixmate_bam" "$coord_sorted_bam"

else
  echo "Dedup BAM already exists, skipping deduplication."
fi

if [[ ! -f "${dedup_bam}.bai" && ! -f "${dedup_bam%.bam}.bai" ]]; then
  samtools index "$dedup_bam"
fi

echo "[1/4] Calling mismatches"

python3 "$scripts_dir/01_call_mismatches.py" \
  --bam "$dedup_bam" \
  --fasta "$ref_fasta" \
  --sample "$sample" \
  --outdir "$out/01_mismatch_calls" \
  --mapq "$mapq"

echo "[2/4] Filtering C>T mismatch events"

python3 "$scripts_dir/02_filter_mismatch_events.py" \
  --mismatches "$out/01_mismatch_calls/${sample}_mismatches.tsv" \
  --sample "$sample" \
  --outdir "$out/02_filtered_events" \
  --min-len "$min_len" \
  --max-len "$max_len" \
  --from3-min "$from3_min" \
  --from3-max "$from3_max" \
  --changes "$changes"

echo "[3/4] Plotting C>T mismatch summaries"

python3 "$scripts_dir/03_plot_mismatch_summary.py" \
  --sample "$sample" \
  --indir "$out/02_filtered_events" \
  --outdir "$out/03_plots" \
  --min-len "$min_len" \
  --max-len "$max_len" \
  --from3-min "$from3_min" \
  --from3-max "$from3_max" \
  --changes "$changes"

echo "[4/4] Plot C>T mismatch position barplots"

python3 "$scripts_dir/04_plot_mismatch_position_bars.py" \
  --events "${out}/02_filtered_events/${sample}_singleMismatch_C_to_T_${from3_min}to${from3_max}nt_from3prime_${min_len}to${max_len}mers.tsv" \
  --sample "$sample" \
  --outdir "$out/04_position_barplots" \
  --min-len 26 \
  --max-len 30 \
  --changes "C>T" \
  --position 5prime \
  --percent

echo "DONE: $sample"
echo "Outputs: $out"