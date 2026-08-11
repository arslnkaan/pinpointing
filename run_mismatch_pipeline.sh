#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=24g
#SBATCH -t 12:00:00
#SBATCH --job-name=mismatch_dedup_all
#SBATCH --output=logs/mismatch_dedup_all_%A_%a.out
#SBATCH --error=logs/mismatch_dedup_all_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu
#SBATCH --array=0-1

set -euo pipefail

module load samtools
module load python/3.9

mkdir -p logs

# ============================================================
# SAMPLE LIST
# Add/remove samples here.
# Array index 0-7 matches 8 samples below.
# ============================================================

samples=(
  "Cumulus-4NQO-2h-r2"
  "Cumulus-4NQO-30m"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# ============================================================
# CONFIG
# ============================================================

ref_fasta="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"
scripts_dir="/work/users/a/r/arslank/mismatch_pipeline"

input_bam="/work/users/a/r/arslank/${sample}_atl_output/genome/${sample}.sorted.bam"

out="/work/users/a/r/arslank/${sample}_mismatch_pipeline"
dedup_dir="${out}/00_dedup_bam"

mkdir -p "$out" "$dedup_dir"

mapq=20
min_len=20
max_len=30
from3_min=6
from3_max=13
changes="G>T,G>A,G>C"

dedup_bam="${dedup_dir}/${sample}.dedup.bam"

echo "============================================================"
echo "Sample: $sample"
echo "Input BAM: $input_bam"
echo "Dedup BAM: $dedup_bam"
echo "Output: $out"
echo "============================================================"

if [[ ! -f "$input_bam" ]]; then
  echo "ERROR: input BAM does not exist:"
  echo "$input_bam"
  exit 1
fi

# ============================================================
# 00. Deduplicate BAM
# This removes coordinate/PCR duplicates using samtools markdup.
# ============================================================

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

# ============================================================
# 01. Call mismatches from deduplicated BAM
# ============================================================

echo "[1/4] Calling mismatches"

python3 "$scripts_dir/01_call_mismatches.py" \
  --bam "$dedup_bam" \
  --fasta "$ref_fasta" \
  --sample "$sample" \
  --outdir "$out/01_mismatch_calls" \
  --mapq "$mapq"

# ============================================================
# 02. Filter mismatch events
# ============================================================

echo "[2/4] Filtering mismatch events"

python3 "$scripts_dir/02_filter_mismatch_events.py" \
  --mismatches "$out/01_mismatch_calls/${sample}_mismatches.tsv" \
  --sample "$sample" \
  --outdir "$out/02_filtered_events" \
  --min-len "$min_len" \
  --max-len "$max_len" \
  --from3-min "$from3_min" \
  --from3-max "$from3_max" \
  --changes "$changes"

# ============================================================
# 03. Plot summary figures
# ============================================================

echo "[3/4] Plotting mismatch summaries"

python3 "$scripts_dir/03_plot_mismatch_summary.py" \
  --sample "$sample" \
  --indir "$out/02_filtered_events" \
  --outdir "$out/03_plots" \
  --min-len "$min_len" \
  --max-len "$max_len" \
  --from3-min "$from3_min" \
  --from3-max "$from3_max" \
  --changes "$changes"

# ============================================================
# 04. Plot position barplots
# ============================================================

echo "[4/4] Plot mismatch position barplots"

python3 "$scripts_dir/04_plot_mismatch_position_bars.py" \
  --events "${out}/02_filtered_events/${sample}_singleMismatch_G_to_T_${from3_min}to${from3_max}nt_from3prime_${min_len}to${max_len}mers.tsv,${out}/02_filtered_events/${sample}_singleMismatch_G_to_C_${from3_min}to${from3_max}nt_from3prime_${min_len}to${max_len}mers.tsv,${out}/02_filtered_events/${sample}_singleMismatch_G_to_A_${from3_min}to${from3_max}nt_from3prime_${min_len}to${max_len}mers.tsv" \
  --sample "$sample" \
  --outdir "$out/04_position_barplots" \
  --min-len 26 \
  --max-len 30 \
  --changes "G>T,G>C,G>A" \
  --position 5prime \
  --percent

echo "DONE: $sample"
echo "Outputs: $out"