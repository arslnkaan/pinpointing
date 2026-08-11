#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --mem=24g
#SBATCH -t 12:00:00
#SBATCH --job-name=mismatch_dedup_2reps_allpos
#SBATCH --output=logs/mismatch_dedup_2reps_allpos_%A_%a.out
#SBATCH --error=logs/mismatch_dedup_2reps_allpos_%A_%a.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu
#SBATCH --array=0-4

set -euo pipefail

module load samtools
module load python/3.9

mkdir -p logs

# ============================================================
# SAMPLE LIST
# Array index 0-1 matches 2 samples below.
# Replace these with the two repeats you want.
# ============================================================

samples=(
  "Ap13NHFWT2hNQO5N_R1_001"
  "Ap13NHFWT2hNQO30T_R1_001"
  "Ap13NHFDDB22hNQO_R1_001"
  "Ap13NHFCSB2hNQO_R1_001"
  "Ap13NHFXPC2hNQO_R1_001"
)

sample="${samples[$SLURM_ARRAY_TASK_ID]}"

# ============================================================
# CONFIG
# ============================================================

ref_fasta="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"
scripts_dir="/work/users/a/r/arslank/mismatch_pipeline"

input_bam="/work/users/a/r/arslank/${sample}_atl_output/genome/${sample}.sorted.bam"

out="/work/users/a/r/arslank/${sample}_mismatch_pipeline_all_mismatches"
dedup_dir="${out}/00_dedup_bam"

mkdir -p "$out" "$dedup_dir"

mapq=20
min_len=20
max_len=30

# Look for mismatches at every allowed position from the 3' end.
# For 20-30 nt reads, this effectively covers the whole read.
from3_min=1
from3_max=30

changes="A>C,A>G,A>T,C>A,C>G,C>T,G>A,G>C,G>T,T>A,T>C,T>G"

dedup_bam="${dedup_dir}/${sample}.dedup.bam"

echo "============================================================"
echo "Sample: $sample"
echo "Input BAM: $input_bam"
echo "Dedup BAM: $dedup_bam"
echo "Output: $out"
echo "MAPQ: $mapq"
echo "Read length: ${min_len}-${max_len}"
echo "Mismatch window from 3' end: ${from3_min}-${from3_max}"
echo "Changes: $changes"
echo "============================================================"

if [[ ! -f "$input_bam" ]]; then
  echo "ERROR: input BAM does not exist:"
  echo "$input_bam"
  exit 1
fi

# ============================================================
# 00. Deduplicate BAM
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

echo "[4/4] Plot mismatch position barplots"

event_files=""
for chg in A_to_C A_to_G A_to_T C_to_A C_to_G C_to_T G_to_A G_to_C G_to_T T_to_A T_to_C T_to_G; do
  f="${out}/02_filtered_events/${sample}_singleMismatch_${chg}_${from3_min}to${from3_max}nt_from3prime_${min_len}to${max_len}mers.tsv"
  if [[ -f "$f" ]]; then
    if [[ -z "$event_files" ]]; then
      event_files="$f"
    else
      event_files="${event_files},${f}"
    fi
  else
    echo "WARNING missing event file: $f"
  fi
done

python3 "$scripts_dir/04_plot_mismatch_position_bars.py" \
  --events "$event_files" \
  --sample "$sample" \
  --outdir "$out/04_position_barplots" \
  --min-len 26 \
  --max-len 30 \
  --changes "A>C,A>G,A>T,C>A,C>G,C>T,G>A,G>C,G>T,T>A,T>C,T>G" \
  --position 5prime \
  --percent

echo "DONE: $sample"
echo "Outputs: $out"