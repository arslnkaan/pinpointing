#!/bin/bash
#SBATCH --job-name=damageseq
#SBATCH --output=damageseq_%j.out
#SBATCH --error=damageseq_%j.err
#SBATCH --time=08:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G

set -euo pipefail

module load bedtools

sample="NHF1-4NQO-30m-r1"
OUTDIR="results/damage_sites"
GENOME_FA="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

BEDPE="${OUTDIR}/${sample}.dedup.bedpe"

if [[ ! -s "$BEDPE" ]]; then
  echo "ERROR: missing or empty input BEDPE: $BEDPE"
  exit 1
fi

awk 'BEGIN{OFS="\t"}{
  strand=$9

  if (strand == "+") {
    s=$2-3
    e=$2+1
  } else if (strand == "-") {
    s=$3-1
    e=$3+3
  } else {
    next
  }

  if (s >= 0) {
    print $1,s,e,$7,$8,strand
  }
}' "$BEDPE" \
  > "${OUTDIR}/${sample}.readstart_minus3_to_0.bed"

echo "Window BED lines:"
wc -l "${OUTDIR}/${sample}.readstart_minus3_to_0.bed"

bedtools getfasta \
  -fi "$GENOME_FA" \
  -bed "${OUTDIR}/${sample}.readstart_minus3_to_0.bed" \
  -s \
  -tab \
  -name \
  > "${OUTDIR}/${sample}.readstart_minus3_to_0.sequence.tsv"

echo "Sequence lines:"
wc -l "${OUTDIR}/${sample}.readstart_minus3_to_0.sequence.tsv"