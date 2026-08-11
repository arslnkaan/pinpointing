#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --mem=32g
#SBATCH -t 08:00:00
#SBATCH --job-name=damage_context
#SBATCH --output=logs/damage_context_%j.out
#SBATCH --error=logs/damage_context_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=arslank@email.unc.edu

set -euo pipefail

mkdir -p logs

module purge
module load python/3.9
module load bedtools
module load samtools

BASE="/work/users/a/r/arslank"

PYTHON_SCRIPT="${BASE}/uv/03_extract_damageseq_contexts.py"

REFERENCE="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

CT_TC_INPUT="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CT_TC_Ccenter_3nt.sequence.tsv"

CC_INPUT="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CC_Ccenter_3nt.sequence.tsv"

OUTDIR="${BASE}/damseq/NHF1_CPD_0h_r1_results/damage_context_denominators"

mkdir -p "$OUTDIR"

python "$PYTHON_SCRIPT" \
  --ct-tc "$CT_TC_INPUT" \
  --cc "$CC_INPUT" \
  --fasta "$REFERENCE" \
  --outdir "$OUTDIR" \
  --cc-center-mode auto

echo
echo "Done."
echo "Output directory: ${OUTDIR}"