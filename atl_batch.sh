#!/usr/bin/env bash
#SBATCH -p general
#SBATCH -N 1
#SBATCH -t 12:00:00
#SBATCH --mem=12g
#SBATCH -n 4
#SBATCH --mail-type=ALL
#SBATCH --mail-user=arslank@email.unc.edu
#SBATCH --output=/work/users/a/r/arslank/logs/atl_batch_%j_4NQO.out
#SBATCH --error=/work/users/a/r/arslank/logs/atl_batch_%j_4NQO.err

set -euo pipefail

SAMPLE_LIST="/work/users/a/r/arslank/samples_atl.txt"

THREADS=8
MAPQ=20

UMI_LEN=10
ADAPTER3="TGGAATTCTCGGGTGCCAAGG"
OVERLAP_ADAPTER=8
POLYA="A{10}"
OVERLAP_A=10

SPIKE_CORE="ATTACGCACTACTGGTCA"
OVERLAP_SPIKE=10
SPIKE_TOTAL_LEN=28
SPIKE_LEFT=5
SPIKE_RIGHT=5
TARGET_SPIKE_UMIS=""

GENOME_DIR="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta"
BOWTIE2_IND="/proj/seq/data/hg38_UCSC/Sequence/Bowtie2Index/genome"
FEATURE_BED="/work/users/a/r/arslank/beds/hg38_5kb_kaan.bed"

LEN_MIN=24
LEN_MAX=30

MONOMER_PY="/work/users/a/r/arslank/scripts/monomer_dimer_batch.py"
FIGURES_R="/work/users/a/r/arslank/scripts/figures_batch.R"

mkdir -p /work/users/a/r/arslank/logs

module purge
module load anaconda
source "$(conda info --base)/etc/profile.d/conda.sh"
conda activate umi_env

module load bowtie2
module load samtools
module load bedtools
module load ucsctools

for t in cutadapt bowtie2 samtools bedtools bedGraphToBigWig python; do
  command -v "$t" >/dev/null 2>&1 || { echo "ERROR: $t not found"; exit 1; }
done

echo "=== BATCH START $(date) ==="

while IFS=$'\t' read -r sample FASTQ HAS_UMI; do

  [[ -z "${sample:-}" ]] && continue
  [[ "$sample" =~ ^# ]] && continue

  echo "======================================"
  echo "Running sample: $sample"
  echo "FASTQ: $FASTQ"
  echo "HAS_UMI: $HAS_UMI"
  echo "======================================"

  OUT="/work/users/a/r/arslank/${sample}_atl_output"
  FIGOUT="/work/users/a/r/arslank/${sample}_figs"

  mkdir -p "$OUT"/{qc,tmp,spike,genome,monomer,spike_qc}
  mkdir -p "$FIGOUT"

  # 1. 3' adapter trim
  echo "[1/13] 3' adapter trim"
  cutadapt -a "$ADAPTER3" -O "$OVERLAP_ADAPTER" -q 30 -j "$THREADS" \
    --discard-untrimmed -m 10 \
    -o "$OUT/${sample}.no3p.fastq" \
    "$FASTQ" \
    > "$OUT/qc/${sample}.cutadapt_3p.log"

  # 2. UMI extraction
  PRE_POLYA="$OUT/${sample}.no3p.fastq"

  if [[ "$HAS_UMI" == "yes" ]]; then
    echo "[2/13] Extracting UMI"

    python - <<PY
import gzip, sys

inp  = "$OUT/${sample}.no3p.fastq"
outp = "$OUT/${sample}.umi.fastq"
k    = int("$UMI_LEN")

def opn(path, mode):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)

n=0; bad=0; short=0
with opn(inp,"rt") as fin, opn(outp,"wt") as fout:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline().rstrip("\\n")
        p = fin.readline()
        q = fin.readline().rstrip("\\n")

        h0 = h.rstrip("\\n")
        if not h0.startswith("@"):
            bad += 1
            continue

        first = h0.split()[0][1:]

        if len(s) != len(q):
            sys.exit("Seq/qual mismatch: " + h0)

        if len(s) <= k:
            short += 1
            continue

        umi = s[-k:]
        s2  = s[:-k]
        q2  = q[:-k]

        fout.write(f"@{first}_{umi}\\n{s2}\\n+\\n{q2}\\n")
        n += 1

print(f"Wrote {n} reads; bad_headers={bad}; too_short={short}", file=sys.stderr)
PY

    PRE_POLYA="$OUT/${sample}.umi.fastq"
  else
    echo "[2/13] No UMI extraction"
  fi

  # 3. polyA trim
  echo "[3/13] polyA trim"
  cutadapt -a "$POLYA" -O "$OVERLAP_A" -q 30 -j "$THREADS" \
    --discard-untrimmed -m 10 \
    -o "$OUT/${sample}.clean.fastq" \
    "$PRE_POLYA" \
    > "$OUT/qc/${sample}.cutadapt_polyA.log"

  # 4. Split spike/nonspike
  echo "[4/13] split spike/nonspike"
  cutadapt -b "$SPIKE_CORE" -O "$OVERLAP_SPIKE" -j "$THREADS" --action=none \
    -o "$OUT/spike/${sample}.spike.fastq" \
    --untrimmed-output "$OUT/spike/${sample}.nonspike.fastq" \
    "$OUT/${sample}.clean.fastq" \
    > "$OUT/qc/${sample}.split_spike.log"

  # 5. Spike strict filter
  echo "[5/13] spike strict filter"

  SPIKE_IN="$OUT/spike/${sample}.spike.fastq"
  SPIKE_STRICT="$OUT/spike/${sample}.spike.strict.fastq"
  COUNTS_TXT="$OUT/qc/${sample}.spike.strict_counts.txt"

  if [[ ! -s "$SPIKE_IN" ]]; then
    : > "$SPIKE_STRICT"
    echo "0 0 0" > "$COUNTS_TXT"
  else
    python - <<PY > "$COUNTS_TXT"
import re

inp = "$SPIKE_IN"
outp = "$SPIKE_STRICT"
core = "$SPIKE_CORE"
L = int("$SPIKE_LEFT")
R = int("$SPIKE_RIGHT")
TOTAL = int("$SPIKE_TOTAL_LEN")

pat = re.compile(rf"^[ACGT]{{{L}}}{re.escape(core)}[ACGT]{{{R}}}$")

total_reads = 0
strict_reads = 0
uniq = set()

with open(inp, "rt") as fin, open(outp, "wt") as fout:
    while True:
        h = fin.readline()
        if not h:
            break
        s = fin.readline().rstrip("\\n")
        p = fin.readline()
        q = fin.readline().rstrip("\\n")

        total_reads += 1

        if len(s) != len(q):
            continue
        if len(s) != TOTAL:
            continue
        if not pat.match(s):
            continue

        fout.write(h)
        fout.write(s + "\\n")
        fout.write(p)
        fout.write(q + "\\n")

        strict_reads += 1
        uniq.add(s[:L] + "_" + s[-R:])

print(total_reads, strict_reads, len(uniq))
PY
  fi

  read -r SPIKE_TOTAL_READS SPIKE_STRICT_READS SPIKE_UNIQ < "$COUNTS_TXT" || true

  echo -e "sample\tspike_total_reads\tspike_strict_reads\tspike_strict_unique_5p5p" \
    > "$OUT/spike/${sample}.spike.strict_stats.tsv"
  echo -e "${sample}\t${SPIKE_TOTAL_READS}\t${SPIKE_STRICT_READS}\t${SPIKE_UNIQ}" \
    >> "$OUT/spike/${sample}.spike.strict_stats.tsv"

  SCALE="1.0"
  if [[ -n "${TARGET_SPIKE_UMIS}" ]]; then
    SCALE="$(python -c "t=float('$TARGET_SPIKE_UMIS'); s=float('$SPIKE_UNIQ'); print('NaN' if s==0 else f'{t/s:.12f}')")"
  fi

  echo -e "sample\tspike_strict_unique_5p5p\ttarget_spike_umis\tscale" \
    > "$OUT/spike/${sample}.spike_scale.tsv"
  echo -e "${sample}\t${SPIKE_UNIQ}\t${TARGET_SPIKE_UMIS:-NA}\t${SCALE}" \
    >> "$OUT/spike/${sample}.spike_scale.tsv"

  # 6. Spike QC
  echo "[6/13] spike QC"

  python - <<PY
import os, re
from collections import Counter
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sample = "$sample"
spike_fastq  = "$SPIKE_IN"
strict_fastq = "$SPIKE_STRICT"
outdir = "$OUT/spike_qc"
os.makedirs(outdir, exist_ok=True)

core = "$SPIKE_CORE"
L = int("$SPIKE_LEFT")
R = int("$SPIKE_RIGHT")
TOTAL = int("$SPIKE_TOTAL_LEN")
pat = re.compile(rf"^[ACGT]{{{L}}}{re.escape(core)}[ACGT]{{{R}}}$")

def fastq_seqs(path):
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return
    with open(path, "rt") as f:
        while True:
            h = f.readline()
            if not h:
                break
            s = f.readline().rstrip("\\n")
            f.readline()
            f.readline()
            yield s

len_counts = Counter()
total_core = 0

for s in fastq_seqs(spike_fastq):
    total_core += 1
    len_counts[len(s)] += 1

with open(os.path.join(outdir, f"{sample}.spike_core.lengths.tsv"), "w") as out:
    out.write("length\\tcount\\n")
    for k in sorted(len_counts):
        out.write(f"{k}\\t{len_counts[k]}\\n")

fail_len = 0
fail_pat = 0
pass_strict = 0

for s in fastq_seqs(spike_fastq):
    if len(s) != TOTAL:
        fail_len += 1
    elif not pat.match(s):
        fail_pat += 1
    else:
        pass_strict += 1

with open(os.path.join(outdir, f"{sample}.spike_strict.fail_reasons.tsv"), "w") as out:
    out.write("category\\tcount\\n")
    out.write(f"pass_strict\\t{pass_strict}\\n")
    out.write(f"fail_length\\t{fail_len}\\n")
    out.write(f"fail_pattern\\t{fail_pat}\\n")

umi_counts = Counter()
strict_reads = 0

for s in fastq_seqs(strict_fastq):
    strict_reads += 1
    umi_counts[s[:L] + "_" + s[-R:]] += 1

with open(os.path.join(outdir, f"{sample}.spike_strict.umi_counts.tsv"), "w") as out:
    out.write("umi_left_right\\treads\\n")
    for k, v in umi_counts.most_common():
        out.write(f"{k}\\t{v}\\n")

dup = Counter()
for v in umi_counts.values():
    dup[min(v, 50)] += 1

with open(os.path.join(outdir, f"{sample}.spike_strict.duplication_dist.tsv"), "w") as out:
    out.write("reads_per_umi\\tumis\\n")
    for k in range(1, 50):
        out.write(f"{k}\\t{dup.get(k,0)}\\n")
    out.write(f"50+\\t{dup.get(50,0)}\\n")

bases = ["A", "C", "G", "T", "N"]
mat = [Counter() for _ in range(TOTAL)]

for s in fastq_seqs(strict_fastq):
    if len(s) != TOTAL:
        continue
    for i, ch in enumerate(s):
        mat[i][ch if ch in bases else "N"] += 1

with open(os.path.join(outdir, f"{sample}.spike_strict.basecomp.tsv"), "w") as out:
    out.write("pos\\tA\\tC\\tG\\tT\\tN\\ttotal\\n")
    for i in range(TOTAL):
        tot = sum(mat[i].values())
        out.write(f"{i+1}\\t{mat[i].get('A',0)}\\t{mat[i].get('C',0)}\\t{mat[i].get('G',0)}\\t{mat[i].get('T',0)}\\t{mat[i].get('N',0)}\\t{tot}\\n")

xs = sorted(len_counts)
ys = [len_counts[x] for x in xs]
plt.figure()
plt.bar(xs, ys)
plt.xlabel("Spike-core read length (nt)")
plt.ylabel("Reads")
plt.title(f"{sample} spike-core length distribution")
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{sample}.spike_core.lengths.png"), dpi=200)
plt.close()

cats = ["pass_strict", "fail_length", "fail_pattern"]
vals = [pass_strict, fail_len, fail_pat]
plt.figure()
plt.bar(cats, vals)
plt.ylabel("Reads")
plt.title(f"{sample} strict filtering outcomes")
plt.xticks(rotation=20, ha="right")
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{sample}.spike_strict.fail_reasons.png"), dpi=200)
plt.close()

xs = list(range(1, 50)) + [50]
ys = [dup.get(i,0) for i in range(1,50)] + [dup.get(50,0)]
plt.figure()
plt.bar(xs, ys)
plt.xlabel("Reads per unique spike molecule; 50=50+")
plt.ylabel("Unique spike molecules")
plt.title(f"{sample} spike duplication")
plt.tight_layout()
plt.savefig(os.path.join(outdir, f"{sample}.spike_strict.duplication_dist.png"), dpi=200)
plt.close()

summary = os.path.join(outdir, f"{sample}.spike_qc_summary.tsv")
uniq = len(umi_counts)
rpu = strict_reads / uniq if uniq else float("nan")

with open(summary, "w") as out:
    out.write("sample\\tspike_core_reads\\tpass_strict\\tfail_length\\tfail_pattern\\tstrict_unique_umis\\tstrict_reads\\treads_per_unique\\n")
    out.write(f"{sample}\\t{total_core}\\t{pass_strict}\\t{fail_len}\\t{fail_pat}\\t{uniq}\\t{strict_reads}\\t{rpu:.4f}\\n")
PY

  # 7. Align
  echo "[7/13] align nonspike"

  bowtie2 -x "$BOWTIE2_IND" -U "$OUT/spike/${sample}.nonspike.fastq" \
    --very-sensitive -p "$THREADS" \
    2> "$OUT/qc/${sample}.bt2.log" \
  | samtools view -b -F 0x904 -q "$MAPQ" - \
  | samtools sort -@ "$THREADS" -o "$OUT/genome/${sample}.sorted.bam" -

  samtools index "$OUT/genome/${sample}.sorted.bam"

  # 8. Dedup
  echo "[8/13] dedup"

  bedtools bamtobed -i "$OUT/genome/${sample}.sorted.bam" \
    > "$OUT/genome/${sample}.bed"

  if [[ "$HAS_UMI" == "yes" ]]; then
    LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 "$OUT/genome/${sample}.bed" \
    | awk -v k="$UMI_LEN" 'BEGIN{OFS="\t"}
        {
          name=$4
          if(length(name) < k) next
          umi=substr(name, length(name)-k+1, k)
          key=$1 FS $2 FS $3 FS $6 FS umi
          if(!(key in seen)){
            seen[key]=1
            print $0
          }
        }' \
    > "$OUT/genome/${sample}.dedup.bed"
  else
    awk 'BEGIN{OFS="\t"}{len=$3-$2; print $0, len}' "$OUT/genome/${sample}.bed" \
    | LC_ALL=C sort -k1,1 -k2,2n -k3,3n -k6,6 -k7,7 \
    | awk 'BEGIN{OFS="\t"}{key=$1 FS $2 FS $3 FS $6 FS $7; if(!(key in seen)){seen[key]=1; print $1,$2,$3,$4,$5,$6}}' \
    > "$OUT/genome/${sample}.dedup.bed"
  fi

  DEDUP_N=$(wc -l < "$OUT/genome/${sample}.dedup.bed")
  echo "$DEDUP_N" > "$OUT/${sample}.dedup_molecules.txt"
  echo "dedup_molecules = $DEDUP_N"

  # 9. Length distribution
  echo "[9/13] length distribution"

  awk '{print $3-$2}' "$OUT/genome/${sample}.dedup.bed" \
  | sort -n \
  | uniq -c \
  | awk 'BEGIN{OFS="\t"; print "Length","Count"}{print $2,$1}' \
  > "$OUT/${sample}_read_length_distribution.txt"

  # 10. BigWigs
  echo "[10/13] bedGraph/bigWig"

  awk '$6=="+"' "$OUT/genome/${sample}.dedup.bed" > "$OUT/genome/${sample}.dedup.plus.bed"
  awk '$6=="-"' "$OUT/genome/${sample}.dedup.bed" > "$OUT/genome/${sample}.dedup.minus.bed"

  SCALE_POS="$SCALE"
  SCALE_NEG=$(python -c "print(f'{-1.0*float(\"$SCALE\"):.12f}')")

  bedtools genomecov -bg -scale "$SCALE_POS" \
    -i "$OUT/genome/${sample}.dedup.plus.bed" \
    -g "$GENOME_DIR/genome.fa.fai" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > "$OUT/${sample}.plus.scaled.bdg"

  bedtools genomecov -bg -scale "$SCALE_NEG" \
    -i "$OUT/genome/${sample}.dedup.minus.bed" \
    -g "$GENOME_DIR/genome.fa.fai" \
  | LC_ALL=C sort -k1,1 -k2,2n \
  > "$OUT/${sample}.minus.scaled.bdg"

  bedGraphToBigWig "$OUT/${sample}.plus.scaled.bdg" \
    "$GENOME_DIR/genome.fa.fai" \
    "$OUT/${sample}.plus.scaled.bw"

  bedGraphToBigWig "$OUT/${sample}.minus.scaled.bdg" \
    "$GENOME_DIR/genome.fa.fai" \
    "$OUT/${sample}.minus.scaled.bw"

  # 11. Monomer BED/FASTA
  echo "[11/13] monomer fasta"

  awk -v lo="$LEN_MIN" -v hi="$LEN_MAX" 'BEGIN{OFS="\t"}{L=$3-$2; if(L>=lo && L<=hi) print $0}' \
    "$OUT/genome/${sample}.dedup.bed" \
  > "$OUT/monomer/${sample}.len${LEN_MIN}_${LEN_MAX}.bed"

  bedtools getfasta \
    -fi "$GENOME_DIR/genome.fa" \
    -bed "$OUT/monomer/${sample}.len${LEN_MIN}_${LEN_MAX}.bed" \
    -fo "$OUT/monomer/${sample}.len${LEN_MIN}_${LEN_MAX}.fa" \
    -s

  # 12. TS/NTS counts using length-filtered BED
  echo "[12/13] TS/NTS counts"

  bedtools intersect -c -a "$FEATURE_BED" \
    -b "$OUT/genome/${sample}.dedup.bed" \
    -wa -s -F 0.5 \
  > "$OUT/${sample}_NTScount.txt"

  bedtools intersect -c -a "$FEATURE_BED" \
    -b "$OUT/genome/${sample}.dedup.bed" \
    -wa -S -F 0.5 \
  > "$OUT/${sample}_TScount.txt"
  
  module purge
  module load python/3.12.4
  # 13. Monomer + figures
  echo "[13/13] monomer analysis and figures"

  python "$MONOMER_PY" \
    "$OUT/monomer/${sample}.len${LEN_MIN}_${LEN_MAX}.fa" \
    --length-min "$LEN_MIN" \
    --length-max "$LEN_MAX"

  module purge
  module load r
  
  Rscript "$FIGURES_R" "$sample" "$OUT" "$FIGOUT" "$LEN_MIN" "$LEN_MAX"

  module purge
  module load anaconda
  source "$(conda info --base)/etc/profile.d/conda.sh"
  conda activate umi_env
  module load bowtie2
  module load samtools
  module load bedtools
  module load ucsctools

  echo "DONE sample: $sample"

done < "$SAMPLE_LIST"

echo "=== BATCH DONE $(date) ==="