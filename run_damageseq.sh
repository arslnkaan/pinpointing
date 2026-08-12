#!/bin/bash
#SBATCH --job-name=damageseq
#SBATCH --output=damageseq_%j.out
#SBATCH --error=damageseq_%j.err
#SBATCH --time=48:00:00
#SBATCH --cpus-per-task=12
#SBATCH --mem=64G

set -euo pipefail

module load cutadapt
module load bowtie2
module load samtools
module load bedtools
module load ucsctools

sample="$1"
R1="$2"
R2="$3"

THREADS=12
ADAPTER="GACTGGTTCCAATTGAAAGTGCTCTTCCGATCT"

BOWTIE2_INDEX="/proj/seq/data/hg38_UCSC/Sequence/Bowtie2Index/genome"
GENOME_FA="/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"
CHROMSIZES="/work/users/a/r/arslank/genome/hg38.chrom.sizes"

OUTDIR="${sample}-results"
mkdir -p ${OUTDIR}/{trimmed,bam,dedup,damage_sites,bigwig,qc,tmp}

echo "Processing ${sample}"
echo "R1: ${R1}"
echo "R2: ${R2}"
echo "Threads: ${THREADS}"

# ============================================================
# 1. DISCARD reads containing adapter
# ============================================================

# If trimmed FASTQs already exist, skip cutadapt
if [[ ! -s ${OUTDIR}/trimmed/${sample}_R1.trimmed.fastq || ! -s ${OUTDIR}/trimmed/${sample}_R2.trimmed.fastq ]]; then

  cutadapt \
    -b ${ADAPTER} \
    -B ${ADAPTER} \
    --discard-trimmed \
    -q 20 \
    -m 20 \
    -o ${OUTDIR}/trimmed/${sample}_R1.trimmed.fastq \
    -p ${OUTDIR}/trimmed/${sample}_R2.trimmed.fastq \
    ${R1} ${R2} \
    > ${OUTDIR}/qc/${sample}.cutadapt.txt

else
  echo "Trimmed FASTQs already exist. Skipping cutadapt."
fi

# ============================================================
# 2. FASTQ sanity check
# ============================================================

awk 'END{if(NR%4!=0){print "ERROR: R1 FASTQ broken"; exit 1}else{print "R1 records:", NR/4}}' \
  ${OUTDIR}/trimmed/${sample}_R1.trimmed.fastq

awk 'END{if(NR%4!=0){print "ERROR: R2 FASTQ broken"; exit 1}else{print "R2 records:", NR/4}}' \
  ${OUTDIR}/trimmed/${sample}_R2.trimmed.fastq

# ============================================================
# 3A. Bowtie2 alignment -> UNSORTED BAM CHECKPOINT
# ============================================================

if [[ ! -s ${OUTDIR}/bam/${sample}.unsorted.bam ]]; then

  echo "Running Bowtie2 alignment..."

  bowtie2 \
    -q \
    --phred33 \
    --local \
    -p ${THREADS} \
    --seed 123 \
    --no-mixed \
    -x ${BOWTIE2_INDEX} \
    -1 ${OUTDIR}/trimmed/${sample}_R1.trimmed.fastq \
    -2 ${OUTDIR}/trimmed/${sample}_R2.trimmed.fastq \
    2> ${OUTDIR}/qc/${sample}.bowtie2.log \
  | samtools view \
      -@ ${THREADS} \
      -bS - \
      > ${OUTDIR}/bam/${sample}.unsorted.bam

  samtools quickcheck -v ${OUTDIR}/bam/${sample}.unsorted.bam

else
  echo "Unsorted BAM already exists. Skipping Bowtie2."
fi

# ============================================================
# 3B. Coordinate sort BAM
# ============================================================

if [[ ! -s ${OUTDIR}/bam/${sample}.sorted.bam ]]; then

  echo "Sorting BAM..."

  samtools sort \
    -@ ${THREADS} \
    -m 4G \
    -T ${OUTDIR}/tmp/${sample}.coord_sort_tmp \
    -o ${OUTDIR}/bam/${sample}.sorted.bam \
    ${OUTDIR}/bam/${sample}.unsorted.bam

  samtools quickcheck -v ${OUTDIR}/bam/${sample}.sorted.bam
  samtools index ${OUTDIR}/bam/${sample}.sorted.bam

else
  echo "Sorted BAM already exists. Skipping coordinate sort."
fi

# ============================================================
# 4. Deduplicate paired fragments
# ============================================================

if [[ ! -s ${OUTDIR}/dedup/${sample}.dedup.bam ]]; then

  echo "Name sorting BAM..."

  samtools sort -n \
    -@ ${THREADS} \
    -m 4G \
    -T ${OUTDIR}/tmp/${sample}.name_sort_tmp \
    -o ${OUTDIR}/bam/${sample}.namesort.bam \
    ${OUTDIR}/bam/${sample}.sorted.bam

  echo "Running fixmate..."

  samtools fixmate -m \
    ${OUTDIR}/bam/${sample}.namesort.bam \
    ${OUTDIR}/bam/${sample}.fixmate.bam

  echo "Position sorting fixmate BAM..."

  samtools sort \
    -@ ${THREADS} \
    -m 4G \
    -T ${OUTDIR}/tmp/${sample}.pos_sort_tmp \
    -o ${OUTDIR}/bam/${sample}.positionsort.bam \
    ${OUTDIR}/bam/${sample}.fixmate.bam

  echo "Running markdup..."

  samtools markdup -r \
    ${OUTDIR}/bam/${sample}.positionsort.bam \
    ${OUTDIR}/dedup/${sample}.dedup.bam

  samtools quickcheck -v ${OUTDIR}/dedup/${sample}.dedup.bam
  samtools index ${OUTDIR}/dedup/${sample}.dedup.bam

  samtools view -c -F 4 ${OUTDIR}/dedup/${sample}.dedup.bam \
    > ${OUTDIR}/qc/${sample}.dedup_mapped_reads.txt

  rm -f ${OUTDIR}/bam/${sample}.fixmate.bam
  rm -f ${OUTDIR}/bam/${sample}.positionsort.bam

else
  echo "Deduplicated BAM already exists. Skipping deduplication."
fi

# ============================================================
# 5. Convert deduplicated paired-end BAM to BEDPE
# ============================================================

if [[ ! -s ${OUTDIR}/damage_sites/${sample}.dedup.bedpe ]]; then

  bedtools bamtobed \
    -bedpe \
    -i ${OUTDIR}/dedup/${sample}.dedup.bam \
    > ${OUTDIR}/damage_sites/${sample}.dedup.bedpe

else
  echo "BEDPE already exists. Skipping BEDPE conversion."
fi

# ============================================================
# 6. Infer predicted 2-nt damage sites
# ============================================================

if [[ ! -s ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed ]]; then

  python infer_damageseq_sites.py \
    --bedpe ${OUTDIR}/damage_sites/${sample}.dedup.bedpe \
    --out ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed

else
  echo "Predicted damage sites already exist. Skipping inference."
fi

# ============================================================
# 7. Extract sequence at predicted damage sites
# ============================================================

if [[ ! -s ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.sequence.tsv ]]; then

  bedtools getfasta \
    -fi ${GENOME_FA} \
    -bed ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed \
    -s \
    -tab \
    -name \
    > ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.sequence.tsv

else
  echo "Damage-site sequence TSV already exists. Skipping getfasta."
fi

# ============================================================
# 8. Classify damage sequence context
# ============================================================

if [[ ! -s ${OUTDIR}/damage_sites/${sample}.classified_sites.tsv ]]; then

  python classify_damageseq_context.py \
    --sites ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed \
    --seq ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.sequence.tsv \
    --out-prefix ${OUTDIR}/damage_sites/${sample}

else
  echo "Classified damage sites already exist. Skipping classification."
fi

# ============================================================
# 9. Make 25-nt normalized bigWig
# ============================================================

if [[ ! -s ${OUTDIR}/bigwig/hg38.25nt.windows.bed ]]; then

  bedtools makewindows \
    -g ${CHROMSIZES} \
    -w 25 \
    > ${OUTDIR}/bigwig/hg38.25nt.windows.bed

fi

bedtools intersect \
  -a ${OUTDIR}/bigwig/hg38.25nt.windows.bed \
  -b ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed \
  -c \
  > ${OUTDIR}/bigwig/${sample}.25nt.raw_counts.bed

python normalize_25nt_by_chromosome.py \
  --counts ${OUTDIR}/bigwig/${sample}.25nt.raw_counts.bed \
  --sites ${OUTDIR}/damage_sites/${sample}.predicted_damage_sites.bed \
  --out ${OUTDIR}/bigwig/${sample}.25nt.chrom_norm.bedGraph

LC_COLLATE=C sort -k1,1 -k2,2n \
  ${OUTDIR}/bigwig/${sample}.25nt.chrom_norm.bedGraph \
  > ${OUTDIR}/bigwig/${sample}.25nt.chrom_norm.sorted.bedGraph

bedGraphToBigWig \
  ${OUTDIR}/bigwig/${sample}.25nt.chrom_norm.sorted.bedGraph \
  ${CHROMSIZES} \
  ${OUTDIR}/bigwig/${sample}.25nt.chrom_norm.bw

echo "Done: ${sample}"