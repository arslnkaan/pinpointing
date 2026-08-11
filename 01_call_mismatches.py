#!/usr/bin/env python3

import argparse
from pathlib import Path
from typing import Dict, List, Tuple

import pysam
from Bio import SeqIO
from Bio.Seq import Seq


def reverse_complement(seq: str) -> str:
    return str(Seq(seq).reverse_complement()).upper()


def load_reference_fasta(fasta: str) -> Dict[str, str]:
    ref = {}
    for record in SeqIO.parse(fasta, "fasta"):
        ref[record.id] = str(record.seq).upper()
    return ref


def get_aligned_pairs_oriented(read, ref_sequences: Dict[str, str]) -> Tuple[str, str, List[int]]:
    """
    Build query/reference strings only across aligned match/mismatch positions.
    Positions are returned as genomic reference coordinates in the same order as the oriented sequence.

    Important behavior:
    - For forward reads: sequence is read 5'->3' as stored.
    - For reverse-strand reads: query and reference are reverse-complemented, so mismatch
      positions are measured in the read's own 5'->3' orientation.
    """
    chrom = read.reference_name
    if chrom not in ref_sequences:
        return "", "", []

    qseq_raw = read.query_sequence
    if qseq_raw is None:
        return "", "", []
    qseq_raw = qseq_raw.upper()

    query_bases = []
    ref_bases = []
    ref_positions = []

    # matches_only=False includes insertions/deletions. We skip indels here because
    # this pipeline focuses on base substitutions.
    for qpos, rpos in read.get_aligned_pairs(matches_only=False):
        if qpos is None or rpos is None:
            continue
        qb = qseq_raw[qpos].upper()
        rb = ref_sequences[chrom][rpos].upper()
        if qb in "ACGTN" and rb in "ACGTN":
            query_bases.append(qb)
            ref_bases.append(rb)
            ref_positions.append(rpos)

    query_seq = "".join(query_bases)
    ref_seq = "".join(ref_bases)

    if read.is_reverse:
        query_seq = reverse_complement(query_seq)
        ref_seq = reverse_complement(ref_seq)
        ref_positions = list(reversed(ref_positions))

    return query_seq, ref_seq, ref_positions


def main():
    parser = argparse.ArgumentParser(
        description="Call base mismatches from a BAM using read-oriented 5'->3' sequences."
    )
    parser.add_argument("--bam", required=True, help="Input sorted/indexed BAM")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA")
    parser.add_argument("--sample", required=True, help="Sample name prefix")
    parser.add_argument("--outdir", default="mismatch_outputs", help="Output directory")
    parser.add_argument("--mapq", type=int, default=0, help="Minimum MAPQ")
    parser.add_argument("--dedup-by-sequence", action="store_true",
                        help="Remove duplicate read sequences, considering reverse complement as duplicate too. Use only if BAM is not already deduplicated.")
    parser.add_argument("--include-perfect", action="store_true",
                        help="Include perfect-match reads in mismatch TSV. Default writes all reads anyway with Mismatches=None; kept for compatibility.")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mismatches_out = outdir / f"{args.sample}_mismatches.tsv"
    mismatch_bed_out = outdir / f"{args.sample}_mismatch_reads.bed"
    aligned_out = outdir / f"{args.sample}_aligned_sequences.tsv"
    summary_out = outdir / f"{args.sample}_mismatch_calling_summary.tsv"

    ref_sequences = load_reference_fasta(args.fasta)
    bam = pysam.AlignmentFile(args.bam, "rb")

    seen_sequences = set()
    n_total = 0
    n_unmapped = 0
    n_low_mapq = 0
    n_no_ref = 0
    n_duplicate_seq = 0
    n_written = 0
    n_with_mismatch = 0

    with open(mismatches_out, "w") as mf, open(mismatch_bed_out, "w") as bf, open(aligned_out, "w") as af:
        mf.write("Read_ID\tChromosome\tStart\tEnd\tStrand\tRead_Length\tAligned_Sequence\tReference_Sequence\tMismatches\tMismatch_Count\tGenomic_Positions\n")
        bf.write("chrom\tstart\tend\tread_id\tscore\tstrand\tsequence\tmismatch_count\n")
        af.write("Read_ID\tChromosome\tStart\tEnd\tStrand\tAligned_Length\tAligned_Sequence\n")

        for read in bam.fetch(until_eof=True):
            n_total += 1

            if read.is_unmapped:
                n_unmapped += 1
                continue
            if read.mapping_quality < args.mapq:
                n_low_mapq += 1
                continue

            chrom = read.reference_name
            if chrom not in ref_sequences:
                n_no_ref += 1
                continue

            strand = "-" if read.is_reverse else "+"
            start = read.reference_start
            end = read.reference_end

            aligned_seq, ref_seq, ref_positions = get_aligned_pairs_oriented(read, ref_sequences)
            read_length = len(aligned_seq)
            if read_length == 0:
                continue

            if args.dedup_by_sequence:
                rc_aligned = reverse_complement(aligned_seq)
                if aligned_seq in seen_sequences or rc_aligned in seen_sequences:
                    n_duplicate_seq += 1
                    continue
                seen_sequences.add(aligned_seq)
                seen_sequences.add(rc_aligned)

            mismatches = []
            mismatch_genomic_positions = []
            for pos_1based, (qb, rb, gpos) in enumerate(zip(aligned_seq, ref_seq, ref_positions), start=1):
                if qb != rb and qb in "ACGT" and rb in "ACGT":
                    mismatches.append(f"{pos_1based}:{rb}>{qb}")
                    mismatch_genomic_positions.append(str(gpos))  # 0-based genomic coordinate of ref base

            mismatch_info = ",".join(mismatches) if mismatches else "None"
            genomic_pos_info = ",".join(mismatch_genomic_positions) if mismatch_genomic_positions else "None"
            mismatch_count = len(mismatches)

            mf.write(
                f"{read.query_name}\t{chrom}\t{start}\t{end}\t{strand}\t{read_length}\t"
                f"{aligned_seq}\t{ref_seq}\t{mismatch_info}\t{mismatch_count}\t{genomic_pos_info}\n"
            )
            af.write(f"{read.query_name}\t{chrom}\t{start}\t{end}\t{strand}\t{read_length}\t{aligned_seq}\n")
            n_written += 1

            if mismatch_count > 0:
                n_with_mismatch += 1
                bf.write(f"{chrom}\t{start}\t{end}\t{read.query_name}\t.\t{strand}\t{aligned_seq}\t{mismatch_count}\n")

    bam.close()

    with open(summary_out, "w") as sf:
        sf.write("metric\tvalue\n")
        for key, value in [
            ("total_bam_records", n_total),
            ("unmapped_skipped", n_unmapped),
            ("low_mapq_skipped", n_low_mapq),
            ("missing_reference_skipped", n_no_ref),
            ("duplicate_sequence_skipped", n_duplicate_seq),
            ("reads_written", n_written),
            ("reads_with_mismatch", n_with_mismatch),
        ]:
            sf.write(f"{key}\t{value}\n")

    print(f"Done: {args.sample}")
    print(f"Mismatch table: {mismatches_out}")
    print(f"Mismatch-read BED: {mismatch_bed_out}")
    print(f"Summary: {summary_out}")


if __name__ == "__main__":
    main()
