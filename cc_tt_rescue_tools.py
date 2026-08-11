#!/usr/bin/env python3

import argparse
import gzip
import pysam


def open_maybe_gz(path, mode):
    if path.endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def generate_variants(in_fastq, out_fastq, map_tsv):
    """
    For each original read, create one variant per overlapping TT occurrence.
    Only that specific TT is changed to CC.

    Example:
    ACTGGATCTGCATTTTCGT

    variants:
    ACTGGATCTGCACCTTCGT
    ACTGGATCTGCATCCTCGT
    ACTGGATCTGCATTCCCGT
    """

    with open_maybe_gz(in_fastq, "rt") as fin, \
         open_maybe_gz(out_fastq, "wt") as fout, \
         open(map_tsv, "w") as m:

        m.write(
            "variant_id\toriginal_read_id\ttt_start_0based\ttt_end_0based\t"
            "tt_start_1based\toriginal_seq\tvariant_seq\n"
        )

        while True:
            name = fin.readline()
            if not name:
                break

            seq = fin.readline().rstrip("\n")
            plus = fin.readline()
            qual = fin.readline().rstrip("\n")

            original_id = name.strip().split()[0].lstrip("@")
            seq_upper = seq.upper()

            for i in range(len(seq_upper) - 1):
                if seq_upper[i:i + 2] != "TT":
                    continue

                variant_seq = seq[:i] + "CC" + seq[i + 2:]
                variant_id = f"{original_id}|TTpos0={i}"

                fout.write(f"@{variant_id}\n")
                fout.write(variant_seq + "\n")
                fout.write("+\n")
                fout.write(qual + "\n")

                m.write(
                    f"{variant_id}\t{original_id}\t{i}\t{i + 2}\t{i + 1}\t"
                    f"{seq}\t{variant_seq}\n"
                )


def load_variant_map(map_tsv):
    d = {}

    with open(map_tsv, "r") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {x: i for i, x in enumerate(header)}

        for line in f:
            fields = line.rstrip("\n").split("\t")

            variant_id = fields[idx["variant_id"]]
            d[variant_id] = {
                "original_read_id": fields[idx["original_read_id"]],
                "tt_start_0based": int(fields[idx["tt_start_0based"]]),
                "tt_end_0based": int(fields[idx["tt_end_0based"]]),
                "tt_start_1based": int(fields[idx["tt_start_1based"]]),
                "original_seq": fields[idx["original_seq"]],
                "variant_seq": fields[idx["variant_seq"]],
            }

    return d


def revcomp(seq):
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]


def call_events(variant_map_tsv, bam_file, ref_fasta, out_tsv, out_bed):
    variant_info = load_variant_map(variant_map_tsv)

    bam = pysam.AlignmentFile(bam_file, "rb")
    ref = pysam.FastaFile(ref_fasta)

    with open(out_tsv, "w") as tsv, open(out_bed, "w") as bed:
        tsv.write(
            "variant_id\toriginal_read_id\tchrom\tstart_0based\tend_0based\t"
            "event_strand\tread_strand\tread_TT_start_1based\tread_TT_end_1based\t"
            "original_read_dinuc\tvariant_read_dinuc\tgenome_dinuc\t"
            "genome_trinuc_leftC\tgenome_trinuc_rightC\t"
            "genome_tetranuc_leftC_plus1\tgenome_tetranuc_minus1_rightC\t"
            "event\n"
        )

        seen = set()

        for read in bam.fetch(until_eof=True):
            if read.is_unmapped:
                continue

            variant_id = read.query_name

            if variant_id not in variant_info:
                continue

            info = variant_info[variant_id]
            q0 = info["tt_start_0based"]
            q1 = q0 + 1

            q_to_r = {
                q: r
                for q, r in read.get_aligned_pairs(matches_only=True)
                if q is not None and r is not None
            }

            if q0 not in q_to_r or q1 not in q_to_r:
                continue

            r0 = q_to_r[q0]
            r1 = q_to_r[q1]

            if abs(r1 - r0) != 1:
                continue

            chrom = read.reference_name
            start = min(r0, r1)
            end = max(r0, r1) + 1

            genome_dinuc = ref.fetch(chrom, start, end).upper()
            read_strand = "-" if read.is_reverse else "+"

            if genome_dinuc == "CC":
                event_strand = "+"
            elif genome_dinuc == "GG":
                event_strand = "-"
            else:
                continue

            key = (variant_id, chrom, start, end, event_strand)
            if key in seen:
                continue
            seen.add(key)

            original_dinuc = info["original_seq"][q0:q0 + 2].upper()
            variant_dinuc = info["variant_seq"][q0:q0 + 2].upper()

            if original_dinuc != "TT":
                continue

            if variant_dinuc != "CC":
                continue

            # Contexts are reported in the event/cytosine-strand orientation.
            # For + strand: genome has CC.
            # For - strand: genome has GG, so reverse-complement the local window.
            chrom_len = ref.get_reference_length(chrom)

            tri_left = "NA"
            tri_right = "NA"
            tetra_left_plus1 = "NA"
            tetra_minus1_right = "NA"

            if event_strand == "+":
                if start - 1 >= 0:
                    tri_left = ref.fetch(chrom, start - 1, start + 2).upper()
                if end + 1 <= chrom_len:
                    tri_right = ref.fetch(chrom, start, end + 1).upper()
                if start - 1 >= 0 and end + 1 <= chrom_len:
                    tetra_left_plus1 = ref.fetch(chrom, start - 1, end + 1).upper()
                if start - 2 >= 0:
                    tetra_minus1_right = ref.fetch(chrom, start - 2, end).upper()

            else:
                # Reference genome has GG. Convert local plus-genome context
                # to the minus-strand CC-oriented context.
                if end + 1 <= chrom_len:
                    tri_left = revcomp(ref.fetch(chrom, start, end + 1).upper())
                if start - 1 >= 0:
                    tri_right = revcomp(ref.fetch(chrom, start - 1, end).upper())
                if start - 1 >= 0 and end + 1 <= chrom_len:
                    tetra_left_plus1 = revcomp(ref.fetch(chrom, start - 1, end + 1).upper())
                if end + 2 <= chrom_len:
                    tetra_minus1_right = revcomp(ref.fetch(chrom, start, end + 2).upper())

            tsv.write(
                f"{variant_id}\t"
                f"{info['original_read_id']}\t"
                f"{chrom}\t{start}\t{end}\t"
                f"{event_strand}\t{read_strand}\t"
                f"{q0 + 1}\t{q1 + 1}\t"
                f"{original_dinuc}\t{variant_dinuc}\t{genome_dinuc}\t"
                f"{tri_left}\t{tri_right}\t"
                f"{tetra_left_plus1}\t{tetra_minus1_right}\t"
                f"CC>TT\n"
            )

            bed.write(
                f"{chrom}\t{start}\t{end}\t{info['original_read_id']}|TTpos0={q0}\t0\t{event_strand}\n"
            )


def main():
    parser = argparse.ArgumentParser()

    sub = parser.add_subparsers(dest="command", required=True)

    p1 = sub.add_parser("make-variants")
    p1.add_argument("--in-fastq", required=True)
    p1.add_argument("--out-fastq", required=True)
    p1.add_argument("--map-tsv", required=True)

    p2 = sub.add_parser("call-events")
    p2.add_argument("--variant-map", required=True)
    p2.add_argument("--bam", required=True)
    p2.add_argument("--ref-fasta", required=True)
    p2.add_argument("--out-tsv", required=True)
    p2.add_argument("--out-bed", required=True)

    args = parser.parse_args()

    if args.command == "make-variants":
        generate_variants(args.in_fastq, args.out_fastq, args.map_tsv)

    elif args.command == "call-events":
        call_events(
            args.variant_map,
            args.bam,
            args.ref_fasta,
            args.out_tsv,
            args.out_bed
        )


if __name__ == "__main__":
    main()