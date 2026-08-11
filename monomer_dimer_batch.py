#!/usr/bin/env python3

import argparse
from collections import defaultdict
from pathlib import Path
from Bio import SeqIO


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("fasta", help="Input FASTA from bedtools getfasta")
    parser.add_argument("--length-min", type=int, default=24)
    parser.add_argument("--length-max", type=int, default=30)
    args = parser.parse_args()

    fasta = Path(args.fasta)
    desired_lengths = list(range(args.length_min, args.length_max + 1))

    mono = {
        length: {
            pos: {"A": 0, "T": 0, "G": 0, "C": 0}
            for pos in range(1, length + 1)
        }
        for length in desired_lengths
    }

    dinuc = {
        length: {
            pos: defaultdict(int)
            for pos in range(1, length)
        }
        for length in desired_lengths
    }

    n_by_length = defaultdict(int)

    for entry in SeqIO.parse(str(fasta), "fasta"):
        seq = str(entry.seq).upper()

        if "N" in seq:
            continue

        length = len(seq)

        if length not in desired_lengths:
            continue

        n_by_length[length] += 1

        for i, nt in enumerate(seq, start=1):
            if nt in mono[length][i]:
                mono[length][i][nt] += 1

        for i in range(length - 1):
            d = seq[i:i+2]
            if set(d).issubset({"A", "T", "G", "C"}):
                dinuc[length][i + 1][d] += 1

    prefix = str(fasta).replace(".fa", "")

    mono_out = prefix + "_monomer_R_df.txt"
    with open(mono_out, "w") as out:
        out.write("Length\tPosition\tBase\tFrequency\n")
        for length in desired_lengths:
            for pos in range(1, length + 1):
                total = sum(mono[length][pos].values())
                for base in ["A", "T", "G", "C"]:
                    val = mono[length][pos][base] / total if total else 0
                    out.write(f"{length}\t{pos}\t{base}\t{val}\n")

    dinuc_out = prefix + "_dinucleotide_R_df.txt"
    all_dinucs = [
        "AA", "AT", "AG", "AC",
        "TA", "TT", "TG", "TC",
        "GA", "GT", "GG", "GC",
        "CA", "CT", "CG", "CC"
    ]

    with open(dinuc_out, "w") as out:
        out.write("Length\tPosition\tDinucleotide\tFrequency\n")
        for length in desired_lengths:
            for pos in range(1, length):
                total = sum(dinuc[length][pos].values())
                for d in all_dinucs:
                    val = dinuc[length][pos][d] / total if total else 0
                    out.write(f"{length}\t{pos}\t{d}\t{val}\n")

    summary_out = prefix + "_sequence_counts_by_length.txt"
    with open(summary_out, "w") as out:
        out.write("Length\tSequences\n")
        for length in desired_lengths:
            out.write(f"{length}\t{n_by_length[length]}\n")


if __name__ == "__main__":
    main()