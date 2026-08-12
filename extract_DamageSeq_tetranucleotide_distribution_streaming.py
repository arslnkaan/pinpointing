#!/usr/bin/env python3

import argparse
import csv
import subprocess
import sys
from collections import Counter
from pathlib import Path


VALID_BASES = set("ACGT")
VALID_DIPYRIMIDINES = {"CC", "CT", "TC", "TT"}


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Stream Damage-seq BED events, extract positions "
            "-3,-2,-1,0, and calculate N[CPD]N distributions "
            "without loading the full dataset into memory."
        )
    )

    parser.add_argument(
        "--input-bed",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--reference",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--strand-column",
        type=int,
        default=6,
        help="One-based strand column.",
    )

    parser.add_argument(
        "--dinucleotide-column",
        type=int,
        default=7,
        help=(
            "One-based dipyrimidine column. "
            "Use 0 to infer from extracted sequence."
        ),
    )

    return parser.parse_args()


def normalize_strand(value):

    value = str(value).strip().lower()

    if value in {
        "+",
        "plus",
        "forward",
        "f",
        "1",
        "+1",
    }:
        return "+"

    if value in {
        "-",
        "minus",
        "reverse",
        "r",
        "-1",
    }:
        return "-"

    return None


def normalize_dinucleotide(value):

    value = str(value).strip().upper().replace(
        "U",
        "T",
    )

    if value in VALID_DIPYRIMIDINES:
        return value

    return None


def read_chromosome_sizes(reference):

    fai = Path(
        f"{reference}.fai"
    )

    if not fai.exists():

        subprocess.run(
            [
                "samtools",
                "faidx",
                str(reference),
            ],
            check=True,
        )

    chromosome_sizes = {}

    with fai.open() as handle:

        for line in handle:

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) >= 2:

                chromosome_sizes[
                    fields[0]
                ] = int(fields[1])

    return chromosome_sizes


def prepare_streaming_windows(
    input_bed,
    output_bed,
    chromosome_sizes,
    strand_column,
    dinucleotide_column,
):

    strand_index = strand_column - 1

    dinucleotide_index = None

    if dinucleotide_column > 0:

        dinucleotide_index = (
            dinucleotide_column - 1
        )

    qc = Counter()
    selected_index = 0

    with (
        input_bed.open() as input_handle,
        output_bed.open("w") as output_handle,
    ):

        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for line_number, line in enumerate(
            input_handle,
            start=1,
        ):

            if (
                not line.strip()
                or line.startswith("#")
            ):
                continue

            qc["input_rows"] += 1

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) < 3:

                qc[
                    "invalid_column_count"
                ] += 1

                continue

            try:

                chrom = fields[0].strip()
                damage_start = int(fields[1])
                damage_end = int(fields[2])

            except ValueError:

                qc[
                    "header_or_invalid_coordinates"
                ] += 1

                continue

            if damage_end - damage_start != 2:

                qc[
                    "not_two_base_damage_interval"
                ] += 1

                continue

            if len(fields) <= strand_index:

                qc[
                    "missing_strand_column"
                ] += 1

                continue

            strand = normalize_strand(
                fields[strand_index]
            )

            if strand is None:

                qc["invalid_strand"] += 1
                continue

            reported_dinucleotide = "."

            if dinucleotide_index is not None:

                if len(fields) <= dinucleotide_index:

                    qc[
                        "missing_dinucleotide_column"
                    ] += 1

                    continue

                reported_dinucleotide = (
                    normalize_dinucleotide(
                        fields[
                            dinucleotide_index
                        ]
                    )
                )

                if reported_dinucleotide is None:

                    qc[
                        "invalid_reported_dinucleotide"
                    ] += 1

                    continue

            if chrom not in chromosome_sizes:

                qc[
                    "chromosome_not_in_reference"
                ] += 1

                continue

            # Original interval contains positions -2 and -1.
            # Add one base on each side to obtain:
            #
            # -3, -2, -1, 0
            #
            # bedtools getfasta -s reverse-complements
            # minus-strand intervals.

            window_start = damage_start - 1
            window_end = damage_end + 1

            if (
                window_start < 0
                or window_end >
                chromosome_sizes[chrom]
            ):

                qc[
                    "window_outside_chromosome"
                ] += 1

                continue

            selected_index += 1

            event_id = (
                f"damage_{selected_index:012d}"
            )

            # BED10:
            # 1 chrom
            # 2 window start
            # 3 window end
            # 4 event ID
            # 5 score
            # 6 strand
            # 7 original damage start
            # 8 original damage end
            # 9 reported dipyrimidine
            # 10 original line number

            writer.writerow(
                [
                    chrom,
                    window_start,
                    window_end,
                    event_id,
                    0,
                    strand,
                    damage_start,
                    damage_end,
                    reported_dinucleotide,
                    line_number,
                ]
            )

            qc["selected_windows"] += 1

    return qc


def extract_and_count(
    reference,
    window_bed,
    sequence_output,
    event_output,
):

    command = [
        "bedtools",
        "getfasta",
        "-fi",
        str(reference),
        "-bed",
        str(window_bed),
        "-s",
        "-bedOut",
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1024 * 1024,
    )

    if process.stdout is None:

        raise RuntimeError(
            "Could not read bedtools output."
        )

    counts = Counter()
    qc = Counter()

    with (
        sequence_output.open("w") as sequence_handle,
        event_output.open("w") as event_handle,
    ):

        sequence_writer = csv.writer(
            sequence_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        event_writer = csv.writer(
            event_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        sequence_writer.writerow(
            [
                "chrom",
                "window_start",
                "window_end",
                "event_id",
                "score",
                "strand",
                "damage_start",
                "damage_end",
                "reported_dinucleotide",
                "input_line_number",
                "sequence",
            ]
        )

        event_writer.writerow(
            [
                "event_id",
                "chrom",
                "damage_start",
                "damage_end",
                "strand",
                "position_minus3",
                "position_minus2",
                "position_minus1",
                "position_0",
                "tetranucleotide",
                "dipyrimidine",
                "reported_dinucleotide",
            ]
        )

        for output_line_number, line in enumerate(
            process.stdout,
            start=1,
        ):

            fields = line.rstrip(
                "\n"
            ).split("\t")

            qc["bedtools_rows"] += 1

            # Original BED10 plus sequence.
            if len(fields) != 11:

                process.kill()

                raise RuntimeError(
                    "Unexpected bedtools output at line "
                    f"{output_line_number}: "
                    f"{len(fields)} columns"
                )

            sequence_writer.writerow(
                fields
            )

            (
                chrom,
                window_start,
                window_end,
                event_id,
                score,
                strand,
                damage_start,
                damage_end,
                reported_dinucleotide,
                input_line_number,
                sequence,
            ) = fields

            sequence = sequence.upper()

            if len(sequence) != 4:

                qc[
                    "sequence_not_four_bases"
                ] += 1

                continue

            if not set(sequence) <= VALID_BASES:

                qc[
                    "non_acgt_sequence"
                ] += 1

                continue

            extracted_dinucleotide = (
                sequence[1:3]
            )

            if (
                extracted_dinucleotide
                not in VALID_DIPYRIMIDINES
            ):

                qc[
                    "central_bases_not_dipyrimidine"
                ] += 1

                continue

            if (
                reported_dinucleotide != "."
                and reported_dinucleotide
                != extracted_dinucleotide
            ):

                qc[
                    "reported_dinucleotide_mismatch"
                ] += 1

                continue

            if reported_dinucleotide != ".":

                qc[
                    "reported_dinucleotide_matches"
                ] += 1

            counts[sequence] += 1

            event_writer.writerow(
                [
                    event_id,
                    chrom,
                    damage_start,
                    damage_end,
                    strand,
                    sequence[0],
                    sequence[1],
                    sequence[2],
                    sequence[3],
                    sequence,
                    extracted_dinucleotide,
                    reported_dinucleotide,
                ]
            )

            qc[
                "valid_tetranucleotides"
            ] += 1

    stderr_text = ""

    if process.stderr is not None:

        stderr_text = process.stderr.read()

    return_code = process.wait()

    if return_code != 0:

        raise RuntimeError(
            "bedtools getfasta failed:\n"
            + stderr_text
        )

    return counts, qc


def tetranucleotide_order():

    return [
        left + dimer + right
        for dimer in (
            "CC",
            "CT",
            "TC",
            "TT",
        )
        for left in "ACGT"
        for right in "ACGT"
    ]


def write_distribution(
    counts,
    output_path,
):

    total = sum(
        counts.values()
    )

    with output_path.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "tetranucleotide",
                "position_minus3",
                "position_minus2",
                "position_minus1",
                "position_0",
                "dipyrimidine",
                "count",
                "damage_total",
                "damage_percentage",
                "percentage_within_dipyrimidine",
            ]
        )

        dimer_totals = Counter()

        for tetranucleotide, count in (
            counts.items()
        ):

            dimer_totals[
                tetranucleotide[1:3]
            ] += count

        for tetranucleotide in (
            tetranucleotide_order()
        ):

            count = counts[
                tetranucleotide
            ]

            dimer = tetranucleotide[
                1:3
            ]

            damage_percentage = (
                100 * count / total
                if total > 0
                else 0
            )

            percentage_within_dimer = (
                100
                * count
                / dimer_totals[dimer]
                if dimer_totals[dimer] > 0
                else 0
            )

            writer.writerow(
                [
                    tetranucleotide,
                    tetranucleotide[0],
                    tetranucleotide[1],
                    tetranucleotide[2],
                    tetranucleotide[3],
                    dimer,
                    count,
                    total,
                    damage_percentage,
                    percentage_within_dimer,
                ]
            )


def write_qc(
    preparation_qc,
    extraction_qc,
    output_path,
):

    with output_path.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "metric",
                "value",
            ]
        )

        for metric, value in sorted(
            preparation_qc.items()
        ):

            writer.writerow(
                [
                    f"preparation_{metric}",
                    value,
                ]
            )

        for metric, value in sorted(
            extraction_qc.items()
        ):

            writer.writerow(
                [
                    f"extraction_{metric}",
                    value,
                ]
            )


def main():

    args = parse_args()

    for input_path in (
        args.input_bed,
        args.reference,
    ):

        if not input_path.exists():

            raise FileNotFoundError(
                f"Missing input: {input_path}"
            )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    chromosome_sizes = (
        read_chromosome_sizes(
            args.reference
        )
    )

    window_bed = (
        args.outdir
        / "DamageSeq.minus3_to_0.windows.hg38.bed10"
    )

    sequence_output = (
        args.outdir
        / "DamageSeq.minus3_to_0.sequences.tsv"
    )

    event_output = (
        args.outdir
        / "DamageSeq.tetranucleotide.events.tsv"
    )

    distribution_output = (
        args.outdir
        / "DamageSeq.tetranucleotide.distribution.tsv"
    )

    qc_output = (
        args.outdir
        / "DamageSeq.tetranucleotide.QC.tsv"
    )

    preparation_qc = prepare_streaming_windows(
        args.input_bed,
        window_bed,
        chromosome_sizes,
        args.strand_column,
        args.dinucleotide_column,
    )

    counts, extraction_qc = (
        extract_and_count(
            args.reference,
            window_bed,
            sequence_output,
            event_output,
        )
    )

    write_distribution(
        counts,
        distribution_output,
    )

    write_qc(
        preparation_qc,
        extraction_qc,
        qc_output,
    )

    print()
    print(
        "Damage-seq tetranucleotide extraction complete"
    )
    print(
        "=============================================="
    )
    print(
        f"Distribution: {distribution_output}"
    )
    print(
        f"Event table: {event_output}"
    )
    print(
        f"QC: {qc_output}"
    )


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)