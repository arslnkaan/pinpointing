#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path


VALID_BASES = set("ACGT")
VALID_DIPYRIMIDINES = {
    "CC",
    "CT",
    "TC",
    "TT",
}

TIMEPOINT_MAP = {
    "30m": "0.5h",
    "30min": "0.5h",
    "0.5h": "0.5h",
    "2h": "2h",
    "4h": "4h",
    "8h": "8h",
}

TIMEPOINT_ORDER = (
    "0.5h",
    "2h",
    "4h",
    "8h",
)


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Extract strand-oriented N[CPD]N tetranucleotides "
            "from a BED file of 2-bp dipyrimidine intervals."
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
        "--timepoint-column",
        type=int,
        default=0,
        help=(
            "One-based BED column containing timepoint. "
            "Use 0 for automatic inference from all columns."
        ),
    )

    parser.add_argument(
        "--strand-column",
        type=int,
        default=6,
        help="One-based strand column. Default: 6.",
    )

    parser.add_argument(
        "--dinucleotide-column",
        type=int,
        default=0,
        help=(
            "Optional one-based column containing the predicted "
            "dipyrimidine. Use 0 to infer it from extracted sequence."
        ),
    )

    return parser.parse_args()


def normalize_timepoint(value):

    value = str(value).strip().lower()
    value = value.replace(" ", "")

    if value in TIMEPOINT_MAP:
        return TIMEPOINT_MAP[value]

    if "30m" in value:
        return "0.5h"

    for timepoint in (
        "2h",
        "4h",
        "8h",
    ):
        if timepoint in value:
            return timepoint

    return None


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

    match = re.fullmatch(
        r"(?:CPD[=:])?(CC|CT|TC|TT)",
        value,
    )

    if match is not None:
        return match.group(1)

    return None


def infer_timepoint(fields):

    for value in fields:
        timepoint = normalize_timepoint(value)

        if timepoint is not None:
            return timepoint

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

            chromosome_sizes[
                fields[0]
            ] = int(fields[1])

    return chromosome_sizes


def prepare_windows(
    input_bed,
    output_bed,
    metadata_file,
    chromosome_sizes,
    timepoint_column,
    strand_column,
    dinucleotide_column,
):

    qc = Counter()

    strand_index = strand_column - 1

    timepoint_index = (
        timepoint_column - 1
        if timepoint_column > 0
        else None
    )

    dinucleotide_index = (
        dinucleotide_column - 1
        if dinucleotide_column > 0
        else None
    )

    with (
        input_bed.open() as input_handle,
        output_bed.open("w") as bed_handle,
        metadata_file.open("w") as metadata_handle,
    ):

        bed_writer = csv.writer(
            bed_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        metadata_columns = [
            "event_id",
            "timepoint",
            "chrom",
            "cpd_start",
            "cpd_end",
            "strand",
            "reported_dinucleotide",
        ]

        metadata_writer = csv.DictWriter(
            metadata_handle,
            fieldnames=metadata_columns,
            delimiter="\t",
            lineterminator="\n",
        )

        metadata_writer.writeheader()

        selected_index = 0

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

                qc["invalid_column_count"] += 1
                continue

            try:

                chrom = fields[0]
                start = int(fields[1])
                end = int(fields[2])

            except ValueError:

                # Allows an optional header.
                qc["header_or_invalid_coordinates"] += 1
                continue

            if end - start != 2:

                qc["not_two_base_interval"] += 1
                continue

            if len(fields) <= strand_index:

                qc["missing_strand_column"] += 1
                continue

            strand = normalize_strand(
                fields[strand_index]
            )

            if strand is None:

                qc["invalid_strand"] += 1
                continue

            if timepoint_index is not None:

                if len(fields) <= timepoint_index:

                    qc["missing_timepoint_column"] += 1
                    continue

                timepoint = normalize_timepoint(
                    fields[timepoint_index]
                )

            else:

                timepoint = infer_timepoint(
                    fields
                )

            if timepoint is None:

                qc["unrecognized_timepoint"] += 1
                timepoint = "unknown"

            reported_dinucleotide = ""

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

                qc["chromosome_not_in_reference"] += 1
                continue

            window_start = start - 1
            window_end = end + 1

            if (
                window_start < 0
                or window_end >
                chromosome_sizes[chrom]
            ):

                qc["window_outside_chromosome"] += 1
                continue

            selected_index += 1

            event_id = (
                f"event_{selected_index:012d}"
            )

            bed_writer.writerow(
                [
                    chrom,
                    window_start,
                    window_end,
                    event_id,
                    0,
                    strand,
                ]
            )

            metadata_writer.writerow(
                {
                    "event_id": event_id,
                    "timepoint": timepoint,
                    "chrom": chrom,
                    "cpd_start": start,
                    "cpd_end": end,
                    "strand": strand,
                    "reported_dinucleotide": (
                        reported_dinucleotide
                    ),
                }
            )

            qc["selected_windows"] += 1

    return qc


def extract_sequences(
    reference,
    window_bed,
    output_file,
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

    with output_file.open("w") as output_handle:

        process = subprocess.run(
            command,
            stdout=output_handle,
            stderr=subprocess.PIPE,
            text=True,
        )

    if process.returncode != 0:

        raise RuntimeError(
            "bedtools getfasta failed:\n"
            + process.stderr
        )


def read_metadata(path):

    metadata = {}

    with path.open() as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        for row in reader:

            metadata[
                row["event_id"]
            ] = row

    return metadata


def process_sequences(
    sequence_file,
    metadata,
    event_output,
):

    qc = Counter()

    overall_counts = Counter()
    timepoint_counts = defaultdict(
        Counter
    )

    with (
        sequence_file.open() as input_handle,
        event_output.open("w") as output_handle,
    ):

        columns = [
            "event_id",
            "timepoint",
            "chrom",
            "cpd_start",
            "cpd_end",
            "strand",
            "tetranucleotide",
            "left_flank",
            "dipyrimidine",
            "right_flank",
            "reported_dinucleotide",
            "reported_matches_extracted",
        ]

        writer = csv.DictWriter(
            output_handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()

        for line_number, line in enumerate(
            input_handle,
            start=1,
        ):

            fields = line.rstrip(
                "\n"
            ).split("\t")

            if len(fields) != 7:

                raise RuntimeError(
                    "Unexpected bedtools output at line "
                    f"{line_number}: {len(fields)} columns"
                )

            event_id = fields[3]
            sequence = fields[6].upper()

            qc["sequence_rows"] += 1

            if event_id not in metadata:

                raise RuntimeError(
                    f"Unknown event ID: {event_id}"
                )

            if len(sequence) != 4:

                qc["not_four_nucleotides"] += 1
                continue

            if not set(sequence) <= VALID_BASES:

                qc["non_acgt_sequence"] += 1
                continue

            dipyrimidine = sequence[1:3]

            if dipyrimidine not in VALID_DIPYRIMIDINES:

                qc[
                    "central_bases_not_dipyrimidine"
                ] += 1

                continue

            row_metadata = metadata[
                event_id
            ]

            reported = row_metadata[
                "reported_dinucleotide"
            ]

            reported_matches = ""

            if reported:

                reported_matches = (
                    reported == dipyrimidine
                )

                if reported_matches:

                    qc[
                        "reported_dinucleotide_matches"
                    ] += 1

                else:

                    qc[
                        "reported_dinucleotide_mismatch"
                    ] += 1

                    continue

            timepoint = row_metadata[
                "timepoint"
            ]

            overall_counts[
                sequence
            ] += 1

            timepoint_counts[
                timepoint
            ][sequence] += 1

            writer.writerow(
                {
                    "event_id": event_id,
                    "timepoint": timepoint,
                    "chrom": row_metadata[
                        "chrom"
                    ],
                    "cpd_start": row_metadata[
                        "cpd_start"
                    ],
                    "cpd_end": row_metadata[
                        "cpd_end"
                    ],
                    "strand": row_metadata[
                        "strand"
                    ],
                    "tetranucleotide": sequence,
                    "left_flank": sequence[0],
                    "dipyrimidine": dipyrimidine,
                    "right_flank": sequence[3],
                    "reported_dinucleotide": (
                        reported
                    ),
                    "reported_matches_extracted": (
                        reported_matches
                    ),
                }
            )

            qc["valid_tetranucleotides"] += 1

    return (
        overall_counts,
        timepoint_counts,
        qc,
    )


def all_possible_tetranucleotides():

    tetranucleotides = []

    for left_base in "ACGT":

        for dimer in (
            "CC",
            "CT",
            "TC",
            "TT",
        ):

            for right_base in "ACGT":

                tetranucleotides.append(
                    left_base
                    + dimer
                    + right_base
                )

    return tetranucleotides


def write_distributions(
    overall_counts,
    timepoint_counts,
    overall_output,
    timepoint_output,
):

    tetranucleotide_order = (
        all_possible_tetranucleotides()
    )

    overall_total = sum(
        overall_counts.values()
    )

    with overall_output.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "tetranucleotide",
                "left_flank",
                "dipyrimidine",
                "right_flank",
                "count",
                "percentage",
            ]
        )

        for tetranucleotide in (
            tetranucleotide_order
        ):

            count = overall_counts[
                tetranucleotide
            ]

            percentage = (
                100 * count / overall_total
                if overall_total > 0
                else 0
            )

            writer.writerow(
                [
                    tetranucleotide,
                    tetranucleotide[0],
                    tetranucleotide[1:3],
                    tetranucleotide[3],
                    count,
                    percentage,
                ]
            )

    with timepoint_output.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "timepoint",
                "tetranucleotide",
                "left_flank",
                "dipyrimidine",
                "right_flank",
                "count",
                "timepoint_total",
                "percentage_within_timepoint",
            ]
        )

        observed_timepoints = list(
            TIMEPOINT_ORDER
        )

        for value in sorted(
            timepoint_counts
        ):

            if value not in observed_timepoints:

                observed_timepoints.append(
                    value
                )

        for timepoint in observed_timepoints:

            total = sum(
                timepoint_counts[
                    timepoint
                ].values()
            )

            for tetranucleotide in (
                tetranucleotide_order
            ):

                count = timepoint_counts[
                    timepoint
                ][tetranucleotide]

                percentage = (
                    100 * count / total
                    if total > 0
                    else 0
                )

                writer.writerow(
                    [
                        timepoint,
                        tetranucleotide,
                        tetranucleotide[0],
                        tetranucleotide[1:3],
                        tetranucleotide[3],
                        count,
                        total,
                        percentage,
                    ]
                )


def write_qc(
    output_path,
    preparation_qc,
    sequence_qc,
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
            sequence_qc.items()
        ):

            writer.writerow(
                [
                    f"sequence_{metric}",
                    value,
                ]
            )


def main():

    args = parse_args()

    for path in (
        args.input_bed,
        args.reference,
    ):

        if not path.exists():

            raise FileNotFoundError(
                f"Missing input: {path}"
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
        / "UV_tetranucleotide_windows.hg38.bed6"
    )

    metadata_file = (
        args.outdir
        / "UV_tetranucleotide_windows.metadata.tsv"
    )

    sequence_file = (
        args.outdir
        / "UV_tetranucleotide_windows.sequence.tsv"
    )

    event_output = (
        args.outdir
        / "UV_tetranucleotide_events.tsv"
    )

    overall_output = (
        args.outdir
        / "UV_tetranucleotide_distribution.overall.tsv"
    )

    timepoint_output = (
        args.outdir
        / "UV_tetranucleotide_distribution.by_timepoint.tsv"
    )

    qc_output = (
        args.outdir
        / "UV_tetranucleotide_distribution.QC.tsv"
    )

    preparation_qc = prepare_windows(
        args.input_bed,
        window_bed,
        metadata_file,
        chromosome_sizes,
        args.timepoint_column,
        args.strand_column,
        args.dinucleotide_column,
    )

    extract_sequences(
        args.reference,
        window_bed,
        sequence_file,
    )

    metadata = read_metadata(
        metadata_file
    )

    (
        overall_counts,
        timepoint_counts,
        sequence_qc,
    ) = process_sequences(
        sequence_file,
        metadata,
        event_output,
    )

    write_distributions(
        overall_counts,
        timepoint_counts,
        overall_output,
        timepoint_output,
    )

    write_qc(
        qc_output,
        preparation_qc,
        sequence_qc,
    )

    print()
    print("Tetranucleotide extraction complete")
    print("===================================")
    print(f"Event table: {event_output}")
    print(f"Overall distribution: {overall_output}")
    print(f"Timepoint distribution: {timepoint_output}")
    print(f"QC: {qc_output}")


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)