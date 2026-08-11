#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from collections import Counter
from pathlib import Path


DIPYRIMIDINES = [
    "CT",
    "TC",
    "TT",
    "CC",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Count CC, CT, TC and TT at one overlapping "
            "dinucleotide position counted from the 3-prime end."
        )
    )

    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="bedtools getfasta -name -tab output.",
    )

    parser.add_argument(
        "--sample",
        required=True,
    )

    parser.add_argument(
        "--timepoint",
        required=True,
    )

    parser.add_argument(
        "--time-h",
        required=True,
        type=float,
    )

    parser.add_argument(
        "--dimer-from-3prime",
        type=int,
        default=7,
        help=(
            "Overlapping dimer number counted from the 3-prime end. "
            "Dimer 1 is the final two bases. "
            "Dimer 7 is bases 8 and 7 from the 3-prime end."
        ),
    )

    parser.add_argument(
        "--events",
        required=True,
        type=Path,
        help="Output event-level dipyrimidine table.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Output count and percentage table.",
    )

    parser.add_argument(
        "--qc",
        required=True,
        type=Path,
        help="Output sample QC table.",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    if not args.input.exists():
        raise FileNotFoundError(
            f"Missing input: {args.input}"
        )

    if args.dimer_from_3prime < 1:
        raise ValueError(
            "--dimer-from-3prime must be >= 1"
        )

    for path in [
        args.events,
        args.summary,
        args.qc,
    ]:
        path.parent.mkdir(
            parents=True,
            exist_ok=True,
        )

    dimer_counts = Counter()
    source_total_counts = Counter()
    source_dipyrimidine_counts = Counter()

    input_rows = 0
    valid_sequence_rows = 0
    invalid_sequence_rows = 0
    too_short_rows = 0
    dipyrimidine_rows = 0

    with open(args.events, "w") as event_handle:
        event_writer = csv.writer(
            event_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        event_writer.writerow([
            "sample",
            "timepoint",
            "time_h",
            "sequence_id",
            "source",
            "read_id",
            "read_length",
            "dimer_number_from_3prime",
            "first_base_position_from_3prime",
            "second_base_position_from_3prime",
            "dinucleotide",
            "sequence",
        ])

        with open(args.input) as input_handle:
            reader = csv.reader(
                input_handle,
                delimiter="\t",
            )

            for fields in reader:
                if len(fields) < 2:
                    continue

                input_rows += 1

                raw_identifier = fields[0]
                sequence = fields[1].strip().upper()

                sequence_id = re.sub(
                    r"::.*$",
                    "",
                    raw_identifier,
                )

                if "|" in sequence_id:
                    source, read_id = sequence_id.split(
                        "|",
                        1,
                    )
                else:
                    source = "unknown"
                    read_id = sequence_id

                if not re.fullmatch(
                    r"[ACGT]+",
                    sequence,
                ):
                    invalid_sequence_rows += 1
                    continue

                sequence_length = len(sequence)

                # Dimer 1 = last two bases.
                # Dimer 7 = bases 8 and 7 from the 3' end.
                minimum_length = (
                    args.dimer_from_3prime + 1
                )

                if sequence_length < minimum_length:
                    too_short_rows += 1
                    continue

                valid_sequence_rows += 1
                source_total_counts[source] += 1

                start_index = (
                    sequence_length -
                    args.dimer_from_3prime -
                    1
                )

                dinucleotide = sequence[
                    start_index:
                    start_index + 2
                ]

                if dinucleotide not in DIPYRIMIDINES:
                    continue

                dipyrimidine_rows += 1
                dimer_counts[dinucleotide] += 1
                source_dipyrimidine_counts[
                    source
                ] += 1

                event_writer.writerow([
                    args.sample,
                    args.timepoint,
                    args.time_h,
                    sequence_id,
                    source,
                    read_id,
                    sequence_length,
                    args.dimer_from_3prime,
                    args.dimer_from_3prime + 1,
                    args.dimer_from_3prime,
                    dinucleotide,
                    sequence,
                ])

    # ========================================================
    # SUMMARY
    # ========================================================

    with open(args.summary, "w") as summary_handle:
        writer = csv.writer(
            summary_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "sample",
            "timepoint",
            "time_h",
            "dimer_number_from_3prime",
            "dinucleotide",
            "count",
            "total_valid_reads",
            "total_dipyrimidine_reads",
            "percent_of_all_reads",
            "percent_of_dipyrimidine_reads",
        ])

        for dinucleotide in DIPYRIMIDINES:
            count = dimer_counts.get(
                dinucleotide,
                0,
            )

            percent_all = (
                100.0 *
                count /
                valid_sequence_rows
                if valid_sequence_rows > 0
                else 0.0
            )

            percent_dipyrimidine = (
                100.0 *
                count /
                dipyrimidine_rows
                if dipyrimidine_rows > 0
                else 0.0
            )

            writer.writerow([
                args.sample,
                args.timepoint,
                args.time_h,
                args.dimer_from_3prime,
                dinucleotide,
                count,
                valid_sequence_rows,
                dipyrimidine_rows,
                percent_all,
                percent_dipyrimidine,
            ])

    # ========================================================
    # QC
    # ========================================================

    with open(args.qc, "w") as qc_handle:
        writer = csv.writer(
            qc_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "sample",
            "timepoint",
            "input_rows",
            "valid_sequence_rows",
            "invalid_sequence_rows",
            "too_short_rows",
            "dipyrimidine_rows",
            "regular_valid_reads",
            "rescued_valid_reads",
            "regular_dipyrimidine_reads",
            "rescued_dipyrimidine_reads",
        ])

        writer.writerow([
            args.sample,
            args.timepoint,
            input_rows,
            valid_sequence_rows,
            invalid_sequence_rows,
            too_short_rows,
            dipyrimidine_rows,
            source_total_counts.get(
                "regular",
                0,
            ),
            source_total_counts.get(
                "rescued",
                0,
            ),
            source_dipyrimidine_counts.get(
                "regular",
                0,
            ),
            source_dipyrimidine_counts.get(
                "rescued",
                0,
            ),
        ])

    print(
        f"{args.sample}: "
        f"{valid_sequence_rows:,} valid reads; "
        f"{dipyrimidine_rows:,} dipyrimidines"
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