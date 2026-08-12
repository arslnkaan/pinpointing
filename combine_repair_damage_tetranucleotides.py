#!/usr/bin/env python3

import argparse
import csv
import math
import sys
from pathlib import Path


TIMEPOINT_ORDER = (
    "0.5h",
    "2h",
    "4h",
    "8h",
)

DIPYRIMIDINE_ORDER = (
    "CC",
    "CT",
    "TC",
    "TT",
)


def parse_args():

    parser = argparse.ArgumentParser(
        description=(
            "Combine UV repair and Damage-seq tetranucleotide "
            "distributions and calculate repair-to-damage ratios."
        )
    )

    parser.add_argument(
        "--repair-overall",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--repair-timepoints",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--damage",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
    )

    return parser.parse_args()


def get_tetranucleotide_order():

    return [
        left_base + dipyrimidine + right_base

        for dipyrimidine in DIPYRIMIDINE_ORDER

        for left_base in (
            "A",
            "C",
            "G",
            "T",
        )

        for right_base in (
            "A",
            "C",
            "G",
            "T",
        )
    ]


def read_table(
    path,
):

    with path.open() as handle:

        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        return list(reader)


def safe_float(value):

    try:
        return float(value)

    except (
        ValueError,
        TypeError,
    ):
        return 0.0


def safe_int(value):

    try:
        return int(float(value))

    except (
        ValueError,
        TypeError,
    ):
        return 0


def main():

    args = parse_args()

    for path in (
        args.repair_overall,
        args.repair_timepoints,
        args.damage,
    ):

        if not path.exists():

            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    repair_overall_rows = read_table(
        args.repair_overall
    )

    repair_timepoint_rows = read_table(
        args.repair_timepoints
    )

    damage_rows = read_table(
        args.damage
    )

    repair_overall = {
        row["tetranucleotide"]: row
        for row in repair_overall_rows
    }

    repair_timepoints = {
        (
            row["timepoint"],
            row["tetranucleotide"],
        ): row
        for row in repair_timepoint_rows
    }

    damage = {
        row["tetranucleotide"]: row
        for row in damage_rows
    }

    tetranucleotide_order = (
        get_tetranucleotide_order()
    )

    # ========================================================
    # COMBINED RAW DISTRIBUTIONS
    # ========================================================

    raw_output = (
        args.outdir
        / "repair_damage.tetranucleotide.raw_distributions.tsv"
    )

    with raw_output.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "dataset",
                "timepoint",
                "tetranucleotide",
                "dipyrimidine",
                "count",
                "total",
                "percentage",
            ]
        )

        for tetranucleotide in tetranucleotide_order:

            repair_row = repair_overall.get(
                tetranucleotide,
                {},
            )

            writer.writerow(
                [
                    "Repair",
                    "Combined",
                    tetranucleotide,
                    tetranucleotide[1:3],
                    safe_int(
                        repair_row.get(
                            "count",
                            0,
                        )
                    ),
                    "",
                    safe_float(
                        repair_row.get(
                            "percentage",
                            0,
                        )
                    ),
                ]
            )

            damage_row = damage.get(
                tetranucleotide,
                {},
            )

            writer.writerow(
                [
                    "Damage-seq",
                    "Damage",
                    tetranucleotide,
                    tetranucleotide[1:3],
                    safe_int(
                        damage_row.get(
                            "count",
                            0,
                        )
                    ),
                    safe_int(
                        damage_row.get(
                            "damage_total",
                            0,
                        )
                    ),
                    safe_float(
                        damage_row.get(
                            "damage_percentage",
                            0,
                        )
                    ),
                ]
            )

        for timepoint in TIMEPOINT_ORDER:

            for tetranucleotide in tetranucleotide_order:

                repair_row = repair_timepoints.get(
                    (
                        timepoint,
                        tetranucleotide,
                    ),
                    {},
                )

                writer.writerow(
                    [
                        "Repair",
                        timepoint,
                        tetranucleotide,
                        tetranucleotide[1:3],
                        safe_int(
                            repair_row.get(
                                "count",
                                0,
                            )
                        ),
                        safe_int(
                            repair_row.get(
                                "timepoint_total",
                                0,
                            )
                        ),
                        safe_float(
                            repair_row.get(
                                "percentage_within_timepoint",
                                0,
                            )
                        ),
                    ]
                )

    # ========================================================
    # NORMALIZED REPAIR / DAMAGE
    # ========================================================

    normalized_output = (
        args.outdir
        / "repair.normalized_to_DamageSeq.tetranucleotides.tsv"
    )

    with normalized_output.open("w") as handle:

        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "timepoint",
                "tetranucleotide",
                "position_minus3",
                "position_minus2",
                "position_minus1",
                "position_0",
                "dipyrimidine",
                "repair_count",
                "repair_total",
                "repair_percentage",
                "damage_count",
                "damage_total",
                "damage_percentage",
                "repair_over_damage",
                "log2_repair_over_damage",
                "normalization_status",
            ]
        )

        for timepoint in TIMEPOINT_ORDER:

            for tetranucleotide in tetranucleotide_order:

                repair_row = repair_timepoints.get(
                    (
                        timepoint,
                        tetranucleotide,
                    ),
                    {},
                )

                damage_row = damage.get(
                    tetranucleotide,
                    {},
                )

                repair_count = safe_int(
                    repair_row.get(
                        "count",
                        0,
                    )
                )

                repair_total = safe_int(
                    repair_row.get(
                        "timepoint_total",
                        0,
                    )
                )

                repair_percentage = safe_float(
                    repair_row.get(
                        "percentage_within_timepoint",
                        0,
                    )
                )

                damage_count = safe_int(
                    damage_row.get(
                        "count",
                        0,
                    )
                )

                damage_total = safe_int(
                    damage_row.get(
                        "damage_total",
                        0,
                    )
                )

                damage_percentage = safe_float(
                    damage_row.get(
                        "damage_percentage",
                        0,
                    )
                )

                if damage_percentage <= 0:

                    repair_over_damage = ""
                    log2_ratio = ""
                    status = "zero_damage"

                else:

                    repair_over_damage = (
                        repair_percentage
                        / damage_percentage
                    )

                    if repair_over_damage > 0:

                        log2_ratio = math.log2(
                            repair_over_damage
                        )

                        status = "supported"

                    else:

                        log2_ratio = ""
                        status = "zero_repair"

                writer.writerow(
                    [
                        timepoint,
                        tetranucleotide,
                        tetranucleotide[0],
                        tetranucleotide[1],
                        tetranucleotide[2],
                        tetranucleotide[3],
                        tetranucleotide[1:3],
                        repair_count,
                        repair_total,
                        repair_percentage,
                        damage_count,
                        damage_total,
                        damage_percentage,
                        repair_over_damage,
                        log2_ratio,
                        status,
                    ]
                )

    print()
    print("Repair/Damage-seq tetranucleotide analysis complete")
    print("===================================================")
    print(f"Raw distributions: {raw_output}")
    print(f"Damage-normalized repair: {normalized_output}")


if __name__ == "__main__":

    try:
        main()

    except Exception as error:

        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )

        sys.exit(1)