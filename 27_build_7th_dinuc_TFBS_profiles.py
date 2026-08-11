#!/usr/bin/env python3

import argparse
import bisect
import csv
import statistics
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Set, Tuple


TIMEPOINTS = [
    "0.5h",
    "2h",
    "4h",
    "8h",
]

DINUCLEOTIDES = [
    "CC",
    "CT",
    "TC",
    "TT",
]

RELATIVE_POSITIONS = list(
    range(-10, 11)
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Create TFBS-centered profiles from the seventh "
            "dinucleotide of final deduplicated UV repair reads."
        )
    )

    parser.add_argument(
        "--tfbs",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--repair",
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


def is_integer(value: str) -> bool:
    try:
        int(value)
        return True
    except (ValueError, TypeError):
        return False


# ============================================================
# TFBS
#
# BED4:
#
# chrom  start  end  TF1,TF2
#
# Exact duplicate TF-site assignments are removed.
# ============================================================

def read_tfbs(
    path: Path,
) -> Tuple[
    Dict[str, List[Tuple[str, int, int, int]]],
    Counter,
]:

    tf_sites: Dict[
        str,
        Set[Tuple[str, int, int, int]],
    ] = defaultdict(set)

    qc = Counter()

    with open(path, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            qc["input_rows"] += 1

            if len(fields) < 4:
                qc["malformed_rows"] += 1
                continue

            if (
                not is_integer(fields[1])
                or not is_integer(fields[2])
            ):
                qc["invalid_coordinates"] += 1
                continue

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            if start < 0 or end <= start:
                qc["invalid_coordinates"] += 1
                continue

            tf_names = [
                value.strip()
                for value in fields[3].split(",")
                if value.strip()
            ]

            if not tf_names:
                qc["missing_tf_name"] += 1
                continue

            # Integer motif-center anchor.
            #
            # For an even-width interval, this selects the
            # right-central genomic base.
            center = (
                start + end
            ) // 2

            for tf_name in tf_names:
                site = (
                    chrom,
                    start,
                    end,
                    center,
                )

                if site in tf_sites[tf_name]:
                    qc[
                        "duplicate_TF_site_assignment_removed"
                    ] += 1

                tf_sites[tf_name].add(site)

    output = {}

    for tf_name, sites in tf_sites.items():
        output[tf_name] = sorted(
            sites,
            key=lambda row: (
                row[0],
                row[1],
                row[2],
            ),
        )

    qc["unique_TFs"] = len(output)

    qc["unique_TF_site_assignments"] = sum(
        len(sites)
        for sites in output.values()
    )

    return output, qc


# ============================================================
# REPAIR EVENTS
#
# Final 13-column seventh-dinucleotide file:
#
#  1 chrom
#  2 dinucleotide start
#  3 dinucleotide end
#  4 unique event ID
#  5 score
#  6 read strand
#  7 sample
#  8 timepoint
#  9 time in hours
# 10 regular or CC>TT_rescue
# 11 read length
# 12 original read ID
# 13 dinucleotide
#
# Each dinucleotide is represented by its leftmost genomic
# base: column 2.
# ============================================================

def read_repair(
    path: Path,
):
    positions = defaultdict(
        lambda: defaultdict(list)
    )

    total_all = Counter()
    total_series = Counter()

    qc = Counter()
    seen_event_ids = set()

    with open(path, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            qc["input_rows"] += 1

            if len(fields) < 13:
                qc["malformed_rows"] += 1
                continue

            chrom = fields[0]
            start_value = fields[1]
            end_value = fields[2]
            event_id = fields[3]
            timepoint = fields[7]
            read_length_value = fields[10]
            dinucleotide = fields[12].upper()

            if event_id == "":
                qc["missing_event_id"] += 1
                continue

            if event_id in seen_event_ids:
                qc["duplicate_event_id_removed"] += 1
                continue

            if (
                not is_integer(start_value)
                or not is_integer(end_value)
                or not is_integer(read_length_value)
            ):
                qc["invalid_numeric_value"] += 1
                continue

            start = int(start_value)
            end = int(end_value)
            read_length = int(read_length_value)

            if end - start != 2:
                qc["non_2bp_dinucleotide"] += 1
                continue

            if not 26 <= read_length <= 30:
                qc["read_outside_26_to_30"] += 1
                continue

            if timepoint not in TIMEPOINTS:
                qc["unexpected_timepoint"] += 1
                continue

            if dinucleotide not in DINUCLEOTIDES:
                qc["non_dipyrimidine"] += 1
                continue

            seen_event_ids.add(event_id)

            # One genomic point per repair event.
            event_position = start

            positions[
                (
                    timepoint,
                    "Combined",
                )
            ][chrom].append(event_position)

            positions[
                (
                    timepoint,
                    dinucleotide,
                )
            ][chrom].append(event_position)

            total_all[timepoint] += 1

            total_series[
                (
                    timepoint,
                    "Combined",
                )
            ] += 1

            total_series[
                (
                    timepoint,
                    dinucleotide,
                )
            ] += 1

            qc["selected_repair_events"] += 1
            qc[f"timepoint_{timepoint}"] += 1
            qc[f"dinucleotide_{dinucleotide}"] += 1

    for dataset in positions.values():
        for chrom in dataset:
            dataset[chrom].sort()

    return (
        positions,
        total_all,
        total_series,
        qc,
    )


# ============================================================
# DAMAGE EVENTS
#
# Existing 10-column 1-bp Damage-seq file:
#
#  1 chrom
#  2 start
#  3 end
#  4 event ID
#  5 score
#  6 strand
#  7 dinucleotide
#  8 sequence
#  9 original window start
# 10 original window end
# ============================================================

def read_damage(
    path: Path,
):
    positions = defaultdict(
        lambda: defaultdict(list)
    )

    total_all = 0
    total_series = Counter()

    qc = Counter()
    seen_event_ids = set()

    with open(path, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            qc["input_rows"] += 1

            if len(fields) < 10:
                qc["malformed_rows"] += 1
                continue

            chrom = fields[0]
            position_value = fields[1]
            event_id = fields[3]
            dinucleotide = fields[6].upper()

            if event_id == "":
                qc["missing_event_id"] += 1
                continue

            if event_id in seen_event_ids:
                qc["duplicate_event_id_removed"] += 1
                continue

            if not is_integer(position_value):
                qc["invalid_position"] += 1
                continue

            if dinucleotide not in DINUCLEOTIDES:
                qc["non_dipyrimidine"] += 1
                continue

            position = int(position_value)

            seen_event_ids.add(event_id)

            positions["Combined"][chrom].append(
                position
            )

            positions[dinucleotide][chrom].append(
                position
            )

            total_all += 1
            total_series["Combined"] += 1
            total_series[dinucleotide] += 1

            qc["selected_damage_events"] += 1
            qc[f"dinucleotide_{dinucleotide}"] += 1

    for dataset in positions.values():
        for chrom in dataset:
            dataset[chrom].sort()

    return (
        positions,
        total_all,
        total_series,
        qc,
    )


# ============================================================
# PROFILE COUNTING
#
# The same genomic event may contribute to two nearby TFBS
# loci. This is expected in a site-centered metaprofile:
# every TFBS locus is treated as one observation.
# ============================================================

def count_profile(
    sites: List[Tuple[str, int, int, int]],
    positions_by_chrom: Dict[str, List[int]],
) -> Tuple[Dict[int, int], int]:

    counts = {
        position: 0
        for position in RELATIVE_POSITIONS
    }

    sites_used = 0

    for chrom, start, end, center in sites:
        window_start = center - 10
        window_end = center + 10

        if window_start < 0:
            continue

        sites_used += 1

        chrom_positions = positions_by_chrom.get(
            chrom
        )

        if not chrom_positions:
            continue

        left_index = bisect.bisect_left(
            chrom_positions,
            window_start,
        )

        right_index = bisect.bisect_right(
            chrom_positions,
            window_end,
        )

        for event_position in chrom_positions[
            left_index:right_index
        ]:
            relative_position = (
                event_position -
                center
            )

            if relative_position in counts:
                counts[relative_position] += 1

    return counts, sites_used


def normalized_signal(
    event_count: int,
    total_events: int,
    sites_used: int,
) -> float:

    if total_events <= 0 or sites_used <= 0:
        return float("nan")

    return (
        event_count
        / total_events
        * 1_000_000
        / sites_used
        * 1_000
    )


def write_qc(
    path: Path,
    sections: Dict[str, Counter],
) -> None:

    with open(path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "section",
                "metric",
                "value",
            ]
        )

        for section, counter in sections.items():
            for metric in sorted(counter):
                writer.writerow(
                    [
                        section,
                        metric,
                        counter[metric],
                    ]
                )


# ============================================================
# MAIN
# ============================================================

def main() -> None:
    args = parse_args()

    for path in [
        args.tfbs,
        args.repair,
        args.damage,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    tf_sites, tfbs_qc = read_tfbs(
        args.tfbs
    )

    (
        repair_positions,
        repair_total_all,
        repair_total_series,
        repair_qc,
    ) = read_repair(
        args.repair
    )

    (
        damage_positions,
        damage_total_all,
        damage_total_series,
        damage_qc,
    ) = read_damage(
        args.damage
    )

    profile_file = (
        args.outdir
        / "TFBS_7th_dinucleotide_minus10_plus10_profiles.tsv"
    )

    summary_file = (
        args.outdir
        / "TFBS_7th_dinucleotide_TF_summary.tsv"
    )

    summary_rows = []

    with open(profile_file, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "TF",
                "n_sites_total",
                "n_sites_used",
                "median_site_width",
                "data_type",
                "series",
                "timepoint",
                "relative_position",
                "event_count",
                "total_all_events",
                "total_series_events",
                "events_per_site",
                "RPM_per_1000_sites_all_library",
                "RPM_per_1000_sites_within_series",
            ]
        )

        sorted_tfs = sorted(
            tf_sites
        )

        for tf_index, tf_name in enumerate(
            sorted_tfs,
            start=1,
        ):
            sites = tf_sites[tf_name]

            widths = [
                end - start
                for chrom, start, end, center
                in sites
            ]

            median_width = statistics.median(
                widths
            )

            valid_sites = [
                site
                for site in sites
                if site[3] - 10 >= 0
            ]

            sites_used = len(valid_sites)

            if sites_used == 0:
                continue

            print(
                f"[{tf_index}/{len(sorted_tfs)}] "
                f"{tf_name}: {sites_used:,} sites"
            )

            summary_rows.append(
                [
                    tf_name,
                    len(sites),
                    sites_used,
                    median_width,
                ]
            )

            # ------------------------------------------------
            # REPAIR PROFILES
            # ------------------------------------------------

            for timepoint in TIMEPOINTS:
                for series in [
                    "Combined",
                    "CC",
                    "CT",
                    "TC",
                    "TT",
                ]:
                    dataset = repair_positions.get(
                        (
                            timepoint,
                            series,
                        ),
                        {},
                    )

                    counts, counted_sites = count_profile(
                        sites=valid_sites,
                        positions_by_chrom=dataset,
                    )

                    if counted_sites != sites_used:
                        raise RuntimeError(
                            "Internal repair site-count mismatch."
                        )

                    total_all_events = repair_total_all.get(
                        timepoint,
                        0,
                    )

                    total_series_events = (
                        repair_total_series.get(
                            (
                                timepoint,
                                series,
                            ),
                            0,
                        )
                    )

                    for relative_position in RELATIVE_POSITIONS:
                        event_count = counts[
                            relative_position
                        ]

                        writer.writerow(
                            [
                                tf_name,
                                len(sites),
                                sites_used,
                                median_width,
                                "repair",
                                series,
                                timepoint,
                                relative_position,
                                event_count,
                                total_all_events,
                                total_series_events,
                                event_count / sites_used,
                                normalized_signal(
                                    event_count,
                                    total_all_events,
                                    sites_used,
                                ),
                                normalized_signal(
                                    event_count,
                                    total_series_events,
                                    sites_used,
                                ),
                            ]
                        )

            # ------------------------------------------------
            # DAMAGE PROFILES
            # ------------------------------------------------

            for series in [
                "Combined",
                "CC",
                "CT",
                "TC",
                "TT",
            ]:
                dataset = damage_positions.get(
                    series,
                    {},
                )

                counts, counted_sites = count_profile(
                    sites=valid_sites,
                    positions_by_chrom=dataset,
                )

                if counted_sites != sites_used:
                    raise RuntimeError(
                        "Internal damage site-count mismatch."
                    )

                total_series_events = (
                    damage_total_series.get(
                        series,
                        0,
                    )
                )

                for relative_position in RELATIVE_POSITIONS:
                    event_count = counts[
                        relative_position
                    ]

                    writer.writerow(
                        [
                            tf_name,
                            len(sites),
                            sites_used,
                            median_width,
                            "damage",
                            series,
                            "0h",
                            relative_position,
                            event_count,
                            damage_total_all,
                            total_series_events,
                            event_count / sites_used,
                            normalized_signal(
                                event_count,
                                damage_total_all,
                                sites_used,
                            ),
                            normalized_signal(
                                event_count,
                                total_series_events,
                                sites_used,
                            ),
                        ]
                    )

    with open(summary_file, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "TF",
                "n_sites_total",
                "n_sites_used",
                "median_site_width",
            ]
        )

        writer.writerows(
            summary_rows
        )

    write_qc(
        args.outdir
        / "TFBS_7th_dinucleotide_profile_QC.tsv",
        {
            "TFBS": tfbs_qc,
            "repair": repair_qc,
            "damage": damage_qc,
        },
    )

    print()
    print("Done.")
    print(profile_file)


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)