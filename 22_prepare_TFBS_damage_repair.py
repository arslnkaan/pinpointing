#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Tuple


DIPYRIMIDINES = {
    "CC",
    "CT",
    "TC",
    "TT",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Prepare active TFBS loci, Damage-seq dipyrimidine "
            "lesions, and pinpointed UV repair coordinates."
        )
    )

    parser.add_argument(
        "--tfbs",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--damage",
        required=True,
        type=Path,
    )

    parser.add_argument(
        "--repair-map",
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
    except ValueError:
        return False


def normalize_sequence(sequence: str) -> str:
    return re.sub(
        r"[^ACGTN]",
        "",
        sequence.upper(),
    )


def find_header_index(
    header: Sequence[str],
    aliases: Sequence[str],
) -> Optional[int]:
    normalized = [
        value.strip().lower()
        for value in header
    ]

    for alias in aliases:
        if alias.lower() in normalized:
            return normalized.index(
                alias.lower()
            )

    return None


def parse_coordinate_name(
    value: str,
) -> Optional[Tuple[str, int, int, str]]:
    """
    Supports names such as:

      chr1:100-104
      chr1:100-104(+)
      chr1:100-104(-)
      chr1:100-104::+
      chr1:100-104::-
    """

    coordinate_match = re.search(
        r"(chr[^:\s|]+):(\d+)-(\d+)",
        value,
    )

    if coordinate_match is None:
        return None

    chrom = coordinate_match.group(1)
    start = int(coordinate_match.group(2))
    end = int(coordinate_match.group(3))

    strand = "."

    strand_match = re.search(
        r"\(([+-])\)",
        value,
    )

    if strand_match is not None:
        strand = strand_match.group(1)
    else:
        strand_match = re.search(
            r"::([+-])",
            value,
        )

        if strand_match is not None:
            strand = strand_match.group(1)

    return (
        chrom,
        start,
        end,
        strand,
    )


def parse_damage_header(
    header: List[str],
) -> Optional[Dict[str, Optional[int]]]:
    chrom_index = find_header_index(
        header,
        [
            "chrom",
            "chromosome",
            "chr",
        ],
    )

    start_index = find_header_index(
        header,
        [
            "start",
            "start_0based",
            "window_start",
            "region_start",
        ],
    )

    end_index = find_header_index(
        header,
        [
            "end",
            "end_0based",
            "window_end",
            "region_end",
        ],
    )

    sequence_index = find_header_index(
        header,
        [
            "sequence",
            "seq",
            "raw_sequence",
            "reference_sequence",
        ],
    )

    strand_index = find_header_index(
        header,
        [
            "strand",
            "read_strand",
        ],
    )

    if (
        chrom_index is None
        or start_index is None
        or end_index is None
        or sequence_index is None
    ):
        return None

    return {
        "chrom": chrom_index,
        "start": start_index,
        "end": end_index,
        "sequence": sequence_index,
        "strand": strand_index,
    }


def parse_header_damage_row(
    fields: List[str],
    indices: Dict[str, Optional[int]],
) -> Optional[Tuple[str, int, int, str, str]]:
    chrom_index = indices["chrom"]
    start_index = indices["start"]
    end_index = indices["end"]
    sequence_index = indices["sequence"]
    strand_index = indices["strand"]

    if (
        chrom_index is None
        or start_index is None
        or end_index is None
        or sequence_index is None
    ):
        return None

    maximum_index = max(
        chrom_index,
        start_index,
        end_index,
        sequence_index,
    )

    if len(fields) <= maximum_index:
        return None

    if (
        not is_integer(fields[start_index])
        or not is_integer(fields[end_index])
    ):
        return None

    chrom = fields[chrom_index]
    start = int(fields[start_index])
    end = int(fields[end_index])

    sequence = normalize_sequence(
        fields[sequence_index]
    )

    strand = "."

    if (
        strand_index is not None
        and len(fields) > strand_index
        and fields[strand_index] in {"+", "-"}
    ):
        strand = fields[strand_index]

    return (
        chrom,
        start,
        end,
        strand,
        sequence,
    )


def parse_generic_damage_row(
    fields: List[str],
) -> Optional[Tuple[str, int, int, str, str]]:
    """
    Supports:

    1. BED-like rows:
       chrom start end ... strand ... sequence

    2. Two-column bedtools getfasta output:
       chr:start-end(strand) sequence
    """

    if (
        len(fields) >= 3
        and is_integer(fields[1])
        and is_integer(fields[2])
    ):
        chrom = fields[0]
        start = int(fields[1])
        end = int(fields[2])

        strand = "."

        for field in fields[3:]:
            if field in {"+", "-"}:
                strand = field
                break

        sequence = ""

        for field in reversed(fields):
            candidate = normalize_sequence(
                field
            )

            if (
                len(candidate) == 4
                and re.fullmatch(
                    r"[ACGTN]{4}",
                    candidate,
                )
                is not None
            ):
                sequence = candidate
                break

        if sequence == "":
            return None

        return (
            chrom,
            start,
            end,
            strand,
            sequence,
        )

    if len(fields) >= 2:
        parsed_name = parse_coordinate_name(
            fields[0]
        )

        if parsed_name is None:
            return None

        chrom, start, end, strand = parsed_name

        sequence = normalize_sequence(
            fields[1]
        )

        if len(sequence) != 4:
            return None

        return (
            chrom,
            start,
            end,
            strand,
            sequence,
        )

    return None


def prepare_tfbs(
    input_path: Path,
    exploded_bed: Path,
    unique_loci_bed: Path,
    summary_path: Path,
    qc_path: Path,
) -> None:
    unique_tf_sites = set()
    unique_loci = set()

    qc = Counter()

    with open(input_path, "r") as handle:
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

            chrom = fields[0]

            if (
                not is_integer(fields[1])
                or not is_integer(fields[2])
            ):
                qc["invalid_coordinates"] += 1
                continue

            start = int(fields[1])
            end = int(fields[2])

            if end <= start:
                qc["invalid_coordinates"] += 1
                continue

            tf_field = fields[3]

            tf_names = [
                value.strip()
                for value in tf_field.split(",")
                if value.strip()
            ]

            if not tf_names:
                qc["missing_tf_names"] += 1
                continue

            unique_loci.add(
                (
                    chrom,
                    start,
                    end,
                )
            )

            for tf_name in tf_names:
                unique_tf_sites.add(
                    (
                        chrom,
                        start,
                        end,
                        tf_name,
                    )
                )

    sorted_tf_sites = sorted(
        unique_tf_sites,
        key=lambda value: (
            value[0],
            value[1],
            value[2],
            value[3],
        ),
    )

    sorted_loci = sorted(
        unique_loci,
        key=lambda value: (
            value[0],
            value[1],
            value[2],
        ),
    )

    with open(exploded_bed, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for row in sorted_tf_sites:
            writer.writerow(row)

    with open(unique_loci_bed, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for row in sorted_loci:
            writer.writerow(row)

    tf_summary = defaultdict(
        lambda: {
            "sites": 0,
            "bp": 0,
        }
    )

    for chrom, start, end, tf_name in sorted_tf_sites:
        tf_summary[tf_name]["sites"] += 1
        tf_summary[tf_name]["bp"] += (
            end - start
        )

    with open(summary_path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow(
            [
                "TF",
                "n_sites",
                "total_bp",
            ]
        )

        for tf_name in sorted(tf_summary):
            writer.writerow(
                [
                    tf_name,
                    tf_summary[tf_name]["sites"],
                    tf_summary[tf_name]["bp"],
                ]
            )

    qc["unique_original_loci"] = len(
        unique_loci
    )

    qc["unique_tf_locus_assignments"] = len(
        unique_tf_sites
    )

    qc["unique_TFs"] = len(
        tf_summary
    )

    with open(qc_path, "w") as handle:
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

        for metric in sorted(qc):
            writer.writerow(
                [
                    metric,
                    qc[metric],
                ]
            )


def prepare_damage(
    input_path: Path,
    output_bed: Path,
    qc_path: Path,
) -> None:
    qc = Counter()

    with open(input_path, "r") as handle:
        raw_lines = [
            line.rstrip("\n")
            for line in handle
            if line.strip()
        ]

    if not raw_lines:
        raise RuntimeError(
            "Damage input is empty."
        )

    first_fields = raw_lines[0].split("\t")

    header_indices = parse_damage_header(
        first_fields
    )

    start_line = 0

    if header_indices is not None:
        start_line = 1
        qc["header_detected"] = 1
    else:
        qc["header_detected"] = 0

    event_number = 0

    with open(output_bed, "w") as output_handle:
        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for line in raw_lines[start_line:]:
            fields = line.split("\t")

            qc["input_rows"] += 1

            parsed = None

            if header_indices is not None:
                parsed = parse_header_damage_row(
                    fields,
                    header_indices,
                )
            else:
                parsed = parse_generic_damage_row(
                    fields
                )

            if parsed is None:
                qc["unparsed_rows"] += 1
                continue

            (
                chrom,
                window_start,
                window_end,
                strand,
                sequence,
            ) = parsed

            if len(sequence) != 4:
                qc["non_4bp_sequences"] += 1
                continue

            if window_end - window_start != 4:
                qc["non_4bp_windows"] += 1
                continue

            dimer = sequence[1:3]

            if dimer not in DIPYRIMIDINES:
                qc["non_dipyrimidine"] += 1
                continue

            # The sequence is ordered:
            #
            # character 1 = -3
            # character 2 = -2
            # character 3 = -1
            # character 4 =  0
            #
            # For either strand, characters 2 and 3 occupy
            # the middle two genomic bases of the 4-bp window.
            lesion_start = window_start + 1
            lesion_end = window_start + 3

            event_number += 1

            event_id = (
                "damage|"
                f"{event_number:012d}"
            )

            writer.writerow(
                [
                    chrom,
                    lesion_start,
                    lesion_end,
                    event_id,
                    0,
                    strand,
                    dimer,
                    sequence,
                    window_start,
                    window_end,
                ]
            )

            qc["selected_damage_events"] += 1
            qc[f"dimer_{dimer}"] += 1

    if qc["selected_damage_events"] == 0:
        raise RuntimeError(
            "No CC/CT/TC/TT damage events were selected."
        )

    with open(qc_path, "w") as handle:
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

        for metric in sorted(qc):
            writer.writerow(
                [
                    metric,
                    qc[metric],
                ]
            )


def parse_repair_event(
    fields: List[str],
) -> Optional[Dict[str, str]]:
    """
    Known C>T event-table layouts.

    Seventeen-column data:
      Change                     = column 11
      Genomic_Position_0based    = column 16

    Sixteen-column data:
      Change                     = column 10
      Genomic_Position_0based    = column 15
    """

    if len(fields) >= 17:
        return {
            "read_id": fields[0],
            "chrom": fields[1],
            "strand": fields[4],
            "change": fields[10],
            "position": fields[15],
        }

    if len(fields) >= 16:
        return {
            "read_id": fields[0],
            "chrom": fields[1],
            "strand": fields[4],
            "change": fields[9],
            "position": fields[14],
        }

    return None


def prepare_repair(
    repair_map_path: Path,
    output_bed: Path,
    qc_path: Path,
) -> None:
    qc_rows = []

    with open(
        repair_map_path,
        "r",
    ) as map_handle, open(
        output_bed,
        "w",
    ) as output_handle:

        map_reader = csv.DictReader(
            map_handle,
            delimiter="\t",
        )

        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        required_map_columns = {
            "sample",
            "timepoint",
            "time_h",
            "replicate",
            "input",
        }

        if map_reader.fieldnames is None:
            raise RuntimeError(
                "Repair map has no header."
            )

        missing = required_map_columns.difference(
            set(map_reader.fieldnames)
        )

        if missing:
            raise RuntimeError(
                "Repair map is missing: "
                + ", ".join(
                    sorted(missing)
                )
            )

        for map_row in map_reader:
            sample = map_row["sample"]
            timepoint = map_row["timepoint"]
            time_h = map_row["time_h"]
            replicate = map_row["replicate"]

            input_path = Path(
                map_row["input"]
            )

            if not input_path.exists():
                raise FileNotFoundError(
                    f"Missing repair input: {input_path}"
                )

            qc = Counter()

            with open(input_path, "r") as input_handle:
                for line in input_handle:
                    fields = line.rstrip("\n").split("\t")

                    if not fields:
                        continue

                    if fields[0] == "Read_ID":
                        continue

                    qc["input_rows"] += 1

                    event = parse_repair_event(
                        fields
                    )

                    if event is None:
                        qc["malformed_rows"] += 1
                        continue

                    if event["change"] != "C>T":
                        qc["non_CtoT_rows"] += 1
                        continue

                    if not is_integer(
                        event["position"]
                    ):
                        qc["invalid_position"] += 1
                        continue

                    position = int(
                        event["position"]
                    )

                    if position < 0:
                        qc["invalid_position"] += 1
                        continue

                    strand = event["strand"]

                    if strand not in {"+", "-"}:
                        qc["invalid_strand"] += 1
                        continue

                    qc["selected_repair_events"] += 1

                    event_id = (
                        f"{sample}|"
                        f"{qc['selected_repair_events']:012d}"
                    )

                    writer.writerow(
                        [
                            event["chrom"],
                            position,
                            position + 1,
                            event_id,
                            0,
                            strand,
                            sample,
                            timepoint,
                            time_h,
                            replicate,
                        ]
                    )

            if qc["selected_repair_events"] == 0:
                raise RuntimeError(
                    f"No C>T repair events selected for {sample}"
                )

            qc_rows.append(
                {
                    "sample": sample,
                    "timepoint": timepoint,
                    "time_h": time_h,
                    "replicate": replicate,
                    "input_rows": qc["input_rows"],
                    "selected_repair_events":
                        qc["selected_repair_events"],
                    "malformed_rows":
                        qc["malformed_rows"],
                    "non_CtoT_rows":
                        qc["non_CtoT_rows"],
                    "invalid_position":
                        qc["invalid_position"],
                    "invalid_strand":
                        qc["invalid_strand"],
                }
            )

    with open(qc_path, "w") as handle:
        fieldnames = [
            "sample",
            "timepoint",
            "time_h",
            "replicate",
            "input_rows",
            "selected_repair_events",
            "malformed_rows",
            "non_CtoT_rows",
            "invalid_position",
            "invalid_strand",
        ]

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writeheader()
        writer.writerows(qc_rows)


def main() -> None:
    args = parse_args()

    for path in [
        args.tfbs,
        args.damage,
        args.repair_map,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    prepare_tfbs(
        input_path=args.tfbs,
        exploded_bed=args.outdir
        / "active_TFBS.exploded.bed",
        unique_loci_bed=args.outdir
        / "active_TFBS.unique_loci.bed",
        summary_path=args.outdir
        / "active_TFBS.TF_summary.tsv",
        qc_path=args.outdir
        / "active_TFBS.QC.tsv",
    )

    prepare_damage(
        input_path=args.damage,
        output_bed=args.outdir
        / "DamageSeq_minus2_minus1_dipyrimidines.bed",
        qc_path=args.outdir
        / "DamageSeq_minus2_minus1_dipyrimidines.QC.tsv",
    )

    prepare_repair(
        repair_map_path=args.repair_map,
        output_bed=args.outdir
        / "UV_repair_pinpointed_1bp.bed",
        qc_path=args.outdir
        / "UV_repair_pinpointed_1bp.QC.tsv",
    )

    print(
        "Prepared TFBS, damage, and repair event files."
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