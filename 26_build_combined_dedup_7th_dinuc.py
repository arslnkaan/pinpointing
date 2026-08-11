#!/usr/bin/env python3

import argparse
import csv
import re
import sys
from collections import Counter
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple


DIPYRIMIDINES = {
    "CC",
    "CT",
    "TC",
    "TT",
}


@dataclass
class ReadRecord:
    chrom: str
    start: int
    end: int
    read_id: str
    score: int
    strand: str
    sample: str
    timepoint: str
    time_h: str
    source: str
    read_length: int
    original_read_id: str
    priority: int


def is_integer(value: str) -> bool:
    try:
        int(value)
        return True
    except (ValueError, TypeError):
        return False


def safe_integer(
    value: str,
    default: int = 0,
) -> int:
    try:
        return int(float(value))
    except (ValueError, TypeError):
        return default


def record_is_better(
    candidate: ReadRecord,
    current: ReadRecord,
) -> bool:
    """
    Lower priority number wins:
      regular ATL = 0
      rescued     = 1

    Within the same source, higher MAPQ/score wins.
    """

    if candidate.priority != current.priority:
        return candidate.priority < current.priority

    if candidate.score != current.score:
        return candidate.score > current.score

    return candidate.read_id < current.read_id


# ============================================================
# REGULAR ATL DEDUP BED
# ============================================================

def read_regular_bed(
    path: Path,
    sample: str,
    timepoint: str,
    time_h: str,
    qc: Counter,
) -> List[ReadRecord]:

    records: List[ReadRecord] = []

    with open(path, "r") as handle:
        for line in handle:
            if not line.strip():
                continue

            if line.startswith("#"):
                continue

            fields = line.rstrip("\n").split("\t")

            qc["regular_input_rows"] += 1

            if len(fields) < 6:
                qc["regular_malformed_rows"] += 1
                continue

            chrom = fields[0]
            start_value = fields[1]
            end_value = fields[2]
            read_id = fields[3]
            score_value = fields[4]
            strand = fields[5]

            if (
                not is_integer(start_value)
                or not is_integer(end_value)
            ):
                qc["regular_invalid_coordinates"] += 1
                continue

            start = int(start_value)
            end = int(end_value)

            if (
                start < 0
                or end <= start
            ):
                qc["regular_invalid_coordinates"] += 1
                continue

            if strand not in {"+", "-"}:
                qc["regular_invalid_strand"] += 1
                continue

            read_length = end - start

            records.append(
                ReadRecord(
                    chrom=chrom,
                    start=start,
                    end=end,
                    read_id=read_id,
                    score=safe_integer(
                        score_value,
                        default=0,
                    ),
                    strand=strand,
                    sample=sample,
                    timepoint=timepoint,
                    time_h=time_h,
                    source="regular",
                    read_length=read_length,
                    original_read_id=read_id,
                    priority=0,
                )
            )

            qc["regular_selected_rows"] += 1

    return records


# ============================================================
# RESCUED CC>TT FULL READS
# ============================================================

def require_columns(
    fieldnames: List[str],
    required: List[str],
    path: Path,
) -> None:

    missing = [
        column
        for column in required
        if column not in fieldnames
    ]

    if missing:
        raise RuntimeError(
            "Missing columns in rescued table:\n"
            + "\n".join(missing)
            + f"\nFile: {path}"
        )


def read_rescued_table(
    path: Path,
    sample: str,
    timepoint: str,
    time_h: str,
    qc: Counter,
) -> List[ReadRecord]:

    records: List[ReadRecord] = []

    with open(path, "r") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise RuntimeError(
                f"Rescued table has no header: {path}"
            )

        required_columns = [
            "original_read_id",
            "event",
            "genome_dinuc",
            "corrected_dinucleotide_observed",
            "bam_strand",
            "alignment_start_0based",
            "alignment_end_0based",
            "mapq",
            "cigar",
            "read_length",
            "reconstruction_status",
        ]

        require_columns(
            reader.fieldnames,
            required_columns,
            path,
        )

        for row in reader:
            qc["rescued_input_rows"] += 1

            if row["event"].strip() != "CC>TT":
                qc["rescued_non_CCTT_event"] += 1
                continue

            if (
                row["reconstruction_status"]
                .strip()
                .lower()
                != "ok"
            ):
                qc["rescued_reconstruction_not_ok"] += 1
                continue

            if (
                row["genome_dinuc"]
                .strip()
                .upper()
                != "CC"
            ):
                qc["rescued_genome_not_CC"] += 1
                continue

            if (
                row["corrected_dinucleotide_observed"]
                .strip()
                .upper()
                != "CC"
            ):
                qc["rescued_corrected_not_CC"] += 1
                continue

            start_value = row[
                "alignment_start_0based"
            ].strip()

            end_value = row[
                "alignment_end_0based"
            ].strip()

            read_length_value = row[
                "read_length"
            ].strip()

            if (
                not is_integer(start_value)
                or not is_integer(end_value)
                or not is_integer(read_length_value)
            ):
                qc["rescued_invalid_coordinates"] += 1
                continue

            start = int(start_value)
            end = int(end_value)
            read_length = int(read_length_value)

            if (
                start < 0
                or end <= start
                or read_length <= 0
            ):
                qc["rescued_invalid_coordinates"] += 1
                continue

            strand = row[
                "bam_strand"
            ].strip()

            if strand not in {"+", "-"}:
                qc["rescued_invalid_strand"] += 1
                continue

            cigar = row[
                "cigar"
            ].strip()

            cigar_match = re.fullmatch(
                r"(\d+)M",
                cigar,
            )

            if cigar_match is None:
                qc["rescued_non_full_length_M_cigar"] += 1
                continue

            cigar_length = int(
                cigar_match.group(1)
            )

            alignment_span = end - start

            if (
                cigar_length != read_length
                or alignment_span != read_length
            ):
                qc["rescued_length_alignment_mismatch"] += 1
                continue

            original_read_id = row[
                "original_read_id"
            ].strip()

            if original_read_id == "":
                qc["rescued_missing_original_read_id"] += 1
                continue

            records.append(
                ReadRecord(
                    chrom=row["chrom"].strip(),
                    start=start,
                    end=end,
                    read_id=original_read_id,
                    score=safe_integer(
                        row["mapq"],
                        default=0,
                    ),
                    strand=strand,
                    sample=sample,
                    timepoint=timepoint,
                    time_h=time_h,
                    source="CC>TT_rescue",
                    read_length=read_length,
                    original_read_id=original_read_id,
                    priority=1,
                )
            )

            qc["rescued_selected_rows"] += 1

    return records


# ============================================================
# FINAL DEDUPLICATION
# ============================================================

def deduplicate_records(
    records: List[ReadRecord],
    qc: Counter,
) -> List[ReadRecord]:

    # --------------------------------------------------------
    # First: one alignment per original read ID.
    #
    # This protects against duplicated rescue rows and against
    # the same original read occurring in both regular and
    # rescued datasets.
    # --------------------------------------------------------

    by_original_read: Dict[
        Tuple[str, str],
        ReadRecord,
    ] = {}

    for record in records:
        key = (
            record.sample,
            record.original_read_id,
        )

        current = by_original_read.get(
            key
        )

        if current is None:
            by_original_read[key] = record
            continue

        qc["duplicate_original_read_ids_removed"] += 1

        if record_is_better(
            record,
            current,
        ):
            by_original_read[key] = record

    read_deduplicated = list(
        by_original_read.values()
    )

    # --------------------------------------------------------
    # Second: fragment-coordinate deduplication.
    #
    # The available regular file is BED6 and has no retained
    # UMI. Therefore the defensible final fragment key is:
    #
    # chrom + start + end + strand
    # --------------------------------------------------------

    by_fragment: Dict[
        Tuple[str, int, int, str],
        ReadRecord,
    ] = {}

    for record in read_deduplicated:
        key = (
            record.chrom,
            record.start,
            record.end,
            record.strand,
        )

        current = by_fragment.get(
            key
        )

        if current is None:
            by_fragment[key] = record
            continue

        qc["duplicate_fragments_removed"] += 1

        if (
            current.source == "regular"
            and record.source == "CC>TT_rescue"
        ):
            qc[
                "rescued_fragments_removed_due_to_regular_priority"
            ] += 1

        if (
            current.source == "CC>TT_rescue"
            and record.source == "regular"
        ):
            qc[
                "rescued_fragments_replaced_by_regular"
            ] += 1

        if record_is_better(
            record,
            current,
        ):
            by_fragment[key] = record

    final_records = sorted(
        by_fragment.values(),
        key=lambda record: (
            record.chrom,
            record.start,
            record.end,
            record.strand,
            record.original_read_id,
        ),
    )

    qc["records_before_deduplication"] = len(
        records
    )

    qc["records_after_read_id_deduplication"] = len(
        read_deduplicated
    )

    qc["records_after_fragment_deduplication"] = len(
        final_records
    )

    qc["final_regular_records"] = sum(
        record.source == "regular"
        for record in final_records
    )

    qc["final_rescued_records"] = sum(
        record.source == "CC>TT_rescue"
        for record in final_records
    )

    return final_records


# ============================================================
# OUTPUT
# ============================================================

def write_combined_bed(
    records: List[ReadRecord],
    path: Path,
) -> None:

    with open(path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for record in records:
            event_id = (
                f"{record.sample}|"
                f"{record.source}|"
                f"{record.original_read_id}"
            )

            writer.writerow(
                [
                    record.chrom,
                    record.start,
                    record.end,
                    event_id,
                    record.score,
                    record.strand,
                    record.sample,
                    record.timepoint,
                    record.time_h,
                    record.source,
                    record.read_length,
                    record.original_read_id,
                ]
            )


def write_seventh_dinucleotide_bed(
    records: List[ReadRecord],
    path: Path,
    qc: Counter,
) -> None:

    with open(path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for record in records:
            if not (
                26 <= record.read_length <= 30
            ):
                qc[
                    "final_reads_outside_26_to_30_nt"
                ] += 1
                continue

            if record.strand == "+":
                dinuc_start = record.end - 8
                dinuc_end = record.end - 6
            else:
                dinuc_start = record.start + 6
                dinuc_end = record.start + 8

            if (
                dinuc_start < record.start
                or dinuc_end > record.end
                or dinuc_end - dinuc_start != 2
            ):
                qc["invalid_seventh_dinucleotide_interval"] += 1
                continue

            event_id = (
                f"{record.sample}|"
                f"{record.source}|"
                f"{record.original_read_id}"
            )

            writer.writerow(
                [
                    record.chrom,
                    dinuc_start,
                    dinuc_end,
                    event_id,
                    record.score,
                    record.strand,
                    record.sample,
                    record.timepoint,
                    record.time_h,
                    record.source,
                    record.read_length,
                    record.original_read_id,
                ]
            )

            qc["seventh_dinucleotide_intervals_written"] += 1
            qc[
                f"seventh_dinucleotide_source_{record.source}"
            ] += 1


def write_qc(
    qc: Counter,
    path: Path,
) -> None:

    with open(path, "w") as handle:
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


# ============================================================
# FILTER BEDTOOLS GETFASTA -BEDOUT RESULTS
# ============================================================

def filter_dinucleotide_sequences(
    input_path: Path,
    output_path: Path,
    qc_path: Path,
) -> None:

    qc = Counter()

    with open(input_path, "r") as input_handle, open(
        output_path,
        "w",
    ) as output_handle:

        writer = csv.writer(
            output_handle,
            delimiter="\t",
            lineterminator="\n",
        )

        for line in input_handle:
            if not line.strip():
                continue

            fields = line.rstrip("\n").split("\t")

            qc["sequence_input_rows"] += 1

            # Original BED has 12 columns.
            # bedtools getfasta -bedOut appends sequence as 13.
            if len(fields) < 13:
                qc["sequence_malformed_rows"] += 1
                continue

            sequence = fields[-1].strip().upper()

            if not re.fullmatch(
                r"[ACGTN]{2}",
                sequence,
            ):
                qc["sequence_not_2bp"] += 1
                continue

            if sequence not in DIPYRIMIDINES:
                qc["sequence_not_dipyrimidine"] += 1
                continue

            writer.writerow(
                fields[:-1]
                + [
                    sequence,
                ]
            )

            qc["selected_dipyrimidine_events"] += 1
            qc[f"dinucleotide_{sequence}"] += 1

            source = fields[9]
            qc[
                f"selected_source_{source}"
            ] += 1

    write_qc(
        qc,
        qc_path,
    )

    if qc["selected_dipyrimidine_events"] == 0:
        raise RuntimeError(
            "No CC, CT, TC, or TT events were selected."
        )


# ============================================================
# COMMAND-LINE INTERFACE
# ============================================================

def build_parser() -> argparse.ArgumentParser:

    parser = argparse.ArgumentParser(
        description=(
            "Combine regular ATL dedup reads with rescued "
            "CC>TT reads, deduplicate again, and extract the "
            "seventh overlapping dinucleotide from the 3' end."
        )
    )

    subparsers = parser.add_subparsers(
        dest="command",
        required=True,
    )

    combine_parser = subparsers.add_parser(
        "combine",
    )

    combine_parser.add_argument(
        "--sample",
        required=True,
    )

    combine_parser.add_argument(
        "--timepoint",
        required=True,
    )

    combine_parser.add_argument(
        "--time-h",
        required=True,
    )

    combine_parser.add_argument(
        "--regular-bed",
        required=True,
        type=Path,
    )

    combine_parser.add_argument(
        "--rescued-table",
        required=True,
        type=Path,
    )

    combine_parser.add_argument(
        "--combined-bed",
        required=True,
        type=Path,
    )

    combine_parser.add_argument(
        "--dinucleotide-bed",
        required=True,
        type=Path,
    )

    combine_parser.add_argument(
        "--qc",
        required=True,
        type=Path,
    )

    filter_parser = subparsers.add_parser(
        "filter-sequences",
    )

    filter_parser.add_argument(
        "--input",
        required=True,
        type=Path,
    )

    filter_parser.add_argument(
        "--output",
        required=True,
        type=Path,
    )

    filter_parser.add_argument(
        "--qc",
        required=True,
        type=Path,
    )

    return parser


def run_combine(
    args: argparse.Namespace,
) -> None:

    for path in [
        args.regular_bed,
        args.rescued_table,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    args.combined_bed.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.dinucleotide_bed.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    args.qc.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    qc = Counter()

    regular_records = read_regular_bed(
        path=args.regular_bed,
        sample=args.sample,
        timepoint=args.timepoint,
        time_h=args.time_h,
        qc=qc,
    )

    rescued_records = read_rescued_table(
        path=args.rescued_table,
        sample=args.sample,
        timepoint=args.timepoint,
        time_h=args.time_h,
        qc=qc,
    )

    combined_records = (
        regular_records
        + rescued_records
    )

    final_records = deduplicate_records(
        records=combined_records,
        qc=qc,
    )

    write_combined_bed(
        records=final_records,
        path=args.combined_bed,
    )

    write_seventh_dinucleotide_bed(
        records=final_records,
        path=args.dinucleotide_bed,
        qc=qc,
    )

    write_qc(
        qc=qc,
        path=args.qc,
    )

    print(
        f"Combined BED: {args.combined_bed}"
    )

    print(
        f"Seventh dinucleotide BED: "
        f"{args.dinucleotide_bed}"
    )


def main() -> None:

    parser = build_parser()
    args = parser.parse_args()

    if args.command == "combine":
        run_combine(args)

    elif args.command == "filter-sequences":
        filter_dinucleotide_sequences(
            input_path=args.input,
            output_path=args.output,
            qc_path=args.qc,
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