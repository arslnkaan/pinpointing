#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
import sys
from pathlib import Path
from typing import Dict, Optional


CIGAR_PATTERN = re.compile(r"(\d+)([MIDNSHP=X])")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Extract accepted CC>TT rescued alignments from a rescue BAM "
            "and write their complete genomic intervals as BED6."
        )
    )

    parser.add_argument(
        "--best-events",
        required=True,
        type=Path,
        help="CCTT_best_per_original_read.tsv file.",
    )

    parser.add_argument(
        "--bam",
        required=True,
        type=Path,
        help="Sorted rescue BAM.",
    )

    parser.add_argument(
        "--output-bed",
        required=True,
        type=Path,
        help="Output BED6 file.",
    )

    parser.add_argument(
        "--summary",
        required=True,
        type=Path,
        help="Output extraction summary.",
    )

    return parser.parse_args()


def find_column(
    header,
    candidates,
    required=True,
) -> Optional[str]:
    lower_map = {
        name.lower(): name
        for name in header
    }

    for candidate in candidates:
        if candidate.lower() in lower_map:
            return lower_map[candidate.lower()]

    if required:
        raise RuntimeError(
            "Could not identify required column. "
            f"Tried: {', '.join(candidates)}. "
            f"Available: {', '.join(header)}"
        )

    return None


def reference_span(cigar: str) -> int:
    """
    CIGAR operations consuming reference:
    M, D, N, =, X
    """

    if cigar == "*":
        return 0

    span = 0

    for length_text, operation in CIGAR_PATTERN.findall(cigar):
        if operation in {"M", "D", "N", "=", "X"}:
            span += int(length_text)

    return span


def main() -> None:
    args = parse_args()

    for path in [
        args.best_events,
        args.bam,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing input: {path}"
            )

    args.output_bed.parent.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # READ ACCEPTED RESCUE IDS
    # ========================================================

    accepted: Dict[str, str] = {}

    with open(args.best_events) as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise RuntimeError(
                "Best-event table has no header."
            )

        variant_column = find_column(
            reader.fieldnames,
            [
                "variant_id",
                "bam_qname",
                "qname",
                "rescued_read_id",
            ],
        )

        original_column = find_column(
            reader.fieldnames,
            [
                "original_read_id",
                "read_id",
                "original_qname",
            ],
        )

        status_column = find_column(
            reader.fieldnames,
            [
                "reconstruction_status",
                "status",
            ],
            required=False,
        )

        rows_total = 0
        rows_accepted = 0

        for row in reader:
            rows_total += 1

            if status_column is not None:
                status = str(
                    row.get(
                        status_column,
                        "",
                    )
                ).strip()

                if status != "ok":
                    continue

            variant_id = str(
                row.get(
                    variant_column,
                    "",
                )
            ).strip()

            original_id = str(
                row.get(
                    original_column,
                    "",
                )
            ).strip()

            if not variant_id or not original_id:
                continue

            accepted[variant_id] = original_id
            rows_accepted += 1

    if not accepted:
        raise RuntimeError(
            "No accepted rescued reads were found."
        )

    # ========================================================
    # STREAM BAM THROUGH SAMTOOLS
    # ========================================================

    command = [
        "samtools",
        "view",
        "-F",
        "0x904",
        str(args.bam),
    ]

    process = subprocess.Popen(
        command,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
    )

    if process.stdout is None:
        raise RuntimeError(
            "Could not read samtools output."
        )

    bam_alignments_scanned = 0
    accepted_alignments_written = 0
    accepted_ids_found = set()

    with open(args.output_bed, "w") as output:
        writer = csv.writer(
            output,
            delimiter="\t",
            lineterminator="\n",
        )

        for line in process.stdout:
            fields = line.rstrip("\n").split("\t")

            if len(fields) < 11:
                continue

            bam_alignments_scanned += 1

            qname = fields[0]

            if qname not in accepted:
                continue

            flag = int(fields[1])
            chrom = fields[2]
            position_1based = int(fields[3])
            mapq = fields[4]
            cigar = fields[5]

            span = reference_span(cigar)

            if span <= 0:
                continue

            start_0based = position_1based - 1
            end_0based = start_0based + span

            strand = (
                "-"
                if flag & 16
                else "+"
            )

            original_id = accepted[qname]

            writer.writerow([
                chrom,
                start_0based,
                end_0based,
                original_id,
                mapq,
                strand,
            ])

            accepted_alignments_written += 1
            accepted_ids_found.add(qname)

    stderr = ""

    if process.stderr is not None:
        stderr = process.stderr.read()

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            "samtools view failed:\n"
            f"{stderr}"
        )

    # ========================================================
    # SUMMARY
    # ========================================================

    with open(args.summary, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "best_table_rows",
            "accepted_best_rows",
            "unique_accepted_variant_ids",
            "bam_alignments_scanned",
            "accepted_alignments_written",
            "accepted_ids_missing_from_bam",
        ])

        writer.writerow([
            rows_total,
            rows_accepted,
            len(accepted),
            bam_alignments_scanned,
            accepted_alignments_written,
            len(accepted) -
            len(accepted_ids_found),
        ])

    print(
        f"Accepted rescued BED rows: "
        f"{accepted_alignments_written:,}"
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