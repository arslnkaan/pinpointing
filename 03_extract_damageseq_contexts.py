#!/usr/bin/env python3

import argparse
import csv
import re
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from typing import Dict, List, Tuple


# Expected record-name structure:
#
# SRR5461431.36652900/1_TC_Ccenter::chr1:31963-31966(+)
# SRR5461431.12345678/1_CT_Ccenter::chr2:100-103(-)
# SRR5461431.12345679/1_CC_Ccenter::chr3:200-203(+)

NAME_PATTERN = re.compile(
    r"^(?P<read_id>.+)_(?P<lesion_type>CT|TC|CC)_Ccenter"
    r"::(?P<chrom>.+):(?P<start>\d+)-(?P<end>\d+)"
    r"\((?P<strand>[+-])\)$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Read headerless Damage-seq sequence tables and generate "
            "CT:NCT, TC:TCN and CC:NCCN context tables."
        )
    )

    parser.add_argument(
        "--ct-tc",
        required=True,
        type=Path,
        help=(
            "Headerless CT/TC sequence table with two columns: "
            "record name and trinucleotide sequence."
        ),
    )

    parser.add_argument(
        "--cc",
        required=True,
        type=Path,
        help=(
            "Headerless CC sequence table with two columns: "
            "record name and C-centered trinucleotide sequence."
        ),
    )

    parser.add_argument(
        "--fasta",
        required=True,
        type=Path,
        help="Reference-genome FASTA used to extend CC to NCCN.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--cc-center-mode",
        choices=["auto", "leftC", "rightC"],
        default="auto",
        help=(
            "How the centered C relates to the CC lesion. "
            "'leftC' means the centered C is the first C, producing NCC. "
            "'rightC' means the centered C is the second C, producing CCN. "
            "'auto' determines this from the trinucleotide and treats CCC "
            "as NCC/leftC. Default: auto."
        ),
    )

    return parser.parse_args()


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )
    return sequence.translate(translation)[::-1]


def parse_record_name(name: str) -> Dict[str, object]:
    match = NAME_PATTERN.match(name)

    if match is None:
        raise ValueError(
            "Could not parse record name: "
            f"{name}"
        )

    result = match.groupdict()

    return {
        "read_id": result["read_id"],
        "lesion_type": result["lesion_type"],
        "chrom": result["chrom"],
        "window_start_0based": int(result["start"]),
        "window_end_0based": int(result["end"]),
        "strand": result["strand"],
    }


def read_sequence_table(
    path: Path,
    expected_lesions: set,
) -> Tuple[List[Dict[str, object]], List[Dict[str, str]]]:

    valid_rows: List[Dict[str, object]] = []
    rejected_rows: List[Dict[str, str]] = []

    with open(path) as handle:
        reader = csv.reader(
            handle,
            delimiter="\t",
        )

        for line_number, fields in enumerate(
            reader,
            start=1,
        ):
            if not fields:
                continue

            if len(fields) < 2:
                rejected_rows.append({
                    "line_number": str(line_number),
                    "record_name": fields[0] if fields else "",
                    "sequence": "",
                    "reason": "fewer_than_two_columns",
                })
                continue

            record_name = fields[0].strip()
            sequence = fields[1].strip().upper()

            try:
                parsed = parse_record_name(
                    record_name
                )
            except ValueError:
                rejected_rows.append({
                    "line_number": str(line_number),
                    "record_name": record_name,
                    "sequence": sequence,
                    "reason": "record_name_not_parsed",
                })
                continue

            lesion_type = str(
                parsed["lesion_type"]
            )

            if lesion_type not in expected_lesions:
                rejected_rows.append({
                    "line_number": str(line_number),
                    "record_name": record_name,
                    "sequence": sequence,
                    "reason": (
                        "unexpected_lesion_type_"
                        f"{lesion_type}"
                    ),
                })
                continue

            expected_length = (
                int(parsed["window_end_0based"]) -
                int(parsed["window_start_0based"])
            )

            if len(sequence) != expected_length:
                rejected_rows.append({
                    "line_number": str(line_number),
                    "record_name": record_name,
                    "sequence": sequence,
                    "reason": (
                        "sequence_length_does_not_match_interval"
                    ),
                })
                continue

            valid_rows.append({
                "line_number": line_number,
                "record_name": record_name,
                "sequence": sequence,
                **parsed,
            })

    return valid_rows, rejected_rows


def process_ct_tc(
    rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:

    output: List[Dict[str, object]] = []

    for index, row in enumerate(rows, start=1):
        lesion_type = str(
            row["lesion_type"]
        )

        sequence = str(
            row["sequence"]
        ).upper()

        if lesion_type == "CT":
            expected = bool(
                re.fullmatch(
                    r"[ACGT]CT",
                    sequence,
                )
            )
            context_class = "NCT"

        elif lesion_type == "TC":
            expected = bool(
                re.fullmatch(
                    r"TC[ACGT]",
                    sequence,
                )
            )
            context_class = "TCN"

        else:
            expected = False
            context_class = ""

        output.append({
            "record_id": f"cttc_{index}",
            "original_record_name": row["record_name"],
            "read_id": row["read_id"],
            "chrom": row["chrom"],
            "window_start_0based": row["window_start_0based"],
            "window_end_0based": row["window_end_0based"],
            "strand": row["strand"],
            "lesion_type": lesion_type,
            "context_class": context_class,
            "trinucleotide": sequence,
            "context_status": (
                "ok"
                if expected
                else "unexpected_reference_context"
            ),
        })

    return output


def determine_cc_orientation(
    sequence: str,
    cc_center_mode: str,
) -> Tuple[str, str]:

    sequence = sequence.upper()

    if cc_center_mode == "leftC":
        return "NCC", "forced_leftC"

    if cc_center_mode == "rightC":
        return "CCN", "forced_rightC"

    # Auto mode
    #
    # GCC, ACC, TCC -> NCC
    # CCA, CCG, CCT -> CCN
    #
    # CCC matches both possibilities. Based on the usual Damage-seq
    # convention in which the called C is followed by the second
    # nucleotide of the dimer, CCC is treated as NCC.

    if sequence == "CCC":
        return "NCC", "auto_CCC_assumed_leftC"

    if re.fullmatch(r"[AGT]CC", sequence):
        return "NCC", "auto_from_NCC"

    if re.fullmatch(r"CC[AGT]", sequence):
        return "CCN", "auto_from_CCN"

    return "", "could_not_locate_CC_in_trinucleotide"


def prepare_cc_windows(
    rows: List[Dict[str, object]],
    cc_center_mode: str,
) -> Tuple[
    List[Dict[str, object]],
    List[Dict[str, object]],
]:

    fetchable: List[Dict[str, object]] = []
    rejected: List[Dict[str, object]] = []

    for index, row in enumerate(rows, start=1):
        sequence = str(
            row["sequence"]
        ).upper()

        orientation, orientation_method = (
            determine_cc_orientation(
                sequence,
                cc_center_mode,
            )
        )

        base_record = {
            "record_id": f"cc_{index}",
            "original_record_name": row["record_name"],
            "read_id": row["read_id"],
            "chrom": row["chrom"],
            "original_window_start_0based":
                row["window_start_0based"],
            "original_window_end_0based":
                row["window_end_0based"],
            "strand": row["strand"],
            "lesion_type": "CC",
            "original_trinucleotide": sequence,
            "CC_orientation": orientation,
            "CC_orientation_method": orientation_method,
        }

        if orientation == "":
            base_record.update({
                "lesion_start_0based": "",
                "lesion_end_0based": "",
                "NCCN_window_start_0based": "",
                "NCCN_window_end_0based": "",
                "NCCN": "",
                "context_status":
                    "could_not_locate_CC_in_trinucleotide",
            })

            rejected.append(base_record)
            continue

        start = int(
            row["window_start_0based"]
        )

        end = int(
            row["window_end_0based"]
        )

        strand = str(
            row["strand"]
        )

        if orientation == "NCC":
            # Oriented trinucleotide:
            #
            # 5' - N C C - 3'
            #
            # Add one nucleotide on the oriented 3' side.

            if strand == "+":
                nccn_start = start
                nccn_end = end + 1

                lesion_start = start + 1
                lesion_end = end

            else:
                nccn_start = start - 1
                nccn_end = end

                lesion_start = start
                lesion_end = end - 1

        else:
            # Oriented trinucleotide:
            #
            # 5' - C C N - 3'
            #
            # Add one nucleotide on the oriented 5' side.

            if strand == "+":
                nccn_start = start - 1
                nccn_end = end

                lesion_start = start
                lesion_end = end - 1

            else:
                nccn_start = start
                nccn_end = end + 1

                lesion_start = start + 1
                lesion_end = end

        if nccn_start < 0:
            base_record.update({
                "lesion_start_0based": lesion_start,
                "lesion_end_0based": lesion_end,
                "NCCN_window_start_0based": nccn_start,
                "NCCN_window_end_0based": nccn_end,
                "NCCN": "",
                "context_status":
                    "NCCN_window_before_chromosome_start",
            })

            rejected.append(base_record)
            continue

        base_record.update({
            "lesion_start_0based": lesion_start,
            "lesion_end_0based": lesion_end,
            "NCCN_window_start_0based": nccn_start,
            "NCCN_window_end_0based": nccn_end,
        })

        fetchable.append(base_record)

    return fetchable, rejected


def fetch_nccn_sequences(
    records: List[Dict[str, object]],
    fasta: Path,
) -> Dict[str, str]:

    if not records:
        return {}

    with tempfile.TemporaryDirectory() as temporary_directory:
        temporary_directory_path = Path(
            temporary_directory
        )

        bed_file = (
            temporary_directory_path /
            "damage_CC_NCCN_windows.bed"
        )

        with open(bed_file, "w") as handle:
            writer = csv.writer(
                handle,
                delimiter="\t",
                lineterminator="\n",
            )

            for row in records:
                writer.writerow([
                    row["chrom"],
                    row["NCCN_window_start_0based"],
                    row["NCCN_window_end_0based"],
                    row["record_id"],
                    0,
                    row["strand"],
                ])

        command = [
            "bedtools",
            "getfasta",
            "-fi",
            str(fasta),
            "-bed",
            str(bed_file),
            "-s",
            "-name",
            "-tab",
        ]

        result = subprocess.run(
            command,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            text=True,
        )

        if result.returncode != 0:
            raise RuntimeError(
                "bedtools getfasta failed:\n"
                f"{result.stderr}"
            )

        sequences: Dict[str, str] = {}

        for line in result.stdout.splitlines():
            fields = line.split("\t")

            if len(fields) != 2:
                continue

            record_id = fields[0].split(
                "::",
                1,
            )[0]

            sequences[record_id] = (
                fields[1]
                .strip()
                .upper()
            )

        return sequences


def finish_cc_records(
    records: List[Dict[str, object]],
    sequences: Dict[str, str],
) -> List[Dict[str, object]]:

    output: List[Dict[str, object]] = []

    for row in records:
        result = dict(row)

        sequence = sequences.get(
            str(row["record_id"]),
            "",
        ).upper()

        result["NCCN"] = sequence

        if re.fullmatch(
            r"[ACGT]CC[ACGT]",
            sequence,
        ):
            result["context_status"] = "ok"
        elif sequence == "":
            result["context_status"] = (
                "sequence_not_returned_by_getfasta"
            )
        else:
            result["context_status"] = (
                "unexpected_reference_context"
            )

        output.append(result)

    return output


def deduplicate_cc_events(
    rows: List[Dict[str, object]],
) -> List[Dict[str, object]]:
    """
    If the same Damage-seq read was represented by both an NCC-centered
    and CCN-centered row, retain only one lesion event.
    """

    best_by_event: Dict[
        Tuple[str, str, object, object, str],
        Dict[str, object]
    ] = {}

    for row in rows:
        event_key = (
            str(row.get("read_id", "")),
            str(row.get("chrom", "")),
            row.get("lesion_start_0based", ""),
            row.get("lesion_end_0based", ""),
            str(row.get("strand", "")),
        )

        existing = best_by_event.get(
            event_key
        )

        if existing is None:
            best_by_event[event_key] = row
            continue

        existing_ok = (
            existing.get("context_status") == "ok"
        )

        current_ok = (
            row.get("context_status") == "ok"
        )

        if current_ok and not existing_ok:
            best_by_event[event_key] = row

    return list(
        best_by_event.values()
    )


def write_tsv(
    path: Path,
    rows: List[Dict[str, object]],
    columns: List[str],
) -> None:

    with open(path, "w") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=columns,
            delimiter="\t",
            lineterminator="\n",
            extrasaction="ignore",
        )

        writer.writeheader()
        writer.writerows(rows)


def write_rejected_input_rows(
    path: Path,
    rejected_rows: List[Dict[str, str]],
) -> None:

    columns = [
        "line_number",
        "record_name",
        "sequence",
        "reason",
    ]

    write_tsv(
        path,
        rejected_rows,
        columns,
    )


def write_ct_tc_counts(
    path: Path,
    rows: List[Dict[str, object]],
) -> None:

    context_order = [
        ("CT", "ACT"),
        ("CT", "CCT"),
        ("CT", "GCT"),
        ("CT", "TCT"),
        ("TC", "TCA"),
        ("TC", "TCC"),
        ("TC", "TCG"),
        ("TC", "TCT"),
    ]

    counts = Counter(
        (
            str(row["lesion_type"]),
            str(row["trinucleotide"]),
        )
        for row in rows
        if row["context_status"] == "ok"
    )

    lesion_totals = Counter(
        str(row["lesion_type"])
        for row in rows
        if row["context_status"] == "ok"
    )

    total_damage = sum(
        lesion_totals.values()
    )

    with open(path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "lesion_type",
            "context_class",
            "trinucleotide",
            "count",
            "lesion_total",
            "total_CT_TC_damage",
            "percent_within_lesion_type",
            "percent_of_all_CT_TC_damage",
        ])

        for lesion_type, context in context_order:
            count = counts.get(
                (lesion_type, context),
                0,
            )

            lesion_total = lesion_totals.get(
                lesion_type,
                0,
            )

            within_lesion_percent = (
                100 * count / lesion_total
                if lesion_total > 0
                else 0
            )

            all_damage_percent = (
                100 * count / total_damage
                if total_damage > 0
                else 0
            )

            writer.writerow([
                lesion_type,
                "NCT" if lesion_type == "CT" else "TCN",
                context,
                count,
                lesion_total,
                total_damage,
                within_lesion_percent,
                all_damage_percent,
            ])


def write_cc_counts(
    path: Path,
    rows: List[Dict[str, object]],
) -> None:

    context_order = [
        f"{left}CC{right}"
        for left in "ACGT"
        for right in "ACGT"
    ]

    counts = Counter(
        str(row["NCCN"])
        for row in rows
        if row["context_status"] == "ok"
    )

    total = sum(
        counts.values()
    )

    with open(path, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "NCCN",
            "count",
            "total_CC_damage",
            "percent",
        ])

        for context in context_order:
            count = counts.get(
                context,
                0,
            )

            percent = (
                100 * count / total
                if total > 0
                else 0
            )

            writer.writerow([
                context,
                count,
                total,
                percent,
            ])


def main() -> None:
    args = parse_args()

    for path in [
        args.ct_tc,
        args.cc,
        args.fasta,
    ]:
        if not path.exists():
            raise FileNotFoundError(
                f"Missing input file: {path}"
            )

        if path.stat().st_size == 0:
            raise RuntimeError(
                f"Input file is empty: {path}"
            )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    # ========================================================
    # CT / TC
    # ========================================================

    ct_tc_input_rows, ct_tc_rejected = (
        read_sequence_table(
            args.ct_tc,
            expected_lesions={"CT", "TC"},
        )
    )

    ct_tc_events = process_ct_tc(
        ct_tc_input_rows
    )

    # ========================================================
    # CC
    # ========================================================

    cc_input_rows, cc_rejected = (
        read_sequence_table(
            args.cc,
            expected_lesions={"CC"},
        )
    )

    cc_fetchable, cc_context_rejected = (
        prepare_cc_windows(
            cc_input_rows,
            args.cc_center_mode,
        )
    )

    cc_sequences = fetch_nccn_sequences(
        cc_fetchable,
        args.fasta,
    )

    cc_events = finish_cc_records(
        cc_fetchable,
        cc_sequences,
    )

    cc_events.extend(
        cc_context_rejected
    )

    cc_events = deduplicate_cc_events(
        cc_events
    )

    # ========================================================
    # OUTPUT PATHS
    # ========================================================

    ct_tc_event_file = (
        args.outdir /
        "DamageSeq_CT_TC_trinucleotide_events.tsv"
    )

    ct_tc_count_file = (
        args.outdir /
        "DamageSeq_CT_TC_trinucleotide_counts.tsv"
    )

    cc_event_file = (
        args.outdir /
        "DamageSeq_CC_NCCN_events.tsv"
    )

    cc_count_file = (
        args.outdir /
        "DamageSeq_CC_NCCN_counts.tsv"
    )

    rejected_ct_tc_file = (
        args.outdir /
        "DamageSeq_CT_TC_rejected_input_rows.tsv"
    )

    rejected_cc_file = (
        args.outdir /
        "DamageSeq_CC_rejected_input_rows.tsv"
    )

    # ========================================================
    # WRITE EVENT TABLES
    # ========================================================

    write_tsv(
        ct_tc_event_file,
        ct_tc_events,
        [
            "record_id",
            "original_record_name",
            "read_id",
            "chrom",
            "window_start_0based",
            "window_end_0based",
            "strand",
            "lesion_type",
            "context_class",
            "trinucleotide",
            "context_status",
        ],
    )

    write_tsv(
        cc_event_file,
        cc_events,
        [
            "record_id",
            "original_record_name",
            "read_id",
            "chrom",
            "original_window_start_0based",
            "original_window_end_0based",
            "strand",
            "lesion_type",
            "original_trinucleotide",
            "CC_orientation",
            "CC_orientation_method",
            "lesion_start_0based",
            "lesion_end_0based",
            "NCCN_window_start_0based",
            "NCCN_window_end_0based",
            "NCCN",
            "context_status",
        ],
    )

    write_rejected_input_rows(
        rejected_ct_tc_file,
        ct_tc_rejected,
    )

    write_rejected_input_rows(
        rejected_cc_file,
        cc_rejected,
    )

    # ========================================================
    # WRITE COUNT TABLES
    # ========================================================

    write_ct_tc_counts(
        ct_tc_count_file,
        ct_tc_events,
    )

    write_cc_counts(
        cc_count_file,
        cc_events,
    )

    # ========================================================
    # SUMMARY
    # ========================================================

    ct_tc_valid = sum(
        row["context_status"] == "ok"
        for row in ct_tc_events
    )

    cc_valid = sum(
        row["context_status"] == "ok"
        for row in cc_events
    )

    summary_file = (
        args.outdir /
        "DamageSeq_context_extraction_summary.tsv"
    )

    with open(summary_file, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "dataset",
            "input_rows_parsed",
            "input_rows_rejected",
            "event_rows_output",
            "valid_context_rows",
        ])

        writer.writerow([
            "CT_TC",
            len(ct_tc_input_rows),
            len(ct_tc_rejected),
            len(ct_tc_events),
            ct_tc_valid,
        ])

        writer.writerow([
            "CC",
            len(cc_input_rows),
            len(cc_rejected),
            len(cc_events),
            cc_valid,
        ])

    print("Done.")
    print(
        f"CT/TC parsed rows: {len(ct_tc_input_rows):,}"
    )
    print(
        f"CT/TC valid contexts: {ct_tc_valid:,}"
    )
    print(
        f"CC parsed rows: {len(cc_input_rows):,}"
    )
    print(
        f"CC valid NCCN contexts: {cc_valid:,}"
    )
    print()
    print(f"CT/TC events: {ct_tc_event_file}")
    print(f"CT/TC counts: {ct_tc_count_file}")
    print(f"CC events: {cc_event_file}")
    print(f"CC counts: {cc_count_file}")
    print(f"Summary: {summary_file}")


if __name__ == "__main__":
    try:
        main()

    except Exception as error:
        print(
            f"ERROR: {error}",
            file=sys.stderr,
        )
        sys.exit(1)