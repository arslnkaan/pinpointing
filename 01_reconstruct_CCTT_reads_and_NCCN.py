#!/usr/bin/env python3

import argparse
import csv
import gzip
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, List, Optional, TextIO, Tuple


REQUIRED_EVENT_COLUMNS = {
    "variant_id",
    "original_read_id",
    "chrom",
    "start_0based",
    "end_0based",
    "event_strand",
    "read_strand",
    "read_TT_start_1based",
    "read_TT_end_1based",
    "original_read_dinuc",
    "variant_read_dinuc",
    "genome_dinuc",
    "genome_trinuc_leftC",
    "genome_trinuc_rightC",
    "genome_tetranuc_leftC_plus1",
    "event",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Reconstruct complete original TT-containing reads from a "
            "TT-to-CC rescue BAM and extract canonical NCCN contexts."
        )
    )

    parser.add_argument(
        "--sample",
        required=True,
        help="Sample name.",
    )

    parser.add_argument(
        "--events",
        required=True,
        type=Path,
        help="Accepted CC>TT tandem-event TSV.",
    )

    parser.add_argument(
        "--rescue-bam",
        required=True,
        type=Path,
        help="Sorted BAM containing corrected TT-to-CC rescue alignments.",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        type=Path,
        help="Output directory.",
    )

    parser.add_argument(
        "--mapq",
        type=int,
        default=10,
        help="Minimum alignment MAPQ. Default: 10.",
    )

    return parser.parse_args()


def open_text(path: Path, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)

    return open(path, mode)


def reverse_complement(sequence: str) -> str:
    translation = str.maketrans(
        "ACGTNacgtn",
        "TGCANtgcan",
    )

    return sequence.translate(translation)[::-1]


def cigar_reference_length(cigar: str) -> int:
    if cigar == "*":
        return 0

    reference_length = 0

    for length_text, operation in re.findall(
        r"(\d+)([MIDNSHP=X])",
        cigar,
    ):
        if operation in {"M", "D", "N", "=", "X"}:
            reference_length += int(length_text)

    return reference_length


def parse_sam_tags(
    fields: List[str],
) -> Tuple[Optional[int], Optional[int]]:
    alignment_score = None
    edit_distance = None

    for field in fields:
        if field.startswith("AS:i:"):
            try:
                alignment_score = int(field[5:])
            except ValueError:
                alignment_score = None

        elif field.startswith("NM:i:"):
            try:
                edit_distance = int(field[5:])
            except ValueError:
                edit_distance = None

    return alignment_score, edit_distance


def alignment_rank(record: Dict[str, object]) -> Tuple[int, int, int]:
    mapq = int(record.get("mapq", 0))

    alignment_score = record.get("alignment_score")

    if alignment_score is None:
        alignment_score = -10**9

    edit_distance = record.get("edit_distance")

    if edit_distance is None:
        edit_distance = 10**9

    return (
        mapq,
        int(alignment_score),
        -int(edit_distance),
    )


def load_event_table(
    path: Path,
) -> List[Dict[str, str]]:
    with open_text(path, "rt") as handle:
        reader = csv.DictReader(
            handle,
            delimiter="\t",
        )

        if reader.fieldnames is None:
            raise RuntimeError(
                f"No header was found in {path}"
            )

        missing_columns = (
            REQUIRED_EVENT_COLUMNS -
            set(reader.fieldnames)
        )

        if missing_columns:
            raise RuntimeError(
                "Missing required columns from event table: "
                + ", ".join(sorted(missing_columns))
            )

        events: List[Dict[str, str]] = []

        for row in reader:
            variant_id = row["variant_id"].strip()

            if not variant_id:
                continue

            # Retain only proper TT-to-CC rescue events at genomic CC.
            if row["original_read_dinuc"].upper() != "TT":
                continue

            if row["variant_read_dinuc"].upper() != "CC":
                continue

            if row["genome_dinuc"].upper() != "CC":
                continue

            row["variant_id"] = variant_id
            row["original_read_id"] = (
                row["original_read_id"].strip()
            )

            events.append(row)

    if not events:
        raise RuntimeError(
            f"No valid CC>TT rescue events were found in {path}"
        )

    return events


def load_rescue_alignments(
    bam_path: Path,
    accepted_variant_ids: set,
    minimum_mapq: int,
) -> Dict[str, Dict[str, object]]:
    command = [
        "samtools",
        "view",
        "-F",
        "0x904",
        "-q",
        str(minimum_mapq),
        str(bam_path),
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

    candidates: Dict[
        str,
        List[Dict[str, object]]
    ] = defaultdict(list)

    for line in process.stdout:
        fields = line.rstrip("\n").split("\t")

        if len(fields) < 11:
            continue

        qname = fields[0]

        if qname not in accepted_variant_ids:
            continue

        flag = int(fields[1])
        chrom = fields[2]
        start_1based = int(fields[3])
        mapq = int(fields[4])
        cigar = fields[5]
        sequence = fields[9]
        quality = fields[10]

        if sequence == "*":
            continue

        reference_length = cigar_reference_length(
            cigar
        )

        start_0based = start_1based - 1
        end_0based = start_0based + reference_length

        alignment_score, edit_distance = parse_sam_tags(
            fields[11:]
        )

        candidates[qname].append({
            "variant_id": qname,
            "flag": flag,
            "chrom": chrom,
            "alignment_start_0based": start_0based,
            "alignment_end_0based": end_0based,
            "alignment_start_1based": start_1based,
            "mapq": mapq,
            "cigar": cigar,
            "bam_sequence": sequence,
            "bam_quality": quality,
            "alignment_score": alignment_score,
            "edit_distance": edit_distance,
            "bam_strand": "-" if flag & 16 else "+",
        })

    stderr_text = ""

    if process.stderr is not None:
        stderr_text = process.stderr.read()

    return_code = process.wait()

    if return_code != 0:
        raise RuntimeError(
            "samtools view failed:\n"
            f"{stderr_text}"
        )

    selected: Dict[str, Dict[str, object]] = {}

    for variant_id, records in candidates.items():
        selected[variant_id] = max(
            records,
            key=alignment_rank,
        )

    return selected


def derive_nccn(
    event_row: Dict[str, str],
) -> Tuple[str, str]:
    """
    Return:
      canonical NCCN context
      extraction method
    """

    direct = (
        event_row
        .get("genome_tetranuc_leftC_plus1", "")
        .upper()
        .strip()
    )

    left_trinuc = (
        event_row
        .get("genome_trinuc_leftC", "")
        .upper()
        .strip()
    )

    right_trinuc = (
        event_row
        .get("genome_trinuc_rightC", "")
        .upper()
        .strip()
    )

    candidates: List[Tuple[str, str]] = []

    if direct:
        candidates.append(
            (direct, "genome_tetranuc_leftC_plus1")
        )

    if (
        len(left_trinuc) == 3 and
        len(right_trinuc) == 3
    ):
        reconstructed = (
            left_trinuc +
            right_trinuc[-1]
        )

        candidates.append(
            (
                reconstructed,
                "left_trinuc_plus_right_flank",
            )
        )

    for context, method in candidates:
        if re.fullmatch(r"[ACGT]CC[ACGT]", context):
            return context, method

        context_rc = reverse_complement(context)

        if re.fullmatch(r"[ACGT]CC[ACGT]", context_rc):
            return (
                context_rc,
                method + "_reverse_complemented",
            )

    return "", "invalid_NCCN"


def reconstruct_original_read(
    event_row: Dict[str, str],
    alignment: Dict[str, object],
) -> Dict[str, object]:
    bam_sequence = str(
        alignment["bam_sequence"]
    ).upper()

    bam_quality = str(
        alignment["bam_quality"]
    )

    flag = int(
        alignment["flag"]
    )

    # SAM/BAM stores reverse-strand SEQ reverse-complemented.
    if flag & 16:
        read_oriented_sequence = reverse_complement(
            bam_sequence
        )

        read_oriented_quality = bam_quality[::-1]

        orientation_method = (
            "reverse_complemented_BAM_sequence"
        )

    else:
        read_oriented_sequence = bam_sequence
        read_oriented_quality = bam_quality
        orientation_method = "BAM_sequence_as_stored"

    try:
        tt_start_1based = int(
            event_row["read_TT_start_1based"]
        )

        tt_end_1based = int(
            event_row["read_TT_end_1based"]
        )

    except ValueError:
        return {
            "reconstruction_status": (
                "invalid_read_TT_coordinates"
            ),
            "orientation_method": orientation_method,
            "rescued_read_sequence": (
                read_oriented_sequence
            ),
            "rescued_read_quality": (
                read_oriented_quality
            ),
            "corrected_dinucleotide_observed": "",
            "reconstructed_original_read_sequence": "",
            "reconstructed_original_read_quality": "",
        }

    start0 = tt_start_1based - 1
    end0 = tt_end_1based

    if (
        start0 < 0 or
        end0 > len(read_oriented_sequence) or
        start0 >= end0
    ):
        return {
            "reconstruction_status": (
                "read_TT_coordinates_out_of_range"
            ),
            "orientation_method": orientation_method,
            "rescued_read_sequence": (
                read_oriented_sequence
            ),
            "rescued_read_quality": (
                read_oriented_quality
            ),
            "corrected_dinucleotide_observed": "",
            "reconstructed_original_read_sequence": "",
            "reconstructed_original_read_quality": "",
        }

    observed_dinucleotide = (
        read_oriented_sequence[start0:end0]
    )

    expected_variant = (
        event_row["variant_read_dinuc"]
        .upper()
    )

    original_dinucleotide = (
        event_row["original_read_dinuc"]
        .upper()
    )

    if observed_dinucleotide != expected_variant:
        # Fallback in case the BAM sequence is already stored in
        # original read orientation by a nonstandard conversion.
        fallback_observed = bam_sequence[start0:end0]

        if fallback_observed == expected_variant:
            read_oriented_sequence = bam_sequence
            read_oriented_quality = bam_quality

            observed_dinucleotide = fallback_observed

            orientation_method = (
                "BAM_sequence_as_stored_fallback"
            )

        else:
            return {
                "reconstruction_status": (
                    "corrected_CC_not_found_at_event_position"
                ),
                "orientation_method": orientation_method,
                "rescued_read_sequence": (
                    read_oriented_sequence
                ),
                "rescued_read_quality": (
                    read_oriented_quality
                ),
                "corrected_dinucleotide_observed": (
                    observed_dinucleotide
                ),
                "reconstructed_original_read_sequence": "",
                "reconstructed_original_read_quality": "",
            }

    reconstructed_sequence = (
        read_oriented_sequence[:start0] +
        original_dinucleotide +
        read_oriented_sequence[end0:]
    )

    return {
        "reconstruction_status": "ok",
        "orientation_method": orientation_method,
        "rescued_read_sequence": (
            read_oriented_sequence
        ),
        "rescued_read_quality": (
            read_oriented_quality
        ),
        "corrected_dinucleotide_observed": (
            observed_dinucleotide
        ),
        "reconstructed_original_read_sequence": (
            reconstructed_sequence
        ),
        "reconstructed_original_read_quality": (
            read_oriented_quality
        ),
    }


def event_rank(row: Dict[str, object]) -> Tuple[int, int, int, int]:
    status_score = (
        1
        if row["reconstruction_status"] == "ok"
        else 0
    )

    mapq = int(
        row.get("mapq", 0) or 0
    )

    alignment_score = row.get("alignment_score")

    if alignment_score in {None, ""}:
        alignment_score = -10**9

    edit_distance = row.get("edit_distance")

    if edit_distance in {None, ""}:
        edit_distance = 10**9

    return (
        status_score,
        mapq,
        int(alignment_score),
        -int(edit_distance),
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


def main() -> None:
    args = parse_args()

    if not args.events.exists():
        raise FileNotFoundError(
            f"Missing event table: {args.events}"
        )

    if not args.rescue_bam.exists():
        raise FileNotFoundError(
            f"Missing rescued BAM: {args.rescue_bam}"
        )

    args.outdir.mkdir(
        parents=True,
        exist_ok=True,
    )

    events = load_event_table(
        args.events
    )

    accepted_variant_ids = {
        row["variant_id"]
        for row in events
    }

    print(
        f"Accepted event rows: {len(events):,}"
    )

    print(
        "Accepted rescue variant IDs: "
        f"{len(accepted_variant_ids):,}"
    )

    alignments = load_rescue_alignments(
        bam_path=args.rescue_bam,
        accepted_variant_ids=accepted_variant_ids,
        minimum_mapq=args.mapq,
    )

    print(
        "Accepted variants with BAM alignments: "
        f"{len(alignments):,}"
    )

    all_rows: List[Dict[str, object]] = []

    for event_row in events:
        variant_id = event_row["variant_id"]

        output_row: Dict[str, object] = {
            "sample": args.sample,
            **event_row,
        }

        nccn, nccn_method = derive_nccn(
            event_row
        )

        output_row["genome_NCCN"] = nccn
        output_row["NCCN_extraction_method"] = (
            nccn_method
        )

        alignment = alignments.get(
            variant_id
        )

        if alignment is None:
            output_row.update({
                "bam_flag": "",
                "bam_strand": "",
                "alignment_start_0based": "",
                "alignment_end_0based": "",
                "alignment_start_1based": "",
                "mapq": "",
                "cigar": "",
                "alignment_score": "",
                "edit_distance": "",
                "bam_sequence": "",
                "read_length": "",
                "orientation_method": "",
                "rescued_read_sequence": "",
                "corrected_dinucleotide_observed": "",
                "reconstructed_original_read_sequence": "",
                "reconstruction_status": (
                    "variant_missing_from_rescue_BAM"
                ),
            })

            all_rows.append(output_row)
            continue

        reconstruction = reconstruct_original_read(
            event_row,
            alignment,
        )

        output_row.update({
            "bam_flag": alignment["flag"],
            "bam_strand": alignment["bam_strand"],
            "alignment_start_0based": (
                alignment["alignment_start_0based"]
            ),
            "alignment_end_0based": (
                alignment["alignment_end_0based"]
            ),
            "alignment_start_1based": (
                alignment["alignment_start_1based"]
            ),
            "mapq": alignment["mapq"],
            "cigar": alignment["cigar"],
            "alignment_score": (
                alignment["alignment_score"]
            ),
            "edit_distance": (
                alignment["edit_distance"]
            ),
            "bam_sequence": (
                alignment["bam_sequence"]
            ),
            "read_length": len(
                str(
                    reconstruction[
                        "rescued_read_sequence"
                    ]
                )
            ),
            **reconstruction,
        })

        all_rows.append(output_row)

    # Keep one best accepted event per original read for
    # molecule-level counting and plotting.
    rows_by_original: Dict[
        str,
        List[Dict[str, object]]
    ] = defaultdict(list)

    for row in all_rows:
        rows_by_original[
            str(row["original_read_id"])
        ].append(row)

    best_rows: List[Dict[str, object]] = []

    for original_read_id, candidate_rows in rows_by_original.items():
        best_row = max(
            candidate_rows,
            key=event_rank,
        )

        best_row = dict(best_row)
        best_row["n_candidate_events_for_original_read"] = (
            len(candidate_rows)
        )

        best_rows.append(best_row)

    best_rows.sort(
        key=lambda row: (
            str(row.get("chrom", "")),
            int(row.get("start_0based", 0) or 0),
            str(row.get("original_read_id", "")),
        )
    )

    output_columns = [
        "sample",
        "variant_id",
        "original_read_id",
        "n_candidate_events_for_original_read",

        "chrom",
        "start_0based",
        "end_0based",
        "event_strand",
        "read_strand",

        "read_TT_start_1based",
        "read_TT_end_1based",

        "original_read_dinuc",
        "variant_read_dinuc",
        "genome_dinuc",

        "genome_trinuc_leftC",
        "genome_trinuc_rightC",
        "genome_tetranuc_leftC_plus1",
        "genome_tetranuc_minus1_rightC",

        "genome_NCCN",
        "NCCN_extraction_method",
        "event",

        "bam_flag",
        "bam_strand",
        "alignment_start_0based",
        "alignment_end_0based",
        "alignment_start_1based",
        "mapq",
        "cigar",
        "alignment_score",
        "edit_distance",

        "bam_sequence",
        "orientation_method",
        "rescued_read_sequence",
        "corrected_dinucleotide_observed",
        "reconstructed_original_read_sequence",
        "read_length",
        "reconstruction_status",
    ]

    all_events_file = (
        args.outdir /
        f"{args.sample}.CCTT_all_accepted_events.tsv"
    )

    best_events_file = (
        args.outdir /
        f"{args.sample}.CCTT_best_per_original_read.tsv"
    )

    write_tsv(
        all_events_file,
        all_rows,
        output_columns,
    )

    write_tsv(
        best_events_file,
        best_rows,
        output_columns,
    )

    # --------------------------------------------------------
    # Reconstructed FASTA and FASTQ
    # --------------------------------------------------------

    fasta_file = (
        args.outdir /
        f"{args.sample}.CCTT_reconstructed_full_reads.fasta"
    )

    fastq_file = (
        args.outdir /
        f"{args.sample}.CCTT_reconstructed_full_reads.fastq.gz"
    )

    with open(fasta_file, "w") as fasta_handle, gzip.open(
        fastq_file,
        "wt",
    ) as fastq_handle:

        for row in best_rows:
            if row["reconstruction_status"] != "ok":
                continue

            sequence = str(
                row["reconstructed_original_read_sequence"]
            )

            quality = str(
                row.get(
                    "reconstructed_original_read_quality",
                    "",
                )
            )

            if not quality or quality == "*":
                quality = "I" * len(sequence)

            header = (
                f"{row['original_read_id']}"
                f"|sample={args.sample}"
                f"|NCCN={row['genome_NCCN']}"
                f"|event={row['event']}"
            )

            fasta_handle.write(
                f">{header}\n{sequence}\n"
            )

            fastq_handle.write(
                f"@{header}\n"
                f"{sequence}\n"
                f"+\n"
                f"{quality}\n"
            )

    # --------------------------------------------------------
    # NCCN count table
    # --------------------------------------------------------

    valid_rows = [
        row
        for row in best_rows
        if (
            row["reconstruction_status"] == "ok" and
            re.fullmatch(
                r"[ACGT]CC[ACGT]",
                str(row["genome_NCCN"]),
            )
        )
    ]

    context_counts = Counter(
        str(row["genome_NCCN"])
        for row in valid_rows
    )

    context_order = [
        f"{left}CC{right}"
        for left in "ACGT"
        for right in "ACGT"
    ]

    total_valid = sum(
        context_counts.values()
    )

    count_file = (
        args.outdir /
        f"{args.sample}.CCTT_NCCN_counts.tsv"
    )

    with open(count_file, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "sample",
            "NCCN",
            "count",
            "percent",
        ])

        for context in context_order:
            count = context_counts.get(
                context,
                0,
            )

            percent = (
                100.0 * count / total_valid
                if total_valid > 0
                else 0.0
            )

            writer.writerow([
                args.sample,
                context,
                count,
                percent,
            ])

    # --------------------------------------------------------
    # Summary
    # --------------------------------------------------------

    status_counts = Counter(
        str(row["reconstruction_status"])
        for row in best_rows
    )

    summary_file = (
        args.outdir /
        f"{args.sample}.CCTT_reconstruction_summary.tsv"
    )

    with open(summary_file, "w") as handle:
        writer = csv.writer(
            handle,
            delimiter="\t",
            lineterminator="\n",
        )

        writer.writerow([
            "sample",
            "accepted_event_rows",
            "accepted_variant_ids",
            "original_reads",
            "reconstructed_ok",
            "valid_NCCN",
            "missing_from_BAM",
            "CC_not_found_at_expected_position",
        ])

        writer.writerow([
            args.sample,
            len(events),
            len(accepted_variant_ids),
            len(best_rows),
            status_counts.get("ok", 0),
            len(valid_rows),
            status_counts.get(
                "variant_missing_from_rescue_BAM",
                0,
            ),
            status_counts.get(
                "corrected_CC_not_found_at_event_position",
                0,
            ),
        ])

    print()
    print("Done.")
    print(f"All accepted events: {all_events_file}")
    print(f"Best event per original read: {best_events_file}")
    print(f"Reconstructed FASTA: {fasta_file}")
    print(f"Reconstructed FASTQ: {fastq_file}")
    print(f"NCCN count table: {count_file}")
    print(f"Summary: {summary_file}")
    print()
    print(
        "Successfully reconstructed original reads: "
        f"{status_counts.get('ok', 0):,}"
    )
    print(
        "Valid NCCN contexts: "
        f"{len(valid_rows):,}"
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