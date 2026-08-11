#!/usr/bin/env python3

import argparse
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import pandas as pd

COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def reverse_complement(seq: str) -> str:
    return str(seq).upper().translate(COMPLEMENT)[::-1]


def parse_changes(changes: str) -> List[str]:
    return [x.strip().upper().replace("U", "T") for x in changes.split(",") if x.strip()]


def change_label(change: str) -> str:
    return change.replace(">", "_to_")


def parse_mismatch_string(s: str) -> List[Tuple[int, str, str, str]]:
    """
    Returns list of: (position_1based, ref, alt, change)
    Input example: 18:G>T,22:G>A
    """
    if pd.isna(s):
        return []
    s = str(s).strip()
    if s == "" or s.lower() == "none":
        return []

    out = []
    for token in s.replace(";", ",").split(","):
        token = token.strip()
        if not token or ":" not in token or ">" not in token:
            continue
        pos_part, change = token.split(":", 1)
        try:
            pos = int(pos_part)
        except ValueError:
            continue
        ref, alt = change.upper().split(">", 1)
        if len(ref) == 1 and len(alt) == 1:
            out.append((pos, ref, alt, f"{ref}>{alt}"))
    return out


def read_mismatch_table(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", dtype=str)

    # Handle older 4-column no-header files if needed.
    expected = {"Read_ID", "Aligned_Sequence", "Reference_Sequence", "Mismatches"}
    if not expected.issubset(set(df.columns)):
        df = pd.read_csv(path, sep="\t", header=None, dtype=str)
        if df.shape[1] == 4:
            df.columns = ["Read_ID", "Aligned_Sequence", "Reference_Sequence", "Mismatches"]
            df["Chromosome"] = "NA"
            df["Start"] = "NA"
            df["End"] = "NA"
            df["Strand"] = "."
            df["Read_Length"] = df["Aligned_Sequence"].str.len().astype(str)
        elif df.shape[1] >= 9:
            df = df.iloc[:, :9].copy()
            df.columns = ["Read_ID", "Chromosome", "Start", "End", "Strand", "Read_Length",
                          "Aligned_Sequence", "Reference_Sequence", "Mismatches"]
        else:
            raise ValueError(f"Could not recognize mismatch table format for {path}")

    for col in ["Read_ID", "Chromosome", "Start", "End", "Strand", "Read_Length",
                "Aligned_Sequence", "Reference_Sequence", "Mismatches"]:
        if col not in df.columns:
            df[col] = "NA"

    df["Aligned_Sequence"] = df["Aligned_Sequence"].astype(str).str.upper()
    df["Reference_Sequence"] = df["Reference_Sequence"].astype(str).str.upper()
    df["Read_Length"] = pd.to_numeric(df["Read_Length"], errors="coerce")
    missing_len = df["Read_Length"].isna()
    df.loc[missing_len, "Read_Length"] = df.loc[missing_len, "Aligned_Sequence"].str.len()
    df["Read_Length"] = df["Read_Length"].astype(int)

    return df


def summarize_size_distribution(df: pd.DataFrame, outpath: Path) -> None:
    size_dist = df["Read_Length"].value_counts().sort_index().reset_index()
    size_dist.columns = ["Read_Length", "Read_Count"]
    size_dist.to_csv(outpath, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(
        description="Filter mismatch events by read length, single-mismatch status, change type, and distance from 3' end."
    )
    parser.add_argument("--mismatches", required=True, help="Mismatch TSV from 01_call_mismatches.py")
    parser.add_argument("--sample", required=True, help="Sample name prefix")
    parser.add_argument("--outdir", default="position_filtered_outputs", help="Output directory")
    parser.add_argument("--min-len", type=int, default=20)
    parser.add_argument("--max-len", type=int, default=30)
    parser.add_argument("--from3-min", type=int, default=6)
    parser.add_argument("--from3-max", type=int, default=13)
    parser.add_argument("--changes", default="G>T,G>A,G>C", help="Comma-separated changes to output, e.g. G>T,G>A,G>C,C>T")
    parser.add_argument("--single-mismatch-only", action="store_true", default=True,
                        help="Keep only reads with exactly one mismatch. Default: true.")
    parser.add_argument("--allow-multi-mismatch", action="store_true",
                        help="Allow reads with multiple mismatches; matching events are output individually.")
    parser.add_argument("--rc-context", action="store_true", default=True,
                        help="Also output reverse-complemented trinucleotide context. Default: true.")
    args = parser.parse_args()

    single_only = not args.allow_multi_mismatch
    changes = parse_changes(args.changes)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    df = read_mismatch_table(args.mismatches)

    # Basic read-length filtering.
    len_df = df[df["Read_Length"].between(args.min_len, args.max_len)].copy()
    len_df.to_csv(outdir / f"{args.sample}_mismatches_{args.min_len}to{args.max_len}mers.tsv", sep="\t", index=False)
    summarize_size_distribution(len_df, outdir / f"{args.sample}_{args.min_len}to{args.max_len}mers_read_length_distribution.tsv")

    events = []
    for _, row in len_df.iterrows():
        parsed = parse_mismatch_string(row["Mismatches"])
        if single_only and len(parsed) != 1:
            continue
        if len(parsed) == 0:
            continue

        read_len = int(row["Read_Length"])
        ref_seq = str(row["Reference_Sequence"]).upper()
        aligned_seq = str(row["Aligned_Sequence"]).upper()

        genomic_positions = []
        if "Genomic_Positions" in row and pd.notna(row["Genomic_Positions"]):
            gp = str(row["Genomic_Positions"])
            if gp.lower() != "none":
                genomic_positions = [x for x in gp.split(",") if x != ""]

        for event_index, (pos, ref, alt, change) in enumerate(parsed):
            if change not in changes:
                continue

            from_3prime = read_len - pos + 1
            if not (args.from3_min <= from_3prime <= args.from3_max):
                continue

            pos0 = pos - 1
            if pos0 <= 0 or pos0 >= len(ref_seq) - 1:
                trinuc = "NA"
                trinuc_rc = "NA"
            else:
                trinuc = ref_seq[pos0 - 1:pos0 + 2].upper()
                trinuc_rc = reverse_complement(trinuc)

            gpos = genomic_positions[event_index] if event_index < len(genomic_positions) else "NA"

            events.append({
                "Read_ID": row["Read_ID"],
                "Chromosome": row["Chromosome"],
                "Start": row["Start"],
                "End": row["End"],
                "Strand": row["Strand"],
                "Read_Length": read_len,
                "Position_from_5prime": pos,
                "Position_from_3prime": from_3prime,
                "Reference_Base": ref,
                "Alt_Base": alt,
                "Change": change,
                "Trinucleotide_Context": trinuc,
                "Trinucleotide_Context_RC": trinuc_rc,
                "Aligned_Sequence": aligned_seq,
                "Reference_Sequence": ref_seq,
                "Genomic_Position_0based": gpos,
                "Original_Mismatches": row["Mismatches"],
            })

    events_df = pd.DataFrame(events)
    all_events_out = outdir / f"{args.sample}_singleMismatch_events_{args.from3_min}to{args.from3_max}nt_from3prime_{args.min_len}to{args.max_len}mers.tsv"
    if len(events_df) == 0:
        # Write an empty table with expected columns.
        events_df = pd.DataFrame(columns=[
            "Read_ID", "Chromosome", "Start", "End", "Strand", "Read_Length",
            "Position_from_5prime", "Position_from_3prime", "Reference_Base", "Alt_Base", "Change",
            "Trinucleotide_Context", "Trinucleotide_Context_RC", "Aligned_Sequence", "Reference_Sequence",
            "Genomic_Position_0based", "Original_Mismatches"
        ])
    events_df.to_csv(all_events_out, sep="\t", index=False)

    # Per-change outputs and summaries.
    summary_rows = []
    for change in changes:
        lab = change_label(change)
        sub = events_df[events_df["Change"] == change].copy()
        prefix = f"{args.sample}_singleMismatch_{lab}_{args.from3_min}to{args.from3_max}nt_from3prime_{args.min_len}to{args.max_len}mers"
        sub.to_csv(outdir / f"{prefix}.tsv", sep="\t", index=False)

        # BED of reads carrying the event.
        if len(sub) > 0 and not (sub["Chromosome"].astype(str) == "NA").all():
            bed = sub[["Chromosome", "Start", "End", "Read_ID", "Change", "Strand"]].copy()
            bed.to_csv(outdir / f"{args.sample}_{lab}_reads.bed", sep="\t", index=False, header=False)

        # Position percentages by read length and 3' position.
        pos_counts = (
            sub.groupby(["Read_Length", "Position_from_3prime"])
            .size()
            .reset_index(name="Event_Count")
        )
        denom = (
            len_df.groupby("Read_Length")
            .size()
            .reset_index(name="Total_Reads_In_Length")
        )
        pos_pct = pos_counts.merge(denom, on="Read_Length", how="left")
        if len(pos_pct) > 0:
            pos_pct["Percentage_of_Reads"] = 100 * pos_pct["Event_Count"] / pos_pct["Total_Reads_In_Length"]
        pos_pct.to_csv(outdir / f"{prefix}_position_percentages.tsv", sep="\t", index=False)

        # Trinucleotide context percentage. Uses RC context for G>X to show C>N damaged-base space.
        if len(sub) > 0:
            ctx_col = "Trinucleotide_Context_RC" if args.rc_context else "Trinucleotide_Context"
            ctx_counts = sub[ctx_col].value_counts().sort_index().reset_index()
            ctx_counts.columns = [ctx_col, "Count"]
            total = ctx_counts["Count"].sum()
            ctx_counts["Percentage"] = 100 * ctx_counts["Count"] / total if total > 0 else 0
        else:
            ctx_col = "Trinucleotide_Context_RC" if args.rc_context else "Trinucleotide_Context"
            ctx_counts = pd.DataFrame(columns=[ctx_col, "Count", "Percentage"])
        ctx_counts.to_csv(outdir / f"{prefix}_trinucleotide_percentages.csv", index=False)

        summary_rows.append({
            "Sample": args.sample,
            "Change": change,
            "Read_Length_Min": args.min_len,
            "Read_Length_Max": args.max_len,
            "From3_Min": args.from3_min,
            "From3_Max": args.from3_max,
            "Total_Reads_Length_Filtered": len(len_df),
            "Events_Passing_Filters": len(sub),
            "Unique_Reads_Passing_Filters": sub["Read_ID"].nunique() if len(sub) > 0 else 0,
        })

    pd.DataFrame(summary_rows).to_csv(outdir / f"{args.sample}_mismatch_filter_summary.tsv", sep="\t", index=False)

    print(f"Done: {args.sample}")
    print(f"Filtered events: {all_events_out}")
    print(f"Summary: {outdir / (args.sample + '_mismatch_filter_summary.tsv')}")


if __name__ == "__main__":
    main()
