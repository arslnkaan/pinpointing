#!/usr/bin/env python3

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd

DEFAULT_CONTEXT_ORDER = [
    "ACA", "ACC", "ACG", "ACT",
    "CCA", "CCC", "CCG", "CCT",
    "GCA", "GCC", "GCG", "GCT",
    "TCA", "TCC", "TCG", "TCT",
]


def change_label(change: str) -> str:
    return change.replace(">", "_to_")


def parse_changes(changes: str):
    return [x.strip().upper() for x in changes.split(",") if x.strip()]


def plot_size_distribution(path: Path, sample: str, outdir: Path):
    if not path.exists():
        return
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return

    plt.figure(figsize=(6, 4))
    plt.bar(df["Read_Length"].astype(str), df["Read_Count"])
    plt.xlabel("Read length")
    plt.ylabel("Read count")
    plt.title(f"{sample}: read length distribution")
    plt.tight_layout()
    plt.savefig(outdir / f"{sample}_read_length_distribution.png", dpi=300)
    plt.savefig(outdir / f"{sample}_read_length_distribution.pdf")
    plt.close()


def plot_contexts(path: Path, sample: str, change: str, outdir: Path, context_order):
    if not path.exists():
        print(f"Missing context file, skipping: {path}")
        return
    df = pd.read_csv(path)
    if df.empty:
        return

    ctx_col = "Trinucleotide_Context_RC" if "Trinucleotide_Context_RC" in df.columns else "Trinucleotide_Context"
    df = df.rename(columns={ctx_col: "context"})
    df["context"] = df["context"].astype(str).str.upper()

    # Use the canonical C-centered order when possible; otherwise use observed order.
    observed = set(df["context"])
    if any(c in observed for c in context_order):
        full = pd.DataFrame({"context": context_order})
        df = full.merge(df, on="context", how="left").fillna({"Count": 0, "Percentage": 0})
    else:
        df = df.sort_values("Percentage", ascending=False)

    plt.figure(figsize=(9.5, 5))
    plt.bar(df["context"], df["Percentage"], edgecolor="black", linewidth=0.5)
    plt.xlabel("Trinucleotide context")
    plt.ylabel("Percentage of events")
    plt.title(f"{sample}: {change} trinucleotide context")
    plt.xticks(rotation=45, ha="right")
    plt.tight_layout()
    lab = change_label(change)
    plt.savefig(outdir / f"{sample}_{lab}_trinucleotide_percentages.png", dpi=300)
    plt.savefig(outdir / f"{sample}_{lab}_trinucleotide_percentages.pdf")
    plt.close()


def plot_position_percentages(path: Path, sample: str, change: str, outdir: Path):
    if not path.exists():
        print(f"Missing position file, skipping: {path}")
        return
    df = pd.read_csv(path, sep="\t")
    if df.empty:
        return

    # Plot each read length separately on same axis.
    plt.figure(figsize=(7, 4.5))
    for read_len, sub in df.groupby("Read_Length"):
        sub = sub.sort_values("Position_from_3prime")
        plt.plot(
            sub["Position_from_3prime"],
            sub["Percentage_of_Reads"],
            marker="o",
            linewidth=1,
            label=f"{int(read_len)} nt"
        )
    plt.xlabel("Mismatch position from 3' end")
    plt.ylabel("Reads with event (%)")
    plt.title(f"{sample}: {change} by position from 3' end")
    plt.legend(title="Read length", fontsize=8, frameon=False)
    plt.tight_layout()
    lab = change_label(change)
    plt.savefig(outdir / f"{sample}_{lab}_position_from3_percentages.png", dpi=300)
    plt.savefig(outdir / f"{sample}_{lab}_position_from3_percentages.pdf")
    plt.close()


def main():
    parser = argparse.ArgumentParser(description="Plot outputs from 02_filter_mismatch_events.py")
    parser.add_argument("--sample", required=True)
    parser.add_argument("--indir", required=True, help="Directory containing filtered mismatch outputs")
    parser.add_argument("--outdir", default=None, help="Plot output directory. Default: indir/plots")
    parser.add_argument("--min-len", type=int, default=20)
    parser.add_argument("--max-len", type=int, default=30)
    parser.add_argument("--from3-min", type=int, default=6)
    parser.add_argument("--from3-max", type=int, default=13)
    parser.add_argument("--changes", default="G>T,G>A,G>C")
    args = parser.parse_args()

    indir = Path(args.indir)
    outdir = Path(args.outdir) if args.outdir else indir / "plots"
    outdir.mkdir(parents=True, exist_ok=True)

    size_path = indir / f"{args.sample}_{args.min_len}to{args.max_len}mers_read_length_distribution.tsv"
    plot_size_distribution(size_path, args.sample, outdir)

    for change in parse_changes(args.changes):
        lab = change_label(change)
        prefix = f"{args.sample}_singleMismatch_{lab}_{args.from3_min}to{args.from3_max}nt_from3prime_{args.min_len}to{args.max_len}mers"
        plot_contexts(indir / f"{prefix}_trinucleotide_percentages.csv", args.sample, change, outdir, DEFAULT_CONTEXT_ORDER)
        plot_position_percentages(indir / f"{prefix}_position_percentages.tsv", args.sample, change, outdir)

    print(f"Done. Plots written to: {outdir}")


if __name__ == "__main__":
    main()
