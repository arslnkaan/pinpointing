#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "sans-serif"
matplotlib.rcParams["font.sans-serif"] = [
    "Helvetica",
    "Arial",
    "Liberation Sans",
    "DejaVu Sans"
]

import matplotlib.pyplot as plt


def read_event_file(path):
    df = pd.read_csv(path, sep="\t")

    required = {"Read_Length", "Position_from_5prime", "Change"}
    missing = required - set(df.columns)

    if missing:
        raise ValueError(f"Missing columns in {path}: {missing}")

    return df.rename(columns={
        "Read_Length": "read_length",
        "Position_from_5prime": "position",
        "Change": "change"
    })


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-len", type=int, default=20)
    parser.add_argument("--max-len", type=int, default=30)
    parser.add_argument("--changes", default="G>T,G>A,G>C")
    parser.add_argument("--percent", action="store_true")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    changes = [x.strip() for x in args.changes.split(",")]
    files = [x.strip() for x in args.events.split(",") if x.strip()]

    df = pd.concat([read_event_file(f) for f in files], ignore_index=True)

    df["read_length"] = pd.to_numeric(df["read_length"], errors="coerce")
    df["position"] = pd.to_numeric(df["position"], errors="coerce")
    df["change"] = df["change"].astype(str)

    df = df.dropna(subset=["read_length", "position", "change"]).copy()
    df["read_length"] = df["read_length"].astype(int)
    df["position"] = df["position"].astype(int)

    df = df[
        (df["read_length"] >= args.min_len) &
        (df["read_length"] <= args.max_len) &
        (df["position"] >= 1) &
        (df["position"] <= df["read_length"]) &
        (df["change"].isin(changes))
    ].copy()

    if df.empty:
        raise ValueError("No events left after filtering.")

    colors = {
        "G>T": "#E69F00",
        "G>A": "#009E73",
        "G>C": "#0072B2",
        "C>T": "darkred",
        "C>A": "#56B4E9",
        "C>G": "#CC79A7"
    }

    lengths = list(range(args.min_len, args.max_len + 1))

    tables = {}
    global_max = 0

    for read_len in lengths:
        sub = df[df["read_length"] == read_len]
        positions = list(range(1, read_len + 1))

        tab = (
            sub.groupby(["position", "change"])
            .size()
            .reset_index(name="count")
            .pivot(index="position", columns="change", values="count")
            .reindex(index=positions, columns=changes)
            .fillna(0)
        )

        if args.percent:
            denom = tab.values.sum()
            tab = tab / denom * 100 if denom > 0 else tab

        tables[read_len] = tab
        global_max = max(global_max, tab.sum(axis=1).max())

    if global_max <= 0:
        global_max = 1

    unit_width = 0.205
    left_margin = 1.15
    right_margin = 0.35
    fig_width = left_margin + args.max_len * unit_width + right_margin
    fig_height = len(lengths) * 1.28 + 2.35

    fig = plt.figure(figsize=(fig_width, fig_height))

    plot_left = left_margin / fig_width
    plot_top = 0.75
    plot_height_total = 0.68
    row_gap = plot_height_total / len(lengths)
    row_height = row_gap * 0.62

    axes = []

    for i, read_len in enumerate(lengths):
        bottom = plot_top - i * row_gap
        width = read_len * unit_width / fig_width

        ax = fig.add_axes([plot_left, bottom, width, row_height])
        axes.append(ax)

        tab = tables[read_len]
        positions = list(range(1, read_len + 1))

        running_bottom = [0] * len(positions)

        for change in changes:
            vals = tab[change].values

            ax.bar(
                positions,
                vals,
                bottom=running_bottom,
                width=0.85,
                color=colors.get(change, "grey"),
                edgecolor="black",
                linewidth=0.25,
                label=change
            )

            running_bottom = [b + v for b, v in zip(running_bottom, vals)]

        ax.set_xlim(0.5, read_len + 0.5)
        ax.set_ylim(0, global_max * 1.15)

        ax.set_xticks(list(range(5, read_len + 1, 5)))

        ax.tick_params(axis="x", labelsize=11, width=1.0, length=4)
        ax.tick_params(axis="y", labelsize=10, width=1.0, length=3.5)

        ax.grid(axis="y", linestyle="--", linewidth=0.45, alpha=0.35)

        for spine in ["top", "right"]:
            ax.spines[spine].set_visible(False)

        ax.spines["left"].set_linewidth(1.0)
        ax.spines["bottom"].set_linewidth(1.0)

    for ax in axes[:-1]:
        ax.set_xticklabels([])

    axes[-1].set_xlabel(
        "Misinsertion position",
        fontsize=16
    )

    fig.text(
        0.07,
        0.50,
        "%" if args.percent else "Count",
        rotation=90,
        va="center",
        ha="center",
        fontsize=16
    )

    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles=handles,
        labels=labels,
        title="Misinsertion",
        ncol=len(changes),
        loc="upper center",
        bbox_to_anchor=(0.50, 0.955),
        frameon=False,
        fontsize=14,
        title_fontsize=14,
        handlelength=2.2,
        columnspacing=2.0
    )

    fig.suptitle(
        f"{args.sample}",
        fontsize=18,
        y=0.985
    )

    suffix = "percent" if args.percent else "counts"

    pdf = outdir / (
        f"{args.sample}_misinsertion_position_5prime_aligned_"
        f"{args.min_len}to{args.max_len}_{suffix}.pdf"
    )

    png = outdir / (
        f"{args.sample}_misinsertion_position_5prime_aligned_"
        f"{args.min_len}to{args.max_len}_{suffix}.png"
    )

    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, dpi=600, bbox_inches="tight")
    plt.close()

    print(f"Wrote: {pdf}")
    print(f"Wrote: {png}")


if __name__ == "__main__":
    main()