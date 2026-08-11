#!/usr/bin/env python3

import argparse
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt


def infer_change_from_filename(path):
    name = Path(path).name
    if "G_to_T" in name:
        return "G>T"
    if "G_to_A" in name:
        return "G>A"
    if "G_to_C" in name:
        return "G>C"
    if "C_to_T" in name:
        return "C>T"
    if "C_to_A" in name:
        return "C>A"
    if "C_to_G" in name:
        return "C>G"
    return None


def read_event_file(path):
    df = pd.read_csv(path, sep="\t")

    if "change" in df.columns:
        df["change"] = df["change"].astype(str)
    elif "Change" in df.columns:
        df["change"] = df["Change"].astype(str)
    elif {"Reference_Base", "Alt_Base"}.issubset(df.columns):
        df["change"] = df["Reference_Base"].astype(str) + ">" + df["Alt_Base"].astype(str)
    else:
        inferred = infer_change_from_filename(path)
        if inferred is None:
            raise ValueError(f"Could not infer mutation class from: {path}")
        df["change"] = inferred

    if "read_length" not in df.columns:
        if "Read_Length" in df.columns:
            df["read_length"] = df["Read_Length"]
        else:
            raise ValueError(f"Missing read length column in: {path}")

    if "position_from_5prime" not in df.columns:
        if "Position_from_5prime" in df.columns:
            df["position_from_5prime"] = df["Position_from_5prime"]
        else:
            raise ValueError(f"Missing 5prime position column in: {path}")

    if "position_from_3prime" not in df.columns:
        if "Position_from_3prime" in df.columns:
            df["position_from_3prime"] = df["Position_from_3prime"]

    return df


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--events", required=True)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--min-len", type=int, default=26)
    parser.add_argument("--max-len", type=int, default=30)
    parser.add_argument("--changes", default="G>T,G>A,G>C")
    parser.add_argument("--percent", action="store_true")
    parser.add_argument("--position", choices=["5prime", "3prime"], default="5prime")
    args = parser.parse_args()

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    changes = [x.strip() for x in args.changes.split(",")]
    event_files = [x.strip() for x in args.events.split(",")]

    dfs = []
    for event_file in event_files:
        if not Path(event_file).exists():
            raise FileNotFoundError(f"Missing file: {event_file}")
        dfs.append(read_event_file(event_file))

    df = pd.concat(dfs, ignore_index=True)

    pos_col = "position_from_5prime" if args.position == "5prime" else "position_from_3prime"

    df["read_length"] = pd.to_numeric(df["read_length"], errors="coerce")
    df[pos_col] = pd.to_numeric(df[pos_col], errors="coerce")

    df = df.dropna(subset=["read_length", pos_col, "change"]).copy()
    df["read_length"] = df["read_length"].astype(int)
    df[pos_col] = df[pos_col].astype(int)

    df = df[
        (df["read_length"] >= args.min_len) &
        (df["read_length"] <= args.max_len) &
        (df["change"].isin(changes))
    ].copy()

    if df.empty:
        raise ValueError("No mismatch events left after filtering.")

    colors = {
        "G>T": "#f39c12",
        "G>A": "#006D2C",
        "G>C": "#0072B2",
        "C>T": "#D55E00",
        "C>A": "#56B4E9",
        "C>G": "#009E73"
    }

    lengths = list(range(args.min_len, args.max_len + 1))

    unit_width = 0.28
    left_margin = 1.25
    right_margin = 0.55
    row_height = 1.28

    fig_width = left_margin + args.max_len * unit_width + right_margin
    fig_height = len(lengths) * row_height + 1.55

    fig = plt.figure(figsize=(fig_width, fig_height))

    axes = []

    # Extra top room for title + legend
    plot_top = 0.72
    plot_height_total = 0.58
    row_gap_factor = plot_height_total / len(lengths)
    row_height_frac = row_gap_factor * 0.72

    for i, read_len in enumerate(lengths):
        bottom = plot_top - i * row_gap_factor
        height = row_height_frac
        width = read_len * unit_width / fig_width
        left = left_margin / fig_width

        ax = fig.add_axes([left, bottom, width, height])
        axes.append(ax)

        sub = df[df["read_length"] == read_len].copy()
        positions = list(range(1, read_len + 1))

        count_table = (
            sub.groupby([pos_col, "change"])
            .size()
            .reset_index(name="count")
            .pivot(index=pos_col, columns="change", values="count")
            .reindex(index=positions, columns=changes)
            .fillna(0)
        )

        if args.percent:
            denom = count_table.values.sum()
            plot_table = count_table / denom * 100 if denom > 0 else count_table
            ylabel = "%"
        else:
            plot_table = count_table
            ylabel = "Count"

        running_bottom = [0] * len(positions)

        for change in changes:
            vals = plot_table[change].values
            ax.bar(
                positions,
                vals,
                bottom=running_bottom,
                width=0.85,
                color=colors.get(change, None),
                edgecolor="black",
                linewidth=0.25,
                label=change
            )
            running_bottom = [b + v for b, v in zip(running_bottom, vals)]

        ax.set_xlim(0.5, read_len + 0.5)
        ax.set_xticks(range(1, read_len + 1))
        ax.tick_params(axis="x", labelsize=6)
        ax.tick_params(axis="y", labelsize=7)
        ax.grid(axis="y", linestyle="--", linewidth=0.4, alpha=0.45)
        ax.set_ylabel(ylabel, fontsize=8)

        ax.text(
            -0.08,
            0.5,
            f"{read_len}-mer",
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=9,
            fontweight="bold"
        )

    handles, labels = axes[0].get_legend_handles_labels()

    fig.legend(
        handles=handles,
        labels=labels,
        title="Mismatch",
        ncol=len(changes),
        loc="upper center",
        bbox_to_anchor=(0.58, 0.92),
        frameon=False,
        fontsize=9,
        title_fontsize=10
    )

    axes[-1].set_xlabel(
        "Mismatch position from 5′ end" if args.position == "5prime"
        else "Mismatch position from 3′ end",
        fontsize=10
    )

    fig.suptitle(
        f"{args.sample}: mismatch position distribution",
        fontsize=13,
        y=0.985
    )

    suffix = "percent" if args.percent else "counts"

    pdf = outdir / (
        f"{args.sample}_mismatch_position_scaled_length_from_{args.position}_"
        f"{args.min_len}to{args.max_len}_{suffix}.pdf"
    )
    png = outdir / (
        f"{args.sample}_mismatch_position_scaled_length_from_{args.position}_"
        f"{args.min_len}to{args.max_len}_{suffix}.png"
    )

    plt.savefig(pdf, bbox_inches="tight")
    plt.savefig(png, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"Wrote: {pdf}")
    print(f"Wrote: {png}")


if __name__ == "__main__":
    main()