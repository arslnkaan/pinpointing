#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(ggplot2)
  library(glue)
  library(dplyr)
  library(readr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 5) {
  stop("Usage: Rscript figures_batch.R <sample> <input_dir> <output_dir> <len_min> <len_max>")
}

sample     <- args[1]
input_dir  <- args[2]
output_dir <- args[3]
len_min    <- as.integer(args[4])
len_max    <- as.integer(args[5])

dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

base_theme <- theme_classic() +
  theme(
    axis.title = element_text(size = 14),
    axis.text  = element_text(size = 12),
    legend.title = element_text(size = 13),
    legend.text  = element_text(size = 11),
    plot.title = element_text(hjust = 0.5, size = 16, margin = margin(b = 10)),
    plot.margin = margin(t = 20, r = 10, b = 10, l = 10),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

save_png_pdf <- function(plot, stem, w, h) {
  ggsave(glue("{stem}.png"), plot = plot, width = w, height = h, dpi = 300, path = output_dir)
  ggsave(glue("{stem}.pdf"), plot = plot, width = w, height = h, dpi = 300, path = output_dir, device = cairo_pdf)
}

base_cols <- c(
  "C" = "dodgerblue4",
  "T" = "orange",
  "A" = "green4",
  "G" = "darkorchid4"
)

# ============================================================
# 1. Read length distribution
# ============================================================

rl_file <- file.path(input_dir, glue("{sample}_read_length_distribution.txt"))

df_rl <- read_tsv(rl_file, show_col_types = FALSE)

p_len <- ggplot(df_rl, aes(x = Length, y = 100 * Count / sum(Count))) +
  geom_col(fill = "#4DBBD5") +
  geom_vline(xintercept = 26, linetype = "dashed", color = "black", linewidth = 0.8) +
  labs(
    y = "Frequency (%)",
    x = "Length (nt)",
    title = glue("{sample} read length distribution")
  ) +
  base_theme

save_png_pdf(p_len, glue("{sample}_read_length"), 12, 8)

# ============================================================
# 2. Monomer position frequency
# ============================================================

mono_file <- file.path(
  input_dir,
  "monomer",
  glue("{sample}.len{len_min}_{len_max}_monomer_R_df.txt")
)

df_mono <- read_tsv(mono_file, show_col_types = FALSE) %>%
  mutate(
    Percent = 100 * Frequency
  )

plot_monomer <- function(base_order, tag) {
  df_plot <- df_mono %>%
    filter(Base %in% base_order) %>%
    mutate(Base = factor(Base, levels = base_order))

  p <- ggplot(df_plot, aes(x = Position, y = Percent, fill = Base)) +
    geom_col(width = 0.95) +
    facet_grid(Length ~ ., scales = "free_x", space = "free_x") +
    scale_fill_manual(values = base_cols[base_order]) +
    labs(
      title = glue("{sample} base distribution: {paste(base_order, collapse = ', ')}"),
      y = "Nucleotide frequency (%)",
      x = "Position"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 11)
    )

  save_png_pdf(p, glue("{sample}_monomer_{tag}"), 12, 20)
}

# ggplot stacks bottom to top in factor level order
plot_monomer(c("G", "A", "C", "T"), "bottom_G_A_C_T")
plot_monomer(c("T", "C", "A", "G"), "bottom_T_C_A_G")

# ============================================================
# 3. Dinucleotide frequency
# ============================================================

dinuc_file <- file.path(
  input_dir,
  "monomer",
  glue("{sample}.len{len_min}_{len_max}_dinucleotide_R_df.txt")
)

df_dinuc <- read_tsv(dinuc_file, show_col_types = FALSE) %>%
  mutate(
    Percent = 100 * Frequency
  )

dinuc_cols <- c(
  "GG" = "#5F0A87",
  "GT" = "#6C63FF",
  "GA" = "#A07BEF",
  "AG" = "#C084FC",
  "TT" = "#C96C1A",
  "CT" = "#E59A3A",
  "TC" = "#F3BE73",
  "CC" = "#F8D9A8"
)

plot_dinuc <- function(dinuc_order, tag) {
  df_plot <- df_dinuc %>%
    filter(Dinucleotide %in% dinuc_order) %>%
    mutate(Dinucleotide = factor(Dinucleotide, levels = dinuc_order))

  p <- ggplot(df_plot, aes(x = Position, y = Percent, fill = Dinucleotide)) +
    geom_col(width = 0.95) +
    facet_grid(Length ~ ., scales = "free_x", space = "free_x") +
    scale_fill_manual(values = dinuc_cols[dinuc_order]) +
    labs(
      title = glue("{sample} dinucleotide distribution: {paste(dinuc_order, collapse = ', ')}"),
      y = "Dinucleotide frequency (%)",
      x = "Position"
    ) +
    base_theme +
    theme(
      axis.text.x = element_text(size = 7),
      strip.text.y = element_text(size = 11)
    )

  save_png_pdf(p, glue("{sample}_dinucleotide_{tag}"), 12, 20)
}

plot_dinuc(c("GG", "GT", "GA", "AG"), "bottom_GG_GT_GA_AG")
plot_dinuc(c("TT", "CT", "TC", "CC"), "bottom_TT_CT_TC_CC")

# ============================================================
# 4. TS/NTS repair profiles
# ============================================================

ts_file <- file.path(input_dir, glue("{sample}_TScount.txt"))
nts_file <- file.path(input_dir, glue("{sample}_NTScount.txt"))
read_count_file <- file.path(input_dir, glue("{sample}.dedup_molecules.txt"))

total_reads <- as.numeric(readLines(read_count_file, n = 1))

load_profile <- function(file, strand_label) {
  read.table(
    file,
    header = FALSE,
    col.names = c("chrom", "start", "end", "bin_num", "name", "strand", "gene_num", "reads")
  ) %>%
    mutate(
      width = end - start,
      RPKM = (reads * 1e9) / (total_reads * width)
    ) %>%
    group_by(bin_num) %>%
    summarise(
      mean = mean(RPKM, na.rm = TRUE),
      se = sd(RPKM, na.rm = TRUE) / sqrt(n()),
      .groups = "drop"
    ) %>%
    mutate(strand = strand_label)
}

ts_summary <- load_profile(ts_file, "TS")
nts_summary <- load_profile(nts_file, "NTS")

combined_summary <- bind_rows(ts_summary, nts_summary)

p_repair <- ggplot(combined_summary, aes(x = bin_num, y = mean, color = strand, fill = strand)) +
  geom_line(linewidth = 1) +
  geom_ribbon(aes(ymin = mean - se, ymax = mean + se), alpha = 0.25, linetype = 0) +
  scale_color_manual(values = c("TS" = "goldenrod2", "NTS" = "dodgerblue3")) +
  scale_fill_manual(values  = c("TS" = "goldenrod2", "NTS" = "dodgerblue3")) +
  scale_x_continuous(
    breaks = c(0, 25, 75, 125, 150),
    labels = c("-2kb", "TSS", "50%", "TES", "2kb")
  ) +
  geom_vline(xintercept = c(25, 125), linetype = "dashed", color = "black", linewidth = 0.5) +
  labs(
    title = glue("{sample} TS vs NTS"),
    y = "Average repair signal (RPKM)",
    x = NULL
  ) +
  base_theme

save_png_pdf(p_repair, glue("{sample}_repair_TS_NTS"), 12, 8)

message("Done. Plots saved to: ", output_dir)