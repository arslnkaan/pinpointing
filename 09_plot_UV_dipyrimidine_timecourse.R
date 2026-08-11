#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# SETTINGS
# ============================================================

indir <- paste0(
  "/work/users/a/r/arslank/uvv/",
  "CSB_UV_dipyrimidine_7th_from_3prime"
)

infile <- file.path(
  indir,
  "UV_dipyrimidine_7th_from_3prime_summary.tsv"
)

outdir <- file.path(
  indir,
  "plots"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

dimer_order <- c(
  "TT",
  "CT",
  "TC",
  "CC"
)

dimer_colors <- c(
  "CT" = "#E59A3A",
  "TC" = "#8E24AA",
  "TT" = "#4F759B",
  "CC" = "#C62828"
)

# ============================================================
# READ
# ============================================================

dat <- read_tsv(
  infile,
  show_col_types = FALSE
) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    dinucleotide = factor(
      dinucleotide,
      levels = dimer_order
    )
  ) %>%
  arrange(
    time_h,
    dinucleotide
  )

write_tsv(
  dat,
  file.path(
    indir,
    "UV_dipyrimidine_7th_from_3prime_plot_table.tsv"
  )
)

# ============================================================
# PLOT 1
# PERCENTAGE OF ALL COMBINED REPAIR READS
# ============================================================

p_all <- ggplot(
  dat,
  aes(
    x = time_h,
    y = percent_of_all_reads,
    color = dinucleotide,
    group = dinucleotide
  )
) +
  geom_line(
    linewidth = 0.9
  ) +
  geom_point(
    size = 2.5
  ) +
  scale_color_manual(
    values = dimer_colors,
    breaks = dimer_order,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    ),
    labels = c(
      "0.5",
      "2",
      "4",
      "8"
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.08
      )
    )
  ) +
  labs(
    title = "Dipyrimidines in UV repair products",
    subtitle = "7th overlapping dinucleotide from the 3′ end",
    x = "Repair time (h)",
    y = "Contribution to all repair reads (%)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    legend.position = "top",

    plot.title = element_text(
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.35
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "UV_dipyrimidine_percent_of_all_reads_timecourse.pdf"
  ),
  plot = p_all,
  width = 6.5,
  height = 4.6,
  units = "in"
)

ggsave(
  filename = file.path(
    outdir,
    "UV_dipyrimidine_percent_of_all_reads_timecourse.png"
  ),
  plot = p_all,
  width = 6.5,
  height = 4.6,
  units = "in",
  dpi = 600
)

# ============================================================
# PLOT 2
# COMPOSITION WITHIN DIPYRIMIDINE READS
# ============================================================

p_composition <- ggplot(
  dat,
  aes(
    x = timepoint,
    y = percent_of_dipyrimidine_reads,
    fill = dinucleotide
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = dimer_colors,
    breaks = dimer_order,
    name = NULL
  ) +
  scale_y_continuous(
    limits = c(
      0,
      100
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  labs(
    title = "Dipyrimidine composition of UV repair products",
    subtitle = "CC + CT + TC + TT sum to 100% at each time point",
    x = "Repair time",
    y = "Dipyrimidine reads (%)"
  ) +
  theme_classic(
    base_size = 12
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    legend.position = "top",

    plot.title = element_text(
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.35
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "WT_UV_dipyrimidine_composition_timecourse.pdf"
  ),
  plot = p_composition,
  width = 5.8,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(
    outdir,
    "WT_UV_dipyrimidine_composition_timecourse.png"
  ),
  plot = p_composition,
  width = 5.8,
  height = 4.8,
  units = "in",
  dpi = 600
)

cat("\nDone.\n")
cat("Plot directory:\n")
cat(outdir, "\n")