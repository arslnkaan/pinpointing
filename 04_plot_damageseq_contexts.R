#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

indir <- paste0(
  "/work/users/a/r/arslank/damseq/",
  "NHF1_CPD_0h_r1_results/",
  "damage_context_denominators"
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

ct_tc_file <- file.path(
  indir,
  "DamageSeq_CT_TC_trinucleotide_counts.tsv"
)

cc_file <- file.path(
  indir,
  "DamageSeq_CC_NCCN_counts.tsv"
)

# ============================================================
# ORDERS
# ============================================================

ct_tc_order <- c(
  "CT:ACT",
  "CT:CCT",
  "CT:GCT",
  "CT:TCT",
  "TC:TCA",
  "TC:TCC",
  "TC:TCG",
  "TC:TCT"
)

nccn_order <- c(
  "ACCA", "ACCC", "ACCG", "ACCT",
  "CCCA", "CCCC", "CCCG", "CCCT",
  "GCCA", "GCCC", "GCCG", "GCCT",
  "TCCA", "TCCC", "TCCG", "TCCT"
)

# ============================================================
# READ DATA
# ============================================================

ct_tc <- read_tsv(
  ct_tc_file,
  show_col_types = FALSE
) %>%
  mutate(
    context_key = paste(
      lesion_type,
      trinucleotide,
      sep = ":"
    ),
    context_key = factor(
      context_key,
      levels = ct_tc_order
    )
  )

cc <- read_tsv(
  cc_file,
  show_col_types = FALSE
) %>%
  mutate(
    NCCN = factor(
      NCCN,
      levels = nccn_order
    ),
    flank_5prime = str_sub(
      as.character(NCCN),
      1,
      1
    )
  )

# ============================================================
# CT / TC PLOT
# ============================================================

ct_tc_colors <- c(
  "CT" = "#C62828",
  "TC" = "#8E24AA"
)

p_ct_tc <- ggplot(
  ct_tc,
  aes(
    x = context_key,
    y = percent_within_lesion_type,
    fill = lesion_type
  )
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = ct_tc_colors,
    name = NULL
  ) +
  labs(
    title = "Damage-seq CT and TC contexts",
    x = NULL,
    y = "Damage within lesion class (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "plain"
    ),
    legend.position = "top",
    panel.border = element_rect(
      fill = NA,
      color = "grey60",
      linewidth = 0.35
    )
  )

# ============================================================
# CC NCCN PLOT
# ============================================================

nccn_colors <- c(
  "A" = "#4F759B",
  "C" = "#D98C00",
  "G" = "#1B7837",
  "T" = "#B2182B"
)

p_cc <- ggplot(
  cc,
  aes(
    x = NCCN,
    y = percent,
    fill = flank_5prime
  )
) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = nccn_colors,
    guide = "none"
  ) +
  labs(
    title = "Damage-seq CC tetranucleotide contexts",
    x = "NCCN context",
    y = "CC damage (%)"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8,
      color = "black"
    ),
    axis.text.y = element_text(
      color = "black"
    ),
    plot.title = element_text(
      hjust = 0.5,
      face = "plain"
    ),
    panel.border = element_rect(
      fill = NA,
      color = "grey60",
      linewidth = 0.35
    )
  )

combined <- p_ct_tc / p_cc +
  plot_layout(
    heights = c(1, 1.1)
  )

# ============================================================
# SAVE
# ============================================================

ggsave(
  filename = "DamageSeq_CT_TC_and_CC_contexts.pdf",
  plot = combined,
  path = outdir,
  width = 7.2,
  height = 7.5,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = "DamageSeq_CT_TC_and_CC_contexts.png",
  plot = combined,
  path = outdir,
  width = 7.2,
  height = 7.5,
  units = "in",
  dpi = 600
)

cat("Done.\n")
cat("Output directory:", outdir, "\n")