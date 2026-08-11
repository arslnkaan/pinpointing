#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

base_dir <- "/work/users/a/r/arslank"

repair_file <- file.path(
  base_dir,
  "UV_CtoT_trinucleotide_timecourse",
  "UV_CtoT_trinucleotide_counts_percentages.tsv"
)

damage_file <- file.path(
  base_dir,
  "damseq",
  "NHF1_CPD_0h_r1_results",
  "damage_CT_TC_contexts",
  "DamageSeq_NCT_TCN_counts.tsv"
)

outdir <- file.path(
  base_dir,
  "UV_NCT_TCN_repair_damage"
)

plotdir <- file.path(
  outdir,
  "plots"
)

dir.create(
  plotdir,
  recursive = TRUE,
  showWarnings = FALSE
)

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

# Exactly four contexts per lesion class.
nct_order <- c(
  "ACT",
  "CCT",
  "GCT",
  "TCT"
)

tcn_order <- c(
  "TCA",
  "TCC",
  "TCG",
  "TCT"
)

nct_colors <- c(
  "ACT" = "#4F759B",
  "CCT" = "#D98C00",
  "GCT" = "#1B7837",
  "TCT" = "#B2182B"
)

tcn_colors <- c(
  "TCA" = "#4F759B",
  "TCC" = "#D98C00",
  "TCG" = "#1B7837",
  "TCT" = "#B2182B"
)

# ============================================================
# CHECK INPUTS
# ============================================================

for (file in c(repair_file, damage_file)) {
  if (!file.exists(file)) {
    stop(
      paste0(
        "Missing input file: ",
        file
      )
    )
  }
}

# ============================================================
# HELPERS
# ============================================================

safe_fraction <- function(
  numerator,
  denominator
) {
  if_else(
    is.finite(denominator) &
      denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

safe_ratio <- function(
  numerator,
  denominator
) {
  if_else(
    is.finite(denominator) &
      denominator > 0,
    numerator / denominator,
    NA_real_
  )
}

safe_log2 <- function(x) {
  case_when(
    is.na(x) ~ NA_real_,
    x > 0 ~ log2(x),
    x == 0 ~ -Inf,
    TRUE ~ NA_real_
  )
}

finite_for_plot <- function(x) {
  if_else(
    is.finite(x),
    x,
    NA_real_
  )
}

# ============================================================
# READ REPAIR DATA
# ============================================================

repair_raw <- read_tsv(
  repair_file,
  show_col_types = FALSE
) %>%
  transmute(
    sample,
    timepoint = as.character(timepoint),
    time_h = as.numeric(time_h),
    replicate,
    trinucleotide = str_to_upper(
      as.character(trinucleotide)
    ),
    repair_count = as.numeric(count)
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate,
    trinucleotide
  ) %>%
  summarise(
    repair_count = sum(
      repair_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

sample_design <- repair_raw %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate
  )

all_repair_totals <- repair_raw %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  summarise(
    total_all_CtoT_repair = sum(
      repair_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

# ============================================================
# BUILD NCT REPAIR TABLE
#
# Fraction denominator:
# ACT + CCT + GCT + TCT
# ============================================================

nct_repair <- crossing(
  sample_design,
  trinucleotide = nct_order
) %>%
  left_join(
    repair_raw,
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate",
      "trinucleotide"
    )
  ) %>%
  mutate(
    repair_count = replace_na(
      repair_count,
      0
    )
  ) %>%
  left_join(
    all_repair_totals,
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate"
    )
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_NCT_repair = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_fraction_within_class =
      safe_fraction(
        repair_count,
        total_NCT_repair
      ),

    repair_percent_within_class =
      100 *
      repair_fraction_within_class,

    repair_fraction_of_all_CtoT =
      safe_fraction(
        repair_count,
        total_all_CtoT_repair
      ),

    repair_percent_of_all_CtoT =
      100 *
      repair_fraction_of_all_CtoT
  ) %>%
  ungroup() %>%
  mutate(
    lesion_class = "NCT",
    orientation_note = if_else(
      trinucleotide == "TCT",
      paste(
        "TCT repair orientation is ambiguous;",
        "the same repair series is shown in NCT and TCN"
      ),
      NA_character_
    ),
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),
    trinucleotide = factor(
      trinucleotide,
      levels = nct_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

# ============================================================
# BUILD TCN REPAIR TABLE
#
# Fraction denominator:
# TCA + TCC + TCG + TCT
# ============================================================

tcn_repair <- crossing(
  sample_design,
  trinucleotide = tcn_order
) %>%
  left_join(
    repair_raw,
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate",
      "trinucleotide"
    )
  ) %>%
  mutate(
    repair_count = replace_na(
      repair_count,
      0
    )
  ) %>%
  left_join(
    all_repair_totals,
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate"
    )
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_TCN_repair = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_fraction_within_class =
      safe_fraction(
        repair_count,
        total_TCN_repair
      ),

    repair_percent_within_class =
      100 *
      repair_fraction_within_class,

    repair_fraction_of_all_CtoT =
      safe_fraction(
        repair_count,
        total_all_CtoT_repair
      ),

    repair_percent_of_all_CtoT =
      100 *
      repair_fraction_of_all_CtoT
  ) %>%
  ungroup() %>%
  mutate(
    lesion_class = "TCN",
    orientation_note = if_else(
      trinucleotide == "TCT",
      paste(
        "TCT repair orientation is ambiguous;",
        "the same repair series is shown in NCT and TCN"
      ),
      NA_character_
    ),
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),
    trinucleotide = factor(
      trinucleotide,
      levels = tcn_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

# ============================================================
# READ DAMAGE DATA
# ============================================================

damage_raw <- read_tsv(
  damage_file,
  show_col_types = FALSE
) %>%
  transmute(
    lesion_type = str_to_upper(
      as.character(lesion_type)
    ),
    trinucleotide = str_to_upper(
      as.character(trinucleotide)
    ),
    damage_count = as.numeric(count)
  )

# ============================================================
# NCT DAMAGE AT 0 h
#
# CT lesions:
# ACT, CCT, GCT, TCT
# ============================================================

nct_damage <- tibble(
  trinucleotide = nct_order
) %>%
  left_join(
    damage_raw %>%
      filter(
        lesion_type == "CT",
        trinucleotide %in%
          nct_order
      ) %>%
      group_by(
        trinucleotide
      ) %>%
      summarise(
        damage_count = sum(
          damage_count,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "trinucleotide"
  ) %>%
  mutate(
    damage_count = replace_na(
      damage_count,
      0
    ),

    total_NCT_damage = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_fraction_within_class =
      safe_fraction(
        damage_count,
        total_NCT_damage
      ),

    damage_percent_within_class =
      100 *
      damage_fraction_within_class,

    lesion_class = "NCT",
    damage_timepoint = "0h",

    trinucleotide = factor(
      trinucleotide,
      levels = nct_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

# ============================================================
# TCN DAMAGE AT 0 h
#
# TC lesions:
# TCA, TCC, TCG, TCT
# ============================================================

tcn_damage <- tibble(
  trinucleotide = tcn_order
) %>%
  left_join(
    damage_raw %>%
      filter(
        lesion_type == "TC",
        trinucleotide %in%
          tcn_order
      ) %>%
      group_by(
        trinucleotide
      ) %>%
      summarise(
        damage_count = sum(
          damage_count,
          na.rm = TRUE
        ),
        .groups = "drop"
      ),
    by = "trinucleotide"
  ) %>%
  mutate(
    damage_count = replace_na(
      damage_count,
      0
    ),

    total_TCN_damage = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_fraction_within_class =
      safe_fraction(
        damage_count,
        total_TCN_damage
      ),

    damage_percent_within_class =
      100 *
      damage_fraction_within_class,

    lesion_class = "TCN",
    damage_timepoint = "0h",

    trinucleotide = factor(
      trinucleotide,
      levels = tcn_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

# ============================================================
# NORMALIZE NCT REPAIR TO NCT DAMAGE
#
# log2[
#   (repair context / all NCT repair contexts)
#   /
#   (damage context / all NCT damage contexts)
# ]
# ============================================================

nct_normalized <- nct_repair %>%
  left_join(
    nct_damage %>%
      transmute(
        trinucleotide =
          as.character(trinucleotide),
        damage_count,
        total_NCT_damage,
        damage_fraction_within_class,
        damage_percent_within_class
      ),
    by = c(
      "trinucleotide"
    )
  ) %>%
  mutate(
    repair_over_damage =
      safe_ratio(
        repair_fraction_within_class,
        damage_fraction_within_class
      ),

    log2_repair_over_damage =
      safe_log2(
        repair_over_damage
      ),

    log2_for_plot =
      finite_for_plot(
        log2_repair_over_damage
      ),

    normalization_status = case_when(
      is.na(damage_fraction_within_class) ~
        "missing_damage",

      damage_fraction_within_class <= 0 ~
        "zero_damage",

      repair_count == 0 ~
        "zero_repair",

      TRUE ~
        "normalized"
    ),

    trinucleotide = factor(
      as.character(trinucleotide),
      levels = nct_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

# ============================================================
# NORMALIZE TCN REPAIR TO TCN DAMAGE
# ============================================================

tcn_normalized <- tcn_repair %>%
  left_join(
    tcn_damage %>%
      transmute(
        trinucleotide =
          as.character(trinucleotide),
        damage_count,
        total_TCN_damage,
        damage_fraction_within_class,
        damage_percent_within_class
      ),
    by = c(
      "trinucleotide"
    )
  ) %>%
  mutate(
    repair_over_damage =
      safe_ratio(
        repair_fraction_within_class,
        damage_fraction_within_class
      ),

    log2_repair_over_damage =
      safe_log2(
        repair_over_damage
      ),

    log2_for_plot =
      finite_for_plot(
        log2_repair_over_damage
      ),

    normalization_status = case_when(
      is.na(damage_fraction_within_class) ~
        "missing_damage",

      damage_fraction_within_class <= 0 ~
        "zero_damage",

      repair_count == 0 ~
        "zero_repair",

      TRUE ~
        "normalized"
    ),

    trinucleotide = factor(
      as.character(trinucleotide),
      levels = tcn_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

# ============================================================
# WRITE TABLES
# ============================================================

write_tsv(
  nct_repair,
  file.path(
    outdir,
    "NCT_repair_over_time.tsv"
  )
)

write_tsv(
  tcn_repair,
  file.path(
    outdir,
    "TCN_repair_over_time.tsv"
  )
)

write_tsv(
  nct_damage,
  file.path(
    outdir,
    "NCT_damage_0h.tsv"
  )
)

write_tsv(
  tcn_damage,
  file.path(
    outdir,
    "TCN_damage_0h.tsv"
  )
)

write_tsv(
  nct_normalized,
  file.path(
    outdir,
    "NCT_log2_repair_over_damage.tsv"
  )
)

write_tsv(
  tcn_normalized,
  file.path(
    outdir,
    "TCN_log2_repair_over_damage.tsv"
  )
)

# ============================================================
# SHARED PLOT LIMITS
# ============================================================

raw_y_max <- max(
  c(
    nct_repair$repair_percent_within_class,
    tcn_repair$repair_percent_within_class
  ),
  na.rm = TRUE
)

if (!is.finite(raw_y_max) || raw_y_max <= 0) {
  raw_y_max <- 100
}

raw_y_max <- raw_y_max * 1.08

damage_y_max <- max(
  c(
    nct_damage$damage_percent_within_class,
    tcn_damage$damage_percent_within_class
  ),
  na.rm = TRUE
)

if (!is.finite(damage_y_max) ||
    damage_y_max <= 0) {
  damage_y_max <- 100
}

damage_y_max <- damage_y_max * 1.08

normalized_values <- c(
  nct_normalized$log2_for_plot,
  tcn_normalized$log2_for_plot
)

normalized_values <- normalized_values[
  is.finite(normalized_values)
]

if (length(normalized_values) > 0) {
  normalized_limit <- max(
    abs(normalized_values),
    na.rm = TRUE
  )
} else {
  normalized_limit <- 1
}

if (!is.finite(normalized_limit) ||
    normalized_limit <= 0) {
  normalized_limit <- 1
}

normalized_limit <- normalized_limit * 1.10

# ============================================================
# COMMON THEME
# ============================================================

common_theme <- theme_classic(
  base_size = 11
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

# ============================================================
# PLOT 1: NCT REPAIR OVER TIME
# ============================================================

p_nct_repair <- ggplot(
  nct_repair,
  aes(
    x = time_h,
    y = repair_percent_within_class,
    color = trinucleotide,
    group = trinucleotide
  )
) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 2.3
  ) +
  scale_color_manual(
    values = nct_colors,
    breaks = nct_order,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      raw_y_max
    )
  ) +
  labs(
    title = "NCT repair over time",
    subtitle = "ACT + CCT + GCT + TCT = 100%",
    x = "Repair time (h)",
    y = "NCT repair composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 2: TCN REPAIR OVER TIME
# ============================================================

p_tcn_repair <- ggplot(
  tcn_repair,
  aes(
    x = time_h,
    y = repair_percent_within_class,
    color = trinucleotide,
    group = trinucleotide
  )
) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 2.3
  ) +
  scale_color_manual(
    values = tcn_colors,
    breaks = tcn_order,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  coord_cartesian(
    ylim = c(
      0,
      raw_y_max
    )
  ) +
  labs(
    title = "TCN repair over time",
    subtitle = "TCA + TCC + TCG + TCT = 100%",
    x = "Repair time (h)",
    y = "TCN repair composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 3: NCT DAMAGE AT 0 h
# ============================================================

p_nct_damage <- ggplot(
  nct_damage,
  aes(
    x = trinucleotide,
    y = damage_percent_within_class,
    fill = trinucleotide
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = nct_colors,
    breaks = nct_order,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = nct_order,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      0,
      damage_y_max
    )
  ) +
  labs(
    title = "NCT Damage-seq at 0 h",
    subtitle = "ACT + CCT + GCT + TCT = 100%",
    x = NULL,
    y = "NCT damage composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 4: TCN DAMAGE AT 0 h
# ============================================================

p_tcn_damage <- ggplot(
  tcn_damage,
  aes(
    x = trinucleotide,
    y = damage_percent_within_class,
    fill = trinucleotide
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = tcn_colors,
    breaks = tcn_order,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = tcn_order,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      0,
      damage_y_max
    )
  ) +
  labs(
    title = "TCN Damage-seq at 0 h",
    subtitle = "TCA + TCC + TCG + TCT = 100%",
    x = NULL,
    y = "TCN damage composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 5: NCT NORMALIZED TO NCT DAMAGE
# ============================================================

p_nct_normalized <- ggplot(
  nct_normalized,
  aes(
    x = time_h,
    y = log2_for_plot,
    color = trinucleotide,
    group = trinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.3,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = nct_colors,
    breaks = nct_order,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  coord_cartesian(
    ylim = c(
      -normalized_limit,
      normalized_limit
    )
  ) +
  labs(
    title = "NCT repair normalized to NCT damage",
    x = "Repair time (h)",
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  common_theme

# ============================================================
# PLOT 6: TCN NORMALIZED TO TCN DAMAGE
# ============================================================

p_tcn_normalized <- ggplot(
  tcn_normalized,
  aes(
    x = time_h,
    y = log2_for_plot,
    color = trinucleotide,
    group = trinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.85,
    na.rm = TRUE
  ) +
  geom_point(
    size = 2.3,
    na.rm = TRUE
  ) +
  scale_color_manual(
    values = tcn_colors,
    breaks = tcn_order,
    name = NULL
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  coord_cartesian(
    ylim = c(
      -normalized_limit,
      normalized_limit
    )
  ) +
  labs(
    title = "TCN repair normalized to TCN damage",
    x = "Repair time (h)",
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  common_theme

# ============================================================
# COMBINED SIX-PANEL FIGURE
# ============================================================

combined_plot <- (
  p_nct_repair |
    p_tcn_repair
) / (
  p_nct_damage |
    p_tcn_damage
) / (
  p_nct_normalized |
    p_tcn_normalized
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plotdir,
    "NCT_TCN_repair_damage_and_normalization.pdf"
  ),
  plot = combined_plot,
  width = 11,
  height = 13,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "NCT_TCN_repair_damage_and_normalization.png"
  ),
  plot = combined_plot,
  width = 11,
  height = 13,
  units = "in",
  dpi = 600
)

# ============================================================
# SAVE INDIVIDUAL FIGURES
# ============================================================

ggsave(
  file.path(
    plotdir,
    "NCT_repair_over_time.pdf"
  ),
  p_nct_repair,
  width = 6.2,
  height = 4.5
)

ggsave(
  file.path(
    plotdir,
    "TCN_repair_over_time.pdf"
  ),
  p_tcn_repair,
  width = 6.2,
  height = 4.5
)

ggsave(
  file.path(
    plotdir,
    "NCT_damage_0h.pdf"
  ),
  p_nct_damage,
  width = 5.2,
  height = 4.3
)

ggsave(
  file.path(
    plotdir,
    "TCN_damage_0h.pdf"
  ),
  p_tcn_damage,
  width = 5.2,
  height = 4.3
)

ggsave(
  file.path(
    plotdir,
    "NCT_log2_repair_over_damage.pdf"
  ),
  p_nct_normalized,
  width = 6.2,
  height = 4.5
)

ggsave(
  file.path(
    plotdir,
    "TCN_log2_repair_over_damage.pdf"
  ),
  p_tcn_normalized,
  width = 6.2,
  height = 4.5
)

# ============================================================
# QC
# ============================================================

qc <- bind_rows(
  nct_repair %>%
    group_by(
      sample,
      timepoint,
      time_h
    ) %>%
    summarise(
      lesion_class = "NCT repair",
      fraction_sum = sum(
        repair_fraction_within_class,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  tcn_repair %>%
    group_by(
      sample,
      timepoint,
      time_h
    ) %>%
    summarise(
      lesion_class = "TCN repair",
      fraction_sum = sum(
        repair_fraction_within_class,
        na.rm = TRUE
      ),
      .groups = "drop"
    ),

  nct_damage %>%
    summarise(
      sample = "DamageSeq",
      timepoint = "0h",
      time_h = 0,
      lesion_class = "NCT damage",
      fraction_sum = sum(
        damage_fraction_within_class,
        na.rm = TRUE
      )
    ),

  tcn_damage %>%
    summarise(
      sample = "DamageSeq",
      timepoint = "0h",
      time_h = 0,
      lesion_class = "TCN damage",
      fraction_sum = sum(
        damage_fraction_within_class,
        na.rm = TRUE
      )
    )
)

write_tsv(
  qc,
  file.path(
    outdir,
    "NCT_TCN_fraction_QC.tsv"
  )
)

cat("\nDone.\n\n")
cat("Combined figure:\n")
cat(
  file.path(
    plotdir,
    "NCT_TCN_repair_damage_and_normalization.pdf"
  ),
  "\n\n"
)

cat("Output tables:\n")
cat(outdir, "\n")