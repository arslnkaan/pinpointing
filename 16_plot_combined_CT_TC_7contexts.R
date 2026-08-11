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
  "CSB_UV_CtoT_trinucleotide_timecourse",
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
  "CSB_UV_combined_CT_TC_7contexts"
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

# Seven unique CT/TC-associated contexts.
#
# TCT occurs only once in repair.
#
# For Damage-seq:
# TCT damage = CT:TCT + TC:TCT.
context_order <- c(
  "ACT",
  "CCT",
  "GCT",
  "TCA",
  "TCC",
  "TCG",
  "TCT"
)

context_colors <- c(
  "ACT" = "#4F759B",
  "CCT" = "#E69F00",
  "GCT" = "#009E73",
  "TCA" = "#D55E00",
  "TCC" = "#CC79A7",
  "TCG" = "#56B4E9",
  "TCT" = "#7A7A7A"
)

# ============================================================
# CHECK INPUTS
# ============================================================

for (input_file in c(repair_file, damage_file)) {
  if (!file.exists(input_file)) {
    stop(
      paste0(
        "Missing input file: ",
        input_file
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

repair_input <- read_tsv(
  repair_file,
  show_col_types = FALSE
)

# Add replicate if absent.
if (!"replicate" %in% names(repair_input)) {
  repair_input <- repair_input %>%
    mutate(
      replicate = "R1"
    )
}

repair_raw <- repair_input %>%
  transmute(
    sample = as.character(sample),
    timepoint = as.character(timepoint),
    time_h = as.numeric(time_h),
    replicate = as.character(replicate),

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

# ============================================================
# COMBINED SEVEN-CONTEXT REPAIR
#
# ACT + CCT + GCT + TCA + TCC + TCG + TCT = 100%
# ============================================================

repair_combined <- crossing(
  sample_design,
  trinucleotide = context_order
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
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_combined_repair = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_fraction = safe_fraction(
      repair_count,
      total_combined_repair
    ),

    repair_percent =
      100 * repair_fraction
  ) %>%
  ungroup() %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    trinucleotide = factor(
      trinucleotide,
      levels = context_order
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
# COMBINED SEVEN-CONTEXT DAMAGE
#
# ACT, CCT, GCT:
#   obtained from CT lesions.
#
# TCA, TCC, TCG:
#   obtained from TC lesions.
#
# TCT:
#   CT:TCT + TC:TCT.
#
# The seven final contexts sum to 100%.
# ============================================================

damage_selected <- damage_raw %>%
  filter(
    (
      lesion_type == "CT" &
        trinucleotide %in%
          c(
            "ACT",
            "CCT",
            "GCT",
            "TCT"
          )
    ) |
      (
        lesion_type == "TC" &
          trinucleotide %in%
            c(
              "TCA",
              "TCC",
              "TCG",
              "TCT"
            )
      )
  ) %>%
  group_by(
    trinucleotide
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_components = paste(
      paste0(
        lesion_type,
        ":",
        trinucleotide
      ),
      collapse = "+"
    ),

    .groups = "drop"
  )

damage_combined <- tibble(
  trinucleotide = context_order
) %>%
  left_join(
    damage_selected,
    by = "trinucleotide"
  ) %>%
  mutate(
    damage_count = replace_na(
      damage_count,
      0
    ),

    damage_components = case_when(
      !is.na(damage_components) ~
        damage_components,

      trinucleotide %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~
        paste0(
          "CT:",
          trinucleotide
        ),

      trinucleotide %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~
        paste0(
          "TC:",
          trinucleotide
        ),

      trinucleotide == "TCT" ~
        "CT:TCT+TC:TCT",

      TRUE ~
        NA_character_
    ),

    total_combined_damage = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_fraction = safe_fraction(
      damage_count,
      total_combined_damage
    ),

    damage_percent =
      100 * damage_fraction,

    damage_timepoint = "0h",

    trinucleotide = factor(
      trinucleotide,
      levels = context_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

# ============================================================
# NORMALIZE REPAIR TO DAMAGE
#
# log2[
#   (context repair / total seven-context repair)
#   /
#   (context damage / total seven-context damage)
# ]
# ============================================================

normalized <- repair_combined %>%
  left_join(
    damage_combined %>%
      transmute(
        trinucleotide =
          as.character(trinucleotide),

        damage_count,
        damage_components,
        total_combined_damage,
        damage_fraction,
        damage_percent
      ),
    by = "trinucleotide"
  ) %>%
  mutate(
    repair_over_damage = safe_ratio(
      repair_fraction,
      damage_fraction
    ),

    log2_repair_over_damage = safe_log2(
      repair_over_damage
    ),

    log2_for_plot = finite_for_plot(
      log2_repair_over_damage
    ),

    normalization_status = case_when(
      is.na(damage_fraction) ~
        "missing_damage",

      damage_fraction <= 0 ~
        "zero_damage",

      repair_count == 0 ~
        "zero_repair",

      TRUE ~
        "normalized"
    ),

    trinucleotide = factor(
      as.character(trinucleotide),
      levels = context_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

# ============================================================
# WRITE MAIN TABLES
# ============================================================

write_tsv(
  repair_combined,
  file.path(
    outdir,
    "combined_7contexts_repair_over_time.tsv"
  )
)

write_tsv(
  damage_combined,
  file.path(
    outdir,
    "combined_7contexts_damage_0h.tsv"
  )
)

write_tsv(
  normalized,
  file.path(
    outdir,
    "combined_7contexts_log2_repair_over_damage.tsv"
  )
)

# ============================================================
# EARLY/LATE TABLE
#
# Early = 0.5 h
# Late  = 8 h
#
# Positive log2(early/late):
#   early-enriched repair.
#
# Negative log2(early/late):
#   late-enriched repair.
#
# Since the same 0 h Damage-seq denominator is used for early
# and late repair, damage cancels:
#
# log2[
#   (repair early / damage) /
#   (repair late / damage)
# ]
#
# =
#
# log2(repair early / repair late)
# ============================================================

early_late_source <- normalized %>%
  mutate(
    timepoint = as.character(timepoint),
    trinucleotide = as.character(trinucleotide)
  ) %>%
  group_by(
    trinucleotide,
    timepoint
  ) %>%
  summarise(
    repair_count = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_fraction = mean(
      repair_fraction,
      na.rm = TRUE
    ),

    repair_percent = mean(
      repair_percent,
      na.rm = TRUE
    ),

    log2_repair_over_damage = mean(
      log2_repair_over_damage,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

early_late <- early_late_source %>%
  select(
    trinucleotide,
    timepoint,
    repair_count,
    repair_fraction,
    repair_percent,
    log2_repair_over_damage
  ) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = c(
      repair_count,
      repair_fraction,
      repair_percent,
      log2_repair_over_damage
    ),
    names_sep = "_"
  ) %>%
  mutate(
    # --------------------------------------------------------
    # Raw early-to-late ratio
    # --------------------------------------------------------

    repair_early_over_late = safe_ratio(
      `repair_fraction_0.5h`,
      repair_fraction_8h
    ),

    repair_log2_early_over_late = safe_log2(
      repair_early_over_late
    ),

    # --------------------------------------------------------
    # Damage-normalized early-to-late ratio
    #
    # This should equal repair_log2_early_over_late.
    # --------------------------------------------------------

    normalized_log2_early_over_late =
      `log2_repair_over_damage_0.5h` -
      log2_repair_over_damage_8h,

    # --------------------------------------------------------
    # Opposite direction: late minus early
    # --------------------------------------------------------

    normalized_change_8h_minus_0.5h =
      log2_repair_over_damage_8h -
      `log2_repair_over_damage_0.5h`,

    # --------------------------------------------------------
    # Percentage-point change
    # --------------------------------------------------------

    repair_percent_change_8h_minus_0.5h =
      repair_percent_8h -
      `repair_percent_0.5h`,

    repair_percent_change_0.5h_minus_8h =
      `repair_percent_0.5h` -
      repair_percent_8h,

    # --------------------------------------------------------
    # Plotting column
    #
    # Infinite values remain in the TSV but are not plotted.
    # --------------------------------------------------------

    repair_log2_early_over_late_for_plot =
      finite_for_plot(
        repair_log2_early_over_late
      ),

    early_late_status = case_when(
      `repair_fraction_0.5h` == 0 &
        repair_fraction_8h == 0 ~
        "zero_early_and_late",

      `repair_fraction_0.5h` == 0 ~
        "zero_early",

      repair_fraction_8h == 0 ~
        "zero_late",

      is.na(
        repair_log2_early_over_late
      ) ~
        "not_calculable",

      repair_log2_early_over_late > 0 ~
        "early_enriched",

      repair_log2_early_over_late < 0 ~
        "late_enriched",

      TRUE ~
        "unchanged"
    ),

    # These should be approximately identical.
    raw_minus_normalized_early_late =
      repair_log2_early_over_late -
      normalized_log2_early_over_late
  ) %>%
  mutate(
    trinucleotide = factor(
      trinucleotide,
      levels = context_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

write_tsv(
  early_late,
  file.path(
    outdir,
    "combined_7contexts_early_late_changes.tsv"
  )
)

# ============================================================
# PLOT LIMITS
# ============================================================

raw_y_max <- max(
  repair_combined$repair_percent,
  na.rm = TRUE
)

if (!is.finite(raw_y_max) ||
    raw_y_max <= 0) {
  raw_y_max <- 100
}

raw_y_max <- raw_y_max * 1.08

damage_y_max <- max(
  damage_combined$damage_percent,
  na.rm = TRUE
)

if (!is.finite(damage_y_max) ||
    damage_y_max <= 0) {
  damage_y_max <- 100
}

damage_y_max <- damage_y_max * 1.08

normalized_values <-
  normalized$log2_for_plot

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

normalized_limit <-
  normalized_limit * 1.10

early_late_values <-
  early_late$repair_log2_early_over_late_for_plot

early_late_values <- early_late_values[
  is.finite(early_late_values)
]

if (length(early_late_values) > 0) {
  early_late_limit <- max(
    abs(early_late_values),
    na.rm = TRUE
  )
} else {
  early_late_limit <- 1
}

if (!is.finite(early_late_limit) ||
    early_late_limit <= 0) {
  early_late_limit <- 1
}

early_late_limit <-
  early_late_limit * 1.10

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
# PLOT 1
# COMBINED SEVEN-CONTEXT REPAIR OVER TIME
# ============================================================

p_repair <- ggplot(
  repair_combined,
  aes(
    x = time_h,
    y = repair_percent,
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
    values = context_colors,
    breaks = context_order,
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
    title = "CT and TC repair contexts",
    x = "Repair time (h)",
    y = "Repair composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 2
# COMBINED SEVEN-CONTEXT DAMAGE AT 0 h
# ============================================================

p_damage <- ggplot(
  damage_combined,
  aes(
    x = trinucleotide,
    y = damage_percent,
    fill = trinucleotide
  )
) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.2
  ) +
  scale_fill_manual(
    values = context_colors,
    breaks = context_order,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = context_order,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      0,
      damage_y_max
    )
  ) +
  labs(
    title = "CT and TC Damage-seq at 0 h",
    x = NULL,
    y = "Damage composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 3
# COMBINED SEVEN-CONTEXT NORMALIZATION
# ============================================================

p_normalized <- ggplot(
  normalized,
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
    values = context_colors,
    breaks = context_order,
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
    title = "CT and TC repair normalized to damage",
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
# PLOT 4
# LOG2 EARLY/LATE REPAIR
#
# Early = 0.5 h
# Late  = 8 h
#
# Positive:
#   early-enriched
#
# Negative:
#   late-enriched
# ============================================================

p_early_late <- ggplot(
  early_late,
  aes(
    x = trinucleotide,
    y = repair_log2_early_over_late_for_plot,
    fill = trinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_col(
    width = 0.72,
    color = "black",
    linewidth = 0.2,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = context_colors,
    breaks = context_order,
    guide = "none"
  ) +
  scale_x_discrete(
    limits = context_order,
    drop = FALSE
  ) +
  coord_cartesian(
    ylim = c(
      -early_late_limit,
      early_late_limit
    )
  ) +
  labs(
    title = "Early versus late repair",
    x = NULL,
    y = expression(
      log[2](
        repair["0.5 h"] /
          repair["8 h"]
      )
    )
  ) +
  common_theme

# ============================================================
# COMBINED FOUR-PANEL FIGURE
# ============================================================

combined_plot <- (
  p_repair /
    p_damage /
    p_normalized /
    p_early_late
) +
  plot_layout(
    guides = "collect",
    heights = c(
      1,
      0.9,
      1,
      0.9
    )
  ) &
  theme(
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plotdir,
    "repair_damage_normalization_early_late.pdf"
  ),
  plot = combined_plot,
  width = 8,
  height = 15,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "repair_damage_normalization_early_late.png"
  ),
  plot = combined_plot,
  width = 8,
  height = 15,
  units = "in",
  dpi = 600
)

# ============================================================
# SAVE INDIVIDUAL PLOTS
# ============================================================

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_repair_over_time.pdf"
  ),
  plot = p_repair,
  width = 7,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_repair_over_time.png"
  ),
  plot = p_repair,
  width = 7,
  height = 4.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_damage_0h.pdf"
  ),
  plot = p_damage,
  width = 6,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_damage_0h.png"
  ),
  plot = p_damage,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_log2_repair_over_damage.pdf"
  ),
  plot = p_normalized,
  width = 7,
  height = 4.8,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_log2_repair_over_damage.png"
  ),
  plot = p_normalized,
  width = 7,
  height = 4.8,
  units = "in",
  dpi = 600
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_log2_early_over_late.pdf"
  ),
  plot = p_early_late,
  width = 6,
  height = 4.5,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "CT_TC_log2_early_over_late.png"
  ),
  plot = p_early_late,
  width = 6,
  height = 4.5,
  units = "in",
  dpi = 600
)

# ============================================================
# QC
# ============================================================

repair_qc <- repair_combined %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  summarise(
    repair_fraction_sum = sum(
      repair_fraction,
      na.rm = TRUE
    ),

    repair_percent_sum = sum(
      repair_percent,
      na.rm = TRUE
    ),

    total_combined_repair = first(
      total_combined_repair
    ),

    .groups = "drop"
  )

damage_qc <- damage_combined %>%
  summarise(
    damage_fraction_sum = sum(
      damage_fraction,
      na.rm = TRUE
    ),

    damage_percent_sum = sum(
      damage_percent,
      na.rm = TRUE
    ),

    total_combined_damage = first(
      total_combined_damage
    )
  )

early_late_qc <- early_late %>%
  transmute(
    trinucleotide,
    repair_log2_early_over_late,
    normalized_log2_early_over_late,
    raw_minus_normalized_early_late,
    early_late_status
  )

write_tsv(
  repair_qc,
  file.path(
    outdir,
    "combined_7contexts_repair_QC.tsv"
  )
)

write_tsv(
  damage_qc,
  file.path(
    outdir,
    "combined_7contexts_damage_QC.tsv"
  )
)

write_tsv(
  early_late_qc,
  file.path(
    outdir,
    "combined_7contexts_early_late_QC.tsv"
  )
)

# ============================================================
# FINISHED
# ============================================================

cat("\nDone.\n\n")

cat("Repair table:\n")
cat(
  file.path(
    outdir,
    "combined_7contexts_repair_over_time.tsv"
  ),
  "\n\n"
)

cat("Damage table:\n")
cat(
  file.path(
    outdir,
    "combined_7contexts_damage_0h.tsv"
  ),
  "\n\n"
)

cat("Normalized table:\n")
cat(
  file.path(
    outdir,
    "combined_7contexts_log2_repair_over_damage.tsv"
  ),
  "\n\n"
)

cat("Early/late table:\n")
cat(
  file.path(
    outdir,
    "combined_7contexts_early_late_changes.tsv"
  ),
  "\n\n"
)

cat("Combined figure:\n")
cat(
  file.path(
    plotdir,
    "repair_damage_normalization_early_late.pdf"
  ),
  "\n\n"
)

cat("Early/late figure:\n")
cat(
  file.path(
    plotdir,
    "CT_TC_log2_early_over_late.pdf"
  ),
  "\n"
)