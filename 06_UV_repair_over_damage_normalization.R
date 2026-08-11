#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

# ============================================================
# SETTINGS
# ============================================================

base_dir <- "/work/users/a/r/arslank"

outdir <- file.path(
  base_dir,
  "CSB_UV_repair_over_damage_normalization"
)

plotdir <- file.path(
  outdir,
  "plots"
)

dir.create(
  plotdir,
  showWarnings = FALSE,
  recursive = TRUE
)

uv_red <- "#C62828"

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

# Full SBS order displayed in the C>T plots.
sbs_order <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

# Only these contexts can be matched to the current CT/TC
# Damage-seq denominator.
ctot_normalizable_contexts <- c(
  "ACT",
  "CCT",
  "GCT",
  "TCA",
  "TCC",
  "TCG",
  "TCT"
)

nccn_order <- c(
  "ACCA", "ACCC", "ACCG", "ACCT",
  "CCCA", "CCCC", "CCCG", "CCCT",
  "GCCA", "GCCC", "GCCG", "GCCT",
  "TCCA", "TCCC", "TCCG", "TCCT"
)

# ============================================================
# INPUT FILES
# ============================================================

ctot_repair_file <- file.path(
  base_dir,
  "CSB_UV_CtoT_trinucleotide_timecourse",
  "UV_CtoT_trinucleotide_counts_percentages.tsv"
)

ct_tc_damage_file <- file.path(
  base_dir,
  "damseq",
  "NHF1_CPD_0h_r1_results",
  "damage_context_denominators",
  "DamageSeq_CT_TC_trinucleotide_counts.tsv"
)

cc_damage_file <- file.path(
  base_dir,
  "damseq",
  "NHF1_CPD_0h_r1_results",
  "damage_CC_NCCN_from_read_starts",
  "NHF1_CPD_0h.minus3_to_0_CC_NCCN_counts.tsv"
)

cctt_samples <- tribble(
  ~sample,            ~timepoint, ~time_h, ~replicate,
  "CSB-UVCPD-30m",   "0.5h",        0.5, "R1",
  "CSB-UVCPD-2h",    "2h",          2.0, "R1",
  "CSB-UVCPD-4h",    "4h",          4.0, "R1",
  "CSB-UVCPD-8h",    "8h",          8.0, "R1"
) %>%
  mutate(
    file = file.path(
      base_dir,
      paste0(
        sample,
        "_CC_to_TT_rescue"
      ),
      "full_reads_and_NCCN",
      paste0(
        sample,
        ".CCTT_best_per_original_read.tsv"
      )
    )
  )

# ============================================================
# CHECK INPUTS
# ============================================================

input_files <- c(
  ctot_repair_file,
  ct_tc_damage_file,
  cc_damage_file,
  cctt_samples$file
)

missing_files <- input_files[
  !file.exists(input_files)
]

if (length(missing_files) > 0) {
  stop(
    paste(
      "Missing input files:",
      paste(
        missing_files,
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}

# ============================================================
# HELPERS
# ============================================================

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

find_first_column <- function(
  dat,
  candidates,
  label
) {
  found <- candidates[
    candidates %in% colnames(dat)
  ]

  if (length(found) == 0) {
    stop(
      paste0(
        "Could not identify ",
        label,
        " column. Available columns: ",
        paste(
          colnames(dat),
          collapse = ", "
        )
      )
    )
  }

  found[[1]]
}

# ============================================================
# PART 1
# C>T REPAIR
# ============================================================

message(
  "Reading C>T repair table: ",
  ctot_repair_file
)

ctot_repair_raw <- read_tsv(
  ctot_repair_file,
  show_col_types = FALSE
) %>%
  mutate(
    trinucleotide = str_to_upper(
      as.character(trinucleotide)
    ),
    count = as.numeric(count),
    timepoint = as.character(timepoint)
  )

required_ctot_columns <- c(
  "sample",
  "timepoint",
  "time_h",
  "replicate",
  "trinucleotide",
  "count"
)

missing_ctot_columns <- setdiff(
  required_ctot_columns,
  colnames(ctot_repair_raw)
)

if (length(missing_ctot_columns) > 0) {
  stop(
    paste0(
      "C>T repair table is missing: ",
      paste(
        missing_ctot_columns,
        collapse = ", "
      )
    )
  )
}

sample_design <- ctot_repair_raw %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate
  )

# Keep only contexts that can be matched to CT/TC damage
# when calculating the repair denominator.
ctot_repair_used <- ctot_repair_raw %>%
  filter(
    trinucleotide %in%
      ctot_normalizable_contexts
  ) %>%
  mutate(
    damage_orientation = case_when(
      trinucleotide %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ "CT",

      trinucleotide %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ "TC",

      trinucleotide == "TCT" ~
        "CT+TC",

      TRUE ~ NA_character_
    )
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_CtoT_CPD_compatible = sum(
      count,
      na.rm = TRUE
    ),

    repair_fraction = if_else(
      total_CtoT_CPD_compatible > 0,
      count /
        total_CtoT_CPD_compatible,
      NA_real_
    )
  ) %>%
  ungroup()

# ============================================================
# CT/TC DAMAGE
# ============================================================

message(
  "Reading CT/TC Damage-seq table: ",
  ct_tc_damage_file
)

ct_tc_damage_raw <- read_tsv(
  ct_tc_damage_file,
  show_col_types = FALSE
)

ct_tc_count_column <- find_first_column(
  ct_tc_damage_raw,
  candidates = c(
    "count",
    "weighted_count"
  ),
  label = "CT/TC damage count"
)

ct_tc_damage <- ct_tc_damage_raw %>%
  transmute(
    lesion_type = str_to_upper(
      as.character(lesion_type)
    ),

    trinucleotide = str_to_upper(
      as.character(trinucleotide)
    ),

    damage_count = as.numeric(
      .data[[ct_tc_count_column]]
    )
  ) %>%
  filter(
    lesion_type %in%
      c(
        "CT",
        "TC"
      )
  ) %>%
  group_by(
    lesion_type,
    trinucleotide
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  )

total_ct_tc_damage <- sum(
  ct_tc_damage$damage_count,
  na.rm = TRUE
)

if (total_ct_tc_damage <= 0) {
  stop(
    "The total CT/TC Damage-seq count is zero."
  )
}

ct_damage_normalization <- ct_tc_damage %>%
  filter(
    lesion_type == "CT",
    trinucleotide %in%
      c(
        "ACT",
        "CCT",
        "GCT"
      )
  ) %>%
  transmute(
    trinucleotide,
    damage_orientation = "CT",
    damage_count
  )

tc_damage_normalization <- ct_tc_damage %>%
  filter(
    lesion_type == "TC",
    trinucleotide %in%
      c(
        "TCA",
        "TCC",
        "TCG"
      )
  ) %>%
  transmute(
    trinucleotide,
    damage_orientation = "TC",
    damage_count
  )

tct_damage_normalization <- ct_tc_damage %>%
  filter(
    lesion_type %in%
      c(
        "CT",
        "TC"
      ),
    trinucleotide == "TCT"
  ) %>%
  summarise(
    trinucleotide = "TCT",
    damage_orientation = "CT+TC",
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    )
  )

ctot_damage_denominators <- bind_rows(
  ct_damage_normalization,
  tc_damage_normalization,
  tct_damage_normalization
) %>%
  complete(
    trinucleotide =
      ctot_normalizable_contexts,
    fill = list(
      damage_count = 0
    )
  ) %>%
  mutate(
    damage_orientation = case_when(
      trinucleotide %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ "CT",

      trinucleotide %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ "TC",

      trinucleotide == "TCT" ~
        "CT+TC"
    ),

    damage_fraction =
      damage_count /
      total_ct_tc_damage
  )

write_tsv(
  ctot_damage_denominators,
  file.path(
    outdir,
    "DamageSeq_CT_TC_normalization_denominators.tsv"
  )
)

# ============================================================
# NORMALIZE C>T REPAIR
# ============================================================

ctot_normalized_used <- ctot_repair_used %>%
  left_join(
    ctot_damage_denominators,
    by = c(
      "trinucleotide",
      "damage_orientation"
    )
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

    repair_percent =
      100 * repair_fraction,

    damage_percent =
      100 * damage_fraction,

    normalization_status = case_when(
      is.na(damage_count) ~
        "missing_damage_context",

      damage_count <= 0 ~
        "zero_damage",

      count == 0 ~
        "zero_repair",

      TRUE ~
        "normalized"
    )
  )

# Add all 16 SBS contexts to every time point.
# The nine contexts not matched to CT/TC damage remain NA.
ctot_normalized <- crossing(
  sample_design,
  trinucleotide = sbs_order
) %>%
  left_join(
    ctot_repair_raw %>%
      select(
        sample,
        timepoint,
        time_h,
        replicate,
        trinucleotide,
        count
      ),
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate",
      "trinucleotide"
    )
  ) %>%
  left_join(
    ctot_normalized_used %>%
      select(
        sample,
        timepoint,
        time_h,
        replicate,
        trinucleotide,
        damage_orientation,
        total_CtoT_CPD_compatible,
        repair_fraction,
        damage_count,
        damage_fraction,
        repair_over_damage,
        log2_repair_over_damage,
        log2_for_plot,
        repair_percent,
        damage_percent,
        normalization_status
      ),
    by = c(
      "sample",
      "timepoint",
      "time_h",
      "replicate",
      "trinucleotide"
    )
  ) %>%
  mutate(
    count = replace_na(
      count,
      0
    ),

    normalization_status = case_when(
      !trinucleotide %in%
        ctot_normalizable_contexts ~
        "not_CT_TC_normalizable",

      TRUE ~
        normalization_status
    ),

    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    trinucleotide = factor(
      trinucleotide,
      levels = sbs_order
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

write_tsv(
  ctot_normalized,
  file.path(
    outdir,
    "UV_CtoT_log2_repair_over_damage.tsv"
  )
)

# Retain the old filename as well so downstream paths do not break.
write_tsv(
  ctot_normalized,
  file.path(
    outdir,
    "UV_CtoT_repair_over_damage.tsv"
  )
)

# ============================================================
# C>T LOG2 CHANGE TABLE
# ============================================================

ctot_change <- ctot_normalized %>%
  filter(
    as.character(trinucleotide) %in%
      ctot_normalizable_contexts
  ) %>%
  mutate(
    timepoint = as.character(timepoint),
    trinucleotide = as.character(trinucleotide)
  ) %>%
  select(
    trinucleotide,
    timepoint,
    log2_repair_over_damage
  ) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = log2_repair_over_damage,
    names_prefix = "log2_repair_over_damage_"
  ) %>%
  mutate(
    log2_change_8h_minus_0.5h =
      log2_repair_over_damage_8h -
      `log2_repair_over_damage_0.5h`
  ) %>%
  mutate(
    trinucleotide = factor(
      trinucleotide,
      levels = sbs_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

write_tsv(
  ctot_change,
  file.path(
    outdir,
    "UV_CtoT_log2_repair_over_damage_changes.tsv"
  )
)

# ============================================================
# PART 2
# CC>TT REPAIR
# ============================================================

read_cctt_sample <- function(
  sample,
  timepoint,
  time_h,
  replicate,
  file
) {
  message(
    "Reading CC>TT repair table: ",
    file
  )

  dat <- read_tsv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )

  context_column <- find_first_column(
    dat,
    candidates = c(
      "genome_NCCN",
      "NCCN",
      "genome_tetranucleotide"
    ),
    label = "CC>TT NCCN"
  )

  if (!"reconstruction_status" %in%
      colnames(dat)) {
    stop(
      paste0(
        "Missing reconstruction_status in ",
        file
      )
    )
  }

  dat %>%
    transmute(
      sample = sample,
      timepoint = timepoint,
      time_h = time_h,
      replicate = replicate,

      reconstruction_status =
        as.character(
          reconstruction_status
        ),

      NCCN = str_to_upper(
        as.character(
          .data[[context_column]]
        )
      )
    ) %>%
    filter(
      reconstruction_status == "ok",
      NCCN %in% nccn_order
    )
}

cctt_events <- pmap_dfr(
  cctt_samples,
  read_cctt_sample
)

cctt_repair_counts <- cctt_events %>%
  count(
    sample,
    timepoint,
    time_h,
    replicate,
    NCCN,
    name = "repair_count"
  ) %>%
  complete(
    nesting(
      sample,
      timepoint,
      time_h,
      replicate
    ),
    NCCN = nccn_order,
    fill = list(
      repair_count = 0L
    )
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_CCTT_repair = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_fraction = if_else(
      total_CCTT_repair > 0,
      repair_count /
        total_CCTT_repair,
      NA_real_
    )
  ) %>%
  ungroup()

# ============================================================
# EXACT NCCN DAMAGE
# ============================================================

message(
  "Reading CC Damage-seq table: ",
  cc_damage_file
)

cc_damage_raw <- read_tsv(
  cc_damage_file,
  show_col_types = FALSE
)

cc_count_column <- find_first_column(
  cc_damage_raw,
  candidates = c(
    "count",
    "weighted_count"
  ),
  label = "CC damage count"
)

cc_damage <- cc_damage_raw %>%
  transmute(
    NCCN = str_to_upper(
      as.character(NCCN)
    ),

    damage_count = as.numeric(
      .data[[cc_count_column]]
    )
  ) %>%
  filter(
    NCCN %in% nccn_order
  ) %>%
  group_by(
    NCCN
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),
    .groups = "drop"
  ) %>%
  complete(
    NCCN = nccn_order,
    fill = list(
      damage_count = 0
    )
  )

total_cc_damage <- sum(
  cc_damage$damage_count,
  na.rm = TRUE
)

if (total_cc_damage <= 0) {
  stop(
    "The total exact NCCN Damage-seq count is zero."
  )
}

cc_damage <- cc_damage %>%
  mutate(
    total_CC_damage =
      total_cc_damage,

    damage_fraction =
      damage_count /
      total_CC_damage
  )

write_tsv(
  cc_damage,
  file.path(
    outdir,
    "DamageSeq_CC_NCCN_normalization_denominators.tsv"
  )
)

# ============================================================
# NORMALIZE CC>TT
# ============================================================

cctt_normalized <- cctt_repair_counts %>%
  left_join(
    cc_damage,
    by = "NCCN"
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

    repair_percent =
      100 * repair_fraction,

    damage_percent =
      100 * damage_fraction,

    normalization_status = case_when(
      is.na(damage_count) ~
        "missing_damage_context",

      damage_count <= 0 ~
        "zero_damage",

      repair_count == 0 ~
        "zero_repair",

      TRUE ~
        "normalized"
    ),

    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    NCCN = factor(
      NCCN,
      levels = nccn_order
    )
  ) %>%
  arrange(
    time_h,
    NCCN
  )

write_tsv(
  cctt_normalized,
  file.path(
    outdir,
    "UV_CCTT_NCCN_log2_repair_over_damage.tsv"
  )
)

write_tsv(
  cctt_normalized,
  file.path(
    outdir,
    "UV_CCTT_NCCN_repair_over_damage.tsv"
  )
)

# ============================================================
# QC SUMMARY
# ============================================================

ctot_qc <- ctot_normalized_used %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate,
    total_CtoT_CPD_compatible
  ) %>%
  transmute(
    analysis = "C>T",
    sample,
    timepoint,
    time_h,
    replicate,
    total_repair_events =
      total_CtoT_CPD_compatible,
    total_damage_events =
      total_ct_tc_damage
  )

cctt_qc <- cctt_normalized %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate,
    total_CCTT_repair
  ) %>%
  transmute(
    analysis = "CC>TT",
    sample,
    timepoint,
    time_h,
    replicate,
    total_repair_events =
      total_CCTT_repair,
    total_damage_events =
      total_cc_damage
  )

normalization_qc <- bind_rows(
  ctot_qc,
  cctt_qc
) %>%
  arrange(
    analysis,
    time_h
  )

write_tsv(
  normalization_qc,
  file.path(
    outdir,
    "UV_repair_over_damage_QC_summary.tsv"
  )
)

# ============================================================
# COMMON PLOT THEME
# ============================================================

common_theme <- theme_classic(
  base_size = 11
) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey60",
      linewidth = 0.35
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.35
    ),

    plot.title = element_text(
      hjust = 0.5
    )
  )

# ============================================================
# PLOT 1
# C>T LOG2 REPAIR/DAMAGE — FULL SBS ORDER
# ============================================================

p_ctot <- ggplot(
  ctot_normalized,
  aes(
    x = trinucleotide,
    y = log2_for_plot
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_col(
    fill = uv_red,
    width = 0.8,
    color = "black",
    linewidth = 0.15,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 2
  ) +
  scale_x_discrete(
    limits = sbs_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.06,
        0.08
      )
    )
  ) +
  labs(
    title = "UV C>T repair normalized to CT/TC damage",
    subtitle = paste(
      "Full SBS order;",
      "blank contexts are not CT/TC-normalizable"
    ),
    x = NULL,
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  common_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

ggsave(
  filename = file.path(
    plotdir,
    "UV_CtoT_log2_repair_over_damage_SBS_order.pdf"
  ),
  plot = p_ctot,
  width = 8.0,
  height = 6.3,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "UV_CtoT_log2_repair_over_damage_SBS_order.png"
  ),
  plot = p_ctot,
  width = 8.0,
  height = 6.3,
  units = "in",
  dpi = 600
)

# ============================================================
# PLOT 2
# C>T LOG2 TRAJECTORIES
# ============================================================

p_ctot_trajectory <- ctot_normalized %>%
  filter(
    as.character(trinucleotide) %in%
      ctot_normalizable_contexts
  ) %>%
  ggplot(
    aes(
      x = time_h,
      y = log2_for_plot
    )
  ) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey50"
  ) +
  geom_line(
    color = uv_red,
    linewidth = 0.7,
    na.rm = TRUE
  ) +
  geom_point(
    color = uv_red,
    size = 1.8,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~trinucleotide,
    ncol = 4,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  labs(
    title = "Damage-normalized C>T trajectories",
    x = "Repair time (h)",
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey65"
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey70",
      linewidth = 0.3
    ),

    plot.title = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    plotdir,
    "UV_CtoT_log2_repair_over_damage_trajectories.pdf"
  ),
  plot = p_ctot_trajectory,
  width = 7.2,
  height = 5.8,
  units = "in"
)

# ============================================================
# PLOT 3
# CC>TT LOG2 REPAIR/DAMAGE
# ============================================================

p_cctt <- ggplot(
  cctt_normalized,
  aes(
    x = NCCN,
    y = log2_for_plot
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.45,
    color = "grey35"
  ) +
  geom_col(
    fill = uv_red,
    width = 0.8,
    color = "black",
    linewidth = 0.15,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 2
  ) +
  scale_x_discrete(
    limits = nccn_order,
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0.06,
        0.08
      )
    )
  ) +
  labs(
    title = "UV CC>TT repair normalized to CC damage",
    x = "NCCN context",
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  common_theme +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 8
    )
  )

ggsave(
  filename = file.path(
    plotdir,
    "UV_CCTT_NCCN_log2_repair_over_damage.pdf"
  ),
  plot = p_cctt,
  width = 7.5,
  height = 6.5,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "UV_CCTT_NCCN_log2_repair_over_damage.png"
  ),
  plot = p_cctt,
  width = 7.5,
  height = 6.5,
  units = "in",
  dpi = 600
)

# ============================================================
# PLOT 4
# CC>TT LOG2 TRAJECTORIES
# ============================================================

p_cctt_trajectory <- ggplot(
  cctt_normalized,
  aes(
    x = time_h,
    y = log2_for_plot
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey50"
  ) +
  geom_line(
    color = uv_red,
    linewidth = 0.65,
    na.rm = TRUE
  ) +
  geom_point(
    color = uv_red,
    size = 1.5,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~NCCN,
    ncol = 4,
    drop = FALSE
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    )
  ) +
  labs(
    title = "Damage-normalized CC>TT trajectories",
    x = "Repair time (h)",
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey65"
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey70",
      linewidth = 0.3
    ),

    plot.title = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    plotdir,
    "UV_CCTT_NCCN_log2_repair_over_damage_trajectories.pdf"
  ),
  plot = p_cctt_trajectory,
  width = 7.5,
  height = 7.2,
  units = "in"
)

# ============================================================
# FINISHED
# ============================================================

cat("\nDone.\n\n")

cat("C>T log2-normalized table:\n")
cat(
  file.path(
    outdir,
    "UV_CtoT_log2_repair_over_damage.tsv"
  ),
  "\n\n"
)

cat("CC>TT log2-normalized table:\n")
cat(
  file.path(
    outdir,
    "UV_CCTT_NCCN_log2_repair_over_damage.tsv"
  ),
  "\n\n"
)

cat("Plots:\n")
cat(
  plotdir,
  "\n"
)