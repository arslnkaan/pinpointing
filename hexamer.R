#!/usr/bin/env Rscript

# ============================================================
# WT UV-CPD HEXANUCLEOTIDE DAMAGE-SEQ NORMALIZATION
#
# This script begins AFTER the memory-safe streaming extraction.
#
# It reads:
#
#   1. WT XR-seq hexamer counts and percentages
#
#   2. The streamed Damage-seq hexamer count table:
#
#      DamageSeq_0h_hexanucleotide_counts_percent_
#      from_damage_tetranucleotide_windows.tsv
#
# It does NOT:
#
#   - read the 38-million-row Damage-seq BED
#   - run bedtools
#   - generate a large FASTA
#   - load Damage-seq sequences into Biostrings
#
# Normalization:
#
#   repair_over_damage =
#       WT XR-seq hexamer percentage
#       --------------------------------
#       Damage-seq hexamer percentage
#
# Interpretation:
#
#   ratio = 1:
#       repair abundance matches starting damage abundance
#
#   ratio > 1:
#       repair enrichment relative to starting damage abundance
#
#   ratio < 1:
#       repair depletion relative to starting damage abundance
#
# This remains a relative repair index, not the literal
# percentage of genomic lesions removed.
# ============================================================

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# PATHS
# ============================================================

base_dir <- "/work/users/a/r/arslank"

input_dir <- file.path(
  base_dir,
  "WT_UV_CPD_hexanucleotide_timecourse"
)

outdir <- file.path(
  input_dir,
  "damage_seq_normalized_visualizations"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

repair_replicate_file <- file.path(
  input_dir,
  "WT_hexanucleotide_counts_percent_per_replicate.tsv"
)

repair_mean_sem_file <- file.path(
  input_dir,
  "WT_hexanucleotide_mean_sem.tsv"
)

damage_hexamer_file <- file.path(
  input_dir,
  paste0(
    "DamageSeq_0h_hexanucleotide_counts_percent_",
    "from_damage_tetranucleotide_windows.tsv"
  )
)

damage_qc_file <- file.path(
  input_dir,
  paste0(
    "DamageSeq_0h_hexanucleotide_extraction_QC_",
    "from_damage_tetranucleotide_windows.tsv"
  )
)

damage_cpd_summary_file <- file.path(
  input_dir,
  "DamageSeq_0h_hexanucleotide_CPD_summary.tsv"
)

# ============================================================
# PARAMETERS
# ============================================================

# Minimum Damage-seq count required for a context to be used
# as a repair/damage denominator.
min_damage_count <- 20L

# Minimum total repair reads across all eight WT samples.
min_total_repair_count <- 20L

# Minimum number of distinct nonzero repair timepoints.
min_nonzero_timepoints <- 2L

top_n_per_timepoint <- 25L
top_n_dynamic <- 24L
top_n_timing <- 25L
top_n_labels <- 12L

# Used only for the optional pseudocount log2 ratio.
# The primary repair_over_damage measurement has no pseudocount.
ratio_pseudocount <- 1e-6

# ============================================================
# COLORS
# ============================================================

wt_color <- "darkred"

cpd_colors <- c(
  "CC" = "#984EA3",
  "CT" = "#377EB8",
  "TC" = "#4DAF4A",
  "TT" = "#FF7F00"
)

regression_colors <- c(
  "Early-enriched"    = "#2166AC",
  "Late-enriched"     = "#B2182B",
  "Flat"              = "grey65",
  "Insufficient data" = "grey85"
)

# ============================================================
# ORDERS
# ============================================================

time_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

replicate_order <- c(
  "R1",
  "R2"
)

bases <- c(
  "A",
  "C",
  "G",
  "T"
)

cpd_order <- c(
  "CC",
  "CT",
  "TC",
  "TT"
)

dinucleotide_order <- expand_grid(
  base_1 = bases,
  base_2 = bases
) %>%
  transmute(
    dinucleotide = paste0(
      base_1,
      base_2
    )
  ) %>%
  pull(
    dinucleotide
  )

hexamer_key <- expand_grid(
  left_flank = dinucleotide_order,
  cpd = cpd_order,
  right_flank = dinucleotide_order
) %>%
  mutate(
    hexamer = paste0(
      left_flank,
      cpd,
      right_flank
    )
  ) %>%
  select(
    hexamer,
    left_flank,
    cpd,
    right_flank
  )

hexamer_order <- hexamer_key$hexamer

if (
  length(hexamer_order) != 1024 ||
  n_distinct(hexamer_order) != 1024
) {
  stop(
    "Failed to generate the expected 1,024 hexanucleotides."
  )
}

# ============================================================
# HELPER FUNCTIONS
# ============================================================

check_required_columns <- function(
  data,
  required_columns,
  table_name
) {

  missing_columns <- setdiff(
    required_columns,
    names(data)
  )

  if (length(missing_columns) > 0) {
    stop(
      table_name,
      " is missing required columns:\n  ",
      paste(
        missing_columns,
        collapse = ", "
      ),
      "\n\nAvailable columns:\n  ",
      paste(
        names(data),
        collapse = ", "
      )
    )
  }
}

safe_finite_max <- function(
  x,
  fallback = 1
) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0) {
    return(fallback)
  }

  value <- max(x)

  if (
    !is.finite(value) ||
    value <= 0
  ) {
    return(fallback)
  }

  value
}

safe_finite_quantile <- function(
  x,
  probability,
  fallback = 1
) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) == 0) {
    return(fallback)
  }

  value <- quantile(
    x,
    probs = probability,
    na.rm = TRUE,
    names = FALSE
  )

  if (
    !is.finite(value) ||
    value <= 0
  ) {
    return(fallback)
  }

  value
}

remove_facet_suffix <- function(x) {

  str_remove(
    x,
    "___.*$"
  )
}

save_plot <- function(
  plot_object,
  filename_stub,
  width,
  height,
  dpi = 400
) {

  ggsave(
    filename = file.path(
      outdir,
      paste0(
        filename_stub,
        ".pdf"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    useDingbats = FALSE
  )

  ggsave(
    filename = file.path(
      outdir,
      paste0(
        filename_stub,
        ".png"
      )
    ),
    plot = plot_object,
    width = width,
    height = height,
    units = "in",
    dpi = dpi,
    bg = "white"
  )
}

# ============================================================
# CHECK INPUT FILES
# ============================================================

required_files <- c(
  repair_replicate_file,
  repair_mean_sem_file,
  damage_hexamer_file
)

missing_files <- required_files[
  !file.exists(required_files)
]

if (length(missing_files) > 0) {

  message(
    "\nMissing required input files:"
  )

  walk(
    missing_files,
    ~ message(
      "  ",
      .x
    )
  )

  stop(
    paste0(
      "\nThe streamed Damage-seq table must already exist.\n",
      "Run the memory-safe Damage-seq extraction script first."
    )
  )
}

# ============================================================
# READ WT REPAIR REPLICATE TABLE
# ============================================================

repair_replicates_raw <- read_tsv(
  repair_replicate_file,
  show_col_types = FALSE,
  progress = FALSE
)

check_required_columns(
  repair_replicates_raw,
  c(
    "sample_id",
    "timepoint",
    "time_h",
    "replicate",
    "hexamer",
    "left_flank",
    "cpd",
    "right_flank",
    "count",
    "percent"
  ),
  "WT repair replicate table"
)

repair_replicates <- repair_replicates_raw %>%
  transmute(
    sample_id = as.character(
      sample_id
    ),

    timepoint = factor(
      as.character(timepoint),
      levels = time_order
    ),

    time_h = as.numeric(
      time_h
    ),

    replicate = factor(
      as.character(replicate),
      levels = replicate_order
    ),

    hexamer = str_to_upper(
      str_trim(
        as.character(hexamer)
      )
    ),

    left_flank = factor(
      as.character(left_flank),
      levels = dinucleotide_order
    ),

    cpd = factor(
      as.character(cpd),
      levels = cpd_order
    ),

    right_flank = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    ),

    count = as.numeric(
      count
    ),

    percent = as.numeric(
      percent
    )
  ) %>%
  filter(
    hexamer %in% hexamer_order,
    is.finite(time_h),
    is.finite(count)
  ) %>%
  arrange(
    time_h,
    replicate,
    cpd,
    left_flank,
    right_flank
  )

# ============================================================
# READ WT REPAIR MEAN/SEM TABLE
#
# This table is checked for consistency, although normalized
# mean and SEM values are recalculated from replicate-level data.
# ============================================================

repair_mean_sem_raw <- read_tsv(
  repair_mean_sem_file,
  show_col_types = FALSE,
  progress = FALSE
)

check_required_columns(
  repair_mean_sem_raw,
  c(
    "timepoint",
    "time_h",
    "hexamer",
    "mean_percent",
    "sem_percent",
    "total_count"
  ),
  "WT repair mean/SEM table"
)

repair_mean_sem <- repair_mean_sem_raw %>%
  mutate(
    timepoint = factor(
      as.character(timepoint),
      levels = time_order
    ),

    time_h = as.numeric(
      time_h
    ),

    hexamer = str_to_upper(
      str_trim(
        as.character(hexamer)
      )
    ),

    mean_percent = as.numeric(
      mean_percent
    ),

    sem_percent = as.numeric(
      sem_percent
    ),

    total_count = as.numeric(
      total_count
    )
  ) %>%
  filter(
    hexamer %in% hexamer_order
  )

# ============================================================
# READ STREAMED DAMAGE-SEQ HEXAMER TABLE
# ============================================================

damage_input_raw <- read_tsv(
  damage_hexamer_file,
  show_col_types = FALSE,
  progress = FALSE
)

check_required_columns(
  damage_input_raw,
  c(
    "hexamer",
    "damage_count"
  ),
  "Damage-seq streamed hexamer table"
)

if (
  !"damage_window_count" %in%
    names(damage_input_raw)
) {
  damage_input_raw <- damage_input_raw %>%
    mutate(
      damage_window_count =
        damage_count
    )
}

damage_observed <- damage_input_raw %>%
  transmute(
    hexamer = str_to_upper(
      str_trim(
        as.character(hexamer)
      )
    ),

    damage_count = as.numeric(
      damage_count
    ),

    damage_window_count = as.numeric(
      damage_window_count
    )
  ) %>%
  filter(
    hexamer %in% hexamer_order,
    is.finite(damage_count),
    damage_count >= 0
  ) %>%
  group_by(
    hexamer
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_window_count = sum(
      damage_window_count,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

# ============================================================
# COMPLETE ALL 1,024 DAMAGE-SEQ CONTEXTS
#
# Damage percentages are recalculated directly from counts.
# ============================================================

damage_counts <- hexamer_key %>%
  left_join(
    damage_observed,
    by = "hexamer"
  ) %>%
  mutate(
    damage_count = replace_na(
      damage_count,
      0
    ),

    damage_window_count = replace_na(
      damage_window_count,
      0
    ),

    total_damage_sites = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_percent = if_else(
      total_damage_sites > 0,

      100 *
        damage_count /
        total_damage_sites,

      NA_real_
    ),

    eligible_damage_denominator = (
      damage_count >=
        min_damage_count &
        is.finite(
          damage_percent
        ) &
        damage_percent > 0
    ),

    left_flank = factor(
      left_flank,
      levels = dinucleotide_order
    ),

    cpd = factor(
      cpd,
      levels = cpd_order
    ),

    right_flank = factor(
      right_flank,
      levels = dinucleotide_order
    )
  ) %>%
  arrange(
    cpd,
    left_flank,
    right_flank
  )

if (nrow(damage_counts) != 1024) {
  stop(
    "The completed Damage-seq table does not contain ",
    "1,024 contexts."
  )
}

total_damage_sites <- sum(
  damage_counts$damage_count,
  na.rm = TRUE
)

if (
  !is.finite(total_damage_sites) ||
  total_damage_sites <= 0
) {
  stop(
    "The total Damage-seq hexamer count is zero or invalid."
  )
}

message(
  "\nDamage-seq total sites: ",
  format(
    total_damage_sites,
    big.mark = ","
  )
)

message(
  "Eligible Damage-seq denominators: ",
  sum(
    damage_counts$
      eligible_damage_denominator,
    na.rm = TRUE
  ),
  " / 1024"
)

# ============================================================
# INPUT VALIDATION SUMMARY
# ============================================================

input_validation <- tibble(
  metric = c(
    "WT replicate table rows",
    "WT replicate-table hexamers",
    "WT repair samples",
    "Damage-seq total sites",
    "Damage-seq positive hexamers",
    paste0(
      "Damage-seq hexamers with count >= ",
      min_damage_count
    ),
    "Damage percentage sum"
  ),

  value = c(
    nrow(
      repair_replicates
    ),

    n_distinct(
      repair_replicates$hexamer
    ),

    n_distinct(
      repair_replicates$sample_id
    ),

    total_damage_sites,

    sum(
      damage_counts$damage_count > 0,
      na.rm = TRUE
    ),

    sum(
      damage_counts$
        eligible_damage_denominator,
      na.rm = TRUE
    ),

    sum(
      damage_counts$damage_percent,
      na.rm = TRUE
    )
  )
)

write_tsv(
  input_validation,
  file.path(
    outdir,
    "WT_damage_normalization_input_validation.tsv"
  )
)

message(
  "\nInput validation:"
)

print(
  input_validation,
  n = Inf
)

# ============================================================
# DAMAGE-SEQ CENTRAL CPD DISTRIBUTION
# ============================================================

damage_cpd_summary <- damage_counts %>%
  group_by(
    cpd
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    total_damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),

    damage_percent = if_else(
      total_damage_count > 0,

      100 *
        damage_count /
        total_damage_count,

      NA_real_
    ),

    cpd = factor(
      as.character(cpd),
      levels = cpd_order
    )
  ) %>%
  arrange(
    cpd
  )

write_tsv(
  damage_cpd_summary,
  file.path(
    outdir,
    "DamageSeq_0h_CPD_family_distribution.tsv"
  )
)

# ============================================================
# REPAIR-OVER-DAMAGE PER WT REPLICATE
# ============================================================

normalized_replicates <- repair_replicates %>%
  left_join(
    damage_counts %>%
      mutate(
        hexamer = as.character(
          hexamer
        )
      ) %>%
      select(
        hexamer,
        damage_count,
        damage_window_count,
        damage_percent,
        eligible_damage_denominator
      ),
    by = "hexamer"
  ) %>%
  mutate(
    repair_over_damage = if_else(
      eligible_damage_denominator,

      percent /
        damage_percent,

      NA_real_
    ),

    log2_repair_over_damage = case_when(
      !eligible_damage_denominator ~
        NA_real_,

      repair_over_damage > 0 ~
        log2(
          repair_over_damage
        ),

      TRUE ~
        NA_real_
    ),

    log2_repair_over_damage_pseudocount =
      if_else(
        eligible_damage_denominator,

        log2(
          (
            percent +
              ratio_pseudocount
          ) /
            (
              damage_percent +
                ratio_pseudocount
            )
        ),

        NA_real_
      )
  ) %>%
  arrange(
    time_h,
    replicate,
    cpd,
    left_flank,
    right_flank
  )

write_tsv(
  normalized_replicates,
  file.path(
    outdir,
    "WT_hexanucleotide_repair_over_damage_per_replicate.tsv"
  )
)

# ============================================================
# MEAN ± SEM ACROSS WT R1 AND R2
#
# The Damage-seq denominator is fixed.
# SEM reflects WT XR-seq replicate variation only.
# ============================================================

normalized_mean_sem <- normalized_replicates %>%
  group_by(
    timepoint,
    time_h,
    hexamer,
    left_flank,
    cpd,
    right_flank,
    damage_count,
    damage_window_count,
    damage_percent,
    eligible_damage_denominator
  ) %>%
  summarise(
    mean_repair_percent = if (
      all(
        is.na(percent)
      )
    ) {
      NA_real_
    } else {
      mean(
        percent,
        na.rm = TRUE
      )
    },

    mean_repair_over_damage = if (
      all(
        is.na(
          repair_over_damage
        )
      )
    ) {
      NA_real_
    } else {
      mean(
        repair_over_damage,
        na.rm = TRUE
      )
    },

    sd_repair_over_damage = if (
      sum(
        !is.na(
          repair_over_damage
        )
      ) > 1
    ) {
      sd(
        repair_over_damage,
        na.rm = TRUE
      )
    } else {
      0
    },

    n_reps = sum(
      !is.na(
        repair_over_damage
      )
    ),

    sem_repair_over_damage = if (
      sum(
        !is.na(
          repair_over_damage
        )
      ) > 1
    ) {
      sd(
        repair_over_damage,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            !is.na(
              repair_over_damage
            )
          )
        )
    } else {
      0
    },

    mean_log2_repair_over_damage = if (
      all(
        is.na(
          log2_repair_over_damage
        )
      )
    ) {
      NA_real_
    } else {
      mean(
        log2_repair_over_damage,
        na.rm = TRUE
      )
    },

    total_repair_count = sum(
      count,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  mutate(
    timepoint = factor(
      as.character(timepoint),
      levels = time_order
    ),

    left_flank = factor(
      as.character(left_flank),
      levels = dinucleotide_order
    ),

    cpd = factor(
      as.character(cpd),
      levels = cpd_order
    ),

    right_flank = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    )
  ) %>%
  arrange(
    time_h,
    cpd,
    left_flank,
    right_flank
  )

write_tsv(
  normalized_mean_sem,
  file.path(
    outdir,
    "WT_hexanucleotide_repair_over_damage_mean_sem.tsv"
  )
)

# ============================================================
# REPAIR CONTEXT TOTALS
# ============================================================

repair_context_totals <- normalized_replicates %>%
  group_by(
    hexamer
  ) %>%
  summarise(
    total_repair_count_all_samples = sum(
      count,
      na.rm = TRUE
    ),

    n_nonzero_repair_samples = sum(
      count > 0,
      na.rm = TRUE
    ),

    n_nonzero_repair_timepoints = n_distinct(
      time_h[
        count > 0
      ]
    ),

    .groups = "drop"
  )

# ============================================================
# DAMAGE-NORMALIZED REGRESSION
#
# Model:
#
#   repair_over_damage ~ time_h
#
# The model is fitted separately for each hexamer using:
#
#   4 timepoints × 2 WT biological replicates
# ============================================================

fit_normalized_regression <- function(data) {

  damage_count_value <- dplyr::first(
    data$damage_count
  )

  total_repair_count <- sum(
    data$count,
    na.rm = TRUE
  )

  nonzero_timepoints <- n_distinct(
    data$time_h[
      data$count > 0
    ]
  )

  eligible <- (
    is.finite(
      damage_count_value
    ) &&
      damage_count_value >=
        min_damage_count &&
      total_repair_count >=
        min_total_repair_count &&
      nonzero_timepoints >=
        min_nonzero_timepoints
  )

  regression_data <- data %>%
    filter(
      is.finite(
        repair_over_damage
      ),
      is.finite(
        time_h
      )
    )

  if (
    !eligible ||
    nrow(regression_data) < 3 ||
    n_distinct(
      regression_data$time_h
    ) < 3
  ) {

    return(
      tibble(
        eligible_for_normalized_regression = FALSE,
        normalized_intercept = NA_real_,
        normalized_slope_per_hour = NA_real_,
        normalized_slope_se = NA_real_,
        normalized_slope_t = NA_real_,
        normalized_slope_p = NA_real_,
        normalized_slope_ci_low = NA_real_,
        normalized_slope_ci_high = NA_real_,
        normalized_r_squared = NA_real_,
        normalized_adjusted_r_squared = NA_real_,
        normalized_n_observations = nrow(
          regression_data
        )
      )
    )
  }

  model <- lm(
    repair_over_damage ~ time_h,
    data = regression_data
  )

  model_summary <- summary(
    model
  )

  coefficient_table <- model_summary$coefficients

  confidence_interval <- tryCatch(
    confint(
      model,
      "time_h",
      level = 0.95
    ),

    error = function(e) {

      matrix(
        c(
          NA_real_,
          NA_real_
        ),
        nrow = 1
      )
    }
  )

  tibble(
    eligible_for_normalized_regression = TRUE,

    normalized_intercept = unname(
      coef(model)[
        "(Intercept)"
      ]
    ),

    normalized_slope_per_hour =
      coefficient_table[
        "time_h",
        "Estimate"
      ],

    normalized_slope_se =
      coefficient_table[
        "time_h",
        "Std. Error"
      ],

    normalized_slope_t =
      coefficient_table[
        "time_h",
        "t value"
      ],

    normalized_slope_p =
      coefficient_table[
        "time_h",
        "Pr(>|t|)"
      ],

    normalized_slope_ci_low =
      confidence_interval[
        1,
        1
      ],

    normalized_slope_ci_high =
      confidence_interval[
        1,
        2
      ],

    normalized_r_squared =
      model_summary$r.squared,

    normalized_adjusted_r_squared =
      model_summary$adj.r.squared,

    normalized_n_observations = nrow(
      regression_data
    )
  )
}

normalized_regression <- normalized_replicates %>%
  mutate(
    hexamer = as.character(
      hexamer
    )
  ) %>%
  group_by(
    hexamer
  ) %>%
  group_modify(
    ~ fit_normalized_regression(.x)
  ) %>%
  ungroup()

normalized_regression$normalized_fdr <- NA_real_

valid_regression_rows <- which(
  normalized_regression$
    eligible_for_normalized_regression &
    is.finite(
      normalized_regression$
        normalized_slope_p
    )
)

if (length(valid_regression_rows) > 0) {

  normalized_regression$
    normalized_fdr[
      valid_regression_rows
    ] <- p.adjust(
      normalized_regression$
        normalized_slope_p[
          valid_regression_rows
        ],
      method = "BH"
    )
}

normalized_regression <- normalized_regression %>%
  mutate(
    normalized_regression_timing = case_when(
      !eligible_for_normalized_regression ~
        "Insufficient data",

      normalized_slope_per_hour < 0 ~
        "Early-enriched",

      normalized_slope_per_hour > 0 ~
        "Late-enriched",

      TRUE ~
        "Flat"
    ),

    normalized_significance = case_when(
      is.na(normalized_fdr) ~ "",

      normalized_fdr < 0.001 ~ "***",

      normalized_fdr < 0.01 ~ "**",

      normalized_fdr < 0.05 ~ "*",

      normalized_fdr < 0.10 ~ "\u00b7",

      TRUE ~ ""
    )
  )

write_tsv(
  normalized_regression,
  file.path(
    outdir,
    "WT_hexanucleotide_damage_normalized_regression.tsv"
  )
)

# ============================================================
# DAMAGE-NORMALIZED WEIGHTED REPAIR TIME
#
# weighted mean time =
#
#   sum(time × mean repair/damage index)
#   ------------------------------------
#   sum(mean repair/damage index)
# ============================================================

normalized_timing <- normalized_mean_sem %>%
  mutate(
    hexamer = as.character(
      hexamer
    ),

    timepoint = as.character(
      timepoint
    )
  ) %>%
  group_by(
    hexamer
  ) %>%
  arrange(
    time_h,
    .by_group = TRUE
  ) %>%
  summarise(
    damage_count = dplyr::first(
      damage_count
    ),

    damage_window_count = dplyr::first(
      damage_window_count
    ),

    damage_percent = dplyr::first(
      damage_percent
    ),

    total_normalized_timecourse_signal = sum(
      mean_repair_over_damage,
      na.rm = TRUE
    ),

    damage_normalized_weighted_mean_time_h =
      if (
        sum(
          mean_repair_over_damage,
          na.rm = TRUE
        ) > 0
      ) {

        sum(
          time_h *
            mean_repair_over_damage,
          na.rm = TRUE
        ) /
          sum(
            mean_repair_over_damage,
            na.rm = TRUE
          )

      } else {

        NA_real_
      },

    peak_normalized_time_h = if (
      all(
        !is.finite(
          mean_repair_over_damage
        )
      )
    ) {

      NA_real_

    } else {

      time_h[
        which.max(
          replace(
            mean_repair_over_damage,
            !is.finite(
              mean_repair_over_damage
            ),
            -Inf
          )
        )
      ][1]
    },

    peak_normalized_timepoint = if (
      all(
        !is.finite(
          mean_repair_over_damage
        )
      )
    ) {

      NA_character_

    } else {

      timepoint[
        which.max(
          replace(
            mean_repair_over_damage,
            !is.finite(
              mean_repair_over_damage
            ),
            -Inf
          )
        )
      ][1]
    },

    normalized_mean_0_5h =
      mean_repair_over_damage[
        time_h == 0.5
      ][1],

    normalized_mean_2h =
      mean_repair_over_damage[
        time_h == 2
      ][1],

    normalized_mean_4h =
      mean_repair_over_damage[
        time_h == 4
      ][1],

    normalized_mean_8h =
      mean_repair_over_damage[
        time_h == 8
      ][1],

    normalized_initial_index_per_hour =
      mean_repair_over_damage[
        time_h == 0.5
      ][1] / 0.5,

    normalized_endpoint_change =
      mean_repair_over_damage[
        time_h == 8
      ][1] -
        mean_repair_over_damage[
          time_h == 0.5
        ][1],

    normalized_early_signal = sum(
      mean_repair_over_damage[
        time_h %in% c(
          0.5,
          2
        )
      ],
      na.rm = TRUE
    ),

    normalized_late_signal = sum(
      mean_repair_over_damage[
        time_h %in% c(
          4,
          8
        )
      ],
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  left_join(
    repair_context_totals,
    by = "hexamer"
  ) %>%
  left_join(
    normalized_regression,
    by = "hexamer"
  ) %>%
  left_join(
    hexamer_key,
    by = "hexamer"
  ) %>%
  mutate(
    eligible_for_normalized_timing = (
      damage_count >=
        min_damage_count &
        total_repair_count_all_samples >=
          min_total_repair_count &
        n_nonzero_repair_timepoints >=
          min_nonzero_timepoints &
        is.finite(
          damage_normalized_weighted_mean_time_h
        )
    )
  )

eligible_normalized_times <- normalized_timing %>%
  filter(
    eligible_for_normalized_timing
  ) %>%
  pull(
    damage_normalized_weighted_mean_time_h
  )

if (length(eligible_normalized_times) < 4) {
  stop(
    "Fewer than four contexts passed the damage-normalized ",
    "timing thresholds.\n",
    "Consider reducing min_damage_count or ",
    "min_total_repair_count."
  )
}

normalized_q25 <- quantile(
  eligible_normalized_times,
  probs = 0.25,
  na.rm = TRUE,
  names = FALSE
)

normalized_median <- median(
  eligible_normalized_times,
  na.rm = TRUE
)

normalized_q75 <- quantile(
  eligible_normalized_times,
  probs = 0.75,
  na.rm = TRUE,
  names = FALSE
)

normalized_timing <- normalized_timing %>%
  mutate(
    normalized_earliest_rank = min_rank(
      if_else(
        eligible_for_normalized_timing,
        damage_normalized_weighted_mean_time_h,
        NA_real_
      )
    ),

    normalized_latest_rank = min_rank(
      if_else(
        eligible_for_normalized_timing,
        -damage_normalized_weighted_mean_time_h,
        NA_real_
      )
    ),

    normalized_timing_class = case_when(
      !eligible_for_normalized_timing ~
        "Insufficient data",

      damage_normalized_weighted_mean_time_h <=
        normalized_q25 ~
        "Earliest repaired",

      damage_normalized_weighted_mean_time_h <=
        normalized_median ~
        "Early-intermediate",

      damage_normalized_weighted_mean_time_h <
        normalized_q75 ~
        "Late-intermediate",

      damage_normalized_weighted_mean_time_h >=
        normalized_q75 ~
        "Latest repaired",

      TRUE ~
        "Insufficient data"
    ),

    normalized_timing_class = factor(
      normalized_timing_class,
      levels = c(
        "Earliest repaired",
        "Early-intermediate",
        "Late-intermediate",
        "Latest repaired",
        "Insufficient data"
      )
    ),

    left_flank = factor(
      left_flank,
      levels = dinucleotide_order
    ),

    cpd = factor(
      cpd,
      levels = cpd_order
    ),

    right_flank = factor(
      right_flank,
      levels = dinucleotide_order
    )
  ) %>%
  arrange(
    damage_normalized_weighted_mean_time_h
  )

write_tsv(
  normalized_timing,
  file.path(
    outdir,
    "WT_hexanucleotide_damage_normalized_timing_metrics.tsv"
  )
)

# ============================================================
# EARLIEST AND LATEST TABLES
# ============================================================

normalized_earliest <- normalized_timing %>%
  filter(
    eligible_for_normalized_timing
  ) %>%
  arrange(
    damage_normalized_weighted_mean_time_h,
    desc(
      normalized_initial_index_per_hour
    )
  )

normalized_latest <- normalized_timing %>%
  filter(
    eligible_for_normalized_timing
  ) %>%
  arrange(
    desc(
      damage_normalized_weighted_mean_time_h
    ),
    normalized_initial_index_per_hour
  )

write_tsv(
  normalized_earliest %>%
    slice_head(
      n = top_n_timing
    ),
  file.path(
    outdir,
    paste0(
      "WT_top",
      top_n_timing,
      "_damage_normalized_earliest_hexamers.tsv"
    )
  )
)

write_tsv(
  normalized_latest %>%
    slice_head(
      n = top_n_timing
    ),
  file.path(
    outdir,
    paste0(
      "WT_top",
      top_n_timing,
      "_damage_normalized_latest_hexamers.tsv"
    )
  )
)

exact_earliest <- normalized_timing %>%
  filter(
    normalized_earliest_rank == 1
  )

exact_latest <- normalized_timing %>%
  filter(
    normalized_latest_rank == 1
  )

write_tsv(
  exact_earliest,
  file.path(
    outdir,
    "WT_damage_normalized_exact_earliest_hexamer.tsv"
  )
)

write_tsv(
  exact_latest,
  file.path(
    outdir,
    "WT_damage_normalized_exact_latest_hexamer.tsv"
  )
)

message(
  "\nExact damage-normalized earliest hexamer:"
)

print(
  exact_earliest %>%
    select(
      hexamer,
      cpd,
      damage_count,
      damage_percent,
      damage_normalized_weighted_mean_time_h,
      peak_normalized_timepoint,
      normalized_initial_index_per_hour,
      normalized_slope_per_hour,
      normalized_fdr
    ),
  n = Inf
)

message(
  "\nExact damage-normalized latest hexamer:"
)

print(
  exact_latest %>%
    select(
      hexamer,
      cpd,
      damage_count,
      damage_percent,
      damage_normalized_weighted_mean_time_h,
      peak_normalized_timepoint,
      normalized_initial_index_per_hour,
      normalized_slope_per_hour,
      normalized_fdr
    ),
  n = Inf
)

# ============================================================
# CENTRAL CPD REPAIR-OVER-DAMAGE
#
# Repair and damage are first aggregated to the CPD-family
# level. The family-level ratio is then calculated.
# ============================================================

repair_cpd_replicates <- repair_replicates %>%
  group_by(
    sample_id,
    timepoint,
    time_h,
    replicate,
    cpd
  ) %>%
  summarise(
    repair_count = sum(
      count,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%
  group_by(
    sample_id,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_repair_count = sum(
      repair_count,
      na.rm = TRUE
    ),

    repair_cpd_percent = case_when(
      total_repair_count > 0 ~
        100 *
          repair_count /
          total_repair_count,

      TRUE ~
        NA_real_
    )
  ) %>%
  ungroup() %>%
  left_join(
    damage_cpd_summary %>%
      mutate(
        cpd = factor(
          as.character(cpd),
          levels = cpd_order
        )
      ) %>%
      select(
        cpd,
        damage_cpd_count =
          damage_count,
        damage_cpd_percent =
          damage_percent
      ),
    by = "cpd"
  ) %>%
  mutate(
    cpd_repair_over_damage = case_when(
      is.finite(
        damage_cpd_percent
      ) &
        damage_cpd_percent > 0 &
        is.finite(
          repair_cpd_percent
        ) ~
        repair_cpd_percent /
          damage_cpd_percent,

      TRUE ~
        NA_real_
    ),

    log2_cpd_repair_over_damage = case_when(
      is.finite(
        cpd_repair_over_damage
      ) &
        cpd_repair_over_damage > 0 ~
        log2(
          cpd_repair_over_damage
        ),

      TRUE ~
        NA_real_
    )
  ) %>%
  arrange(
    time_h,
    replicate,
    cpd
  )

cpd_normalized_mean_sem <- repair_cpd_replicates %>%
  group_by(
    timepoint,
    time_h,
    cpd,
    damage_cpd_count,
    damage_cpd_percent
  ) %>%
  summarise(
    mean_repair_cpd_percent = if (
      all(
        is.na(
          repair_cpd_percent
        )
      )
    ) {
      NA_real_
    } else {
      mean(
        repair_cpd_percent,
        na.rm = TRUE
      )
    },

    sem_repair_cpd_percent = if (
      sum(
        !is.na(
          repair_cpd_percent
        )
      ) > 1
    ) {
      sd(
        repair_cpd_percent,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            !is.na(
              repair_cpd_percent
            )
          )
        )
    } else {
      0
    },

    mean_cpd_repair_over_damage = if (
      all(
        is.na(
          cpd_repair_over_damage
        )
      )
    ) {
      NA_real_
    } else {
      mean(
        cpd_repair_over_damage,
        na.rm = TRUE
      )
    },

    sd_cpd_repair_over_damage = if (
      sum(
        !is.na(
          cpd_repair_over_damage
        )
      ) > 1
    ) {
      sd(
        cpd_repair_over_damage,
        na.rm = TRUE
      )
    } else {
      0
    },

    n_reps = sum(
      !is.na(
        cpd_repair_over_damage
      )
    ),

    sem_cpd_repair_over_damage = if (
      sum(
        !is.na(
          cpd_repair_over_damage
        )
      ) > 1
    ) {
      sd(
        cpd_repair_over_damage,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            !is.na(
              cpd_repair_over_damage
            )
          )
        )
    } else {
      0
    },

    mean_log2_cpd_repair_over_damage = if (
      all(
        is.na(
          log2_cpd_repair_over_damage
        )
      )
    ) {
      NA_real_
    } else {
      mean(
        log2_cpd_repair_over_damage,
        na.rm = TRUE
      )
    },

    .groups = "drop"
  ) %>%
  mutate(
    timepoint = factor(
      as.character(timepoint),
      levels = time_order
    ),

    cpd = factor(
      as.character(cpd),
      levels = cpd_order
    )
  ) %>%
  arrange(
    time_h,
    cpd
  )

write_tsv(
  repair_cpd_replicates,
  file.path(
    outdir,
    "WT_CPD_family_repair_over_damage_per_replicate.tsv"
  )
)

write_tsv(
  cpd_normalized_mean_sem,
  file.path(
    outdir,
    "WT_CPD_family_repair_over_damage_mean_sem.tsv"
  )
)

message(
  "\nDamage-seq central CPD distribution:"
)

print(
  damage_cpd_summary,
  n = Inf
)

message(
  "\nWT central CPD repair-over-damage:"
)

print(
  cpd_normalized_mean_sem,
  n = Inf
)

# ============================================================
# 1. DAMAGE-SEQ HEXAMER DISTRIBUTION HEATMAP
# ============================================================

damage_heatmap_df <- damage_counts %>%
  mutate(
    left_flank_plot = factor(
      as.character(left_flank),
      levels = rev(
        dinucleotide_order
      )
    ),

    right_flank_plot = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    )
  )

p_damage_heatmap <- ggplot(
  damage_heatmap_df,
  aes(
    x = right_flank_plot,
    y = left_flank_plot,
    fill = damage_percent
  )
) +
  geom_tile() +
  facet_wrap(
    ~cpd,
    ncol = 2
  ) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",
    name = "Damage-seq\nsites (%)"
  ) +
  labs(
    x = "Two 3′ flanking bases",
    y = "Two 5′ flanking bases",
    title = "Damage-seq 0-hour hexanucleotide distribution"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    )
  )

save_plot(
  p_damage_heatmap,
  "DamageSeq_0h_hexanucleotide_distribution_heatmap",
  width = 8,
  height = 7
)

# ============================================================
# 2. ALL-TIMEPOINT REPAIR/DAMAGE HEATMAP
#
# Values above the 99th percentile are capped for display only.
# The output tables retain the uncapped values.
# ============================================================

normalized_heatmap_df <- normalized_mean_sem %>%
  mutate(
    left_flank_plot = factor(
      as.character(left_flank),
      levels = rev(
        dinucleotide_order
      )
    ),

    right_flank_plot = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    )
  )

normalized_heatmap_max <- safe_finite_quantile(
  normalized_heatmap_df$
    mean_repair_over_damage,
  probability = 0.99,
  fallback = 1
)

normalized_heatmap_df <- normalized_heatmap_df %>%
  mutate(
    repair_over_damage_plot = if_else(
      is.finite(
        mean_repair_over_damage
      ),

      pmin(
        mean_repair_over_damage,
        normalized_heatmap_max
      ),

      NA_real_
    )
  )

p_normalized_heatmaps <- ggplot(
  normalized_heatmap_df,
  aes(
    x = right_flank_plot,
    y = left_flank_plot,
    fill = repair_over_damage_plot
  )
) +
  geom_tile() +
  facet_grid(
    rows = vars(cpd),
    cols = vars(timepoint)
  ) +
  scale_fill_viridis_c(
    option = "magma",
    trans = "sqrt",

    limits = c(
      0,
      normalized_heatmap_max
    ),

    na.value = "grey90",

    name = paste0(
      "Repair / damage\n",
      "(99th-percentile cap)"
    )
  ) +
  labs(
    x = "Two 3′ flanking bases",
    y = "Two 5′ flanking bases",

    title = paste0(
      "WT UV-CPD hexanucleotide repair normalized ",
      "to Damage-seq"
    ),

    subtitle = paste0(
      "WT repair percentage / Damage-seq percentage"
    )
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5.5,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 5.5,
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    )
  )

save_plot(
  p_normalized_heatmaps,
  paste0(
    "WT_hexanucleotide_repair_over_damage_",
    "all_timepoints_heatmap"
  ),
  width = 13,
  height = 12
)

# ============================================================
# 3. TOP NORMALIZED HEXAMERS PER TIMEPOINT
# ============================================================

top_normalized_df <- normalized_mean_sem %>%
  filter(
    eligible_damage_denominator,
    is.finite(
      mean_repair_over_damage
    )
  ) %>%
  group_by(
    timepoint
  ) %>%
  slice_max(
    order_by = mean_repair_over_damage,
    n = top_n_per_timepoint,
    with_ties = FALSE
  ) %>%
  ungroup() %>%
  mutate(
    plot_id = paste0(
      hexamer,
      "___",
      as.character(
        timepoint
      )
    )
  ) %>%
  arrange(
    timepoint,
    mean_repair_over_damage
  ) %>%
  mutate(
    plot_id = factor(
      plot_id,
      levels = unique(
        plot_id
      )
    )
  )

write_tsv(
  top_normalized_df,
  file.path(
    outdir,
    paste0(
      "WT_top",
      top_n_per_timepoint,
      "_repair_over_damage_hexamers_per_timepoint.tsv"
    )
  )
)

p_top_normalized <- ggplot(
  top_normalized_df,
  aes(
    x = mean_repair_over_damage,
    y = plot_id
  )
) +
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "grey35",
    linewidth = 0.35
  ) +
  geom_col(
    width = 0.72,
    fill = wt_color,
    color = "black",
    linewidth = 0.18
  ) +
  geom_segment(
    aes(
      x = pmax(
        mean_repair_over_damage -
          sem_repair_over_damage,
        0
      ),

      xend =
        mean_repair_over_damage +
        sem_repair_over_damage,

      y = plot_id,
      yend = plot_id
    ),
    linewidth = 0.3,
    color = "black"
  ) +
  facet_wrap(
    ~timepoint,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = remove_facet_suffix
  ) +
  scale_x_continuous(
    expand = expansion(
      mult = c(
        0,
        0.08
      )
    )
  ) +
  labs(
    x = "WT repair / Damage-seq index",
    y = NULL,

    title = paste0(
      "Top ",
      top_n_per_timepoint,
      " damage-normalized WT hexamers"
    ),

    subtitle = paste0(
      "Dashed line at 1: repair frequency equals ",
      "starting damage frequency"
    )
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      size = 7,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 5.7,
      family = "mono",
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    axis.line.y = element_blank(),

    axis.ticks.y = element_blank(),

    panel.spacing.x = grid::unit(
      0.8,
      "lines"
    ),

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    )
  )

save_plot(
  p_top_normalized,
  paste0(
    "WT_top",
    top_n_per_timepoint,
    "_repair_over_damage_barplots"
  ),
  width = 13,
  height = 7
)

# ============================================================
# 4. MOST DYNAMIC NORMALIZED TRAJECTORIES
# ============================================================

dynamic_normalized <- normalized_timing %>%
  filter(
    eligible_for_normalized_regression,
    is.finite(
      normalized_slope_per_hour
    )
  ) %>%
  slice_max(
    order_by = abs(
      normalized_slope_per_hour
    ),
    n = top_n_dynamic,
    with_ties = FALSE
  ) %>%
  arrange(
    desc(
      abs(
        normalized_slope_per_hour
      )
    )
  )

dynamic_hexamers <- dynamic_normalized$hexamer

dynamic_normalized_mean <- normalized_mean_sem %>%
  filter(
    hexamer %in%
      dynamic_hexamers
  ) %>%
  mutate(
    hexamer = factor(
      hexamer,
      levels = dynamic_hexamers
    )
  )

dynamic_normalized_replicates <- normalized_replicates %>%
  filter(
    hexamer %in%
      dynamic_hexamers
  ) %>%
  mutate(
    hexamer = factor(
      hexamer,
      levels = dynamic_hexamers
    )
  )

dynamic_labels <- dynamic_normalized %>%
  transmute(
    hexamer = factor(
      hexamer,
      levels = dynamic_hexamers
    ),

    x = 0.65,
    y = Inf,

    label = paste0(
      "slope = ",
      sprintf(
        "%.4f",
        normalized_slope_per_hour
      ),
      "/h\nFDR = ",
      if_else(
        is.finite(
          normalized_fdr
        ),
        formatC(
          normalized_fdr,
          format = "g",
          digits = 2
        ),
        "NA"
      )
    )
  )

p_dynamic_normalized <- ggplot() +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.3,
    color = "grey55"
  ) +
  geom_line(
    data = dynamic_normalized_mean,
    aes(
      x = time_h,
      y = mean_repair_over_damage,
      group = hexamer
    ),
    color = wt_color,
    linewidth = 0.65
  ) +
  geom_errorbar(
    data = dynamic_normalized_mean,
    aes(
      x = time_h,

      ymin = pmax(
        mean_repair_over_damage -
          sem_repair_over_damage,
        0
      ),

      ymax =
        mean_repair_over_damage +
        sem_repair_over_damage
    ),
    width = 0.12,
    linewidth = 0.28,
    color = "black"
  ) +
  geom_point(
    data = dynamic_normalized_replicates,
    aes(
      x = time_h,
      y = repair_over_damage
    ),
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 1,
    size = 1.25,
    stroke = 0.3,
    color = "grey25"
  ) +
  geom_point(
    data = dynamic_normalized_mean,
    aes(
      x = time_h,
      y = mean_repair_over_damage
    ),
    shape = 21,
    size = 2,
    stroke = 0.3,
    fill = wt_color,
    color = "black"
  ) +
  geom_text(
    data = dynamic_labels,
    aes(
      x = x,
      y = y,
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 0,
    vjust = 1.05,
    size = 2.2,
    lineheight = 0.9
  ) +
  facet_wrap(
    ~hexamer,
    ncol = 6,
    scales = "free_y"
  ) +
  scale_x_continuous(
    breaks = c(
      0.5,
      2,
      4,
      8
    ),
    limits = c(
      0.35,
      8.15
    )
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.20
      )
    )
  ) +
  labs(
    x = "Time after UV (hours)",
    y = "WT repair / Damage-seq index",

    title = paste0(
      "Top ",
      top_n_dynamic,
      " damage-normalized dynamic hexamers"
    ),

    subtitle = "Selected by absolute damage-normalized slope"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      size = 6,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      family = "mono",
      size = 7
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    )
  )

save_plot(
  p_dynamic_normalized,
  paste0(
    "WT_top",
    top_n_dynamic,
    "_damage_normalized_dynamic_trajectories"
  ),
  width = 12,
  height = 9
)

# ============================================================
# 5. NORMALIZED 0.5H VERSUS 8H SCATTER
# ============================================================

normalized_scatter_df <- normalized_timing %>%
  filter(
    eligible_for_normalized_timing,
    is.finite(
      normalized_mean_0_5h
    ),
    is.finite(
      normalized_mean_8h
    )
  )

scatter_limit <- safe_finite_quantile(
  c(
    normalized_scatter_df$
      normalized_mean_0_5h,

    normalized_scatter_df$
      normalized_mean_8h
  ),
  probability = 0.995,
  fallback = 1
)

normalized_scatter_labels <- normalized_scatter_df %>%
  slice_max(
    order_by = abs(
      normalized_endpoint_change
    ),
    n = top_n_labels,
    with_ties = FALSE
  )

p_normalized_scatter <- ggplot(
  normalized_scatter_df,
  aes(
    x = normalized_mean_0_5h,
    y = normalized_mean_8h,
    color = cpd
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey35"
  ) +
  geom_point(
    alpha = 0.68,
    size = 1.6
  ) +
  geom_text(
    data = normalized_scatter_labels,
    aes(
      label = hexamer
    ),
    nudge_y = scatter_limit * 0.015,
    size = 2.4,
    family = "mono",
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = cpd_colors,
    drop = FALSE
  ) +
  scale_x_sqrt(
    limits = c(
      0,
      scatter_limit * 1.08
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  scale_y_sqrt(
    limits = c(
      0,
      scatter_limit * 1.08
    ),
    expand = expansion(
      mult = c(
        0,
        0
      )
    )
  ) +
  coord_equal() +
  labs(
    x = "Damage-normalized WT 0.5h",
    y = "Damage-normalized WT 8h",

    title = paste0(
      "Damage-normalized WT hexamers: ",
      "0.5h versus 8h"
    ),

    subtitle = paste0(
      "Below diagonal = early-enriched; ",
      "above diagonal = late-enriched"
    ),

    color = "Central CPD"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    ),

    legend.position = "bottom",

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.3
    )
  )

save_plot(
  p_normalized_scatter,
  "WT_damage_normalized_0_5h_vs_8h_scatter",
  width = 6,
  height = 5.5
)

# ============================================================
# 6. NORMALIZED SLOPE VERSUS FDR
# ============================================================

normalized_regression_scatter <- normalized_timing %>%
  filter(
    eligible_for_normalized_regression,
    is.finite(
      normalized_slope_per_hour
    ),
    is.finite(
      normalized_fdr
    )
  ) %>%
  mutate(
    negative_log10_fdr = -log10(
      pmax(
        normalized_fdr,
        1e-300
      )
    )
  )

normalized_regression_labels <-
  normalized_regression_scatter %>%
  arrange(
    normalized_fdr,
    desc(
      abs(
        normalized_slope_per_hour
      )
    )
  ) %>%
  slice_head(
    n = top_n_labels
  )

p_normalized_regression <- ggplot(
  normalized_regression_scatter,
  aes(
    x = normalized_slope_per_hour,
    y = negative_log10_fdr,
    color = normalized_regression_timing
  )
) +
  geom_vline(
    xintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey35"
  ) +
  geom_hline(
    yintercept = -log10(0.05),
    linetype = "dotted",
    linewidth = 0.4,
    color = "grey35"
  ) +
  geom_point(
    alpha = 0.72,
    size = 1.6
  ) +
  geom_text(
    data = normalized_regression_labels,
    aes(
      label = hexamer
    ),
    nudge_y = 0.12,
    size = 2.4,
    family = "mono",
    check_overlap = TRUE,
    show.legend = FALSE
  ) +
  scale_color_manual(
    values = regression_colors,
    drop = FALSE
  ) +
  labs(
    x = "Damage-normalized temporal slope per hour",
    y = expression(-log[10]("BH FDR")),

    title = paste0(
      "Damage-normalized WT hexanucleotide regression"
    ),

    subtitle = "Regression: repair-over-damage index ~ time",

    color = "Timing"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    ),

    legend.position = "bottom",

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.3
    )
  )

save_plot(
  p_normalized_regression,
  "WT_damage_normalized_hexanucleotide_slope_vs_FDR",
  width = 6.2,
  height = 5
)

# ============================================================
# 7. CENTRAL CPD REPAIR/DAMAGE TIME COURSE
# ============================================================

p_cpd_normalized <- ggplot(
  cpd_normalized_mean_sem,
  aes(
    x = time_h,
    y = mean_cpd_repair_over_damage,
    color = cpd,
    shape = cpd,
    group = cpd
  )
) +
  geom_hline(
    yintercept = 1,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.75
  ) +
  geom_errorbar(
    aes(
      ymin = pmax(
        mean_cpd_repair_over_damage -
          sem_cpd_repair_over_damage,
        0
      ),

      ymax =
        mean_cpd_repair_over_damage +
        sem_cpd_repair_over_damage
    ),
    width = 0.12,
    linewidth = 0.3
  ) +
  geom_point(
    size = 2.5
  ) +
  geom_point(
    data = repair_cpd_replicates,
    aes(
      x = time_h,
      y = cpd_repair_over_damage,
      color = cpd
    ),
    inherit.aes = FALSE,
    position = position_jitter(
      width = 0.08,
      height = 0
    ),
    shape = 1,
    size = 1.5,
    stroke = 0.35
  ) +
  scale_color_manual(
    values = cpd_colors,
    breaks = cpd_order,
    drop = FALSE
  ) +
  scale_shape_manual(
    values = c(
      "CC" = 16,
      "CT" = 17,
      "TC" = 15,
      "TT" = 18
    ),
    breaks = cpd_order,
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
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.07
      )
    )
  ) +
  labs(
    x = "Time after UV (hours)",
    y = "CPD-family repair / damage index",

    title = paste0(
      "WT central CPD repair normalized to Damage-seq"
    ),

    subtitle = paste0(
      "Dashed line at 1: repair composition matches ",
      "starting damage composition"
    ),

    color = "Central CPD",
    shape = "Central CPD"
  ) +
  theme_classic(
    base_size = 9
  ) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    plot.title = element_text(
      face = "bold",
      size = 11,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    ),

    legend.position = "bottom",

    panel.border = element_rect(
      fill = NA,
      color = "black",
      linewidth = 0.3
    )
  )

save_plot(
  p_cpd_normalized,
  "WT_CPD_family_repair_over_damage_timecourse",
  width = 5.8,
  height = 4.8
)

# ============================================================
# 8. EARLIEST AND LATEST NORMALIZED LOLLIPOP
# ============================================================

earliest_lollipop <- normalized_earliest %>%
  slice_head(
    n = top_n_timing
  ) %>%
  mutate(
    timing_group = "Earliest",

    plot_id = paste0(
      hexamer,
      "___Earliest"
    )
  )

latest_lollipop <- normalized_latest %>%
  slice_head(
    n = top_n_timing
  ) %>%
  mutate(
    timing_group = "Latest",

    plot_id = paste0(
      hexamer,
      "___Latest"
    )
  )

lollipop_df <- bind_rows(
  earliest_lollipop,
  latest_lollipop
) %>%
  arrange(
    timing_group,
    damage_normalized_weighted_mean_time_h
  ) %>%
  mutate(
    plot_id = factor(
      plot_id,
      levels = unique(
        plot_id
      )
    ),

    timing_group = factor(
      timing_group,
      levels = c(
        "Earliest",
        "Latest"
      )
    )
  )

p_normalized_lollipop <- ggplot(
  lollipop_df,
  aes(
    x = damage_normalized_weighted_mean_time_h,
    y = plot_id
  )
) +
  geom_segment(
    aes(
      x = 0.5,

      xend =
        damage_normalized_weighted_mean_time_h,

      y = plot_id,
      yend = plot_id,
      color = timing_group
    ),
    linewidth = 0.55
  ) +
  geom_point(
    aes(
      color = timing_group
    ),
    size = 2.3
  ) +
  facet_wrap(
    ~timing_group,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_y_discrete(
    labels = remove_facet_suffix
  ) +
  scale_color_manual(
    values = c(
      "Earliest" = "#2166AC",
      "Latest"   = "#B2182B"
    )
  ) +
  scale_x_continuous(
    limits = c(
      0.5,
      8
    ),
    breaks = c(
      0.5,
      2,
      4,
      6,
      8
    )
  ) +
  labs(
    x = "Damage-normalized weighted mean time (hours)",
    y = NULL,

    title = paste0(
      "Top ",
      top_n_timing,
      " damage-normalized earliest and latest hexamers"
    )
  ) +
  guides(
    color = "none"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      size = 7,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      family = "mono",
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    axis.line.y = element_blank(),

    axis.ticks.y = element_blank(),

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    )
  )

save_plot(
  p_normalized_lollipop,
  paste0(
    "WT_top",
    top_n_timing,
    "_damage_normalized_earliest_latest_lollipop"
  ),
  width = 8,
  height = 7
)

# ============================================================
# 9. REPAIR PERCENTAGE VERSUS DAMAGE PERCENTAGE
# ============================================================

repair_damage_scatter <- normalized_mean_sem %>%
  filter(
    eligible_damage_denominator,
    is.finite(
      mean_repair_percent
    ),
    is.finite(
      damage_percent
    )
  )

repair_damage_limit <- safe_finite_max(
  c(
    repair_damage_scatter$
      mean_repair_percent,

    repair_damage_scatter$
      damage_percent
  ),
  fallback = 1
)

p_repair_damage_scatter <- ggplot(
  repair_damage_scatter,
  aes(
    x = damage_percent,
    y = mean_repair_percent,
    color = cpd
  )
) +
  geom_abline(
    slope = 1,
    intercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey40"
  ) +
  geom_point(
    alpha = 0.65,
    size = 1
  ) +
  facet_wrap(
    ~timepoint,
    nrow = 1
  ) +
  scale_color_manual(
    values = cpd_colors,
    drop = FALSE
  ) +
  scale_x_sqrt(
    limits = c(
      0,
      repair_damage_limit * 1.05
    )
  ) +
  scale_y_sqrt(
    limits = c(
      0,
      repair_damage_limit * 1.05
    )
  ) +
  coord_equal() +
  labs(
    x = "Damage-seq hexamer percentage",
    y = "WT repair hexamer percentage",

    title = paste0(
      "WT repair distribution versus starting ",
      "Damage-seq distribution"
    ),

    color = "Central CPD"
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    axis.text = element_text(
      color = "black"
    ),

    legend.position = "bottom",

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    )
  )

save_plot(
  p_repair_damage_scatter,
  "WT_repair_percentage_vs_DamageSeq_percentage",
  width = 11,
  height = 4
)

# ============================================================
# 10. DAMAGE-NORMALIZED TIMING HEATMAP
# ============================================================

normalized_timing_heatmap <- normalized_timing %>%
  mutate(
    left_flank_plot = factor(
      as.character(left_flank),
      levels = rev(
        dinucleotide_order
      )
    ),

    right_flank_plot = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    ),

    weighted_time_plot = if_else(
      eligible_for_normalized_timing,
      damage_normalized_weighted_mean_time_h,
      NA_real_
    )
  )

p_normalized_timing_heatmap <- ggplot(
  normalized_timing_heatmap,
  aes(
    x = right_flank_plot,
    y = left_flank_plot,
    fill = weighted_time_plot
  )
) +
  geom_tile() +
  facet_wrap(
    ~cpd,
    ncol = 2
  ) +
  scale_fill_viridis_c(
    option = "plasma",

    limits = c(
      0.5,
      8
    ),

    breaks = c(
      0.5,
      2,
      4,
      6,
      8
    ),

    na.value = "grey90",

    name = "Weighted mean\ntime (hours)"
  ) +
  labs(
    x = "Two 3′ flanking bases",
    y = "Two 5′ flanking bases",

    title = paste0(
      "Damage-normalized WT hexanucleotide repair timing"
    ),

    subtitle = paste0(
      "Lower weighted time = earlier; ",
      "higher weighted time = later"
    )
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    )
  )

save_plot(
  p_normalized_timing_heatmap,
  "WT_damage_normalized_weighted_timing_heatmap",
  width = 8,
  height = 7
)

# ============================================================
# 11. DAMAGE-NORMALIZED SLOPE HEATMAP
# ============================================================

eligible_slopes <- normalized_timing %>%
  filter(
    eligible_for_normalized_regression,
    is.finite(
      normalized_slope_per_hour
    )
  ) %>%
  pull(
    normalized_slope_per_hour
  )

slope_limit <- safe_finite_max(
  abs(
    eligible_slopes
  ),
  fallback = 1
)

normalized_slope_heatmap <- normalized_timing %>%
  mutate(
    left_flank_plot = factor(
      as.character(left_flank),
      levels = rev(
        dinucleotide_order
      )
    ),

    right_flank_plot = factor(
      as.character(right_flank),
      levels = dinucleotide_order
    ),

    slope_plot = if_else(
      eligible_for_normalized_regression,
      normalized_slope_per_hour,
      NA_real_
    )
  )

p_normalized_slope_heatmap <- ggplot(
  normalized_slope_heatmap,
  aes(
    x = right_flank_plot,
    y = left_flank_plot,
    fill = slope_plot
  )
) +
  geom_tile() +
  facet_wrap(
    ~cpd,
    ncol = 2
  ) +
  scale_fill_gradient2(
    low = "#2166AC",
    mid = "white",
    high = "#B2182B",
    midpoint = 0,

    limits = c(
      -slope_limit,
      slope_limit
    ),

    na.value = "grey90",

    name = "Normalized\nslope/hour"
  ) +
  labs(
    x = "Two 3′ flanking bases",
    y = "Two 5′ flanking bases",

    title = paste0(
      "Damage-normalized WT hexanucleotide temporal slopes"
    ),

    subtitle = paste0(
      "Negative = early-enriched; ",
      "positive = late-enriched"
    )
  ) +
  theme_classic(
    base_size = 8
  ) +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 6,
      color = "black"
    ),

    axis.text.y = element_text(
      size = 6,
      color = "black"
    ),

    strip.text = element_text(
      face = "bold",
      size = 8
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey40",
      linewidth = 0.25
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey40",
      linewidth = 0.25
    ),

    plot.title = element_text(
      face = "bold",
      size = 10,
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      size = 8,
      hjust = 0.5
    )
  )

save_plot(
  p_normalized_slope_heatmap,
  "WT_damage_normalized_temporal_slope_heatmap",
  width = 8,
  height = 7
)

# ============================================================
# 12. COMBINED DAMAGE-NORMALIZED SUMMARY
# ============================================================

combined_normalized_summary <- (
  p_normalized_scatter |
    p_normalized_regression |
    p_cpd_normalized
) /
  (
    p_normalized_lollipop |
      p_repair_damage_scatter
  ) +
  plot_layout(
    heights = c(
      1,
      1.25
    )
  ) +
  plot_annotation(
    tag_levels = "A"
  )

save_plot(
  combined_normalized_summary,
  "WT_hexanucleotide_damage_normalized_summary",
  width = 15,
  height = 11
)

# ============================================================
# ANALYSIS SUMMARY
# ============================================================

analysis_summary <- tibble(
  metric = c(
    "Total possible hexamers",

    "Damage-seq hexamers with positive count",

    paste0(
      "Damage-seq hexamers with count >= ",
      min_damage_count
    ),

    "Hexamers eligible for normalized timing",

    "Hexamers eligible for normalized regression",

    "Normalized regressions FDR < 0.05",

    "Normalized regressions FDR < 0.10"
  ),

  value = c(
    1024,

    sum(
      damage_counts$damage_count > 0,
      na.rm = TRUE
    ),

    sum(
      damage_counts$
        eligible_damage_denominator,
      na.rm = TRUE
    ),

    sum(
      normalized_timing$
        eligible_for_normalized_timing,
      na.rm = TRUE
    ),

    sum(
      normalized_timing$
        eligible_for_normalized_regression,
      na.rm = TRUE
    ),

    sum(
      normalized_timing$
        normalized_fdr < 0.05,
      na.rm = TRUE
    ),

    sum(
      normalized_timing$
        normalized_fdr < 0.10,
      na.rm = TRUE
    )
  )
)

write_tsv(
  analysis_summary,
  file.path(
    outdir,
    "WT_damage_normalized_hexanucleotide_analysis_summary.tsv"
  )
)

message(
  "\nDamage-normalized analysis summary:"
)

print(
  analysis_summary,
  n = Inf
)

# ============================================================
# OUTPUT INDEX
# ============================================================

output_index <- tibble(
  output_type = c(
    "Table",
    "Table",
    "Table",
    "Table",
    "Table",
    "Table",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure",
    "Figure"
  ),

  filename = c(
    "WT_hexanucleotide_repair_over_damage_per_replicate.tsv",

    "WT_hexanucleotide_repair_over_damage_mean_sem.tsv",

    "WT_hexanucleotide_damage_normalized_regression.tsv",

    "WT_hexanucleotide_damage_normalized_timing_metrics.tsv",

    "WT_damage_normalized_exact_earliest_hexamer.tsv",

    "WT_damage_normalized_exact_latest_hexamer.tsv",

    "DamageSeq_0h_hexanucleotide_distribution_heatmap.pdf",

    paste0(
      "WT_hexanucleotide_repair_over_damage_",
      "all_timepoints_heatmap.pdf"
    ),

    paste0(
      "WT_top",
      top_n_per_timepoint,
      "_repair_over_damage_barplots.pdf"
    ),

    paste0(
      "WT_top",
      top_n_dynamic,
      "_damage_normalized_dynamic_trajectories.pdf"
    ),

    "WT_damage_normalized_0_5h_vs_8h_scatter.pdf",

    "WT_damage_normalized_hexanucleotide_slope_vs_FDR.pdf",

    "WT_CPD_family_repair_over_damage_timecourse.pdf",

    paste0(
      "WT_top",
      top_n_timing,
      "_damage_normalized_earliest_latest_lollipop.pdf"
    ),

    "WT_damage_normalized_weighted_timing_heatmap.pdf",

    "WT_damage_normalized_temporal_slope_heatmap.pdf"
  )
)

write_tsv(
  output_index,
  file.path(
    outdir,
    "WT_damage_normalized_output_index.tsv"
  )
)

# ============================================================
# FINAL MESSAGES
# ============================================================

message(
  "\nDamage-seq normalization completed successfully."
)

message(
  "\nDamage-seq table used:"
)

message(
  "  ",
  damage_hexamer_file
)

message(
  "\nDamage-seq total sites:"
)

message(
  "  ",
  format(
    total_damage_sites,
    big.mark = ","
  )
)

message(
  "\nOutput directory:"
)

message(
  "  ",
  normalizePath(
    outdir
  )
)

message(
  "\nPrimary normalization:"
)

message(
  "  repair_over_damage = ",
  "WT repair percentage / Damage-seq percentage"
)

message(
  "\nPrimary timing measurement:"
)

message(
  "  damage_normalized_weighted_mean_time_h"
)

message(
  "\nPrimary temporal regression:"
)

message(
  "  repair_over_damage ~ time_h"
)

message(
  "\nImportant:"
)

message(
  "  Damage-seq provides one fixed 0-hour denominator."
)

message(
  "  SEM reflects WT XR-seq R1/R2 variation only."
)

message(
  "  repair_over_damage is a relative repair index, not the ",
  "literal fraction of genomic lesions removed."
)