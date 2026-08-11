#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(grid)
  library(splines)
})

# ============================================================
# WT UV-CPD TETRANUCLEOTIDE REPAIR-TIMING ANALYSIS
#
# PRIMARY ANALYSIS
#
#   CLR ~ ns(time_hr, df = 2) + replicate
#
# For each of 64 tetranucleotide contexts.
#
# Overall time test:
#
#   Null:
#     CLR ~ replicate
#
#   Time model:
#     CLR ~ ns(time_hr, df = 2) + replicate
#
#   Nested F-test:
#     null versus spline model
#
#   BH correction:
#     across all 64 tetranucleotides
#
# Classification:
#
#   Early repair:
#     BH-FDR <= cutoff
#     AND
#     fitted CLR at 8 h < fitted CLR at 0.5 h
#
#   Late repair:
#     BH-FDR <= cutoff
#     AND
#     fitted CLR at 8 h > fitted CLR at 0.5 h
#
#   No significant trend:
#     BH-FDR > cutoff
#
#
# SECONDARY ANALYSES
#
# 1. Nonlinearity:
#
#      CLR ~ time_hr + replicate
#
#      versus
#
#      CLR ~ ns(time_hr, df = 2) + replicate
#
# 2. Piecewise:
#
#      0.5 -> 4 h slope
#      4   -> 8 h slope
#
# 3. Early-phase sensitivity:
#
#      CLR ~ time_hr + replicate
#      using 0.5, 2, 4 h
#
# 4. All-time linear sensitivity:
#
#      CLR ~ time_hr + replicate
#
#
# COMPOSITIONAL ANALYSIS
#
# CLR is calculated across ALL 64 tetranucleotide contexts
# together within each repair sample.
#
# TT, CT, TC and CC are NOT transformed separately.
#
# Negative CLR change:
#   relative depletion over repair time
#
# Positive CLR change:
#   relative enrichment over repair time
#
#
# FIGURES
#
# 1. Combined TT | CT | TC | CC figure
#
# 2. Individual TT, CT, TC, CC figures
#
#    EXACT INDIVIDUAL DIMENSIONS:
#
#      width  = 60 pt
#      height = 200 pt
#
#    Individual figures use thinner:
#
#      bars
#      CI lines
#      zero line
#      text
#      margins
#
# ============================================================


# ============================================================
# INPUT FILE
# ============================================================

input_file <- paste0(
  "/work/users/a/r/arslank/",
  "UV_tetranucleotide_distribution/",
  "repair_damage_comparison/",
  "plots_updated/",
  "WT_R1_R2_all_64_sequences_percent_clear/",
  "WT_R1_R2_repair_percentages_per_replicate.tsv"
)


# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- file.path(
  dirname(input_file),
  "WT_tetranucleotide_spline_CLR_BH_horizontal_bars_95CI"
)

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)


# Individual manuscript-panel figures go here.
individual_plot_dir <- file.path(
  outdir,
  "individual_CPD_60x200pt"
)

dir.create(
  individual_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


if (!file.exists(input_file)) {

  stop(
    "\nInput file does not exist:\n  ",
    input_file,
    "\n"
  )
}


# ============================================================
# ANALYSIS SETTINGS
# ============================================================

time_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

time_values <- c(
  "0.5h" = 0.5,
  "2h"   = 2,
  "4h"   = 4,
  "8h"   = 8
)

replicate_order <- c(
  "R1",
  "R2"
)

cpd_order <- c(
  "TT",
  "CT",
  "TC",
  "CC"
)

dna_order <- c(
  "A",
  "C",
  "G",
  "T"
)


# ============================================================
# PRIMARY SIGNIFICANCE THRESHOLD
# ============================================================

FDR_CUTOFF <- 0.01


# ============================================================
# NATURAL SPLINE
# ============================================================

SPLINE_DF <- 2


# ============================================================
# PIECEWISE MODEL KNOT
# ============================================================

PIECEWISE_KNOT_H <- 4


# ============================================================
# DESCRIPTIVE PLATEAU THRESHOLD
# ============================================================

PLATEAU_SLOPE_RATIO_CUTOFF <- 0.50


# ============================================================
# CLR ZERO REPLACEMENT
# ============================================================

ZERO_REPLACEMENT_FRACTION <- 0.50


# ============================================================
# CLASSIFICATION
# ============================================================

classification_order <- c(
  "Early repair",
  "No significant trend",
  "Late repair"
)

classification_colors <- c(
  "Early repair" = "green3",
  "No significant trend" = "gray70",
  "Late repair" = "red3"
)

shape_order <- c(
  "Approximately linear",
  "Evidence of nonlinearity",
  "No significant time effect"
)


# ============================================================
# FIGURE DIMENSIONS
# ============================================================

# Combined figure.
COMBINED_EFFECT_WIDTH <- 8.5
COMBINED_EFFECT_HEIGHT <- 4.2


# Individual figures.
#
# 72 pt = 1 inch
#
# Required:
#
#   width  = 60 pt
#   height = 200 pt
#
INDIVIDUAL_EFFECT_WIDTH <- 60 / 72
INDIVIDUAL_EFFECT_HEIGHT <- 200 / 72


# ============================================================
# COMBINED FIGURE GEOMETRY
# ============================================================

COMBINED_BAR_WIDTH <- 0.28
COMBINED_BAR_BORDER_WIDTH <- 0.16

COMBINED_CI_LINE_WIDTH <- 0.26
COMBINED_CI_CAP_WIDTH <- 0.12

COMBINED_ZERO_LINE_WIDTH <- 0.26

COMBINED_STAR_SIZE <- 2.0


# ============================================================
# INDIVIDUAL FIGURE GEOMETRY
#
# Intentionally much thinner.
# ============================================================

INDIVIDUAL_BAR_WIDTH <- 0.18
INDIVIDUAL_BAR_BORDER_WIDTH <- 0.10

INDIVIDUAL_CI_LINE_WIDTH <- 0.16
INDIVIDUAL_CI_CAP_WIDTH <- 0.065

INDIVIDUAL_ZERO_LINE_WIDTH <- 0.16

INDIVIDUAL_STAR_SIZE <- 1.25


# ============================================================
# HELPER:
# STANDARD ERROR
# ============================================================

standard_error <- function(x) {

  x <- x[
    is.finite(x)
  ]

  if (length(x) <= 1) {
    return(0)
  }

  sd(x) / sqrt(length(x))
}


# ============================================================
# HELPER:
# EXTRACT MODEL TERM
# ============================================================

extract_model_term <- function(
  tidy_table,
  term_name,
  column_name
) {

  result <- tidy_table %>%
    filter(
      .data$term == term_name
    ) %>%
    pull(
      all_of(column_name)
    )

  if (length(result) != 1) {
    return(NA_real_)
  }

  as.numeric(
    result[[1]]
  )
}


# ============================================================
# HELPER:
# NESTED MODEL F-TEST
# ============================================================

safe_nested_model_p <- function(
  smaller_model,
  larger_model
) {

  comparison <- tryCatch(
    anova(
      smaller_model,
      larger_model
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(comparison) ||
      nrow(comparison) < 2 ||
      !"Pr(>F)" %in% colnames(comparison)
  ) {
    return(NA_real_)
  }

  as.numeric(
    comparison[2, "Pr(>F)"]
  )
}


# ============================================================
# HELPER:
# MODEL GLANCE VALUE
# ============================================================

safe_glance_value <- function(
  model,
  column_name
) {

  model_glance <- tryCatch(
    broom::glance(model),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(model_glance) ||
      !column_name %in% colnames(model_glance)
  ) {
    return(NA_real_)
  }

  as.numeric(
    model_glance[[column_name]][1]
  )
}


# ============================================================
# HELPER:
# GENERAL LINEAR CONTRAST
# ============================================================

safe_linear_contrast <- function(
  model,
  weights,
  confidence_level = 0.95
) {

  coefficients <- coef(model)

  covariance_matrix <- vcov(model)

  if (
    any(
      !names(weights) %in%
        names(coefficients)
    )
  ) {

    return(
      tibble(
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  contrast_vector <- setNames(
    rep(
      0,
      length(coefficients)
    ),
    names(coefficients)
  )

  contrast_vector[
    names(weights)
  ] <- weights

  if (
    any(
      !is.finite(coefficients)
    )
  ) {

    return(
      tibble(
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  estimate <- sum(
    contrast_vector *
      coefficients
  )

  variance <- as.numeric(
    t(contrast_vector) %*%
      covariance_matrix %*%
      contrast_vector
  )

  if (
    !is.finite(variance) ||
      variance < 0
  ) {

    return(
      tibble(
        estimate = estimate,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  std_error <- sqrt(
    variance
  )

  residual_df <- df.residual(
    model
  )

  if (
    !is.finite(std_error) ||
      std_error <= 0 ||
      residual_df <= 0
  ) {

    return(
      tibble(
        estimate = estimate,
        std_error = std_error,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  statistic <- estimate /
    std_error

  p_value <- 2 *
    pt(
      abs(statistic),
      df = residual_df,
      lower.tail = FALSE
    )

  alpha <- 1 -
    confidence_level

  critical_value <- qt(
    1 -
      alpha / 2,
    df = residual_df
  )

  tibble(
    estimate = estimate,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,
    conf_low = estimate -
      critical_value *
      std_error,
    conf_high = estimate +
      critical_value *
      std_error
  )
}


# ============================================================
# HELPER:
# SPLINE ENDPOINT CONTRAST
#
# Average over R1/R2:
#
#   fitted CLR at 8 h
#             -
#   fitted CLR at 0.5 h
#
# Returns:
#
#   estimate
#   SE
#   t
#   P
#   95% CI
# ============================================================

safe_spline_endpoint_contrast <- function(
  model,
  time_early = 0.5,
  time_late = 8,
  confidence_level = 0.95
) {

  early_data <- tibble(
    time_hr = rep(
      time_early,
      length(replicate_order)
    ),

    replicate = factor(
      replicate_order,
      levels = replicate_order
    )
  )

  late_data <- tibble(
    time_hr = rep(
      time_late,
      length(replicate_order)
    ),

    replicate = factor(
      replicate_order,
      levels = replicate_order
    )
  )

  model_terms <- delete.response(
    terms(model)
  )

  X_early <- tryCatch(
    model.matrix(
      model_terms,
      early_data,
      contrasts.arg = model$contrasts,
      xlev = model$xlevels
    ),
    error = function(e) {
      NULL
    }
  )

  X_late <- tryCatch(
    model.matrix(
      model_terms,
      late_data,
      contrasts.arg = model$contrasts,
      xlev = model$xlevels
    ),
    error = function(e) {
      NULL
    }
  )

  if (
    is.null(X_early) ||
      is.null(X_late)
  ) {

    return(
      tibble(
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  mean_X_early <- colMeans(
    X_early
  )

  mean_X_late <- colMeans(
    X_late
  )

  contrast <- mean_X_late -
    mean_X_early

  coefficients <- coef(
    model
  )

  covariance_matrix <- vcov(
    model
  )

  if (
    !all(
      names(coefficients) %in%
        names(contrast)
    )
  ) {

    return(
      tibble(
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  contrast <- contrast[
    names(coefficients)
  ]

  if (
    any(
      !is.finite(contrast)
    ) ||
      any(
        !is.finite(coefficients)
      )
  ) {

    return(
      tibble(
        estimate = NA_real_,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  estimate <- sum(
    contrast *
      coefficients
  )

  variance <- as.numeric(
    t(contrast) %*%
      covariance_matrix %*%
      contrast
  )

  if (
    !is.finite(variance) ||
      variance < 0
  ) {

    return(
      tibble(
        estimate = estimate,
        std_error = NA_real_,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  std_error <- sqrt(
    variance
  )

  residual_df <- df.residual(
    model
  )

  if (
    !is.finite(std_error) ||
      std_error <= 0 ||
      residual_df <= 0
  ) {

    return(
      tibble(
        estimate = estimate,
        std_error = std_error,
        statistic = NA_real_,
        p_value = NA_real_,
        conf_low = NA_real_,
        conf_high = NA_real_
      )
    )
  }

  statistic <- estimate /
    std_error

  p_value <- 2 *
    pt(
      abs(statistic),
      df = residual_df,
      lower.tail = FALSE
    )

  alpha <- 1 -
    confidence_level

  critical_value <- qt(
    1 -
      alpha / 2,
    df = residual_df
  )

  tibble(
    estimate = estimate,
    std_error = std_error,
    statistic = statistic,
    p_value = p_value,

    conf_low = estimate -
      critical_value *
      std_error,

    conf_high = estimate +
      critical_value *
      std_error
  )
}


# ============================================================
# HELPER:
# MULTIPLICATIVE ZERO REPLACEMENT
# ============================================================

replace_zeros_multiplicative <- function(
  x,
  replacement_fraction = 0.50
) {

  x <- as.numeric(x)

  if (
    any(
      !is.finite(x)
    )
  ) {
    stop(
      "Non-finite value found during compositional transformation."
    )
  }

  if (
    any(
      x < 0
    )
  ) {
    stop(
      "Negative composition component found."
    )
  }

  total <- sum(x)

  if (
    !is.finite(total) ||
      total <= 0
  ) {
    stop(
      "Composition has a non-positive total."
    )
  }

  x <- x / total

  zero_positions <- x <= 0

  if (
    !any(zero_positions)
  ) {
    return(x)
  }

  positive_values <- x[
    !zero_positions
  ]

  if (
    length(positive_values) == 0
  ) {
    stop(
      "Composition contains no positive components."
    )
  }

  n_zeros <- sum(
    zero_positions
  )

  replacement_value <-
    replacement_fraction *
    min(
      positive_values
    )

  maximum_allowed_replacement <-
    0.95 /
    n_zeros

  replacement_value <- min(
    replacement_value,
    maximum_allowed_replacement
  )

  remaining_mass <-
    1 -
    n_zeros *
    replacement_value

  output <- numeric(
    length(x)
  )

  output[
    zero_positions
  ] <- replacement_value

  output[
    !zero_positions
  ] <-
    x[
      !zero_positions
    ] /
    sum(
      x[
        !zero_positions
      ]
    ) *
    remaining_mass

  output /
    sum(output)
}


# ============================================================
# HELPER:
# GET PREDICTION AT SPECIFIC TIME
# ============================================================

get_prediction_at_time <- function(
  prediction_table,
  requested_time
) {

  result <- prediction_table %>%
    filter(
      abs(
        time_hr -
          requested_time
      ) <
        1e-10
    ) %>%
    pull(
      predicted_clr
    )

  if (
    length(result) != 1
  ) {
    return(NA_real_)
  }

  as.numeric(
    result[[1]]
  )
}


# ============================================================
# READ INPUT
# ============================================================

message(
  "\nReading:\n  ",
  input_file
)

raw_data <- read_tsv(
  input_file,
  show_col_types = FALSE,
  progress = FALSE
)


# ============================================================
# REQUIRED COLUMNS
# ============================================================

required_columns <- c(
  "timepoint",
  "replicate",
  "tetranucleotide",
  "count",
  "total_count",
  "percent",
  "cpd",
  "left",
  "right",
  "context_number"
)

missing_columns <- setdiff(
  required_columns,
  names(raw_data)
)

if (
  length(missing_columns) > 0
) {

  stop(
    "\nInput table is missing required columns:\n  ",
    paste(
      missing_columns,
      collapse = ", "
    ),
    "\n"
  )
}


# ============================================================
# CLEAN INPUT
# ============================================================

repair_data <- raw_data %>%
  transmute(

    timepoint = factor(
      as.character(timepoint),
      levels = time_order
    ),

    replicate = factor(
      as.character(replicate),
      levels = replicate_order
    ),

    tetranucleotide = toupper(
      as.character(tetranucleotide)
    ),

    count = suppressWarnings(
      as.numeric(count)
    ),

    total_count = suppressWarnings(
      as.numeric(total_count)
    ),

    percent = suppressWarnings(
      as.numeric(percent)
    ),

    cpd = factor(
      toupper(
        as.character(cpd)
      ),
      levels = cpd_order
    ),

    left = toupper(
      as.character(left)
    ),

    right = toupper(
      as.character(right)
    ),

    context_number = suppressWarnings(
      as.integer(context_number)
    )
  ) %>%

  mutate(
    time_hr = unname(
      time_values[
        as.character(timepoint)
      ]
    )
  ) %>%

  filter(
    !is.na(timepoint),
    !is.na(replicate),
    !is.na(cpd),

    str_detect(
      tetranucleotide,
      "^[ACGT]{4}$"
    ),

    left %in% dna_order,
    right %in% dna_order,

    is.finite(count),
    count >= 0,

    is.finite(total_count),
    total_count > 0,

    is.finite(percent),
    percent >= 0,

    is.finite(time_hr),
    is.finite(context_number)
  ) %>%

  arrange(
    context_number,
    time_hr,
    replicate
  )


if (
  nrow(repair_data) == 0
) {

  stop(
    "No valid repair rows remained after cleaning."
  )
}


# ============================================================
# CONTEXT METADATA
# ============================================================

context_metadata <- repair_data %>%
  select(
    tetranucleotide,
    cpd,
    left,
    right,
    context_number
  ) %>%
  distinct() %>%
  arrange(
    context_number
  )


metadata_conflicts <- context_metadata %>%
  count(
    tetranucleotide,
    name = "n_metadata_rows"
  ) %>%
  filter(
    n_metadata_rows != 1
  )


if (
  nrow(metadata_conflicts) > 0
) {

  print(
    metadata_conflicts,
    n = Inf
  )

  stop(
    "Some tetranucleotides have inconsistent metadata."
  )
}


if (
  nrow(context_metadata) != 64
) {

  stop(
    "\nExpected 64 tetranucleotide contexts, but found ",
    nrow(context_metadata),
    ".\n"
  )
}


# ============================================================
# VERIFY METADATA AGAINST SEQUENCE
# ============================================================

metadata_sequence_qc <- context_metadata %>%
  mutate(

    expected_tetranucleotide = paste0(
      left,
      as.character(cpd),
      right
    ),

    metadata_matches_sequence =
      tetranucleotide ==
      expected_tetranucleotide
  )


if (
  any(
    !metadata_sequence_qc$
      metadata_matches_sequence
  )
) {

  print(
    metadata_sequence_qc %>%
      filter(
        !metadata_matches_sequence
      ),
    n = Inf
  )

  stop(
    "Some tetranucleotide / CPD / flank metadata do not agree."
  )
}


# ============================================================
# CPD CONTEXT QC
# ============================================================

cpd_context_qc <- context_metadata %>%
  count(
    cpd,
    name = "n_contexts"
  ) %>%
  complete(

    cpd = factor(
      cpd_order,
      levels = cpd_order
    ),

    fill = list(
      n_contexts = 0L
    )
  ) %>%
  arrange(
    cpd
  )


message(
  "\nContexts per CPD class:"
)

print(
  cpd_context_qc,
  n = Inf
)


if (
  any(
    cpd_context_qc$n_contexts != 16
  )
) {

  stop(
    "Each CPD class must contain exactly 16 tetranucleotides."
  )
}


# ============================================================
# SAMPLE QC
# ============================================================

sample_qc <- repair_data %>%
  group_by(
    timepoint,
    time_hr,
    replicate
  ) %>%
  summarise(

    n_contexts = n_distinct(
      tetranucleotide
    ),

    total_context_count = sum(
      count,
      na.rm = TRUE
    ),

    reported_total_count = first(
      total_count
    ),

    minimum_total_count = min(
      total_count,
      na.rm = TRUE
    ),

    maximum_total_count = max(
      total_count,
      na.rm = TRUE
    ),

    sum_percent = sum(
      percent,
      na.rm = TRUE
    ),

    minimum_percent = min(
      percent,
      na.rm = TRUE
    ),

    maximum_percent = max(
      percent,
      na.rm = TRUE
    ),

    n_zero_contexts = sum(
      percent <= 0,
      na.rm = TRUE
    ),

    .groups = "drop"
  ) %>%

  mutate(

    has_all_64_contexts =
      n_contexts == 64,

    total_count_is_constant =
      minimum_total_count ==
      maximum_total_count,

    count_matches_reported_total =
      abs(
        total_context_count -
          reported_total_count
      ) <= 1,

    percent_sums_to_100 =
      abs(
        sum_percent -
          100
      ) <= 0.05,

    all_values_nonnegative =
      minimum_percent >= 0
  )


message(
  "\nSample-level QC:"
)

print(
  sample_qc,
  n = Inf
)


bad_sample_qc <- sample_qc %>%
  filter(
    !has_all_64_contexts |
      !total_count_is_constant |
      !count_matches_reported_total |
      !percent_sums_to_100 |
      !all_values_nonnegative
  )


if (
  nrow(bad_sample_qc) > 0
) {

  message(
    "\nSamples failing strict QC:"
  )

  print(
    bad_sample_qc,
    n = Inf
  )

  stop(
    "\nAt least one sample failed QC.\n"
  )
}


# ============================================================
# COMPLETE DESIGN QC
# ============================================================

design_qc <- repair_data %>%
  count(
    tetranucleotide,
    timepoint,
    replicate,
    name = "n_rows"
  )


bad_design <- design_qc %>%
  filter(
    n_rows != 1
  )


expected_rows <-
  64 *
  length(time_order) *
  length(replicate_order)


if (
  nrow(repair_data) != expected_rows ||
    nrow(bad_design) > 0
) {

  print(
    bad_design,
    n = Inf
  )

  stop(
    "\nExpected ",
    expected_rows,
    " rows representing 64 contexts x 4 timepoints x ",
    "2 replicates, but found ",
    nrow(repair_data),
    ".\n"
  )
}


message(
  "\nConfirmed complete design: 64 contexts x 4 timepoints x 2 replicates."
)


# ============================================================
# RAW PERCENTAGE TIMECOURSE SUMMARY
# ============================================================

summary_time <- repair_data %>%
  group_by(
    timepoint,
    time_hr,
    tetranucleotide,
    cpd,
    left,
    right,
    context_number
  ) %>%
  summarise(

    mean_percent = mean(
      percent,
      na.rm = TRUE
    ),

    sd_percent = sd(
      percent,
      na.rm = TRUE
    ),

    n = sum(
      is.finite(percent)
    ),

    sem_percent = standard_error(
      percent
    ),

    .groups = "drop"
  ) %>%

  mutate(

    sd_percent = replace_na(
      sd_percent,
      0
    ),

    sem_percent = replace_na(
      sem_percent,
      0
    )
  )


# ============================================================
# RAW ENDPOINT SUMMARY
# ============================================================

raw_endpoint_summary <- summary_time %>%
  group_by(
    tetranucleotide
  ) %>%
  summarise(

    observed_mean_percent_0_5h =
      mean_percent[
        which(
          abs(
            time_hr - 0.5
          ) < 1e-10
        )
      ][1],

    observed_mean_percent_2h =
      mean_percent[
        which(
          abs(
            time_hr - 2
          ) < 1e-10
        )
      ][1],

    observed_mean_percent_4h =
      mean_percent[
        which(
          abs(
            time_hr - 4
          ) < 1e-10
        )
      ][1],

    observed_mean_percent_8h =
      mean_percent[
        which(
          abs(
            time_hr - 8
          ) < 1e-10
        )
      ][1],

    .groups = "drop"
  ) %>%

  mutate(

    observed_percent_change_0_5_to_8h =
      observed_mean_percent_8h -
      observed_mean_percent_0_5h,

    observed_percent_change_0_5_to_4h =
      observed_mean_percent_4h -
      observed_mean_percent_0_5h,

    observed_percent_change_4_to_8h =
      observed_mean_percent_8h -
      observed_mean_percent_4h
  )


# ============================================================
# CLR TRANSFORMATION
#
# ALL 64 CONTEXTS TOGETHER WITHIN EACH SAMPLE
# ============================================================

dat_clr <- repair_data %>%
  group_by(
    timepoint,
    time_hr,
    replicate
  ) %>%

  arrange(
    context_number,
    .by_group = TRUE
  ) %>%

  mutate(

    proportion =
      percent /
      sum(percent),

    adjusted_proportion =
      replace_zeros_multiplicative(
        proportion,
        replacement_fraction =
          ZERO_REPLACEMENT_FRACTION
      ),

    log_adjusted_proportion =
      log(
        adjusted_proportion
      ),

    clr =
      log_adjusted_proportion -
      mean(
        log_adjusted_proportion
      )
  ) %>%

  ungroup() %>%

  select(
    timepoint,
    time_hr,
    replicate,
    tetranucleotide,
    cpd,
    left,
    right,
    context_number,
    count,
    total_count,
    percent,
    proportion,
    adjusted_proportion,
    clr
  )


# ============================================================
# CLR QC
# ============================================================

clr_qc <- dat_clr %>%
  group_by(
    timepoint,
    time_hr,
    replicate
  ) %>%
  summarise(

    n_contexts = n(),

    n_zero_components_before_replacement =
      sum(
        proportion <= 0
      ),

    minimum_adjusted_proportion =
      min(
        adjusted_proportion
      ),

    sum_adjusted_proportion =
      sum(
        adjusted_proportion
      ),

    sum_clr =
      sum(clr),

    .groups = "drop"
  ) %>%

  mutate(

    zero_replacement_applied =
      n_zero_components_before_replacement > 0,

    adjusted_composition_sums_to_one =
      abs(
        sum_adjusted_proportion - 1
      ) < 1e-10,

    clr_sums_to_zero =
      abs(
        sum_clr
      ) < 1e-10
  )


message(
  "\nCLR transformation QC:"
)

print(
  clr_qc,
  n = Inf
)


if (
  any(
    !clr_qc$
      adjusted_composition_sums_to_one
  ) ||
    any(
      !clr_qc$
        clr_sums_to_zero
    )
) {

  stop(
    "CLR transformation QC failed."
  )
}


# ============================================================
# MODEL FITTING FUNCTION
# ============================================================

fit_context_models <- function(
  context_data
) {

  context_data <- context_data %>%
    mutate(
      replicate = factor(
        replicate,
        levels = replicate_order
      )
    )


  model_null <- lm(
    clr ~ replicate,
    data = context_data
  )


  model_linear <- lm(
    clr ~
      time_hr +
      replicate,
    data = context_data
  )


  model_spline <- lm(
    clr ~
      splines::ns(
        time_hr,
        df = SPLINE_DF
      ) +
      replicate,
    data = context_data
  )


  piecewise_data <- context_data %>%
    mutate(

      early_time =
        pmin(
          time_hr,
          PIECEWISE_KNOT_H
        ) -
        min(time_hr),

      late_time =
        pmax(
          time_hr -
            PIECEWISE_KNOT_H,
          0
        )
    )


  model_piecewise <- lm(
    clr ~
      early_time +
      late_time +
      replicate,
    data = piecewise_data
  )


  early_phase_data <- context_data %>%
    filter(
      time_hr <=
        PIECEWISE_KNOT_H
    )


  model_early_phase <- lm(
    clr ~
      time_hr +
      replicate,
    data = early_phase_data
  )


  list(
    null = model_null,
    linear = model_linear,
    spline = model_spline,
    piecewise = model_piecewise,
    early_phase = model_early_phase
  )
}


# ============================================================
# SPLINE PREDICTIONS AT OBSERVED TIMES
# ============================================================

predict_spline_at_observed_times <- function(
  model
) {

  prediction_data <- tidyr::expand_grid(

    time_hr = unname(
      time_values
    ),

    replicate = factor(
      replicate_order,
      levels = replicate_order
    )
  )


  prediction_data$
    predicted_clr <- as.numeric(
      predict(
        model,
        newdata = prediction_data
      )
    )


  prediction_data %>%
    group_by(
      time_hr
    ) %>%
    summarise(

      predicted_clr =
        mean(
          predicted_clr
        ),

      .groups = "drop"
    ) %>%
    arrange(
      time_hr
    )
}


# ============================================================
# FIT MODELS TO ALL 64 CONTEXTS
# ============================================================

message(
  "\nFitting CLR spline, linear, piecewise, and early-phase models..."
)


model_table <- dat_clr %>%
  group_by(
    tetranucleotide
  ) %>%
  nest() %>%

  mutate(

    fitted_models = map(
      data,
      fit_context_models
    ),

    model_null = map(
      fitted_models,
      "null"
    ),

    model_linear = map(
      fitted_models,
      "linear"
    ),

    model_spline = map(
      fitted_models,
      "spline"
    ),

    model_piecewise = map(
      fitted_models,
      "piecewise"
    ),

    model_early_phase = map(
      fitted_models,
      "early_phase"
    ),

    tidy_linear = map(
      model_linear,
      ~broom::tidy(
        .x,
        conf.int = TRUE,
        conf.level = 0.95
      )
    ),

    tidy_piecewise = map(
      model_piecewise,
      ~broom::tidy(
        .x,
        conf.int = TRUE,
        conf.level = 0.95
      )
    ),

    tidy_early_phase = map(
      model_early_phase,
      ~broom::tidy(
        .x,
        conf.int = TRUE,
        conf.level = 0.95
      )
    ),

    predicted_observed_times = map(
      model_spline,
      predict_spline_at_observed_times
    ),

    segment_slope_difference = map(
      model_piecewise,
      ~safe_linear_contrast(
        .x,
        weights = c(
          early_time = -1,
          late_time = 1
        )
      )
    ),

    endpoint_spline_contrast = map(
      model_spline,
      ~safe_spline_endpoint_contrast(
        .x,
        time_early = 0.5,
        time_late = 8,
        confidence_level = 0.95
      )
    )
  ) %>%

  unnest_wider(
    segment_slope_difference,
    names_sep = "_"
  ) %>%

  unnest_wider(
    endpoint_spline_contrast,
    names_sep = "_"
  )


# ============================================================
# EXTRACT MODEL RESULTS
# ============================================================

model_results <- model_table %>%

  mutate(

    # --------------------------------------------------------
    # PRIMARY OVERALL SPLINE TEST
    # --------------------------------------------------------

    p_time = map2_dbl(
      model_null,
      model_spline,
      safe_nested_model_p
    ),


    # --------------------------------------------------------
    # NONLINEARITY
    # --------------------------------------------------------

    p_curvature = map2_dbl(
      model_linear,
      model_spline,
      safe_nested_model_p
    ),


    # --------------------------------------------------------
    # ALL-TIME LINEAR
    # --------------------------------------------------------

    all_time_linear_slope = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_hr",
        "estimate"
      )
    ),

    all_time_linear_slope_se = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_hr",
        "std.error"
      )
    ),

    all_time_linear_p = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_hr",
        "p.value"
      )
    ),

    all_time_linear_CI95_lower = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_hr",
        "conf.low"
      )
    ),

    all_time_linear_CI95_upper = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_hr",
        "conf.high"
      )
    ),


    # --------------------------------------------------------
    # PIECEWISE EARLY
    # --------------------------------------------------------

    early_segment_slope_0_5_to_4h = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "early_time",
        "estimate"
      )
    ),

    early_segment_slope_se = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "early_time",
        "std.error"
      )
    ),

    early_segment_p = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "early_time",
        "p.value"
      )
    ),

    early_segment_CI95_lower = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "early_time",
        "conf.low"
      )
    ),

    early_segment_CI95_upper = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "early_time",
        "conf.high"
      )
    ),


    # --------------------------------------------------------
    # PIECEWISE LATE
    # --------------------------------------------------------

    late_segment_slope_4_to_8h = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "late_time",
        "estimate"
      )
    ),

    late_segment_slope_se = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "late_time",
        "std.error"
      )
    ),

    late_segment_p = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "late_time",
        "p.value"
      )
    ),

    late_segment_CI95_lower = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "late_time",
        "conf.low"
      )
    ),

    late_segment_CI95_upper = map_dbl(
      tidy_piecewise,
      ~extract_model_term(
        .x,
        "late_time",
        "conf.high"
      )
    ),


    # --------------------------------------------------------
    # EARLY PHASE
    # --------------------------------------------------------

    early_phase_linear_slope = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_hr",
        "estimate"
      )
    ),

    early_phase_linear_slope_se = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_hr",
        "std.error"
      )
    ),

    early_phase_linear_p = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_hr",
        "p.value"
      )
    ),

    early_phase_linear_CI95_lower = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_hr",
        "conf.low"
      )
    ),

    early_phase_linear_CI95_upper = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_hr",
        "conf.high"
      )
    ),


    # --------------------------------------------------------
    # SPLINE PREDICTIONS
    # --------------------------------------------------------

    predicted_CLR_0_5h = map_dbl(
      predicted_observed_times,
      ~get_prediction_at_time(
        .x,
        0.5
      )
    ),

    predicted_CLR_2h = map_dbl(
      predicted_observed_times,
      ~get_prediction_at_time(
        .x,
        2
      )
    ),

    predicted_CLR_4h = map_dbl(
      predicted_observed_times,
      ~get_prediction_at_time(
        .x,
        4
      )
    ),

    predicted_CLR_8h = map_dbl(
      predicted_observed_times,
      ~get_prediction_at_time(
        .x,
        8
      )
    ),

    fitted_CLR_change_0_5_to_8h =
      predicted_CLR_8h -
      predicted_CLR_0_5h,

    fitted_CLR_change_0_5_to_4h =
      predicted_CLR_4h -
      predicted_CLR_0_5h,

    fitted_CLR_change_4_to_8h =
      predicted_CLR_8h -
      predicted_CLR_4h,


    # --------------------------------------------------------
    # MODEL FIT
    # --------------------------------------------------------

    spline_r_squared = map_dbl(
      model_spline,
      ~safe_glance_value(
        .x,
        "r.squared"
      )
    ),

    spline_adjusted_r_squared = map_dbl(
      model_spline,
      ~safe_glance_value(
        .x,
        "adj.r.squared"
      )
    ),

    spline_residual_standard_error = map_dbl(
      model_spline,
      ~safe_glance_value(
        .x,
        "sigma"
      )
    ),

    piecewise_r_squared = map_dbl(
      model_piecewise,
      ~safe_glance_value(
        .x,
        "r.squared"
      )
    ),

    piecewise_adjusted_r_squared = map_dbl(
      model_piecewise,
      ~safe_glance_value(
        .x,
        "adj.r.squared"
      )
    )
  ) %>%

  ungroup() %>%

  left_join(
    context_metadata,
    by = "tetranucleotide"
  ) %>%

  mutate(

    # ========================================================
    # BH ADJUSTMENT ACROSS ALL 64 CONTEXTS
    # ========================================================

    padj_time = p.adjust(
      p_time,
      method = "BH"
    ),

    padj_curvature = p.adjust(
      p_curvature,
      method = "BH"
    ),

    padj_all_time_linear = p.adjust(
      all_time_linear_p,
      method = "BH"
    ),

    padj_early_segment = p.adjust(
      early_segment_p,
      method = "BH"
    ),

    padj_late_segment = p.adjust(
      late_segment_p,
      method = "BH"
    ),

    padj_segment_slope_difference = p.adjust(
      segment_slope_difference_p_value,
      method = "BH"
    ),

    padj_early_phase_linear = p.adjust(
      early_phase_linear_p,
      method = "BH"
    ),

    padj_endpoint_contrast = p.adjust(
      endpoint_spline_contrast_p_value,
      method = "BH"
    ),


    # ========================================================
    # PRIMARY CLASSIFICATION
    # ========================================================

    classification = case_when(

      is.finite(padj_time) &
        padj_time <= FDR_CUTOFF &
        fitted_CLR_change_0_5_to_8h < 0 ~
        "Early repair",

      is.finite(padj_time) &
        padj_time <= FDR_CUTOFF &
        fitted_CLR_change_0_5_to_8h > 0 ~
        "Late repair",

      TRUE ~
        "No significant trend"
    ),

    classification = factor(
      classification,
      levels = classification_order
    ),


    trajectory_shape = case_when(

      !is.finite(padj_time) |
        padj_time > FDR_CUTOFF ~
        "No significant time effect",

      is.finite(padj_curvature) &
        padj_curvature <= FDR_CUTOFF ~
        "Evidence of nonlinearity",

      TRUE ~
        "Approximately linear"
    ),

    trajectory_shape = factor(
      trajectory_shape,
      levels = shape_order
    ),


    late_to_early_absolute_slope_ratio =
      if_else(

        is.finite(
          early_segment_slope_0_5_to_4h
        ) &
          abs(
            early_segment_slope_0_5_to_4h
          ) >
          .Machine$double.eps,

        abs(
          late_segment_slope_4_to_8h
        ) /
          abs(
            early_segment_slope_0_5_to_4h
          ),

        NA_real_
      ),


    descriptive_plateau_pattern = case_when(

      is.finite(
        late_to_early_absolute_slope_ratio
      ) &
        late_to_early_absolute_slope_ratio <=
        PLATEAU_SLOPE_RATIO_CUTOFF ~
        "Plateau-like reduction in slope",

      TRUE ~
        "No clear plateau-like reduction"
    ),


    statistical_slope_change = case_when(

      is.finite(
        padj_segment_slope_difference
      ) &
        padj_segment_slope_difference <=
        FDR_CUTOFF ~
        "Significant early-versus-late slope difference",

      TRUE ~
        "No significant early-versus-late slope difference"
    ),


    sig_label = case_when(

      is.finite(padj_time) &
        padj_time <= FDR_CUTOFF ~
        "*",

      TRUE ~
        ""
    ),


    cpd = factor(
      as.character(cpd),
      levels = cpd_order
    )
  ) %>%

  left_join(
    raw_endpoint_summary,
    by = "tetranucleotide"
  )


# ============================================================
# PRIMARY RESULTS
# ============================================================

primary_results <- model_results %>%
  select(

    tetranucleotide,
    cpd,
    left,
    right,
    context_number,

    predicted_CLR_0_5h,
    predicted_CLR_2h,
    predicted_CLR_4h,
    predicted_CLR_8h,

    fitted_CLR_change_0_5_to_8h,
    fitted_CLR_change_0_5_to_4h,
    fitted_CLR_change_4_to_8h,

    endpoint_spline_contrast_estimate,
    endpoint_spline_contrast_std_error,
    endpoint_spline_contrast_statistic,
    endpoint_spline_contrast_p_value,
    padj_endpoint_contrast,
    endpoint_spline_contrast_conf_low,
    endpoint_spline_contrast_conf_high,

    p_time,
    padj_time,

    p_curvature,
    padj_curvature,

    classification,
    trajectory_shape,

    spline_r_squared,
    spline_adjusted_r_squared,
    spline_residual_standard_error,

    observed_mean_percent_0_5h,
    observed_mean_percent_2h,
    observed_mean_percent_4h,
    observed_mean_percent_8h,

    observed_percent_change_0_5_to_8h,
    observed_percent_change_0_5_to_4h,
    observed_percent_change_4_to_8h,

    sig_label
  ) %>%

  arrange(
    cpd,
    fitted_CLR_change_0_5_to_8h
  )


# ============================================================
# PIECEWISE RESULTS
# ============================================================

piecewise_results <- model_results %>%
  select(

    tetranucleotide,
    cpd,
    left,
    right,
    context_number,

    early_segment_slope_0_5_to_4h,
    early_segment_slope_se,
    early_segment_CI95_lower,
    early_segment_CI95_upper,
    early_segment_p,
    padj_early_segment,

    late_segment_slope_4_to_8h,
    late_segment_slope_se,
    late_segment_CI95_lower,
    late_segment_CI95_upper,
    late_segment_p,
    padj_late_segment,

    segment_slope_difference_estimate,
    segment_slope_difference_std_error,
    segment_slope_difference_statistic,
    segment_slope_difference_conf_low,
    segment_slope_difference_conf_high,
    segment_slope_difference_p_value,
    padj_segment_slope_difference,

    late_to_early_absolute_slope_ratio,
    descriptive_plateau_pattern,
    statistical_slope_change,

    piecewise_r_squared,
    piecewise_adjusted_r_squared
  ) %>%

  arrange(
    cpd,
    late_to_early_absolute_slope_ratio
  )


# ============================================================
# EARLY-PHASE RESULTS
# ============================================================

early_phase_results <- model_results %>%
  select(

    tetranucleotide,
    cpd,
    left,
    right,
    context_number,

    early_phase_linear_slope,
    early_phase_linear_slope_se,
    early_phase_linear_CI95_lower,
    early_phase_linear_CI95_upper,
    early_phase_linear_p,
    padj_early_phase_linear
  ) %>%

  arrange(
    cpd,
    early_phase_linear_slope
  )


# ============================================================
# ALL-TIME LINEAR RESULTS
# ============================================================

all_time_linear_results <- model_results %>%
  select(

    tetranucleotide,
    cpd,
    left,
    right,
    context_number,

    all_time_linear_slope,
    all_time_linear_slope_se,
    all_time_linear_CI95_lower,
    all_time_linear_CI95_upper,
    all_time_linear_p,
    padj_all_time_linear
  ) %>%

  arrange(
    cpd,
    all_time_linear_slope
  )


# ============================================================
# EARLY / LATE / NONSIGNIFICANT
# ============================================================

early_repair_contexts <- primary_results %>%
  filter(
    classification ==
      "Early repair"
  ) %>%
  arrange(
    cpd,
    fitted_CLR_change_0_5_to_8h
  )


late_repair_contexts <- primary_results %>%
  filter(
    classification ==
      "Late repair"
  ) %>%
  arrange(
    cpd,
    desc(
      fitted_CLR_change_0_5_to_8h
    )
  )


nonsignificant_contexts <- primary_results %>%
  filter(
    classification ==
      "No significant trend"
  ) %>%
  arrange(
    cpd,
    fitted_CLR_change_0_5_to_8h
  )


# ============================================================
# CLASSIFICATION SUMMARY
# ============================================================

classification_summary <- primary_results %>%
  count(
    cpd,
    classification,
    name = "n_contexts"
  ) %>%
  complete(

    cpd = factor(
      cpd_order,
      levels = cpd_order
    ),

    classification = factor(
      classification_order,
      levels = classification_order
    ),

    fill = list(
      n_contexts = 0L
    )
  ) %>%
  arrange(
    cpd,
    classification
  )


# ============================================================
# TRAJECTORY SHAPE SUMMARY
# ============================================================

shape_summary <- primary_results %>%
  count(
    cpd,
    trajectory_shape,
    name = "n_contexts"
  ) %>%
  complete(

    cpd = factor(
      cpd_order,
      levels = cpd_order
    ),

    trajectory_shape = factor(
      shape_order,
      levels = shape_order
    ),

    fill = list(
      n_contexts = 0L
    )
  ) %>%
  arrange(
    cpd,
    trajectory_shape
  )


# ============================================================
# ENDPOINT CONTRAST QC
# ============================================================

endpoint_contrast_qc <- primary_results %>%
  mutate(

    absolute_difference = abs(
      fitted_CLR_change_0_5_to_8h -
        endpoint_spline_contrast_estimate
    )
  ) %>%
  select(

    tetranucleotide,
    cpd,

    fitted_CLR_change_0_5_to_8h,
    endpoint_spline_contrast_estimate,

    endpoint_spline_contrast_conf_low,
    endpoint_spline_contrast_conf_high,

    absolute_difference
  )


message(
  "\nEndpoint contrast QC:"
)

print(
  endpoint_contrast_qc %>%
    arrange(
      desc(
        absolute_difference
      )
    ) %>%
    slice_head(
      n = 10
    ),
  n = Inf
)


# ============================================================
# WRITE TABLES
# ============================================================

write_tsv(
  repair_data %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "replicate_level_raw_percentages.tsv"
  )
)


write_tsv(
  dat_clr %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "replicate_level_CLR_values.tsv"
  )
)


write_tsv(
  sample_qc %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "sample_percentage_QC.tsv"
  )
)


write_tsv(
  clr_qc %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "CLR_transformation_QC.tsv"
  )
)


write_tsv(
  context_metadata %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "tetranucleotide_context_metadata.tsv"
  )
)


write_tsv(
  cpd_context_qc %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "CPD_context_QC.tsv"
  )
)


write_tsv(
  summary_time %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "raw_percentage_timecourse_mean_SEM.tsv"
  )
)


write_tsv(
  endpoint_contrast_qc %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_endpoint_contrast_QC.tsv"
  )
)


write_tsv(
  primary_results %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_primary_CLR_spline_results.tsv"
  )
)


write_tsv(
  piecewise_results %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_piecewise_4h_results.tsv"
  )
)


write_tsv(
  early_phase_results %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_early_phase_0_5_to_4h_results.tsv"
  )
)


write_tsv(
  all_time_linear_results %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_all_time_linear_results.tsv"
  )
)


write_tsv(
  early_repair_contexts %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_significant_EARLY_repair_tetranucleotides.tsv"
  )
)


write_tsv(
  late_repair_contexts %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_significant_LATE_repair_tetranucleotides.tsv"
  )
)


write_tsv(
  nonsignificant_contexts %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_nonsignificant_tetranucleotides.tsv"
  )
)


write_tsv(
  classification_summary %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_classification_summary.tsv"
  )
)


write_tsv(
  shape_summary %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_trajectory_shape_summary.tsv"
  )
)


# ============================================================
# ONE PRIMARY TABLE PER CPD CLASS
# ============================================================

for (
  cpd_name in cpd_order
) {

  primary_results %>%
    filter(
      cpd == cpd_name
    ) %>%
    arrange(
      fitted_CLR_change_0_5_to_8h
    ) %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ) %>%
    write_tsv(
      file.path(
        outdir,
        paste0(
          "WT_UV_",
          cpd_name,
          "_tetranucleotide_CLR_spline_results.tsv"
        )
      )
    )
}


# ============================================================
# ============================================================
#
# PLOTTING
#
# ============================================================
# ============================================================


# ============================================================
# COMBINED PLOT DATA
#
# Composite factor allows independent sequence labels in each
# CPD facet.
# ============================================================

plot_level_table <- primary_results %>%
  arrange(
    cpd,
    fitted_CLR_change_0_5_to_8h,
    tetranucleotide
  ) %>%
  transmute(
    plot_key = paste(
      as.character(cpd),
      tetranucleotide,
      sep = "___"
    )
  )


# Reverse so strongest negative / early context appears at top
# after coord_flip().
plot_levels <- rev(
  plot_level_table$plot_key
)


effect_plot_dat <- primary_results %>%
  mutate(

    plot_key = paste(
      as.character(cpd),
      tetranucleotide,
      sep = "___"
    ),

    plot_key = factor(
      plot_key,
      levels = plot_levels
    )
  ) %>%
  arrange(
    cpd,
    fitted_CLR_change_0_5_to_8h,
    tetranucleotide
  )


# ============================================================
# CPD-SPECIFIC EFFECT RANGES
# ============================================================

cpd_effect_ranges <- effect_plot_dat %>%
  group_by(
    cpd
  ) %>%
  summarise(

    minimum_effect = min(
      fitted_CLR_change_0_5_to_8h,
      na.rm = TRUE
    ),

    maximum_effect = max(
      fitted_CLR_change_0_5_to_8h,
      na.rm = TRUE
    ),

    minimum_CI95 = min(
      endpoint_spline_contrast_conf_low,
      na.rm = TRUE
    ),

    maximum_CI95 = max(
      endpoint_spline_contrast_conf_high,
      na.rm = TRUE
    ),

    .groups = "drop"
  )


write_tsv(
  cpd_effect_ranges %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "WT_UV_tetranucleotide_CPD_specific_effect_ranges.tsv"
  )
)


message(
  "\nCPD-specific effect ranges:"
)

print(
  cpd_effect_ranges,
  n = Inf
)


# ============================================================
# SIGNIFICANCE STAR POSITIONS
# ============================================================

effect_plot_dat <- effect_plot_dat %>%
  group_by(
    cpd
  ) %>%
  mutate(

    panel_minimum = min(
      endpoint_spline_contrast_conf_low,
      fitted_CLR_change_0_5_to_8h,
      na.rm = TRUE
    ),

    panel_maximum = max(
      endpoint_spline_contrast_conf_high,
      fitted_CLR_change_0_5_to_8h,
      na.rm = TRUE
    ),

    panel_range =
      panel_maximum -
      panel_minimum,

    panel_range = if_else(
      is.finite(panel_range) &
        panel_range > 0,
      panel_range,
      1
    ),

    star_offset =
      0.035 *
      panel_range,

    star_position = case_when(

      sig_label == "" ~
        NA_real_,

      fitted_CLR_change_0_5_to_8h < 0 ~
        endpoint_spline_contrast_conf_low -
        star_offset,

      fitted_CLR_change_0_5_to_8h > 0 ~
        endpoint_spline_contrast_conf_high +
        star_offset,

      TRUE ~
        NA_real_
    )
  ) %>%
  ungroup()


# ============================================================
# COMBINED FOUR-PANEL PLOT
# ============================================================

p_effect_combined <- ggplot(
  effect_plot_dat,
  aes(
    x = plot_key,
    y = fitted_CLR_change_0_5_to_8h,
    fill = classification
  )
) +

  geom_col(
    width = COMBINED_BAR_WIDTH,
    color = "black",
    linewidth = COMBINED_BAR_BORDER_WIDTH
  ) +

  geom_errorbar(
    aes(
      ymin = endpoint_spline_contrast_conf_low,
      ymax = endpoint_spline_contrast_conf_high
    ),
    width = COMBINED_CI_CAP_WIDTH,
    linewidth = COMBINED_CI_LINE_WIDTH,
    color = "black",
    na.rm = TRUE
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = COMBINED_ZERO_LINE_WIDTH,
    color = "black"
  ) +

  geom_text(
    aes(
      y = star_position,
      label = sig_label
    ),
    size = COMBINED_STAR_SIZE,
    fontface = "bold",
    color = "black",
    show.legend = FALSE,
    na.rm = TRUE
  ) +

  facet_wrap(
    facets = vars(cpd),
    nrow = 1,
    scales = "free",
    drop = TRUE
  ) +

  scale_fill_manual(
    values = classification_colors,
    breaks = classification_order,
    limits = classification_order,
    drop = FALSE,
    name = NULL
  ) +

  scale_x_discrete(
    drop = TRUE,

    labels = function(x) {
      sub(
        "^[A-Z]{2}___",
        "",
        x
      )
    },

    expand = expansion(
      add = c(
        0.25,
        0.25
      )
    )
  ) +

  scale_y_continuous(
    breaks = scales::pretty_breaks(
      n = 5
    ),
    expand = expansion(
      mult = c(
        0.10,
        0.10
      )
    )
  ) +

  coord_flip(
    clip = "off"
  ) +

  labs(
    title =
      "WT UV-CPD tetranucleotide repair timing",

    subtitle = paste0(
      "CLR spline model; bars = fitted 8 h - 0.5 h change; ",
      "error bars = 95% CI; BH-FDR <= ",
      FDR_CUTOFF
    ),

    x = NULL,

    y =
      "Spline-fitted CLR change (8 h - 0.5 h)"
  ) +

  guides(
    fill = guide_legend(
      title = NULL,
      nrow = 1,
      byrow = TRUE
    )
  ) +

  theme_classic(
    base_size = 8
  ) +

  theme(

    legend.position = "top",

    legend.justification = "center",

    legend.direction = "horizontal",

    legend.text = element_text(
      size = 6.5,
      color = "black"
    ),

    legend.key.width = unit(
      0.35,
      "cm"
    ),

    legend.key.height = unit(
      0.18,
      "cm"
    ),

    legend.spacing.x = unit(
      0.07,
      "cm"
    ),

    legend.margin = margin(
      t = 0,
      r = 0,
      b = -2,
      l = 0
    ),

    strip.background = element_rect(
      fill = "grey96",
      color = "grey70",
      linewidth = 0.22
    ),

    strip.text.x = element_text(
      size = 8,
      face = "bold",
      color = "black"
    ),

    panel.spacing.x = unit(
      0.65,
      "lines"
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey70",
      linewidth = 0.22
    ),

    axis.text.y = element_text(
      size = 6.0,
      color = "black"
    ),

    axis.text.x = element_text(
      size = 5.9,
      color = "black"
    ),

    axis.title.x = element_text(
      size = 7.2,
      color = "black",
      margin = margin(
        t = 3
      )
    ),

    axis.title.y = element_blank(),

    axis.line = element_line(
      color = "black",
      linewidth = 0.28
    ),

    axis.ticks = element_line(
      color = "black",
      linewidth = 0.25
    ),

    axis.ticks.length = unit(
      1.2,
      "pt"
    ),

    panel.grid = element_blank(),

    plot.title = element_text(
      hjust = 0.5,
      size = 9,
      face = "bold"
    ),

    plot.subtitle = element_text(
      hjust = 0.5,
      size = 6.1,
      color = "black"
    ),

    plot.margin = margin(
      t = 2,
      r = 5,
      b = 3,
      l = 4
    )
  )


# ============================================================
# SAVE COMBINED FIGURE
# ============================================================

combined_pdf <- file.path(
  outdir,
  "WT_UV_tetranucleotide_CLR_spline_effect_ALL_CPD_BARS_95CI.pdf"
)

combined_png <- file.path(
  outdir,
  "WT_UV_tetranucleotide_CLR_spline_effect_ALL_CPD_BARS_95CI.png"
)


ggsave(
  filename = combined_pdf,
  plot = p_effect_combined,
  width = COMBINED_EFFECT_WIDTH,
  height = COMBINED_EFFECT_HEIGHT,
  units = "in",
  device = cairo_pdf
)


ggsave(
  filename = combined_png,
  plot = p_effect_combined,
  width = COMBINED_EFFECT_WIDTH,
  height = COMBINED_EFFECT_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)


# ============================================================
# ============================================================
#
# INDIVIDUAL 60 pt x 200 pt FIGURES
#
# ============================================================
# ============================================================


# ============================================================
# INDIVIDUAL PLOT FUNCTION
#
# This function deliberately rebuilds the factor ordering
# independently for each CPD class.
#
# That avoids the odd spacing / compressed axes caused by
# retaining factor levels from the other three CPD classes.
# ============================================================

make_individual_CPD_plot <- function(
  cpd_name
) {

  individual_dat <- effect_plot_dat %>%
    filter(
      as.character(cpd) ==
        cpd_name
    ) %>%
    arrange(
      fitted_CLR_change_0_5_to_8h,
      tetranucleotide
    )


  # ----------------------------------------------------------
  # EXACT LOCAL ORDER
  #
  # Most negative = top
  # Most positive = bottom
  # ----------------------------------------------------------

  local_levels <- individual_dat %>%
    arrange(
      fitted_CLR_change_0_5_to_8h,
      tetranucleotide
    ) %>%
    pull(
      tetranucleotide
    ) %>%
    as.character()


  individual_dat <- individual_dat %>%
    mutate(

      context_local = factor(
        tetranucleotide,
        levels = rev(
          local_levels
        )
      )
    )


  # ----------------------------------------------------------
  # Recalculate star offset within this individual panel.
  # ----------------------------------------------------------

  local_minimum <- min(
    individual_dat$
      endpoint_spline_contrast_conf_low,
    individual_dat$
      fitted_CLR_change_0_5_to_8h,
    na.rm = TRUE
  )


  local_maximum <- max(
    individual_dat$
      endpoint_spline_contrast_conf_high,
    individual_dat$
      fitted_CLR_change_0_5_to_8h,
    na.rm = TRUE
  )


  local_range <-
    local_maximum -
    local_minimum


  if (
    !is.finite(local_range) ||
      local_range <= 0
  ) {

    local_range <- 1
  }


  local_star_offset <-
    0.030 *
    local_range


  individual_dat <- individual_dat %>%
    mutate(

      local_star_position = case_when(

        sig_label == "" ~
          NA_real_,

        fitted_CLR_change_0_5_to_8h < 0 ~
          endpoint_spline_contrast_conf_low -
          local_star_offset,

        fitted_CLR_change_0_5_to_8h > 0 ~
          endpoint_spline_contrast_conf_high +
          local_star_offset,

        TRUE ~
          NA_real_
      )
    )


  # ----------------------------------------------------------
  # PLOT
  # ----------------------------------------------------------

  ggplot(
    individual_dat,
    aes(
      x = context_local,
      y = fitted_CLR_change_0_5_to_8h,
      fill = classification
    )
  ) +

    # --------------------------------------------------------
    # THIN BARS
    # --------------------------------------------------------

    geom_col(
      width = INDIVIDUAL_BAR_WIDTH,
      color = "black",
      linewidth =
        INDIVIDUAL_BAR_BORDER_WIDTH
    ) +

    # --------------------------------------------------------
    # THIN 95% CI
    # --------------------------------------------------------

    geom_errorbar(
      aes(
        ymin =
          endpoint_spline_contrast_conf_low,
        ymax =
          endpoint_spline_contrast_conf_high
      ),
      width =
        INDIVIDUAL_CI_CAP_WIDTH,
      linewidth =
        INDIVIDUAL_CI_LINE_WIDTH,
      color = "black",
      na.rm = TRUE
    ) +

    # --------------------------------------------------------
    # THIN ZERO LINE
    # --------------------------------------------------------

    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth =
        INDIVIDUAL_ZERO_LINE_WIDTH,
      color = "black"
    ) +

    # --------------------------------------------------------
    # SMALL FDR STAR
    # --------------------------------------------------------

    geom_text(
      aes(
        y =
          local_star_position,
        label =
          sig_label
      ),
      size =
        INDIVIDUAL_STAR_SIZE,
      fontface = "bold",
      color = "black",
      show.legend = FALSE,
      na.rm = TRUE
    ) +

    # --------------------------------------------------------
    # COLORS
    # --------------------------------------------------------

    scale_fill_manual(
      values =
        classification_colors,
      breaks =
        classification_order,
      limits =
        classification_order,
      drop = FALSE,
      name = NULL
    ) +

    # --------------------------------------------------------
    # CONTEXT AXIS
    # --------------------------------------------------------

    scale_x_discrete(
      drop = TRUE,
      expand = expansion(
        add = c(
          0.15,
          0.15
        )
      )
    ) +

    # --------------------------------------------------------
    # CLR AXIS
    #
    # Independent for each CPD.
    # --------------------------------------------------------

    scale_y_continuous(
      breaks =
        scales::pretty_breaks(
          n = 3
        ),
      expand =
        expansion(
          mult = c(
            0.12,
            0.12
          )
        )
    ) +

    coord_flip(
      clip = "off"
    ) +

    labs(
      title =
        cpd_name,

      x =
        NULL,

      # Intentionally short because figure is only 60 pt wide.
      y =
        "CLR change"
    ) +

    guides(
      fill = "none"
    ) +

    theme_classic(
      base_size = 5
    ) +

    theme(

      # ------------------------------------------------------
      # NO LEGEND
      # ------------------------------------------------------

      legend.position =
        "none",

      # ------------------------------------------------------
      # SEQUENCE LABELS
      # ------------------------------------------------------

      axis.text.y = element_text(
        size = 4.4,
        color = "black",
        margin = margin(
          r = 0.4,
          unit = "pt"
        )
      ),

      # ------------------------------------------------------
      # CLR TICK LABELS
      # ------------------------------------------------------

      axis.text.x = element_text(
        size = 4.0,
        color = "black",
        margin = margin(
          t = 0.5,
          unit = "pt"
        )
      ),

      # ------------------------------------------------------
      # AXIS TITLE
      # ------------------------------------------------------

      axis.title.x = element_text(
        size = 4.5,
        color = "black",
        margin = margin(
          t = 1.3,
          unit = "pt"
        )
      ),

      axis.title.y =
        element_blank(),

      # ------------------------------------------------------
      # CPD TITLE
      # ------------------------------------------------------

      plot.title = element_text(
        size = 5.5,
        face = "bold",
        hjust = 0.5,
        margin = margin(
          b = 1,
          unit = "pt"
        )
      ),

      plot.subtitle =
        element_blank(),

      # ------------------------------------------------------
      # THIN AXES
      # ------------------------------------------------------

      axis.line = element_line(
        linewidth = 0.16,
        color = "black"
      ),

      axis.ticks = element_line(
        linewidth = 0.14,
        color = "black"
      ),

      axis.ticks.length = unit(
        0.65,
        "pt"
      ),

      # ------------------------------------------------------
      # NO GRID / BORDER
      # ------------------------------------------------------

      panel.grid =
        element_blank(),

      panel.border =
        element_blank(),

      # ------------------------------------------------------
      # VERY SMALL MARGINS
      # ------------------------------------------------------

      plot.margin = margin(
        t = 1.2,
        r = 1.2,
        b = 1.5,
        l = 1.0,
        unit = "pt"
      )
    )
}

# ============================================================
# ============================================================
#
# INDIVIDUAL CPD FIGURES
#
# EXACT SIZE:
#
#   60 pt wide
#   200 pt high
#
# Y-axis labels:
#
#   RIGHT side
#
# Sequence labels:
#
#   left..right
#
# Examples:
#
#   ACTG  -> A..G   for CT
#   GTTA  -> G..A   for TT
#
# Bars:
#
#   thicker than previous version
#
# ============================================================
# ============================================================


# ============================================================
# INDIVIDUAL OUTPUT DIRECTORY
# ============================================================

individual_plot_dir <- file.path(
  outdir,
  "individual_CPD_60x200pt"
)

dir.create(
  individual_plot_dir,
  showWarnings = FALSE,
  recursive = TRUE
)


# ============================================================
# EXACT PDF DIMENSIONS
#
# PDF point = 1 / 72 inch
# ============================================================

INDIVIDUAL_WIDTH_PT <- 60
INDIVIDUAL_HEIGHT_PT <- 200

INDIVIDUAL_WIDTH_IN <-
  INDIVIDUAL_WIDTH_PT / 72

INDIVIDUAL_HEIGHT_IN <-
  INDIVIDUAL_HEIGHT_PT / 72


# ============================================================
# PNG DIMENSIONS AT 600 DPI
# ============================================================

INDIVIDUAL_WIDTH_PX_600DPI <- round(
  INDIVIDUAL_WIDTH_IN * 600
)

INDIVIDUAL_HEIGHT_PX_600DPI <- round(
  INDIVIDUAL_HEIGHT_IN * 600
)


# ============================================================
# INDIVIDUAL BAR GEOMETRY
#
# Thicker bars
# ============================================================

INDIVIDUAL_BAR_WIDTH <- 0.46

INDIVIDUAL_BAR_BORDER_WIDTH <- 0.12

INDIVIDUAL_CI_LINE_WIDTH <- 0.18

INDIVIDUAL_CI_CAP_WIDTH <- 0.14

INDIVIDUAL_ZERO_LINE_WIDTH <- 0.17

INDIVIDUAL_STAR_SIZE <- 1.35


# ============================================================
# FUNCTION:
# INDIVIDUAL CPD PLOT
# ============================================================

make_individual_CPD_plot <- function(
  cpd_name
) {

  # ----------------------------------------------------------
  # Select one CPD class
  # ----------------------------------------------------------

  individual_dat <- primary_results %>%

    filter(
      as.character(cpd) ==
        cpd_name
    ) %>%

    mutate(

      # ------------------------------------------------------
      # Short display label
      #
      # Instead of:
      #
      #   ACTG
      #
      # show:
      #
      #   A..G
      #
      # The actual tetranucleotide is still retained in the
      # data and output tables.
      # ------------------------------------------------------

      flank_label =
        paste0(
          left,
          "..",
          right
        )
    ) %>%

    arrange(
      fitted_CLR_change_0_5_to_8h,
      tetranucleotide
    )


  # ----------------------------------------------------------
  # QC:
  # There should be exactly 16 flank combinations.
  # ----------------------------------------------------------

  if (
    nrow(individual_dat) != 16
  ) {

    stop(
      "Expected 16 contexts for ",
      cpd_name,
      ", but found ",
      nrow(individual_dat),
      "."
    )
  }


  if (
    n_distinct(
      individual_dat$flank_label
    ) != 16
  ) {

    stop(
      "Flanking-base labels are not unique for ",
      cpd_name,
      "."
    )
  }


  # ----------------------------------------------------------
  # ORDER
  #
  # Most negative CLR change at top.
  # Most positive CLR change at bottom.
  # ----------------------------------------------------------

  local_levels <- individual_dat %>%

    arrange(
      fitted_CLR_change_0_5_to_8h,
      tetranucleotide
    ) %>%

    pull(
      flank_label
    ) %>%

    as.character()


  individual_dat <- individual_dat %>%

    mutate(

      flank_label =
        factor(
          flank_label,
          levels =
            rev(
              local_levels
            )
        )
    )


  # ----------------------------------------------------------
  # LOCAL PANEL RANGE
  # ----------------------------------------------------------

  local_minimum <- min(
    c(
      individual_dat$
        endpoint_spline_contrast_conf_low,

      individual_dat$
        fitted_CLR_change_0_5_to_8h
    ),
    na.rm = TRUE
  )


  local_maximum <- max(
    c(
      individual_dat$
        endpoint_spline_contrast_conf_high,

      individual_dat$
        fitted_CLR_change_0_5_to_8h
    ),
    na.rm = TRUE
  )


  local_range <-
    local_maximum -
    local_minimum


  if (
    !is.finite(local_range) ||
      local_range <= 0
  ) {

    local_range <- 1
  }


  # ----------------------------------------------------------
  # STAR OFFSET
  # ----------------------------------------------------------

  local_star_offset <-
    0.030 *
    local_range


  individual_dat <- individual_dat %>%

    mutate(

      local_star_position = case_when(

        sig_label == "" ~
          NA_real_,

        fitted_CLR_change_0_5_to_8h < 0 ~
          endpoint_spline_contrast_conf_low -
          local_star_offset,

        fitted_CLR_change_0_5_to_8h > 0 ~
          endpoint_spline_contrast_conf_high +
          local_star_offset,

        TRUE ~
          NA_real_
      )
    )


  # ----------------------------------------------------------
  # PLOT
  #
  # IMPORTANT:
  #
  # x = categorical flank label
  # y = numerical CLR effect
  #
  # scale_x_discrete(position = "top")
  #
  # before coord_flip() moves the categorical axis to the
  # RIGHT after coord_flip().
  # ----------------------------------------------------------

  p <- ggplot(
    individual_dat,
    aes(
      x =
        flank_label,
      y =
        fitted_CLR_change_0_5_to_8h,
      fill =
        classification
    )
  ) +


    # ========================================================
    # THICKER BARS
    # ========================================================

    geom_col(
      width =
        INDIVIDUAL_BAR_WIDTH,
      color =
        "black",
      linewidth =
        INDIVIDUAL_BAR_BORDER_WIDTH
    ) +


    # ========================================================
    # 95% CI
    # ========================================================

    geom_errorbar(

      aes(
        ymin =
          endpoint_spline_contrast_conf_low,

        ymax =
          endpoint_spline_contrast_conf_high
      ),

      width =
        INDIVIDUAL_CI_CAP_WIDTH,

      linewidth =
        INDIVIDUAL_CI_LINE_WIDTH,

      color =
        "black",

      na.rm =
        TRUE
    ) +


    # ========================================================
    # ZERO LINE
    #
    # After coord_flip() this becomes vertical.
    # ========================================================

    geom_hline(
      yintercept =
        0,
      linetype =
        "dashed",
      linewidth =
        INDIVIDUAL_ZERO_LINE_WIDTH,
      color =
        "black"
    ) +


    # ========================================================
    # SIGNIFICANCE STAR
    #
    # Primary overall spline test.
    # ========================================================

    geom_text(

      aes(
        y =
          local_star_position,
        label =
          sig_label
      ),

      size =
        INDIVIDUAL_STAR_SIZE,

      fontface =
        "bold",

      color =
        "black",

      show.legend =
        FALSE,

      na.rm =
        TRUE
    ) +


    # ========================================================
    # CLASSIFICATION COLORS
    # ========================================================

    scale_fill_manual(

      values =
        classification_colors,

      breaks =
        classification_order,

      limits =
        classification_order,

      drop =
        FALSE,

      name =
        NULL
    ) +


    # ========================================================
    # CATEGORICAL AXIS
    #
    # position = "top" BEFORE coord_flip()
    #
    # becomes RIGHT-side y-axis AFTER coord_flip().
    # ========================================================

    scale_x_discrete(

      position =
        "top",

      drop =
        TRUE,

      expand =
        expansion(
          add = c(
            0.12,
            0.12
          )
        )
    ) +


    # ========================================================
    # CLR EFFECT SCALE
    # ========================================================

    scale_y_continuous(

      breaks =
        scales::pretty_breaks(
          n = 3
        ),

      expand =
        expansion(
          mult = c(
            0.12,
            0.12
          )
        )
    ) +


    # ========================================================
    # HORIZONTAL BARS
    # ========================================================

    coord_flip(
      clip =
        "off"
    ) +


    # ========================================================
    # LABELS
    # ========================================================

    labs(

      title =
        cpd_name,

      x =
        NULL,

      y =
        "CLR change"
    ) +


    guides(
      fill =
        "none"
    ) +


    # ========================================================
    # COMPACT 60 x 200 pt THEME
    # ========================================================

    theme_classic(
      base_size =
        5
    ) +

    theme(

      # ------------------------------------------------------
      # NO LEGEND
      # ------------------------------------------------------

      legend.position =
        "none",


      # ------------------------------------------------------
      # RIGHT-SIDE N..N LABELS
      #
      # After coord_flip(), these are the categorical
      # y-axis labels.
      # ------------------------------------------------------

      axis.text.y = element_text(
        size =
          4.5,
        color =
          "black",
        margin =
          margin(
            l = 0.8,
            r = 0,
            unit = "pt"
          )
      ),


      # ------------------------------------------------------
      # CLR NUMERIC LABELS
      # ------------------------------------------------------

      axis.text.x = element_text(
        size =
          4.0,
        color =
          "black",
        margin =
          margin(
            t = 0.5,
            unit = "pt"
          )
      ),


      # ------------------------------------------------------
      # CLR AXIS TITLE
      # ------------------------------------------------------

      axis.title.x = element_text(
        size =
          4.4,
        color =
          "black",
        margin =
          margin(
            t = 1.2,
            unit = "pt"
          )
      ),


      axis.title.y =
        element_blank(),


      # ------------------------------------------------------
      # CPD TITLE
      # ------------------------------------------------------

      plot.title = element_text(
        size =
          5.5,
        face =
          "bold",
        hjust =
          0.5,
        margin =
          margin(
            b = 1,
            unit = "pt"
          )
      ),


      plot.subtitle =
        element_blank(),


      # ------------------------------------------------------
      # AXES
      # ------------------------------------------------------

      axis.line = element_line(
        linewidth =
          0.16,
        color =
          "black"
      ),


      axis.ticks = element_line(
        linewidth =
          0.14,
        color =
          "black"
      ),


      axis.ticks.length = unit(
        0.65,
        "pt"
      ),


      # ------------------------------------------------------
      # NO GRID / BORDER
      # ------------------------------------------------------

      panel.grid =
        element_blank(),


      panel.border =
        element_blank(),


      # ------------------------------------------------------
      # TIGHT MARGINS
      #
      # Slightly more room on RIGHT because labels are there.
      # ------------------------------------------------------

      plot.margin = margin(
        t =
          1.2,
        r =
          2.5,
        b =
          1.5,
        l =
          0.8,
        unit =
          "pt"
      )
    )


  return(
    p
  )
}


# ============================================================
# GENERATE ALL FOUR INDIVIDUAL FIGURES
# ============================================================

individual_plot_index <- list()


for (
  cpd_name in
  cpd_order
) {

  message(
    "\nGenerating ",
    cpd_name,
    " individual plot..."
  )


  p_individual <- make_individual_CPD_plot(
    cpd_name
  )


  # ----------------------------------------------------------
  # FILE NAMES
  # ----------------------------------------------------------

  individual_pdf <- file.path(
    individual_plot_dir,
    paste0(
      "WT_UV_",
      cpd_name,
      "_CLR_spline_BARS_flanks_right_60x200pt.pdf"
    )
  )


  individual_png <- file.path(
    individual_plot_dir,
    paste0(
      "WT_UV_",
      cpd_name,
      "_CLR_spline_BARS_flanks_right_60x200pt.png"
    )
  )


  # ==========================================================
  # PDF
  #
  # Explicit Cairo PDF device:
  #
  # EXACT PAGE SIZE = 60 x 200 pt
  # ==========================================================

  grDevices::cairo_pdf(
    filename =
      individual_pdf,
    width =
      INDIVIDUAL_WIDTH_IN,
    height =
      INDIVIDUAL_HEIGHT_IN,
    onefile =
      FALSE
  )


  print(
    p_individual
  )


  grDevices::dev.off()


  # ==========================================================
  # PNG
  #
  # Equivalent physical size at 600 DPI.
  # ==========================================================

  ggsave(
    filename =
      individual_png,

    plot =
      p_individual,

    width =
      INDIVIDUAL_WIDTH_PX_600DPI,

    height =
      INDIVIDUAL_HEIGHT_PX_600DPI,

    units =
      "px",

    dpi =
      600,

    bg =
      "white"
  )


  # ----------------------------------------------------------
  # INDEX
  # ----------------------------------------------------------

  individual_plot_index[[cpd_name]] <- tibble(

    cpd =
      cpd_name,

    width_pt =
      INDIVIDUAL_WIDTH_PT,

    height_pt =
      INDIVIDUAL_HEIGHT_PT,

    width_in =
      INDIVIDUAL_WIDTH_IN,

    height_in =
      INDIVIDUAL_HEIGHT_IN,

    bar_width =
      INDIVIDUAL_BAR_WIDTH,

    labels =
      "left..right",

    categorical_axis_side =
      "right",

    pdf =
      individual_pdf,

    png =
      individual_png
  )
}


# ============================================================
# SAVE INDEX
# ============================================================

individual_plot_index <- bind_rows(
  individual_plot_index
)


write_tsv(
  individual_plot_index,
  file.path(
    individual_plot_dir,
    "WT_UV_individual_CPD_plot_index.tsv"
  )
)


# ============================================================
# PRINT OUTPUTS
# ============================================================

cat(
  "\n============================================================\n"
)

cat(
  "INDIVIDUAL FIGURES COMPLETE\n"
)

cat(
  "============================================================\n"
)


cat(
  "\nDimensions:\n",
  "  60 pt wide\n",
  "  200 pt high\n",
  sep = ""
)


cat(
  "\nLabels:\n",
  "  left..right\n",
  "  e.g. A..G\n",
  sep = ""
)


cat(
  "\nCategorical y-axis:\n",
  "  RIGHT side\n",
  sep = ""
)


cat(
  "\nBar width:\n  ",
  INDIVIDUAL_BAR_WIDTH,
  "\n",
  sep = ""
)


for (
  cpd_name in
  cpd_order
) {

  cat(
    "\n",
    cpd_name,
    ":\n  ",
    file.path(
      individual_plot_dir,
      paste0(
        "WT_UV_",
        cpd_name,
        "_CLR_spline_BARS_flanks_right_60x200pt.pdf"
      )
    ),
    "\n",
    sep = ""
  )
}