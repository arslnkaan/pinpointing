#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(broom)
  library(grid)
  library(splines)
})

# ============================================================
# NHF1 4NQO G>T TRINUCLEOTIDE REPAIR-TIMING ANALYSIS
#
# PRIMARY ANALYSIS
#
#   CLR-transformed contribution ~ ns(time_h, df = 2) + replicate
#
# The two-degree-of-freedom natural spline permits a nonlinear,
# plateau-like trajectory while remaining conservative for only
# four unique timepoints.
#
# Overall time test:
#
#   Null model:
#     CLR ~ replicate
#
#   Time model:
#     CLR ~ ns(time_h, df = 2) + replicate
#
# The models are compared with a nested F-test. The resulting
# p-values are BH-adjusted across the 16 trinucleotides.
#
# Direction:
#
#   Early repair:
#     significant overall time effect and fitted CLR at 8 h
#     is lower than fitted CLR at 0.5 h
#
#   Late repair:
#     significant overall time effect and fitted CLR at 8 h
#     is higher than fitted CLR at 0.5 h
#
#   No significant trend:
#     BH-FDR > 0.05
#
# SECONDARY ANALYSES
#
# 1. Nonlinearity test:
#
#      Linear model versus two-df spline model
#
# 2. Fixed-knot piecewise model:
#
#      CLR ~ early_time + late_time + replicate
#
#    early_time estimates the slope from 0.5 to 4 h.
#    late_time estimates the slope from 4 to 8 h.
#
# 3. Early-phase sensitivity analysis:
#
#      CLR ~ time_h + replicate
#
#    using only 0.5, 2, and 4 h.
#
# IMPORTANT
#
# Percentages are compositional because all 16 contexts sum to
# 100% within each sample. The primary models therefore use a
# centered log-ratio transformation. The results describe
# relative enrichment or depletion among trinucleotides, not an
# absolute lesion-removal rate.
#
# The replicate term assumes R1 and R2 are corresponding
# replicate series across the four timepoints.
# ============================================================

# ============================================================
# OUTPUT DIRECTORY
# ============================================================

outdir <- "4NQO_trinuc_spline_CLR_BH"

dir.create(
  outdir,
  showWarnings = FALSE,
  recursive = TRUE
)

# ============================================================
# SETTINGS
# ============================================================

DEFAULT_CONTEXT_ORDER <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

time_order <- c(
  "30m",
  "2h",
  "4h",
  "8h"
)

time_values <- c(
  "30m" = 0.5,
  "2h" = 2,
  "4h" = 4,
  "8h" = 8
)

replicate_order <- c(
  "R1",
  "R2"
)

FDR_CUTOFF <- 0.01

# The spline has two time-related degrees of freedom:
# one approximately linear component and one curvature component.
SPLINE_DF <- 2

# Fixed knot for the secondary piecewise model.
PIECEWISE_KNOT_H <- 4

# Descriptive threshold for identifying a plateau-like reduction
# in slope. This is not itself a hypothesis test.
PLATEAU_SLOPE_RATIO_CUTOFF <- 0.50

# For CLR transformation, zero components are replaced with:
#
#   ZERO_REPLACEMENT_FRACTION × smallest positive component
#
# within the same sample, followed by renormalization.
ZERO_REPLACEMENT_FRACTION <- 0.50

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

# Figure dimensions.
EFFECT_WIDTH <- 2.8
EFFECT_HEIGHT <- 3.65

RAW_TRAJECTORY_WIDTH <- 6.5
RAW_TRAJECTORY_HEIGHT <- 4.8

CLR_TRAJECTORY_WIDTH <- 6.5
CLR_TRAJECTORY_HEIGHT <- 4.8

HEATMAP_WIDTH <- 3.2
HEATMAP_HEIGHT <- 4.0

# ============================================================
# INPUT FILES
# ============================================================

files <- tribble(
  ~timepoint, ~time_h, ~replicate, ~change, ~file,

  "30m", 0.5, "R1", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-30m-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-30m-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "30m", 0.5, "R2", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-30m-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-30m-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "2h", 2, "R1", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-2h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-2h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "2h", 2, "R2", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-2h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-2h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "4h", 4, "R1", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-4h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-4h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "4h", 4, "R2", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-4h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-4h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "8h", 8, "R1", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-8h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-8h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "8h", 8, "R2", "G>T",
  "/work/users/a/r/arslank/NHF1-4NQO-8h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-8h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv"
)

# ============================================================
# HELPER FUNCTIONS
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

safe_glance_value <- function(
  model,
  column_name
) {

  model_glance <- tryCatch(
    broom::glance(
      model
    ),
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

safe_linear_contrast <- function(
  model,
  weights,
  confidence_level = 0.95
) {

  coefficients <- coef(
    model
  )

  covariance_matrix <- vcov(
    model
  )

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

replace_zeros_multiplicative <- function(
  x,
  replacement_fraction = 0.50
) {

  x <- as.numeric(
    x
  )

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

  total <- sum(
    x
  )

  if (
    !is.finite(total) ||
      total <= 0
  ) {
    stop(
      "Composition has a non-positive total."
    )
  }

  x <- x /
    total

  zero_positions <- x <= 0

  if (
    !any(zero_positions)
  ) {
    return(
      x
    )
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

  replacement_value <- replacement_fraction *
    min(
      positive_values
    )

  # Ensure that zero replacement cannot consume all mass.
  maximum_allowed_replacement <- 0.95 /
    n_zeros

  replacement_value <- min(
    replacement_value,
    maximum_allowed_replacement
  )

  remaining_mass <- 1 -
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
  ] <- x[
    !zero_positions
  ] /
    sum(
      x[
        !zero_positions
      ]
    ) *
    remaining_mass

  output /
    sum(
      output
    )
}

get_prediction_at_time <- function(
  prediction_table,
  requested_time
) {

  result <- prediction_table %>%
    filter(
      abs(
        time_h -
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
# READ ONE TRINUCLEOTIDE TABLE
# ============================================================

read_trinuc <- function(
  file,
  timepoint,
  time_h,
  replicate,
  change
) {

  if (
    !file.exists(file)
  ) {
    stop(
      "\nMissing file:\n  ",
      file,
      "\n"
    )
  }

  input_table <- read_csv(
    file,
    show_col_types = FALSE,
    progress = FALSE
  )

  context_col <- intersect(
    c(
      "Trinucleotide_Context_RC",
      "Trinucleotide_Context",
      "trinucleotide",
      "Trinucleotide",
      "context",
      "Context"
    ),
    colnames(input_table)
  )[1]

  pct_col <- intersect(
    c(
      "Percentage",
      "Percent",
      "percentage",
      "percent"
    ),
    colnames(input_table)
  )[1]

  if (
    is.na(context_col)
  ) {
    stop(
      "\nNo trinucleotide context column found in:\n  ",
      file,
      "\n"
    )
  }

  if (
    is.na(pct_col)
  ) {
    stop(
      "\nNo percentage column found in:\n  ",
      file,
      "\n"
    )
  }

  output <- input_table %>%
    transmute(
      timepoint = timepoint,

      time_h = as.numeric(
        time_h
      ),

      replicate = replicate,

      change = change,

      context = toupper(
        as.character(
          .data[[context_col]]
        )
      ),

      percent = suppressWarnings(
        as.numeric(
          .data[[pct_col]]
        )
      )
    )

  finite_percent <- output$percent[
    is.finite(
      output$percent
    )
  ]

  # Convert proportions to percentages when the largest value
  # is no greater than one.
  if (
    length(finite_percent) > 0 &&
      max(
        finite_percent,
        na.rm = TRUE
      ) <= 1
  ) {

    output <- output %>%
      mutate(
        percent = percent *
          100
      )
  }

  output
}

# ============================================================
# READ ALL FILES
# ============================================================

message(
  "\nReading 4NQO trinucleotide percentage tables..."
)

dat <- pmap_dfr(
  files,
  read_trinuc
) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = time_order
    ),

    replicate = factor(
      replicate,
      levels = replicate_order
    ),

    context = factor(
      context,
      levels = DEFAULT_CONTEXT_ORDER
    )
  ) %>%
  filter(
    !is.na(timepoint),
    !is.na(replicate),
    !is.na(context),
    is.finite(percent),
    is.finite(time_h)
  ) %>%
  arrange(
    timepoint,
    replicate,
    context
  )

if (
  nrow(dat) == 0
) {
  stop(
    "No valid rows remained after reading and cleaning."
  )
}

# ============================================================
# SAMPLE-LEVEL QC
# ============================================================

sample_qc <- dat %>%
  group_by(
    timepoint,
    time_h,
    replicate,
    change
  ) %>%
  summarise(
    n_contexts = n_distinct(
      context
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
    has_all_16_contexts =
      n_contexts ==
      length(
        DEFAULT_CONTEXT_ORDER
      ),

    percent_sums_to_100 =
      abs(
        sum_percent -
          100
      ) <=
      0.5,

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
    !has_all_16_contexts |
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
    "\nAt least one sample failed the context-count, ",
    "percentage-total, or nonnegative-value QC.\n"
  )
}

# ============================================================
# COMPLETE DESIGN QC
# ============================================================

design_qc <- dat %>%
  count(
    timepoint,
    replicate,
    context,
    name = "n_rows"
  )

bad_design <- design_qc %>%
  filter(
    n_rows != 1
  )

expected_rows <- length(
  time_order
) *
  length(
    replicate_order
  ) *
  length(
    DEFAULT_CONTEXT_ORDER
  )

if (
  nrow(dat) != expected_rows ||
    nrow(bad_design) > 0
) {

  print(
    bad_design,
    n = Inf
  )

  stop(
    "\nExpected ",
    expected_rows,
    " rows representing 4 timepoints × 2 replicates × ",
    "16 contexts, but found ",
    nrow(dat),
    ".\n"
  )
}

message(
  "\nConfirmed complete design: 16 contexts × 4 timepoints × ",
  "2 replicates."
)

# ============================================================
# RAW PERCENTAGE SUMMARY
# ============================================================

summary_time <- dat %>%
  group_by(
    timepoint,
    time_h,
    change,
    context
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
      is.finite(
        percent
      )
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

raw_endpoint_summary <- summary_time %>%
  group_by(
    change,
    context
  ) %>%
  summarise(
    observed_mean_percent_0_5h =
      mean_percent[
        which(
          abs(
            time_h -
              0.5
          ) <
            1e-10
        )
      ][1],

    observed_mean_percent_4h =
      mean_percent[
        which(
          abs(
            time_h -
              4
          ) <
            1e-10
        )
      ][1],

    observed_mean_percent_8h =
      mean_percent[
        which(
          abs(
            time_h -
              8
          ) <
            1e-10
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
# CENTERED LOG-RATIO TRANSFORMATION
# ============================================================

dat_clr <- dat %>%
  group_by(
    timepoint,
    time_h,
    replicate,
    change
  ) %>%
  arrange(
    context,
    .by_group = TRUE
  ) %>%
  mutate(
    proportion = percent /
      sum(
        percent
      ),

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
    time_h,
    replicate,
    change,
    context,
    percent,
    proportion,
    adjusted_proportion,
    clr
  )

clr_qc <- dat_clr %>%
  group_by(
    timepoint,
    time_h,
    replicate,
    change
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
      sum(
        clr
      ),

    .groups = "drop"
  ) %>%
  mutate(
    zero_replacement_applied =
      n_zero_components_before_replacement > 0,

    adjusted_composition_sums_to_one =
      abs(
        sum_adjusted_proportion -
          1
      ) <
      1e-10,

    clr_sums_to_zero =
      abs(
        sum_clr
      ) <
      1e-10
  )

message(
  "\nCLR transformation QC:"
)

print(
  clr_qc,
  n = Inf
)

# ============================================================
# MODEL-FITTING FUNCTIONS
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
    clr ~ time_h + replicate,
    data = context_data
  )

  model_spline <- lm(
    clr ~
      splines::ns(
        time_h,
        df = SPLINE_DF
      ) +
      replicate,
    data = context_data
  )

  piecewise_data <- context_data %>%
    mutate(
      early_time =
        pmin(
          time_h,
          PIECEWISE_KNOT_H
        ) -
        min(
          time_h
        ),

      late_time =
        pmax(
          time_h -
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
      time_h <= PIECEWISE_KNOT_H
    )

  model_early_phase <- lm(
    clr ~
      time_h +
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

predict_spline_at_observed_times <- function(
  model
) {

  prediction_data <- tidyr::expand_grid(
    time_h = unname(
      time_values
    ),

    replicate = factor(
      replicate_order,
      levels = replicate_order
    )
  )

  prediction_data$predicted_clr <- as.numeric(
    predict(
      model,
      newdata = prediction_data
    )
  )

  prediction_data %>%
    group_by(
      time_h
    ) %>%
    summarise(
      predicted_clr = mean(
        predicted_clr
      ),

      .groups = "drop"
    ) %>%
    arrange(
      time_h
    )
}

predict_spline_curve <- function(
  model
) {

  prediction_data <- tidyr::expand_grid(
    time_h = seq(
      min(
        unname(
          time_values
        )
      ),
      max(
        unname(
          time_values
        )
      ),
      length.out = 200
    ),

    replicate = factor(
      replicate_order,
      levels = replicate_order
    )
  )

  prediction_data$predicted_clr <- as.numeric(
    predict(
      model,
      newdata = prediction_data
    )
  )

  prediction_data %>%
    group_by(
      time_h
    ) %>%
    summarise(
      predicted_clr = mean(
        predicted_clr
      ),

      .groups = "drop"
    )
}

# ============================================================
# FIT MODELS FOR EACH TRINUCLEOTIDE
# ============================================================

message(
  "\nFitting CLR spline, linear, piecewise, and early-phase models..."
)

model_table <- dat_clr %>%
  group_by(
    change,
    context
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
    )
  ) %>%
  unnest_wider(
    segment_slope_difference,
    names_sep = "_"
  )

# ============================================================
# EXTRACT MODEL RESULTS
# ============================================================

model_results <- model_table %>%
  mutate(
    # --------------------------------------------------------
    # Primary overall spline test
    # --------------------------------------------------------

    p_time = map2_dbl(
      model_null,
      model_spline,
      safe_nested_model_p
    ),

    # --------------------------------------------------------
    # Does the spline improve over a straight line?
    # --------------------------------------------------------

    p_curvature = map2_dbl(
      model_linear,
      model_spline,
      safe_nested_model_p
    ),

    # --------------------------------------------------------
    # Straight-line model using all four timepoints
    # --------------------------------------------------------

    all_time_linear_slope = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_h",
        "estimate"
      )
    ),

    all_time_linear_slope_se = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_h",
        "std.error"
      )
    ),

    all_time_linear_p = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_h",
        "p.value"
      )
    ),

    all_time_linear_CI95_lower = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_h",
        "conf.low"
      )
    ),

    all_time_linear_CI95_upper = map_dbl(
      tidy_linear,
      ~extract_model_term(
        .x,
        "time_h",
        "conf.high"
      )
    ),

    # --------------------------------------------------------
    # Piecewise slopes
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
    # Early-phase linear sensitivity analysis: 0.5, 2, and 4 h
    # --------------------------------------------------------

    early_phase_linear_slope = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_h",
        "estimate"
      )
    ),

    early_phase_linear_slope_se = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_h",
        "std.error"
      )
    ),

    early_phase_linear_p = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_h",
        "p.value"
      )
    ),

    early_phase_linear_CI95_lower = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_h",
        "conf.low"
      )
    ),

    early_phase_linear_CI95_upper = map_dbl(
      tidy_early_phase,
      ~extract_model_term(
        .x,
        "time_h",
        "conf.high"
      )
    ),

    # --------------------------------------------------------
    # Spline-fitted values averaged over R1 and R2
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
    # Model fit statistics
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
  mutate(
    # BH adjustment is performed separately for each family of
    # 16 tests.

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

    classification = case_when(
      is.finite(
        padj_time
      ) &
        padj_time <= FDR_CUTOFF &
        fitted_CLR_change_0_5_to_8h < 0 ~
        "Early repair",

      is.finite(
        padj_time
      ) &
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
      !is.finite(
        padj_time
      ) |
        padj_time > FDR_CUTOFF ~
        "No significant time effect",

      is.finite(
        padj_curvature
      ) &
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

      is.finite(
        padj_time
      ) &
        padj_time <= FDR_CUTOFF ~
        "*",

      TRUE ~
        ""
    )
  ) %>%
  left_join(
    raw_endpoint_summary,
    by = c(
      "change",
      "context"
    )
  )

# ============================================================
# FINAL RESULT TABLES
# ============================================================

primary_results <- model_results %>%
  select(
    change,
    context,
    predicted_CLR_0_5h,
    predicted_CLR_2h,
    predicted_CLR_4h,
    predicted_CLR_8h,
    fitted_CLR_change_0_5_to_8h,
    fitted_CLR_change_0_5_to_4h,
    fitted_CLR_change_4_to_8h,
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
    observed_mean_percent_4h,
    observed_mean_percent_8h,
    observed_percent_change_0_5_to_8h,
    observed_percent_change_0_5_to_4h,
    observed_percent_change_4_to_8h,
    sig_label
  ) %>%
  arrange(
    fitted_CLR_change_0_5_to_8h
  )

piecewise_results <- model_results %>%
  select(
    change,
    context,
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
    late_to_early_absolute_slope_ratio
  )

early_phase_results <- model_results %>%
  select(
    change,
    context,
    early_phase_linear_slope,
    early_phase_linear_slope_se,
    early_phase_linear_CI95_lower,
    early_phase_linear_CI95_upper,
    early_phase_linear_p,
    padj_early_phase_linear
  ) %>%
  arrange(
    early_phase_linear_slope
  )

all_time_linear_results <- model_results %>%
  select(
    change,
    context,
    all_time_linear_slope,
    all_time_linear_slope_se,
    all_time_linear_CI95_lower,
    all_time_linear_CI95_upper,
    all_time_linear_p,
    padj_all_time_linear
  ) %>%
  arrange(
    all_time_linear_slope
  )

early_repair_contexts <- primary_results %>%
  filter(
    classification ==
      "Early repair"
  ) %>%
  arrange(
    fitted_CLR_change_0_5_to_8h
  )

late_repair_contexts <- primary_results %>%
  filter(
    classification ==
      "Late repair"
  ) %>%
  arrange(
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
    fitted_CLR_change_0_5_to_8h
  )

classification_summary <- tibble(
  classification = factor(
    classification_order,
    levels = classification_order
  )
) %>%
  left_join(
    primary_results %>%
      count(
        classification,
        name = "n_contexts"
      ),
    by = "classification"
  ) %>%
  mutate(
    n_contexts = replace_na(
      n_contexts,
      0L
    )
  )

shape_summary <- tibble(
  trajectory_shape = factor(
    shape_order,
    levels = shape_order
  )
) %>%
  left_join(
    primary_results %>%
      count(
        trajectory_shape,
        name = "n_contexts"
      ),
    by = "trajectory_shape"
  ) %>%
  mutate(
    n_contexts = replace_na(
      n_contexts,
      0L
    )
  )

# ============================================================
# WRITE TABLES
# ============================================================

write_tsv(
  dat %>%
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
  primary_results %>%
    mutate(
      across(
        where(is.factor),
        as.character
      )
    ),
  file.path(
    outdir,
    "4NQO_trinuc_primary_CLR_spline_results.tsv"
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
    "4NQO_trinuc_piecewise_4h_results.tsv"
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
    "4NQO_trinuc_early_phase_0_5_to_4h_results.tsv"
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
    "4NQO_trinuc_original_all_time_linear_results.tsv"
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
    "4NQO_significant_early_repair_trinucleotides.tsv"
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
    "4NQO_significant_late_repair_trinucleotides.tsv"
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
    "4NQO_nonsignificant_trinucleotides.tsv"
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
    "4NQO_trinuc_classification_summary.tsv"
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
    "4NQO_trinuc_trajectory_shape_summary.tsv"
  )
)

# ============================================================
# BASE THEME
# ============================================================

base_theme <- theme_classic(
  base_size = 8
) +
  theme(
    axis.text = element_text(
      size = 6.5,
      color = "black"
    ),

    axis.title = element_text(
      size = 7.5,
      color = "black"
    ),

    plot.title = element_text(
      hjust = 0.5,
      size = 8,
      face = "bold"
    ),

    plot.subtitle = element_blank(),

    axis.line = element_line(
      color = "black",
      linewidth = 0.35
    ),

    axis.ticks = element_line(
      color = "black",
      linewidth = 0.30
    ),

    axis.ticks.length = unit(
      1.4,
      "pt"
    ),

    panel.grid = element_blank()
  )

# ============================================================
# RANKED PRIMARY SPLINE EFFECT PLOT
# ============================================================

ranked_context_order <- primary_results %>%
  arrange(
    fitted_CLR_change_0_5_to_8h,
    context
  ) %>%
  pull(
    context
  ) %>%
  as.character()

effect_plot_dat <- primary_results %>%
  mutate(
    context_sorted = factor(
      as.character(
        context
      ),
      levels = ranked_context_order
    )
  ) %>%
  arrange(
    context_sorted
  )

effect_xlim <- max(
  abs(
    effect_plot_dat$fitted_CLR_change_0_5_to_8h
  ),
  na.rm = TRUE
) *
  1.35

if (
  !is.finite(effect_xlim) ||
    effect_xlim <= 0
) {
  effect_xlim <- 1
}

star_offset <- 0.055 *
  effect_xlim

effect_plot_dat <- effect_plot_dat %>%
  mutate(
    star_position = if_else(
      fitted_CLR_change_0_5_to_8h >= 0,

      fitted_CLR_change_0_5_to_8h +
        star_offset,

      fitted_CLR_change_0_5_to_8h -
        star_offset
    )
  )

p_effect <- ggplot(
  effect_plot_dat,
  aes(
    x = context_sorted,
    y = fitted_CLR_change_0_5_to_8h,
    fill = classification
  )
) +

  geom_col(
    width = 0.68,
    color = "black",
    linewidth = 0.22
  ) +

  geom_text(
    aes(
      y = star_position,
      label = sig_label
    ),
    size = 2.4,
    color = "black",
    fontface = "bold",
    show.legend = FALSE
  ) +

  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "black"
  ) +

  coord_flip(
    ylim = c(
      -effect_xlim,
      effect_xlim
    ),
    clip = "off"
  ) +

  scale_fill_manual(
    values = classification_colors,
    breaks = classification_order,
    limits = classification_order,
    drop = FALSE,
    name = NULL
  ) +

  labs(
    title = "Time-dependent shift in 4NQO trinucleotide context",

    x = NULL,

    y = "Spline-fitted CLR change (8 h − 0.5 h)",

    fill = NULL
  ) +

  base_theme +

  theme(
    legend.position = "top",

    legend.text = element_text(
      size = 6.3
    ),

    legend.key.size = unit(
      0.45,
      "lines"
    ),

    legend.margin = margin(
      0,
      0,
      0,
      0
    ),

    legend.box.margin = margin(
      0,
      0,
      0,
      0
    ),

    axis.text.y = element_text(
      size = 6.5,
      color = "black"
    ),

    axis.text.x = element_text(
      size = 6.2,
      color = "black"
    ),

    plot.margin = margin(
      3,
      9,
      4,
      3
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_CLR_spline_effect_ranked_BH_FDR.pdf"
  ),
  plot = p_effect,
  width = EFFECT_WIDTH,
  height = EFFECT_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_CLR_spline_effect_ranked_BH_FDR.png"
  ),
  plot = p_effect,
  width = EFFECT_WIDTH,
  height = EFFECT_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# RAW PERCENTAGE TRAJECTORY PLOT
# ============================================================

global_ymax <- max(
  summary_time$mean_percent +
    summary_time$sem_percent,
  na.rm = TRUE
)

global_ymax <- ceiling(
  global_ymax *
    1.05
)

if (
  !is.finite(global_ymax) ||
    global_ymax <= 0
) {
  global_ymax <- 1
}

p_raw_trajectory <- ggplot(
  summary_time,
  aes(
    x = time_h,
    y = mean_percent
  )
) +

  geom_line(
    linewidth = 0.40,
    color = "black"
  ) +

  geom_point(
    size = 1.0,
    color = "black"
  ) +

  geom_errorbar(
    aes(
      ymin = pmax(
        mean_percent -
          sem_percent,
        0
      ),

      ymax =
        mean_percent +
        sem_percent
    ),

    width = 0.12,
    linewidth = 0.22,
    color = "black"
  ) +

  facet_wrap(
    facets = vars(
      context
    ),
    scales = "fixed",
    ncol = 4
  ) +

  scale_x_continuous(
    breaks = unname(
      time_values
    ),

    labels = names(
      time_values
    )
  ) +

  scale_y_continuous(
    limits = c(
      0,
      global_ymax
    ),

    expand = expansion(
      mult = c(
        0,
        0.03
      )
    )
  ) +

  labs(
    title = "Observed 4NQO trinucleotide trajectories",

    x = "Repair time",

    y = "Repair contribution (%)"
  ) +

  base_theme +

  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.22
    ),

    strip.text = element_text(
      size = 5.8,
      face = "plain"
    ),

    axis.text.x = element_text(
      size = 5.5,
      angle = 45,
      hjust = 1
    ),

    panel.spacing = unit(
      0.20,
      "lines"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_raw_percent_trajectories.pdf"
  ),
  plot = p_raw_trajectory,
  width = RAW_TRAJECTORY_WIDTH,
  height = RAW_TRAJECTORY_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_raw_percent_trajectories.png"
  ),
  plot = p_raw_trajectory,
  width = RAW_TRAJECTORY_WIDTH,
  height = RAW_TRAJECTORY_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# CLR SPLINE FIT PLOT
# ============================================================

spline_curve_dat <- model_table %>%
  transmute(
    change,
    context,
    spline_curve = map(
      model_spline,
      predict_spline_curve
    )
  ) %>%
  unnest(
    spline_curve
  ) %>%
  mutate(
    context = factor(
      context,
      levels = DEFAULT_CONTEXT_ORDER
    )
  )

p_clr_trajectory <- ggplot() +

  geom_line(
    data = spline_curve_dat,
    aes(
      x = time_h,
      y = predicted_clr
    ),
    linewidth = 0.55,
    color = "black"
  ) +

  geom_point(
    data = dat_clr,
    aes(
      x = time_h,
      y = clr,
      shape = replicate
    ),
    size = 1.15,
    stroke = 0.25,
    color = "black"
  ) +

  facet_wrap(
    facets = vars(
      context
    ),
    scales = "fixed",
    ncol = 4
  ) +

  scale_x_continuous(
    breaks = unname(
      time_values
    ),

    labels = names(
      time_values
    )
  ) +

  scale_shape_manual(
    values = c(
      "R1" = 16,
      "R2" = 17
    ),
    breaks = replicate_order,
    name = NULL
  ) +

  labs(
    title = "CLR-transformed spline trajectories",

    x = "Repair time",

    y = "Centered log-ratio"
  ) +

  base_theme +

  theme(
    legend.position = "top",

    legend.text = element_text(
      size = 6.3
    ),

    legend.key.size = unit(
      0.45,
      "lines"
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "black",
      linewidth = 0.22
    ),

    strip.text = element_text(
      size = 5.8,
      face = "plain"
    ),

    axis.text.x = element_text(
      size = 5.5,
      angle = 45,
      hjust = 1
    ),

    panel.spacing = unit(
      0.20,
      "lines"
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_CLR_spline_trajectories.pdf"
  ),
  plot = p_clr_trajectory,
  width = CLR_TRAJECTORY_WIDTH,
  height = CLR_TRAJECTORY_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_CLR_spline_trajectories.png"
  ),
  plot = p_clr_trajectory,
  width = CLR_TRAJECTORY_WIDTH,
  height = CLR_TRAJECTORY_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# RAW PERCENTAGE HEATMAP
# ============================================================

p_heat <- summary_time %>%
  mutate(
    context = factor(
      context,
      levels = DEFAULT_CONTEXT_ORDER
    )
  ) %>%
  ggplot(
    aes(
      x = timepoint,
      y = context,
      fill = mean_percent
    )
  ) +

  geom_tile(
    color = "white",
    linewidth = 0.30
  ) +

  scale_fill_gradient(
    low = "white",
    high = "#08306B"
  ) +

  labs(
    title = "4NQO trinucleotide context over repair time",

    x = NULL,

    y = NULL,

    fill = "%"
  ) +

  base_theme +

  theme(
    legend.position = "right",

    legend.title = element_text(
      size = 6.5
    ),

    legend.text = element_text(
      size = 6
    )
  )

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_raw_percent_heatmap.pdf"
  ),
  plot = p_heat,
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(
    outdir,
    "4NQO_trinuc_raw_percent_heatmap.png"
  ),
  plot = p_heat,
  width = HEATMAP_WIDTH,
  height = HEATMAP_HEIGHT,
  units = "in",
  dpi = 600,
  bg = "white"
)

# ============================================================
# FINAL CONSOLE SUMMARY
# ============================================================

cat(
  "\nDone.\n"
)

cat(
  "\nPrimary model:\n",
  "  CLR contribution ~ ns(time_h, df = ",
  SPLINE_DF,
  ") + replicate\n",
  sep = ""
)

cat(
  "\nPrimary overall time test:\n",
  "  replicate-only model versus spline time model\n",
  sep = ""
)

cat(
  "\nBH-FDR cutoff:\n  ",
  FDR_CUTOFF,
  "\n",
  sep = ""
)

cat(
  "\nClassification summary:\n"
)

print(
  classification_summary,
  n = Inf
)

cat(
  "\nTrajectory-shape summary:\n"
)

print(
  shape_summary,
  n = Inf
)

cat(
  "\nSignificant early-repair contexts:\n"
)

print(
  early_repair_contexts %>%
    select(
      context,
      fitted_CLR_change_0_5_to_8h,
      p_time,
      padj_time,
      p_curvature,
      padj_curvature,
      trajectory_shape,
      observed_percent_change_0_5_to_8h
    ),
  n = Inf
)

cat(
  "\nSignificant late-repair contexts:\n"
)

print(
  late_repair_contexts %>%
    select(
      context,
      fitted_CLR_change_0_5_to_8h,
      p_time,
      padj_time,
      p_curvature,
      padj_curvature,
      trajectory_shape,
      observed_percent_change_0_5_to_8h
    ),
  n = Inf
)

cat(
  "\nContexts with evidence of nonlinearity:\n"
)

print(
  primary_results %>%
    filter(
      trajectory_shape ==
        "Evidence of nonlinearity"
    ) %>%
    select(
      context,
      fitted_CLR_change_0_5_to_8h,
      p_time,
      padj_time,
      p_curvature,
      padj_curvature,
      classification
    ),
  n = Inf
)

cat(
  "\nPrimary results:\n  ",
  file.path(
    outdir,
    "4NQO_trinuc_primary_CLR_spline_results.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nPiecewise plateau analysis:\n  ",
  file.path(
    outdir,
    "4NQO_trinuc_piecewise_4h_results.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nEarly-phase sensitivity analysis:\n  ",
  file.path(
    outdir,
    "4NQO_trinuc_early_phase_0_5_to_4h_results.tsv"
  ),
  "\n",
  sep = ""
)

cat(
  "\nPrimary ranked figure:\n  ",
  file.path(
    outdir,
    "4NQO_trinuc_CLR_spline_effect_ranked_BH_FDR.pdf"
  ),
  "\n",
  sep = ""
)

cat(
  "\nOutput directory:\n  ",
  outdir,
  "\n",
  sep = ""
)