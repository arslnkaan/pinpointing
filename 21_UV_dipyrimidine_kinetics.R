#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

base_dir <- "/work/users/a/r/arslank/uv"

infile <- file.path(
  base_dir,
  "UV_dipyrimidine_7th_from_3prime",
  "UV_dipyrimidine_7th_from_3prime_summary.tsv"
)

outdir <- file.path(
  base_dir,
  "UV_dipyrimidine_kinetics"
)

plotdir <- file.path(
  outdir,
  "plots"
)

individual_plotdir <- file.path(
  plotdir,
  "individual_dipyrimidines"
)

dir.create(
  individual_plotdir,
  recursive = TRUE,
  showWarnings = FALSE
)

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

dimer_order <- c(
  "CC",
  "TC",
  "CT",
  "TT"
)

dimer_colors <- c(
  "CC" = "#D98C00",
  "TC" = "#8E24AA",
  "CT" = "#C62828",
  "TT" = "#4F759B"
)

metric_order <- c(
  "All repair reads",
  "Dipyrimidine composition"
)

# ============================================================
# CHECK INPUT
# ============================================================

if (!file.exists(infile)) {
  stop(
    paste0(
      "Missing input file: ",
      infile
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

value_at_time <- function(
  time_vector,
  value_vector,
  target_time
) {
  selected <- value_vector[
    is.finite(time_vector) &
      abs(
        time_vector -
          target_time
      ) < 1e-8
  ]

  if (length(selected) == 0) {
    return(NA_real_)
  }

  selected[1]
}

linear_slope <- function(
  time_vector,
  value_vector
) {
  keep <- is.finite(time_vector) &
    is.finite(value_vector)

  x <- time_vector[keep]
  y <- value_vector[keep]

  if (length(unique(x)) < 2) {
    return(NA_real_)
  }

  unname(
    coef(
      lm(
        y ~ x
      )
    )[2]
  )
}

linear_r_squared <- function(
  time_vector,
  value_vector
) {
  keep <- is.finite(time_vector) &
    is.finite(value_vector)

  x <- time_vector[keep]
  y <- value_vector[keep]

  if (length(unique(x)) < 2) {
    return(NA_real_)
  }

  summary(
    lm(
      y ~ x
    )
  )$r.squared
}

classify_kinetics <- function(
  time_vector,
  value_vector
) {
  keep <- is.finite(time_vector) &
    is.finite(value_vector)

  x <- time_vector[keep]
  y <- value_vector[keep]

  if (length(y) < 2) {
    return("insufficient_data")
  }

  order_index <- order(x)
  y <- y[order_index]

  differences <- diff(y)

  if (all(differences == 0)) {
    return("stable")
  }

  if (
    all(differences >= 0) &&
      any(differences > 0)
  ) {
    return("non_decreasing")
  }

  if (
    all(differences <= 0) &&
      any(differences < 0)
  ) {
    return("non_increasing")
  }

  "non_monotonic"
}

mean_sem <- function(
  data,
  grouping_columns,
  value_column
) {
  data %>%
    group_by(
      across(
        all_of(
          grouping_columns
        )
      )
    ) %>%
    summarise(
      mean_value = mean(
        .data[[value_column]],
        na.rm = TRUE
      ),

      sem_value = if_else(
        sum(
          is.finite(
            .data[[value_column]]
          )
        ) > 1,

        sd(
          .data[[value_column]],
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                .data[[value_column]]
              )
            )
          ),

        NA_real_
      ),

      n_replicates = sum(
        is.finite(
          .data[[value_column]]
        )
      ),

      .groups = "drop"
    )
}

# ============================================================
# READ DATA
# ============================================================

dat <- read_tsv(
  infile,
  show_col_types = FALSE
)

required_columns <- c(
  "sample",
  "timepoint",
  "time_h",
  "dinucleotide",
  "count",
  "total_valid_reads",
  "total_dipyrimidine_reads",
  "percent_of_all_reads",
  "percent_of_dipyrimidine_reads"
)

missing_columns <- setdiff(
  required_columns,
  names(dat)
)

if (length(missing_columns) > 0) {
  stop(
    paste0(
      "Missing required columns: ",
      paste(
        missing_columns,
        collapse = ", "
      )
    )
  )
}

# Add replicate column if it does not exist.
#
# This allows the same script to work later if replicate-level
# dipyrimidine tables are generated.
if (!"replicate" %in% names(dat)) {
  dat <- dat %>%
    mutate(
      replicate = case_when(
        str_detect(
          sample,
          regex(
            "(^|[-_])r1([-_]|$)",
            ignore_case = TRUE
          )
        ) ~ "R1",

        str_detect(
          sample,
          regex(
            "(^|[-_])r2([-_]|$)",
            ignore_case = TRUE
          )
        ) ~ "R2",

        TRUE ~ "R1"
      )
    )
}

dat <- dat %>%
  transmute(
    sample = as.character(sample),
    replicate = as.character(replicate),
    timepoint = as.character(timepoint),
    time_h = as.numeric(time_h),

    dinucleotide = str_to_upper(
      as.character(dinucleotide)
    ),

    count = as.numeric(count),

    total_valid_reads = as.numeric(
      total_valid_reads
    ),

    total_dipyrimidine_reads = as.numeric(
      total_dipyrimidine_reads
    ),

    percent_of_all_reads = as.numeric(
      percent_of_all_reads
    ),

    percent_of_dipyrimidine_reads = as.numeric(
      percent_of_dipyrimidine_reads
    )
  ) %>%
  filter(
    dinucleotide %in%
      dimer_order
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
    replicate,
    dinucleotide,
    time_h
  )

# ============================================================
# REPLICATE QC
# ============================================================

replicate_qc <- dat %>%
  distinct(
    sample,
    replicate,
    timepoint,
    time_h
  ) %>%
  count(
    timepoint,
    time_h,
    name = "n_replicates"
  ) %>%
  arrange(
    time_h
  )

write_tsv(
  replicate_qc,
  file.path(
    outdir,
    "dipyrimidine_replicate_QC.tsv"
  )
)

# ============================================================
# LONG-FORM METRIC TABLE
#
# Metric 1:
#   percentage of all combined repair reads
#
# Metric 2:
#   percentage within CC + CT + TC + TT
# ============================================================

kinetics_long <- dat %>%
  select(
    sample,
    replicate,
    timepoint,
    time_h,
    dinucleotide,
    count,
    total_valid_reads,
    total_dipyrimidine_reads,
    percent_of_all_reads,
    percent_of_dipyrimidine_reads
  ) %>%
  pivot_longer(
    cols = c(
      percent_of_all_reads,
      percent_of_dipyrimidine_reads
    ),
    names_to = "metric_code",
    values_to = "percent"
  ) %>%
  mutate(
    metric = recode(
      metric_code,

      "percent_of_all_reads" =
        "All repair reads",

      "percent_of_dipyrimidine_reads" =
        "Dipyrimidine composition"
    ),

    metric = factor(
      metric,
      levels = metric_order
    )
  ) %>%
  group_by(
    replicate,
    dinucleotide,
    metric
  ) %>%
  arrange(
    time_h,
    .by_group = TRUE
  ) %>%
  mutate(
    baseline_percent_0.5h =
      value_at_time(
        time_h,
        percent,
        0.5
      ),

    percentage_point_change_vs_0.5h =
      percent -
      baseline_percent_0.5h,

    fold_change_vs_0.5h =
      safe_ratio(
        percent,
        baseline_percent_0.5h
      ),

    percent_change_vs_0.5h =
      if_else(
        is.finite(
          baseline_percent_0.5h
        ) &
          baseline_percent_0.5h > 0,

        100 *
          (
            percent -
              baseline_percent_0.5h
          ) /
          baseline_percent_0.5h,

        NA_real_
      ),

    log2_fold_change_vs_0.5h =
      safe_log2(
        fold_change_vs_0.5h
      ),

    log2_fold_change_for_plot =
      finite_for_plot(
        log2_fold_change_vs_0.5h
      ),

    previous_time_h = lag(
      time_h
    ),

    previous_percent = lag(
      percent
    ),

    interval_hours =
      time_h -
      previous_time_h,

    interval_percentage_point_change =
      percent -
      previous_percent,

    interval_percent_change =
      if_else(
        is.finite(
          previous_percent
        ) &
          previous_percent > 0,

        100 *
          (
            percent -
              previous_percent
          ) /
          previous_percent,

        NA_real_
      ),

    interval_fold_change =
      safe_ratio(
        percent,
        previous_percent
      ),

    interval_log2_fold_change =
      safe_log2(
        interval_fold_change
      ),

    interval_slope_pp_per_hour =
      if_else(
        is.finite(
          interval_hours
        ) &
          interval_hours > 0,

        interval_percentage_point_change /
          interval_hours,

        NA_real_
      ),

    interval = case_when(
      abs(
        previous_time_h -
          0.5
      ) < 1e-8 &
        abs(
          time_h -
            2
        ) < 1e-8 ~
        "0.5–2 h",

      abs(
        previous_time_h -
          2
      ) < 1e-8 &
        abs(
          time_h -
            4
        ) < 1e-8 ~
        "2–4 h",

      abs(
        previous_time_h -
          4
      ) < 1e-8 &
        abs(
          time_h -
            8
        ) < 1e-8 ~
        "4–8 h",

      TRUE ~
        NA_character_
    )
  ) %>%
  ungroup() %>%
  mutate(
    interval = factor(
      interval,
      levels = c(
        "0.5–2 h",
        "2–4 h",
        "4–8 h"
      )
    )
  )

write_tsv(
  kinetics_long,
  file.path(
    outdir,
    "dipyrimidine_kinetics_all_values.tsv"
  )
)

# ============================================================
# TIME-COURSE MEAN ± SEM
# ============================================================

timecourse_summary <- kinetics_long %>%
  group_by(
    metric,
    timepoint,
    time_h,
    dinucleotide
  ) %>%
  summarise(
    mean_percent = mean(
      percent,
      na.rm = TRUE
    ),

    sem_percent = if_else(
      sum(
        is.finite(
          percent
        )
      ) > 1,

      sd(
        percent,
        na.rm = TRUE
      ) /
        sqrt(
          sum(
            is.finite(
              percent
            )
          )
        ),

      NA_real_
    ),

    mean_percentage_point_change =
      mean(
        percentage_point_change_vs_0.5h,
        na.rm = TRUE
      ),

    sem_percentage_point_change =
      if_else(
        sum(
          is.finite(
            percentage_point_change_vs_0.5h
          )
        ) > 1,

        sd(
          percentage_point_change_vs_0.5h,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                percentage_point_change_vs_0.5h
              )
            )
          ),

        NA_real_
      ),

    mean_percent_change =
      mean(
        percent_change_vs_0.5h,
        na.rm = TRUE
      ),

    sem_percent_change =
      if_else(
        sum(
          is.finite(
            percent_change_vs_0.5h
          )
        ) > 1,

        sd(
          percent_change_vs_0.5h,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                percent_change_vs_0.5h
              )
            )
          ),

        NA_real_
      ),

    mean_log2_fold_change =
      mean(
        log2_fold_change_for_plot,
        na.rm = TRUE
      ),

    sem_log2_fold_change =
      if_else(
        sum(
          is.finite(
            log2_fold_change_for_plot
          )
        ) > 1,

        sd(
          log2_fold_change_for_plot,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                log2_fold_change_for_plot
              )
            )
          ),

        NA_real_
      ),

    n_replicates = n_distinct(
      replicate
    ),

    .groups = "drop"
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    dinucleotide = factor(
      dinucleotide,
      levels = dimer_order
    ),

    metric = factor(
      metric,
      levels = metric_order
    )
  ) %>%
  arrange(
    metric,
    dinucleotide,
    time_h
  )

write_tsv(
  timecourse_summary,
  file.path(
    outdir,
    "dipyrimidine_kinetics_timecourse_mean_sem.tsv"
  )
)

# ============================================================
# INTERVAL-SPECIFIC KINETICS
#
# 0.5–2 h
# 2–4 h
# 4–8 h
# ============================================================

interval_kinetics <- kinetics_long %>%
  filter(
    !is.na(interval)
  ) %>%
  group_by(
    metric,
    dinucleotide,
    interval
  ) %>%
  summarise(
    mean_percentage_point_change =
      mean(
        interval_percentage_point_change,
        na.rm = TRUE
      ),

    sem_percentage_point_change =
      if_else(
        sum(
          is.finite(
            interval_percentage_point_change
          )
        ) > 1,

        sd(
          interval_percentage_point_change,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                interval_percentage_point_change
              )
            )
          ),

        NA_real_
      ),

    mean_percent_change =
      mean(
        interval_percent_change,
        na.rm = TRUE
      ),

    mean_log2_fold_change =
      mean(
        finite_for_plot(
          interval_log2_fold_change
        ),
        na.rm = TRUE
      ),

    mean_slope_pp_per_hour =
      mean(
        interval_slope_pp_per_hour,
        na.rm = TRUE
      ),

    sem_slope_pp_per_hour =
      if_else(
        sum(
          is.finite(
            interval_slope_pp_per_hour
          )
        ) > 1,

        sd(
          interval_slope_pp_per_hour,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                interval_slope_pp_per_hour
              )
            )
          ),

        NA_real_
      ),

    n_replicates = n_distinct(
      replicate
    ),

    .groups = "drop"
  ) %>%
  mutate(
    dinucleotide = factor(
      dinucleotide,
      levels = dimer_order
    ),

    metric = factor(
      metric,
      levels = metric_order
    )
  ) %>%
  arrange(
    metric,
    dinucleotide,
    interval
  )

write_tsv(
  interval_kinetics,
  file.path(
    outdir,
    "dipyrimidine_interval_kinetics.tsv"
  )
)

# ============================================================
# EARLY VERSUS LATE
#
# Early = 0.5 h
# Late  = 8 h
#
# Positive log2(early/late):
#   early-enriched
#
# Negative:
#   late-enriched
# ============================================================

early_late_replicates <- kinetics_long %>%
  group_by(
    replicate,
    metric,
    dinucleotide
  ) %>%
  summarise(
    early_percent_0.5h =
      value_at_time(
        time_h,
        percent,
        0.5
      ),

    late_percent_8h =
      value_at_time(
        time_h,
        percent,
        8
      ),

    percentage_point_change_8h_minus_0.5h =
      late_percent_8h -
      early_percent_0.5h,

    percent_change_8h_vs_0.5h =
      if_else(
        is.finite(
          early_percent_0.5h
        ) &
          early_percent_0.5h > 0,

        100 *
          (
            late_percent_8h -
              early_percent_0.5h
          ) /
          early_percent_0.5h,

        NA_real_
      ),

    fold_change_8h_vs_0.5h =
      safe_ratio(
        late_percent_8h,
        early_percent_0.5h
      ),

    log2_fold_change_8h_vs_0.5h =
      safe_log2(
        fold_change_8h_vs_0.5h
      ),

    early_over_late =
      safe_ratio(
        early_percent_0.5h,
        late_percent_8h
      ),

    log2_early_over_late =
      safe_log2(
        early_over_late
      ),

    .groups = "drop"
  ) %>%
  mutate(
    log2_early_over_late_for_plot =
      finite_for_plot(
        log2_early_over_late
      ),

    kinetic_direction = case_when(
      log2_early_over_late > 0 ~
        "early_enriched",

      log2_early_over_late < 0 ~
        "late_enriched",

      log2_early_over_late == 0 ~
        "unchanged",

      TRUE ~
        "not_calculable"
    )
  )

early_late_summary <- early_late_replicates %>%
  group_by(
    metric,
    dinucleotide
  ) %>%
  summarise(
    mean_early_percent =
      mean(
        early_percent_0.5h,
        na.rm = TRUE
      ),

    mean_late_percent =
      mean(
        late_percent_8h,
        na.rm = TRUE
      ),

    mean_percentage_point_change =
      mean(
        percentage_point_change_8h_minus_0.5h,
        na.rm = TRUE
      ),

    mean_percent_change =
      mean(
        percent_change_8h_vs_0.5h,
        na.rm = TRUE
      ),

    mean_log2_fold_change_8h_vs_0.5h =
      mean(
        finite_for_plot(
          log2_fold_change_8h_vs_0.5h
        ),
        na.rm = TRUE
      ),

    mean_log2_early_over_late =
      mean(
        log2_early_over_late_for_plot,
        na.rm = TRUE
      ),

    sem_log2_early_over_late =
      if_else(
        sum(
          is.finite(
            log2_early_over_late_for_plot
          )
        ) > 1,

        sd(
          log2_early_over_late_for_plot,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                log2_early_over_late_for_plot
              )
            )
          ),

        NA_real_
      ),

    n_replicates = n_distinct(
      replicate
    ),

    .groups = "drop"
  ) %>%
  mutate(
    dinucleotide = factor(
      dinucleotide,
      levels = dimer_order
    ),

    metric = factor(
      metric,
      levels = metric_order
    )
  ) %>%
  arrange(
    metric,
    dinucleotide
  )

write_tsv(
  early_late_replicates,
  file.path(
    outdir,
    "dipyrimidine_early_late_replicates.tsv"
  )
)

write_tsv(
  early_late_summary,
  file.path(
    outdir,
    "dipyrimidine_early_late_summary.tsv"
  )
)

# ============================================================
# OVERALL KINETIC FEATURES
# ============================================================

kinetic_features_replicates <- kinetics_long %>%
  group_by(
    replicate,
    metric,
    dinucleotide
  ) %>%
  summarise(
    linear_slope_pp_per_hour =
      linear_slope(
        time_h,
        percent
      ),

    linear_r_squared =
      linear_r_squared(
        time_h,
        percent
      ),

    minimum_percent = min(
      percent,
      na.rm = TRUE
    ),

    maximum_percent = max(
      percent,
      na.rm = TRUE
    ),

    dynamic_range_percentage_points =
      maximum_percent -
      minimum_percent,

    peak_time_h = time_h[
      which.max(
        percent
      )
    ][1],

    minimum_time_h = time_h[
      which.min(
        percent
      )
    ][1],

    kinetic_pattern =
      classify_kinetics(
        time_h,
        percent
      ),

    .groups = "drop"
  )

kinetic_features_summary <-
  kinetic_features_replicates %>%
  group_by(
    metric,
    dinucleotide
  ) %>%
  summarise(
    mean_linear_slope_pp_per_hour =
      mean(
        linear_slope_pp_per_hour,
        na.rm = TRUE
      ),

    sem_linear_slope_pp_per_hour =
      if_else(
        sum(
          is.finite(
            linear_slope_pp_per_hour
          )
        ) > 1,

        sd(
          linear_slope_pp_per_hour,
          na.rm = TRUE
        ) /
          sqrt(
            sum(
              is.finite(
                linear_slope_pp_per_hour
              )
            )
          ),

        NA_real_
      ),

    mean_linear_r_squared =
      mean(
        linear_r_squared,
        na.rm = TRUE
      ),

    mean_dynamic_range_percentage_points =
      mean(
        dynamic_range_percentage_points,
        na.rm = TRUE
      ),

    most_common_kinetic_pattern = names(
      sort(
        table(
          kinetic_pattern
        ),
        decreasing = TRUE
      )
    )[1],

    n_replicates = n_distinct(
      replicate
    ),

    .groups = "drop"
  ) %>%
  mutate(
    metric = factor(
      metric,
      levels = metric_order
    ),

    dinucleotide = factor(
      dinucleotide,
      levels = dimer_order
    )
  ) %>%
  arrange(
    metric,
    dinucleotide
  )

write_tsv(
  kinetic_features_replicates,
  file.path(
    outdir,
    "dipyrimidine_kinetic_features_replicates.tsv"
  )
)

write_tsv(
  kinetic_features_summary,
  file.path(
    outdir,
    "dipyrimidine_kinetic_features_summary.tsv"
  )
)

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
# CONTRIBUTION TO ALL REPAIR READS
# ============================================================

plot_all_reads <- timecourse_summary %>%
  filter(
    metric ==
      "All repair reads"
  )

p_all_reads <- ggplot(
  plot_all_reads,
  aes(
    x = time_h,
    y = mean_percent,
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
  geom_errorbar(
    aes(
      ymin = mean_percent -
        sem_percent,

      ymax = mean_percent +
        sem_percent
    ),
    width = 0.12,
    na.rm = TRUE
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
    title = "Dipyrimidine repair kinetics",
    subtitle = "Contribution to all combined repair reads",
    x = "Repair time (h)",
    y = "All repair reads (%)"
  ) +
  common_theme

# ============================================================
# PLOT 2
# COMPOSITION WITHIN DIPYRIMIDINES
# ============================================================

plot_composition <- timecourse_summary %>%
  filter(
    metric ==
      "Dipyrimidine composition"
  )

p_composition <- ggplot(
  plot_composition,
  aes(
    x = time_h,
    y = mean_percent,
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
  geom_errorbar(
    aes(
      ymin = mean_percent -
        sem_percent,

      ymax = mean_percent +
        sem_percent
    ),
    width = 0.12,
    na.rm = TRUE
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
    )
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
    title = "Dipyrimidine composition kinetics",
    subtitle = "CC + TC + CT + TT = 100% at each time point",
    x = "Repair time (h)",
    y = "Dipyrimidine composition (%)"
  ) +
  common_theme

# ============================================================
# PLOT 3
# PERCENT CHANGE RELATIVE TO 0.5 h
#
# Main requested kinetics plot.
# ============================================================

change_plot_data <- timecourse_summary %>%
  filter(
    metric ==
      "All repair reads"
  )

p_percent_change <- ggplot(
  change_plot_data,
  aes(
    x = time_h,
    y = mean_percent_change,
    color = dinucleotide,
    group = dinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 2.3
  ) +
  geom_errorbar(
    aes(
      ymin = mean_percent_change -
        sem_percent_change,

      ymax = mean_percent_change +
        sem_percent_change
    ),
    width = 0.12,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~dinucleotide,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = dimer_colors,
    guide = "none"
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
    title = "Dipyrimidine percentage change over time",
    subtitle = "Change relative to the 0.5 h contribution to all repair reads",
    x = "Repair time (h)",
    y = "Change from 0.5 h (%)"
  ) +
  common_theme +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "grey65"
    ),

    strip.text = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 4
# LOG2 FOLD CHANGE RELATIVE TO 0.5 h
# ============================================================

p_log2_change <- ggplot(
  change_plot_data,
  aes(
    x = time_h,
    y = mean_log2_fold_change,
    color = dinucleotide,
    group = dinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_line(
    linewidth = 0.85
  ) +
  geom_point(
    size = 2.3
  ) +
  geom_errorbar(
    aes(
      ymin = mean_log2_fold_change -
        sem_log2_fold_change,

      ymax = mean_log2_fold_change +
        sem_log2_fold_change
    ),
    width = 0.12,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~dinucleotide,
    ncol = 2,
    scales = "free_y"
  ) +
  scale_color_manual(
    values = dimer_colors,
    guide = "none"
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
    title = "Dipyrimidine fold-change kinetics",
    subtitle = "Contribution to all repair reads relative to 0.5 h",
    x = "Repair time (h)",
    y = expression(
      log[2](
        percentage[t] /
          percentage["0.5 h"]
      )
    )
  ) +
  common_theme +
  theme(
    strip.background = element_rect(
      fill = "grey95",
      color = "grey65"
    ),

    strip.text = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 5
# INTERVAL-SPECIFIC SLOPES
# ============================================================

slope_plot_data <- interval_kinetics %>%
  filter(
    metric ==
      "All repair reads"
  )

p_interval_slopes <- ggplot(
  slope_plot_data,
  aes(
    x = interval,
    y = mean_slope_pp_per_hour,
    fill = dinucleotide
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
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = mean_slope_pp_per_hour -
        sem_slope_pp_per_hour,

      ymax = mean_slope_pp_per_hour +
        sem_slope_pp_per_hour
    ),
    width = 0.16,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~dinucleotide,
    nrow = 1,
    scales = "free_y"
  ) +
  scale_fill_manual(
    values = dimer_colors,
    guide = "none"
  ) +
  labs(
    title = "Interval-specific dipyrimidine kinetics",
    subtitle = "Contribution to all repair reads",
    x = NULL,
    y = "Slope (percentage points/hour)"
  ) +
  common_theme +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey65"
    ),

    strip.text = element_text(
      face = "bold"
    )
  )

# ============================================================
# PLOT 6
# EARLY/LATE
# ============================================================

early_late_plot_data <- early_late_summary %>%
  filter(
    metric ==
      "All repair reads"
  )

early_late_values <-
  early_late_plot_data$mean_log2_early_over_late

early_late_values <- early_late_values[
  is.finite(
    early_late_values
  )
]

if (length(early_late_values) > 0) {
  early_late_limit <- max(
    abs(
      early_late_values
    ),
    na.rm = TRUE
  )
} else {
  early_late_limit <- 1
}

if (
  !is.finite(
    early_late_limit
  ) ||
    early_late_limit <= 0
) {
  early_late_limit <- 1
}

early_late_limit <-
  early_late_limit * 1.10

p_early_late <- ggplot(
  early_late_plot_data,
  aes(
    x = dinucleotide,
    y = mean_log2_early_over_late,
    fill = dinucleotide
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
    linewidth = 0.2
  ) +
  geom_errorbar(
    aes(
      ymin = mean_log2_early_over_late -
        sem_log2_early_over_late,

      ymax = mean_log2_early_over_late +
        sem_log2_early_over_late
    ),
    width = 0.16,
    na.rm = TRUE
  ) +
  scale_fill_manual(
    values = dimer_colors,
    breaks = dimer_order,
    guide = "none"
  ) +
  coord_cartesian(
    ylim = c(
      -early_late_limit,
      early_late_limit
    )
  ) +
  labs(
    title = "Early versus late dipyrimidine repair",
    subtitle = "Positive values indicate early enrichment",
    x = NULL,
    y = expression(
      log[2](
        percentage["0.5 h"] /
          percentage["8 h"]
      )
    )
  ) +
  common_theme

# ============================================================
# COMBINED FIGURES
# ============================================================

raw_figure <- (
  p_all_reads /
    p_composition
) +
  plot_layout(
    guides = "collect"
  ) &
  theme(
    legend.position = "top"
  )

kinetics_figure <- (
  p_percent_change /
    p_log2_change /
    p_interval_slopes /
    p_early_late
) +
  plot_layout(
    heights = c(
      1,
      1,
      0.9,
      0.8
    )
  )

master_figure <- (
  p_all_reads /
    p_percent_change /
    p_interval_slopes /
    p_early_late
) +
  plot_layout(
    guides = "collect",
    heights = c(
      0.8,
      1.2,
      0.8,
      0.7
    )
  ) &
  theme(
    legend.position = "top"
  )

# ============================================================
# SAVE MAIN FIGURES
# ============================================================

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_raw_timecourses.pdf"
  ),
  plot = raw_figure,
  width = 7.5,
  height = 8,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_kinetics_analysis.pdf"
  ),
  plot = kinetics_figure,
  width = 8,
  height = 15,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_kinetics_master.pdf"
  ),
  plot = master_figure,
  width = 8,
  height = 14,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_kinetics_master.png"
  ),
  plot = master_figure,
  width = 8,
  height = 14,
  units = "in",
  dpi = 600
)

# ============================================================
# SAVE INDIVIDUAL FIGURES
# ============================================================

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_percent_of_all_reads.pdf"
  ),
  plot = p_all_reads,
  width = 6.5,
  height = 4.6,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_percent_change_vs_0.5h.pdf"
  ),
  plot = p_percent_change,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_log2_fold_change_vs_0.5h.pdf"
  ),
  plot = p_log2_change,
  width = 7,
  height = 6,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_interval_slopes.pdf"
  ),
  plot = p_interval_slopes,
  width = 9,
  height = 4.2,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "dipyrimidine_log2_early_over_late.pdf"
  ),
  plot = p_early_late,
  width = 5.8,
  height = 4.5,
  units = "in"
)

# ============================================================
# ONE FIGURE PER DIPYRIMIDINE
# ============================================================

for (
  current_dimer in dimer_order
) {
  raw_data <- timecourse_summary %>%
    filter(
      metric ==
        "All repair reads",

      as.character(
        dinucleotide
      ) ==
        current_dimer
    )

  slope_data <- interval_kinetics %>%
    filter(
      metric ==
        "All repair reads",

      as.character(
        dinucleotide
      ) ==
        current_dimer
    )

  early_late_data <- early_late_summary %>%
    filter(
      metric ==
        "All repair reads",

      as.character(
        dinucleotide
      ) ==
        current_dimer
    )

  current_color <-
    dimer_colors[
      current_dimer
    ]

  p_current_raw <- ggplot(
    raw_data,
    aes(
      x = time_h,
      y = mean_percent
    )
  ) +
    geom_line(
      linewidth = 0.9,
      color = current_color
    ) +
    geom_point(
      size = 2.5,
      color = current_color
    ) +
    geom_errorbar(
      aes(
        ymin = mean_percent -
          sem_percent,

        ymax = mean_percent +
          sem_percent
      ),
      width = 0.12,
      color = current_color,
      na.rm = TRUE
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
      title = paste0(
        current_dimer,
        " repair contribution"
      ),
      x = "Repair time (h)",
      y = "All repair reads (%)"
    ) +
    common_theme +
    theme(
      legend.position = "none"
    )

  p_current_change <- ggplot(
    raw_data,
    aes(
      x = time_h,
      y = mean_percent_change
    )
  ) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey40"
    ) +
    geom_line(
      linewidth = 0.9,
      color = current_color
    ) +
    geom_point(
      size = 2.5,
      color = current_color
    ) +
    geom_errorbar(
      aes(
        ymin = mean_percent_change -
          sem_percent_change,

        ymax = mean_percent_change +
          sem_percent_change
      ),
      width = 0.12,
      color = current_color,
      na.rm = TRUE
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
      title = paste0(
        current_dimer,
        " change from 0.5 h"
      ),
      x = "Repair time (h)",
      y = "Change from 0.5 h (%)"
    ) +
    common_theme +
    theme(
      legend.position = "none"
    )

  p_current_slope <- ggplot(
    slope_data,
    aes(
      x = interval,
      y = mean_slope_pp_per_hour
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
      fill = current_color,
      color = "black",
      linewidth = 0.2
    ) +
    geom_errorbar(
      aes(
        ymin = mean_slope_pp_per_hour -
          sem_slope_pp_per_hour,

        ymax = mean_slope_pp_per_hour +
          sem_slope_pp_per_hour
      ),
      width = 0.16,
      na.rm = TRUE
    ) +
    labs(
      title = paste0(
        current_dimer,
        " interval kinetics"
      ),
      x = NULL,
      y = "Slope (percentage points/hour)"
    ) +
    common_theme +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      ),

      legend.position = "none"
    )

  current_figure <- (
    p_current_raw |
      p_current_change
  ) / p_current_slope +
    plot_layout(
      heights = c(
        1,
        0.8
      )
    )

  ggsave(
    filename = file.path(
      individual_plotdir,
      paste0(
        current_dimer,
        "_kinetics.pdf"
      )
    ),
    plot = current_figure,
    width = 9,
    height = 7,
    units = "in"
  )

  ggsave(
    filename = file.path(
      individual_plotdir,
      paste0(
        current_dimer,
        "_kinetics.png"
      )
    ),
    plot = current_figure,
    width = 9,
    height = 7,
    units = "in",
    dpi = 600
  )
}

# ============================================================
# FINISHED
# ============================================================

cat("\nDone.\n\n")

cat("Main kinetics table:\n")
cat(
  file.path(
    outdir,
    "dipyrimidine_kinetics_all_values.tsv"
  ),
  "\n\n"
)

cat("Interval kinetics:\n")
cat(
  file.path(
    outdir,
    "dipyrimidine_interval_kinetics.tsv"
  ),
  "\n\n"
)

cat("Early/late summary:\n")
cat(
  file.path(
    outdir,
    "dipyrimidine_early_late_summary.tsv"
  ),
  "\n\n"
)

cat("Kinetic feature summary:\n")
cat(
  file.path(
    outdir,
    "dipyrimidine_kinetic_features_summary.tsv"
  ),
  "\n\n"
)

cat("Master figure:\n")
cat(
  file.path(
    plotdir,
    "dipyrimidine_kinetics_master.pdf"
  ),
  "\n\n"
)

cat("Individual dipyrimidine figures:\n")
cat(
  individual_plotdir,
  "\n"
)