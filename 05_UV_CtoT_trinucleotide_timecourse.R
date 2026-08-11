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
  "UV_CtoT_trinucleotide_timecourse"
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

context_order <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

samples <- tribble(
  ~sample,           ~timepoint, ~time_h, ~replicate,
  "NHF1-UVCPD-30m",  "0.5h",        0.5, "R1",
  "NHF1-UVCPD-2h",   "2h",          2.0, "R1",
  "NHF1-UVCPD-4h",   "4h",          4.0, "R1",
  "NHF1-UVCPD-8h",   "8h",          8.0, "R1",
  "NHF1-UVCPD-30m-r2",  "0.5h",        0.5, "R2",
  "NHF1-UVCPD-2h-r2",   "2h",          2.0, "R2",
  "NHF1-UVCPD-4h-r2",   "4h",          4.0, "R2",
  "NHF1-UVCPD-8h-r2",   "8h",          8.0, "R2"
) %>%
  mutate(
    file = file.path(
      base_dir,
      paste0(
        sample,
        "_mismatch_pipeline_CtoT"
      ),
      "02_filtered_events",
      paste0(
        sample,
        "_singleMismatch_C_to_T_6to13nt_from3prime_20to30mers.tsv"
      )
    )
  )

# ============================================================
# CHECK FILES
# ============================================================

missing_files <- samples %>%
  filter(!file.exists(file))

if (nrow(missing_files) > 0) {
  stop(
    paste(
      "Missing input files:",
      paste(
        missing_files$file,
        collapse = "\n"
      ),
      sep = "\n"
    )
  )
}

# ============================================================
# COLUMN DEFINITIONS
# ============================================================

# Actual 16-field data format.
columns_16 <- c(
  "Read_ID",
  "Chromosome",
  "Start",
  "End",
  "Strand",
  "Read_Length",
  "Position_from_5prime",
  "Position_from_3prime",
  "Reference_Base",
  "Change",
  "Trinucleotide_Context",
  "Trinucleotide_Context_RC",
  "Aligned_Sequence",
  "Reference_Sequence",
  "Genomic_Position_0based",
  "Original_Mismatches"
)

# Correct 17-field format, in case another file really contains Alt_Base.
columns_17 <- c(
  "Read_ID",
  "Chromosome",
  "Start",
  "End",
  "Strand",
  "Read_Length",
  "Position_from_5prime",
  "Position_from_3prime",
  "Reference_Base",
  "Alt_Base",
  "Change",
  "Trinucleotide_Context",
  "Trinucleotide_Context_RC",
  "Aligned_Sequence",
  "Reference_Sequence",
  "Genomic_Position_0based",
  "Original_Mismatches"
)

# ============================================================
# READ ONE FILE ROBUSTLY
# ============================================================

read_one <- function(
  sample,
  timepoint,
  time_h,
  replicate,
  file
) {

  message("Reading: ", file)

  # Skip the malformed header and read the actual fields.
  dat <- read_tsv(
    file,
    col_names = FALSE,
    skip = 1,
    show_col_types = FALSE,
    progress = FALSE,
    name_repair = "minimal"
  )

  observed_columns <- ncol(dat)

  message(
    "  Observed fields per row: ",
    observed_columns
  )

  if (observed_columns == 16) {

    colnames(dat) <- columns_16

    dat <- dat %>%
      mutate(
        Alt_Base = str_extract(
          Change,
          "(?<=>)[ACGT]$"
        )
      )

  } else if (observed_columns == 17) {

    colnames(dat) <- columns_17

  } else {

    stop(
      paste0(
        "Unexpected number of fields in ",
        file,
        ": ",
        observed_columns,
        ". Expected 16 or 17."
      )
    )
  }

  dat %>%
    mutate(
      sample = sample,
      timepoint = timepoint,
      time_h = time_h,
      replicate = replicate,
      input_file = file,
      input_row_number = row_number()
    )
}

# ============================================================
# READ ALL SAMPLES
# ============================================================

raw_events <- pmap_dfr(
  samples,
  read_one
)

# ============================================================
# CLEAN AND VALIDATE EVENTS
# ============================================================

classified_events <- raw_events %>%
  mutate(
    Reference_Base = str_to_upper(
      str_trim(
        as.character(Reference_Base)
      )
    ),

    Alt_Base = str_to_upper(
      str_trim(
        as.character(Alt_Base)
      )
    ),

    Change = str_to_upper(
      str_trim(
        as.character(Change)
      )
    ),

    trinucleotide = str_to_upper(
      str_trim(
        as.character(
          Trinucleotide_Context
        )
      )
    ),

    valid_change = (
      Reference_Base == "C" &
      Change == "C>T"
    ),

    valid_context = (
      trinucleotide %in%
        context_order
    ),

    retain_event = (
      valid_change &
      valid_context
    )
  )

# ============================================================
# FILTER SUMMARY
# ============================================================

filter_summary <- classified_events %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  summarise(
    input_rows = n(),

    reference_C_rows = sum(
      Reference_Base == "C",
      na.rm = TRUE
    ),

    change_CtoT_rows = sum(
      Change == "C>T",
      na.rm = TRUE
    ),

    valid_context_rows = sum(
      valid_context,
      na.rm = TRUE
    ),

    retained_rows_before_dedup = sum(
      retain_event,
      na.rm = TRUE
    ),

    .groups = "drop"
  )

# ============================================================
# RETAIN AND DEDUPLICATE EVENTS
# ============================================================

events_before_dedup <- classified_events %>%
  filter(
    retain_event
  )

events <- events_before_dedup %>%
  arrange(
    sample,
    Read_ID,
    Genomic_Position_0based,
    input_row_number
  ) %>%
  distinct(
    sample,
    Read_ID,
    Genomic_Position_0based,
    .keep_all = TRUE
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    trinucleotide = factor(
      trinucleotide,
      levels = context_order
    ),

    damage_class = case_when(
      as.character(trinucleotide) %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ "CT_NCT",

      as.character(trinucleotide) %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ "TC_TCN",

      as.character(trinucleotide) ==
        "TCT" ~ "CT_or_TC_ambiguous",

      TRUE ~ "not_CT_TC_normalizable"
    ),

    damage_context_key = case_when(
      as.character(trinucleotide) %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ paste0(
          "CT:",
          as.character(trinucleotide)
        ),

      as.character(trinucleotide) %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ paste0(
          "TC:",
          as.character(trinucleotide)
        ),

      as.character(trinucleotide) ==
        "TCT" ~ "CT:TCT + TC:TCT",

      TRUE ~ NA_character_
    )
  )

dedup_summary <- events_before_dedup %>%
  count(
    sample,
    timepoint,
    name = "retained_rows_before_dedup_check"
  ) %>%
  left_join(
    events %>%
      count(
        sample,
        timepoint,
        name = "unique_events"
      ),
    by = c(
      "sample",
      "timepoint"
    )
  ) %>%
  mutate(
    duplicate_rows_removed =
      retained_rows_before_dedup_check -
      unique_events
  )

filter_summary <- filter_summary %>%
  left_join(
    dedup_summary,
    by = c(
      "sample",
      "timepoint"
    )
  )

write_tsv(
  filter_summary,
  file.path(
    outdir,
    "UV_CtoT_filter_summary.tsv"
  )
)

# ============================================================
# EVENT-LEVEL OUTPUT
# ============================================================

event_output <- events %>%
  transmute(
    sample,
    timepoint,
    time_h,
    replicate,

    event_id = paste(
      sample,
      Read_ID,
      Genomic_Position_0based,
      sep = "|"
    ),

    Read_ID,
    Chromosome,
    Start,
    End,
    Strand,
    Read_Length,
    Position_from_5prime,
    Position_from_3prime,
    Reference_Base,
    Alt_Base,
    Change,
    trinucleotide,
    trinucleotide_RC =
      Trinucleotide_Context_RC,
    damage_class,
    damage_context_key,
    Aligned_Sequence,
    Reference_Sequence,
    Genomic_Position_0based,
    Original_Mismatches,
    input_file,
    input_row_number
  )

write_tsv(
  event_output,
  file.path(
    outdir,
    "UV_CtoT_trinucleotide_events.tsv"
  )
)

# ============================================================
# COMPLETE COUNT TABLE
# ============================================================

complete_design <- crossing(
  samples %>%
    select(
      sample,
      timepoint,
      time_h,
      replicate
    ),
  trinucleotide = context_order
)

observed_counts <- events %>%
  mutate(
    trinucleotide =
      as.character(trinucleotide),
    timepoint =
      as.character(timepoint)
  ) %>%
  count(
    sample,
    timepoint,
    time_h,
    replicate,
    trinucleotide,
    name = "count"
  )

context_summary <- complete_design %>%
  left_join(
    observed_counts,
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
      0L
    )
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_CtoT_events = sum(
      count
    ),

    percentage = if_else(
      total_CtoT_events > 0,
      100 * count /
        total_CtoT_events,
      NA_real_
    )
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
    ),

    damage_class = case_when(
      as.character(trinucleotide) %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ "CT_NCT",

      as.character(trinucleotide) %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ "TC_TCN",

      as.character(trinucleotide) ==
        "TCT" ~ "CT_or_TC_ambiguous",

      TRUE ~ "not_CT_TC_normalizable"
    ),

    damage_context_key = case_when(
      as.character(trinucleotide) %in%
        c(
          "ACT",
          "CCT",
          "GCT"
        ) ~ paste0(
          "CT:",
          as.character(trinucleotide)
        ),

      as.character(trinucleotide) %in%
        c(
          "TCA",
          "TCC",
          "TCG"
        ) ~ paste0(
          "TC:",
          as.character(trinucleotide)
        ),

      as.character(trinucleotide) ==
        "TCT" ~ "CT:TCT + TC:TCT",

      TRUE ~ NA_character_
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

write_tsv(
  context_summary,
  file.path(
    outdir,
    "UV_CtoT_trinucleotide_counts_percentages.tsv"
  )
)

# ============================================================
# SAMPLE TOTALS
# ============================================================

sample_totals <- context_summary %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate,
    total_CtoT_events
  ) %>%
  arrange(
    time_h
  )

write_tsv(
  sample_totals,
  file.path(
    outdir,
    "UV_CtoT_sample_totals.tsv"
  )
)

# ============================================================
# TIME-COURSE CHANGES
# ============================================================

context_changes <- context_summary %>%
  mutate(
    timepoint =
      as.character(timepoint),
    trinucleotide =
      as.character(trinucleotide)
  ) %>%
  select(
    trinucleotide,
    timepoint,
    percentage
  ) %>%
  pivot_wider(
    names_from = timepoint,
    values_from = percentage,
    names_prefix = "percentage_"
  ) %>%
  mutate(
    change_8h_minus_0.5h =
      percentage_8h -
      `percentage_0.5h`,

    fold_change_8h_vs_0.5h =
      if_else(
        `percentage_0.5h` > 0,
        percentage_8h /
          `percentage_0.5h`,
        NA_real_
      ),

    direction = case_when(
      change_8h_minus_0.5h > 0 ~
        "Increasing",

      change_8h_minus_0.5h < 0 ~
        "Decreasing",

      TRUE ~
        "No change"
    ),

    trinucleotide = factor(
      trinucleotide,
      levels = context_order
    )
  ) %>%
  arrange(
    trinucleotide
  )

write_tsv(
  context_changes,
  file.path(
    outdir,
    "UV_CtoT_trinucleotide_timecourse_changes.tsv"
  )
)

# ============================================================
# PLOT THEME
# ============================================================

common_theme <- theme_classic(
  base_size = 11
) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.text.x = element_text(
      angle = 45,
      hjust = 1
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
# SBS PERCENT PLOT
# ============================================================

p_percent <- ggplot(
  context_summary,
  aes(
    x = trinucleotide,
    y = percentage
  )
) +
  geom_col(
    fill = uv_red,
    width = 0.8,
    color = "black",
    linewidth = 0.15
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 2
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  scale_y_continuous(
    expand = expansion(
      mult = c(
        0,
        0.05
      )
    )
  ) +
  labs(
    title = "UV XR-seq C>T trinucleotide contexts",
    x = NULL,
    y = "C>T events (%)"
  ) +
  common_theme

ggsave(
  filename = file.path(
    plotdir,
    "DDB2_UV_CtoT_SBS_timecourse_percent.pdf"
  ),
  plot = p_percent,
  width = 7.2,
  height = 6.2,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "DDB2_UV_CtoT_SBS_timecourse_percent.png"
  ),
  plot = p_percent,
  width = 7.2,
  height = 6.2,
  units = "in",
  dpi = 600
)

# ============================================================
# CONTEXT TRAJECTORIES
# ============================================================

p_trajectory <- ggplot(
  context_summary,
  aes(
    x = time_h,
    y = percentage
  )
) +
  geom_line(
    color = uv_red,
    linewidth = 0.7
  ) +
  geom_point(
    color = uv_red,
    size = 1.8
  ) +
  facet_wrap(
    ~trinucleotide,
    ncol = 4
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
    title = "UV C>T context trajectories",
    x = "Repair time (h)",
    y = "C>T events (%)"
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
    "CSB_UV_CtoT_context_trajectories.pdf"
  ),
  plot = p_trajectory,
  width = 7.4,
  height = 7.2,
  units = "in"
)

# ============================================================
# CHANGE PLOT
# ============================================================

direction_colors <- c(
  "Increasing" = "#B2182B",
  "Decreasing" = "#2166AC",
  "No change" = "grey70"
)

p_change <- ggplot(
  context_changes,
  aes(
    x = trinucleotide,
    y = change_8h_minus_0.5h,
    fill = direction
  )
) +
  geom_hline(
    yintercept = 0,
    linewidth = 0.4
  ) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.15
  ) +
  scale_fill_manual(
    values = direction_colors,
    name = NULL
  ) +
  scale_x_discrete(
    drop = FALSE
  ) +
  labs(
    title = "Change in UV C>T context frequency",
    subtitle = "8 h minus 0.5 h",
    x = NULL,
    y = "Change in percentage points"
  ) +
  common_theme +
  theme(
    legend.position = "top"
  )

ggsave(
  filename = file.path(
    plotdir,
    "CSB_UV_CtoT_change_8h_vs_0.5h.pdf"
  ),
  plot = p_change,
  width = 7.2,
  height = 4.5,
  units = "in"
)

cat("\nDone.\n\n")

cat("Sample totals:\n")
print(
  sample_totals,
  n = Inf
)

cat("\nFilter summary:\n")
print(
  filter_summary,
  n = Inf
)

cat("\nOutput directory:\n")
cat(outdir, "\n")