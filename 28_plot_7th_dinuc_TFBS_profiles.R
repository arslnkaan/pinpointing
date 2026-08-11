#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
})

# ============================================================
# SETTINGS
# ============================================================

base_dir <- "/work/users/a/r/arslank"

indir <- file.path(
  base_dir,
  "UV_TFBS_7th_dinucleotide_repair",
  "motif_center_profiles"
)

infile <- file.path(
  indir,
  "TFBS_7th_dinucleotide_minus10_plus10_profiles.tsv"
)

outdir <- file.path(
  indir,
  "plots"
)

individual_dir <- file.path(
  outdir,
  "individual_TFs"
)

dir.create(
  individual_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

dinucleotide_order <- c(
  "CC",
  "CT",
  "TC",
  "TT"
)

timepoint_colors <- c(
  "0.5h" = "#C62828",
  "2h"   = "#E69F00",
  "4h"   = "#009E73",
  "8h"   = "#4F759B"
)

dinucleotide_colors <- c(
  "CC" = "#D98C00",
  "CT" = "#C62828",
  "TC" = "#8E24AA",
  "TT" = "#4F759B"
)

# Plot every TF present in active_TFBS.bed.
minimum_sites <- 1

# TRUE adds the Damage-seq profile below the repair profiles.
show_damage_panel <- TRUE

# FALSE avoids producing many large PNG files.
save_individual_png <- FALSE

# Dimer-specific repair profiles should normally use
# all-library normalization so CC, CT, TC and TT abundance
# remains directly comparable.
#
# Options:
#   "all_library"
#   "within_dinucleotide"
dimer_normalization <- "all_library"

# Use fixed scales to preserve magnitude differences between
# the four dipyrimidines.
#
# Options:
#   "fixed"
#   "free_y"
dimer_facet_scales <- "fixed"

# ============================================================
# CHECK INPUT
# ============================================================

if (!file.exists(infile)) {
  stop(
    paste0(
      "Missing input file:\n",
      infile,
      "\n\nRun 27_build_7th_dinuc_TFBS_profiles.py first."
    )
  )
}

# ============================================================
# HELPERS
# ============================================================

sanitize_filename <- function(x) {
  x %>%
    str_replace_all(
      "[^A-Za-z0-9._-]+",
      "_"
    ) %>%
    str_replace_all(
      "^_+|_+$",
      ""
    )
}

position_scale <- function() {
  scale_x_continuous(
    limits = c(-10, 10),
    breaks = seq(-10, 10, by = 5),
    minor_breaks = seq(-10, 10, by = 1),
    expand = expansion(mult = c(0, 0))
  )
}

make_motif_band <- function(median_width) {
  half_width <- median_width / 2

  tibble(
    xmin = max(-10, -half_width),
    xmax = min(10, half_width),
    ymin = -Inf,
    ymax = Inf
  )
}

# ============================================================
# READ PROFILE TABLE
# ============================================================

profiles <- read_tsv(
  infile,
  show_col_types = FALSE
) %>%
  transmute(
    TF = as.character(TF),

    n_sites_total = as.integer(
      n_sites_total
    ),

    n_sites_used = as.integer(
      n_sites_used
    ),

    median_site_width = as.numeric(
      median_site_width
    ),

    data_type = as.character(
      data_type
    ),

    series = as.character(
      series
    ),

    timepoint = as.character(
      timepoint
    ),

    relative_position = as.integer(
      relative_position
    ),

    event_count = as.numeric(
      event_count
    ),

    total_all_events = as.numeric(
      total_all_events
    ),

    total_series_events = as.numeric(
      total_series_events
    ),

    events_per_site = as.numeric(
      events_per_site
    ),

    signal_all_library = as.numeric(
      RPM_per_1000_sites_all_library
    ),

    signal_within_series = as.numeric(
      RPM_per_1000_sites_within_series
    )
  ) %>%
  filter(
    relative_position >= -10,
    relative_position <= 10,
    !is.na(TF),
    TF != ""
  )

if (nrow(profiles) == 0) {
  stop(
    "The profile input contains no usable rows."
  )
}

# ============================================================
# TF SUMMARY
# ============================================================

tf_summary <- profiles %>%
  distinct(
    TF,
    n_sites_total,
    n_sites_used,
    median_site_width
  ) %>%
  filter(
    n_sites_used >= minimum_sites
  ) %>%
  arrange(
    desc(n_sites_used),
    TF
  )

if (nrow(tf_summary) == 0) {
  stop(
    "No TFs passed minimum_sites."
  )
}

# ============================================================
# COMBINED REPAIR PROFILE
#
# CC + CT + TC + TT
#
# Normalized to all selected dipyrimidine repair events at
# each time point and to the number of TFBS loci.
# ============================================================

repair_combined <- profiles %>%
  filter(
    data_type == "repair",
    series == "Combined",
    timepoint %in% timepoint_order,
    TF %in% tf_summary$TF
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    )
  ) %>%
  arrange(
    TF,
    timepoint,
    relative_position
  )

# ============================================================
# REPAIR BY DIPYRIMIDINE
# ============================================================

repair_by_dimer <- profiles %>%
  filter(
    data_type == "repair",
    series %in% dinucleotide_order,
    timepoint %in% timepoint_order,
    TF %in% tf_summary$TF
  ) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),

    series = factor(
      series,
      levels = dinucleotide_order
    ),

    dimer_signal = case_when(
      dimer_normalization == "all_library" ~
        signal_all_library,

      dimer_normalization == "within_dinucleotide" ~
        signal_within_series,

      TRUE ~ NA_real_
    )
  ) %>%
  arrange(
    TF,
    series,
    timepoint,
    relative_position
  )

if (all(is.na(repair_by_dimer$dimer_signal))) {
  stop(
    paste0(
      "Invalid dimer_normalization setting: ",
      dimer_normalization
    )
  )
}

# ============================================================
# DAMAGE-SEQ PROFILES
# ============================================================

damage_combined <- profiles %>%
  filter(
    data_type == "damage",
    series == "Combined",
    TF %in% tf_summary$TF
  ) %>%
  arrange(
    TF,
    relative_position
  )

damage_by_dimer <- profiles %>%
  filter(
    data_type == "damage",
    series %in% dinucleotide_order,
    TF %in% tf_summary$TF
  ) %>%
  mutate(
    series = factor(
      series,
      levels = dinucleotide_order
    )
  ) %>%
  arrange(
    TF,
    series,
    relative_position
  )

# ============================================================
# WRITE PLOTTING TABLES
# ============================================================

write_tsv(
  repair_combined,
  file.path(
    indir,
    "TFBS_combined_repair_profiles.tsv"
  )
)

write_tsv(
  repair_by_dimer,
  file.path(
    indir,
    "TFBS_repair_profiles_by_dinucleotide.tsv"
  )
)

write_tsv(
  damage_combined,
  file.path(
    indir,
    "TFBS_combined_damage_profiles.tsv"
  )
)

write_tsv(
  damage_by_dimer,
  file.path(
    indir,
    "TFBS_damage_profiles_by_dinucleotide.tsv"
  )
)

write_tsv(
  tf_summary,
  file.path(
    indir,
    "TFBS_profile_plotting_summary.tsv"
  )
)

# ============================================================
# PLOT THEME
# ============================================================

profile_theme <- theme_classic(
  base_size = 11
) +
  theme(
    axis.text = element_text(
      color = "black"
    ),

    axis.title = element_text(
      color = "black"
    ),

    axis.ticks = element_line(
      linewidth = 0.3
    ),

    plot.title = element_text(
      face = "bold",
      hjust = 0.5
    ),

    plot.subtitle = element_text(
      hjust = 0.5
    ),

    strip.background = element_rect(
      fill = "grey95",
      color = "grey70",
      linewidth = 0.4
    ),

    strip.text = element_text(
      face = "bold"
    ),

    panel.border = element_rect(
      fill = NA,
      color = "grey65",
      linewidth = 0.35
    ),

    legend.position = "top"
  )

# ============================================================
# FUNCTION: ONE TF FIGURE
# ============================================================

make_tf_plot <- function(current_tf) {

  current_meta <- tf_summary %>%
    filter(
      TF == current_tf
    ) %>%
    slice(1)

  if (nrow(current_meta) != 1) {
    stop(
      paste0(
        "Could not retrieve TF metadata for ",
        current_tf
      )
    )
  }

  motif_band <- make_motif_band(
    current_meta$median_site_width
  )

  current_repair_combined <- repair_combined %>%
    filter(
      TF == current_tf
    )

  current_repair_dimers <- repair_by_dimer %>%
    filter(
      TF == current_tf
    )

  current_damage_combined <- damage_combined %>%
    filter(
      TF == current_tf
    )

  current_damage_dimers <- damage_by_dimer %>%
    filter(
      TF == current_tf
    )

  # ----------------------------------------------------------
  # PANEL 1: COMBINED REPAIR
  # ----------------------------------------------------------

  p_repair_combined <- ggplot(
    current_repair_combined,
    aes(
      x = relative_position,
      y = signal_all_library,
      color = timepoint,
      group = timepoint
    )
  ) +
    geom_rect(
      data = motif_band,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "grey92",
      color = NA
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.4,
      color = "grey30"
    ) +
    geom_line(
      linewidth = 0.9,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.7,
      na.rm = TRUE
    ) +
    scale_color_manual(
      values = timepoint_colors,
      breaks = timepoint_order,
      name = "Repair time",
      drop = FALSE
    ) +
    position_scale() +
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    labs(
      title = "Combined repair profile",
      subtitle = paste0(
        "CC + CT + TC + TT at the seventh ",
        "dinucleotide from the 3′ end"
      ),
      x = "Position relative to motif center (bp)",
      y = paste0(
        "Repair events per million\n",
        "per 1,000 TFBS loci"
      )
    ) +
    profile_theme

  # ----------------------------------------------------------
  # PANEL 2: REPAIR BY DINUCLEOTIDE
  # ----------------------------------------------------------

  dimer_subtitle <- case_when(
    dimer_normalization == "all_library" ~
      paste0(
        "Normalized to all dipyrimidine repair events ",
        "at each time point"
      ),

    dimer_normalization == "within_dinucleotide" ~
      paste0(
        "Each dipyrimidine normalized to its own ",
        "library total"
      ),

    TRUE ~ ""
  )

  p_repair_dimers <- ggplot(
    current_repair_dimers,
    aes(
      x = relative_position,
      y = dimer_signal,
      color = timepoint,
      group = timepoint
    )
  ) +
    geom_rect(
      data = motif_band,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "grey92",
      color = NA
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey30"
    ) +
    geom_line(
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.2,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~series,
      ncol = 2,
      scales = dimer_facet_scales
    ) +
    scale_color_manual(
      values = timepoint_colors,
      breaks = timepoint_order,
      name = "Repair time",
      drop = FALSE
    ) +
    position_scale() +
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    labs(
      title = "Repair profiles by dipyrimidine",
      subtitle = dimer_subtitle,
      x = "Position relative to motif center (bp)",
      y = paste0(
        "Repair events per million\n",
        "per 1,000 TFBS loci"
      )
    ) +
    profile_theme

  # ----------------------------------------------------------
  # PANEL 3: DAMAGE-SEQ
  # ----------------------------------------------------------

  p_damage <- ggplot(
    current_damage_dimers,
    aes(
      x = relative_position,
      y = signal_all_library,
      color = series,
      group = series
    )
  ) +
    geom_rect(
      data = motif_band,
      aes(
        xmin = xmin,
        xmax = xmax,
        ymin = ymin,
        ymax = ymax
      ),
      inherit.aes = FALSE,
      fill = "grey92",
      color = NA
    ) +
    geom_vline(
      xintercept = 0,
      linetype = "dashed",
      linewidth = 0.35,
      color = "grey30"
    ) +
    geom_line(
      linewidth = 0.75,
      na.rm = TRUE
    ) +
    geom_point(
      size = 1.2,
      na.rm = TRUE
    ) +
    facet_wrap(
      ~series,
      ncol = 4,
      scales = dimer_facet_scales
    ) +
    scale_color_manual(
      values = dinucleotide_colors,
      breaks = dinucleotide_order,
      guide = "none",
      drop = FALSE
    ) +
    position_scale() +
    scale_y_continuous(
      expand = expansion(
        mult = c(0, 0.08)
      )
    ) +
    labs(
      title = "Initial Damage-seq profile",
      subtitle = "CC, CT, TC and TT at Damage-seq positions −2/−1",
      x = "Position relative to motif center (bp)",
      y = paste0(
        "Damage events per million\n",
        "per 1,000 TFBS loci"
      )
    ) +
    profile_theme +
    theme(
      legend.position = "none"
    )

  figure_title <- paste0(
    current_tf,
    " motif-centered repair profile"
  )

  figure_subtitle <- paste0(
    format(
      current_meta$n_sites_used,
      big.mark = ","
    ),
    " binding loci; median motif interval width = ",
    round(
      current_meta$median_site_width,
      1
    ),
    " bp; −10 to +10 bp; 1-bp resolution; no smoothing"
  )

  if (
    show_damage_panel &&
      nrow(current_damage_dimers) > 0
  ) {
    combined_figure <- (
      p_repair_combined /
        p_repair_dimers /
        p_damage
    ) +
      plot_layout(
        heights = c(
          0.9,
          1.7,
          0.9
        ),
        guides = "collect"
      )
  } else {
    combined_figure <- (
      p_repair_combined /
        p_repair_dimers
    ) +
      plot_layout(
        heights = c(
          0.9,
          1.7
        ),
        guides = "collect"
      )
  }

  combined_figure +
    plot_annotation(
      title = figure_title,
      subtitle = figure_subtitle,
      theme = theme(
        plot.title = element_text(
          size = 16,
          face = "bold",
          hjust = 0.5
        ),

        plot.subtitle = element_text(
          size = 10,
          hjust = 0.5
        )
      )
    ) &
    theme(
      legend.position = "top"
    )
}

# ============================================================
# TF ORDER
# ============================================================

tf_order <- tf_summary %>%
  arrange(
    desc(n_sites_used),
    TF
  ) %>%
  pull(TF)

# ============================================================
# MULTIPAGE PDF AND INDIVIDUAL FILES
# ============================================================

multipage_pdf <- file.path(
  outdir,
  "all_TFs_7th_dinucleotide_minus10_plus10_profiles.pdf"
)

if (show_damage_panel) {
  figure_height <- 12
} else {
  figure_height <- 9
}

pdf(
  file = multipage_pdf,
  width = 10,
  height = figure_height,
  onefile = TRUE
)

manifest <- vector(
  mode = "list",
  length = length(tf_order)
)

for (tf_index in seq_along(tf_order)) {

  current_tf <- tf_order[
    tf_index
  ]

  message(
    "[",
    tf_index,
    "/",
    length(tf_order),
    "] Plotting ",
    current_tf
  )

  current_plot <- make_tf_plot(
    current_tf
  )

  print(
    current_plot
  )

  safe_tf <- sanitize_filename(
    current_tf
  )

  individual_pdf <- file.path(
    individual_dir,
    paste0(
      safe_tf,
      "_7th_dinucleotide_minus10_plus10.pdf"
    )
  )

  ggsave(
    filename = individual_pdf,
    plot = current_plot,
    width = 10,
    height = figure_height,
    units = "in"
  )

  individual_png <- NA_character_

  if (save_individual_png) {

    individual_png <- file.path(
      individual_dir,
      paste0(
        safe_tf,
        "_7th_dinucleotide_minus10_plus10.png"
      )
    )

    ggsave(
      filename = individual_png,
      plot = current_plot,
      width = 10,
      height = figure_height,
      units = "in",
      dpi = 500
    )
  }

  current_n_sites <- tf_summary %>%
    filter(
      TF == current_tf
    ) %>%
    pull(n_sites_used)

  manifest[[tf_index]] <- tibble(
    TF = current_tf,
    n_sites_used = current_n_sites[1],
    individual_pdf = individual_pdf,
    individual_png = individual_png
  )
}

dev.off()

plot_manifest <- bind_rows(
  manifest
)

write_tsv(
  plot_manifest,
  file.path(
    indir,
    "TFBS_7th_dinucleotide_plot_manifest.tsv"
  )
)

# ============================================================
# STANDALONE CTCF FIGURE
# ============================================================

if ("CTCF" %in% tf_order) {

  ctcf_plot <- make_tf_plot(
    "CTCF"
  )

  ggsave(
    filename = file.path(
      outdir,
      "CTCF_7th_dinucleotide_minus10_plus10.pdf"
    ),
    plot = ctcf_plot,
    width = 10,
    height = figure_height,
    units = "in"
  )

  ggsave(
    filename = file.path(
      outdir,
      "CTCF_7th_dinucleotide_minus10_plus10.png"
    ),
    plot = ctcf_plot,
    width = 10,
    height = figure_height,
    units = "in",
    dpi = 600
  )
}

# ============================================================
# FINISHED
# ============================================================

cat("\nDone.\n\n")

cat("Multipage TF profile PDF:\n")
cat(
  multipage_pdf,
  "\n\n"
)

cat("Individual TF profiles:\n")
cat(
  individual_dir,
  "\n\n"
)

if ("CTCF" %in% tf_order) {
  cat("CTCF profile:\n")
  cat(
    file.path(
      outdir,
      "CTCF_7th_dinucleotide_minus10_plus10.pdf"
    ),
    "\n"
  )
}