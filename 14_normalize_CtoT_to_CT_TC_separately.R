#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
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
  "UV_CtoT_CT_TC_separate_normalization"
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

timepoint_order <- c(
  "0.5h",
  "2h",
  "4h",
  "8h"
)

sbs_order <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

ct_contexts <- c(
  "ACT",
  "CCT",
  "GCT"
)

tc_contexts <- c(
  "TCA",
  "TCC",
  "TCG"
)

ambiguous_context <- "TCT"

ct_color <- "#C62828"
tc_color <- "#8E24AA"
ambiguous_color <- "#4F759B"

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

# ============================================================
# CHECK INPUTS
# ============================================================

for (file in c(repair_file, damage_file)) {
  if (!file.exists(file)) {
    stop(
      paste0(
        "Missing input: ",
        file
      )
    )
  }
}

# ============================================================
# READ REPAIR
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
  )

sample_design <- repair_raw %>%
  distinct(
    sample,
    timepoint,
    time_h,
    replicate
  )

# ============================================================
# READ DAMAGE
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
    damage_count = as.numeric(count),
    lesion_total = as.numeric(lesion_total),
    total_CT_TC_damage = as.numeric(
      total_CT_TC_damage
    )
  )

# ============================================================
# UNAMBIGUOUS CT NORMALIZATION
#
# Repair denominator:
# ACT + CCT + GCT
#
# Damage denominator:
# CT:ACT + CT:CCT + CT:GCT
#
# TCT is excluded here because repair orientation is ambiguous.
# ============================================================

ct_repair <- repair_raw %>%
  filter(
    trinucleotide %in%
      ct_contexts
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_CT_repair = sum(
      repair_count,
      na.rm = TRUE
    ),
    repair_fraction = if_else(
      total_CT_repair > 0,
      repair_count / total_CT_repair,
      NA_real_
    )
  ) %>%
  ungroup()

ct_damage <- damage_raw %>%
  filter(
    lesion_type == "CT",
    trinucleotide %in%
      ct_contexts
  ) %>%
  mutate(
    total_CT_damage_used = sum(
      damage_count,
      na.rm = TRUE
    ),
    damage_fraction = if_else(
      total_CT_damage_used > 0,
      damage_count /
        total_CT_damage_used,
      NA_real_
    )
  )

ct_normalized <- ct_repair %>%
  left_join(
    ct_damage %>%
      select(
        trinucleotide,
        damage_count,
        total_CT_damage_used,
        damage_fraction
      ),
    by = "trinucleotide"
  ) %>%
  mutate(
    lesion_class = "CT",
    repair_over_damage = safe_ratio(
      repair_fraction,
      damage_fraction
    ),
    log2_repair_over_damage = safe_log2(
      repair_over_damage
    ),
    log2_for_plot = finite_for_plot(
      log2_repair_over_damage
    )
  )

# ============================================================
# UNAMBIGUOUS TC NORMALIZATION
#
# Repair denominator:
# TCA + TCC + TCG
#
# Damage denominator:
# TC:TCA + TC:TCC + TC:TCG
# ============================================================

tc_repair <- repair_raw %>%
  filter(
    trinucleotide %in%
      tc_contexts
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_TC_repair = sum(
      repair_count,
      na.rm = TRUE
    ),
    repair_fraction = if_else(
      total_TC_repair > 0,
      repair_count / total_TC_repair,
      NA_real_
    )
  ) %>%
  ungroup()

tc_damage <- damage_raw %>%
  filter(
    lesion_type == "TC",
    trinucleotide %in%
      tc_contexts
  ) %>%
  mutate(
    total_TC_damage_used = sum(
      damage_count,
      na.rm = TRUE
    ),
    damage_fraction = if_else(
      total_TC_damage_used > 0,
      damage_count /
        total_TC_damage_used,
      NA_real_
    )
  )

tc_normalized <- tc_repair %>%
  left_join(
    tc_damage %>%
      select(
        trinucleotide,
        damage_count,
        total_TC_damage_used,
        damage_fraction
      ),
    by = "trinucleotide"
  ) %>%
  mutate(
    lesion_class = "TC",
    repair_over_damage = safe_ratio(
      repair_fraction,
      damage_fraction
    ),
    log2_repair_over_damage = safe_log2(
      repair_over_damage
    ),
    log2_for_plot = finite_for_plot(
      log2_repair_over_damage
    )
  )

# ============================================================
# AMBIGUOUS TCT
#
# Repair:
# all TCT events
#
# Damage:
# CT:TCT + TC:TCT
#
# This is shown separately and is not assigned to CT or TC.
# ============================================================

tct_repair <- repair_raw %>%
  filter(
    trinucleotide ==
      ambiguous_context
  ) %>%
  group_by(
    sample,
    timepoint,
    time_h,
    replicate
  ) %>%
  mutate(
    total_repair_all_CPD_contexts = sum(
      repair_raw$repair_count[
        repair_raw$sample == first(sample) &
        repair_raw$timepoint == first(timepoint) &
        repair_raw$trinucleotide %in%
          c(
            ct_contexts,
            tc_contexts,
            ambiguous_context
          )
      ],
      na.rm = TRUE
    ),
    repair_fraction = if_else(
      total_repair_all_CPD_contexts > 0,
      repair_count /
        total_repair_all_CPD_contexts,
      NA_real_
    )
  ) %>%
  ungroup()

tct_damage_count <- damage_raw %>%
  filter(
    trinucleotide == "TCT",
    lesion_type %in%
      c(
        "CT",
        "TC"
      )
  ) %>%
  summarise(
    damage_count = sum(
      damage_count,
      na.rm = TRUE
    ),
    total_damage = first(
      total_CT_TC_damage
    )
  ) %>%
  mutate(
    damage_fraction = if_else(
      total_damage > 0,
      damage_count /
        total_damage,
      NA_real_
    )
  )

tct_normalized <- tct_repair %>%
  mutate(
    damage_count =
      tct_damage_count$damage_count,
    damage_fraction =
      tct_damage_count$damage_fraction,
    lesion_class = "CT/TC ambiguous",
    repair_over_damage = safe_ratio(
      repair_fraction,
      damage_fraction
    ),
    log2_repair_over_damage = safe_log2(
      repair_over_damage
    ),
    log2_for_plot = finite_for_plot(
      log2_repair_over_damage
    )
  )

# ============================================================
# COMBINE
# ============================================================

normalized_used <- bind_rows(
  ct_normalized,
  tc_normalized,
  tct_normalized
) %>%
  mutate(
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),
    trinucleotide = factor(
      trinucleotide,
      levels = sbs_order
    ),
    lesion_class = factor(
      lesion_class,
      levels = c(
        "CT",
        "TC",
        "CT/TC ambiguous"
      )
    )
  ) %>%
  arrange(
    time_h,
    trinucleotide
  )

write_tsv(
  normalized_used,
  file.path(
    outdir,
    "UV_CtoT_CT_TC_separate_repair_over_damage.tsv"
  )
)

# ============================================================
# FULL SBS PLOT TABLE
#
# Keep all SBS positions visible.
# Non-normalizable positions remain NA.
# ============================================================

plot_table <- crossing(
  sample_design,
  trinucleotide = sbs_order
) %>%
  left_join(
    normalized_used %>%
      mutate(
        timepoint =
          as.character(timepoint),
        trinucleotide =
          as.character(trinucleotide)
      ) %>%
      select(
        sample,
        timepoint,
        time_h,
        replicate,
        trinucleotide,
        lesion_class,
        repair_count,
        damage_count,
        repair_fraction,
        damage_fraction,
        repair_over_damage,
        log2_repair_over_damage,
        log2_for_plot
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
    timepoint = factor(
      timepoint,
      levels = timepoint_order
    ),
    trinucleotide = factor(
      trinucleotide,
      levels = sbs_order
    )
  )

write_tsv(
  plot_table,
  file.path(
    outdir,
    "UV_CtoT_CT_TC_separate_full_SBS_table.tsv"
  )
)

# ============================================================
# COLORS
# ============================================================

class_colors <- c(
  "CT" = ct_color,
  "TC" = tc_color,
  "CT/TC ambiguous" =
    ambiguous_color
)

# ============================================================
# PLOT 1
# FULL SBS ORDER
# ============================================================

p_sbs <- ggplot(
  plot_table,
  aes(
    x = trinucleotide,
    y = log2_for_plot,
    fill = lesion_class
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.4,
    color = "grey40"
  ) +
  geom_col(
    width = 0.8,
    color = "black",
    linewidth = 0.15,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~timepoint,
    ncol = 2
  ) +
  scale_fill_manual(
    values = class_colors,
    name = NULL,
    na.value = "transparent"
  ) +
  scale_x_discrete(
    limits = sbs_order,
    drop = FALSE
  ) +
  labs(
    title = "C>T repair normalized to CT and TC damage separately",
    subtitle = paste(
      "CT and TC denominators are calculated independently;",
      "TCT is shown as orientation-ambiguous"
    ),
    x = NULL,
    y = expression(
      log[2](
        repair~fraction /
          damage~fraction
      )
    )
  ) +
  theme_classic(
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
    legend.position = "top",
    strip.background = element_rect(
      fill = "grey95",
      color = "grey60"
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
    ),
    plot.subtitle = element_text(
      hjust = 0.5
    )
  )

ggsave(
  filename = file.path(
    plotdir,
    "UV_CtoT_CT_TC_separate_log2_repair_over_damage_SBS.pdf"
  ),
  plot = p_sbs,
  width = 8.2,
  height = 6.4,
  units = "in"
)

ggsave(
  filename = file.path(
    plotdir,
    "UV_CtoT_CT_TC_separate_log2_repair_over_damage_SBS.png"
  ),
  plot = p_sbs,
  width = 8.2,
  height = 6.4,
  units = "in",
  dpi = 600
)

# ============================================================
# PLOT 2
# TRAJECTORIES
# ============================================================

p_trajectory <- ggplot(
  normalized_used,
  aes(
    x = time_h,
    y = log2_for_plot,
    color = lesion_class,
    group = trinucleotide
  )
) +
  geom_hline(
    yintercept = 0,
    linetype = "dashed",
    linewidth = 0.35,
    color = "grey50"
  ) +
  geom_line(
    linewidth = 0.75,
    na.rm = TRUE
  ) +
  geom_point(
    size = 1.8,
    na.rm = TRUE
  ) +
  facet_wrap(
    ~trinucleotide,
    ncol = 4,
    drop = FALSE
  ) +
  scale_color_manual(
    values = class_colors,
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
  labs(
    title = "CT- and TC-specific repair trajectories",
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
    legend.position = "top",
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
    "UV_CtoT_CT_TC_separate_log2_trajectories.pdf"
  ),
  plot = p_trajectory,
  width = 7.4,
  height = 5.8,
  units = "in"
)

cat("\nDone.\n\n")
cat("Normalized table:\n")
cat(
  file.path(
    outdir,
    "UV_CtoT_CT_TC_separate_repair_over_damage.tsv"
  ),
  "\n\n"
)
cat("Plots:\n")
cat(plotdir, "\n")