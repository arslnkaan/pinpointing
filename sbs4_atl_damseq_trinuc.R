#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
  library(patchwork)
  library(cowplot)
})

outdir <- "ATL_SBS4_DamageSeq_percent_and_repair_over_damage_compact"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

sbs4_file <- "SBS4_GRCh38_CtoA_only_renormalized_100_percent.csv"
damage_percent_file <- "NHF1-4NQO-30m-r1_minus1G_trinuc_RC_SBS4_percentages.tsv"

damage_files <- tribble(
  ~replicate, ~file,
  "R1", "damseq/NHF1-4NQO-30m-r1-results/damage_sites/NHF1-4NQO-30m-r1.fragment_oriented_minus1G_trinuc_percentages.tsv",
  "R2", "damseq/NHF1-4NQO-30m-r2-results/damage_sites/NHF1-4NQO-30m-r2.fragment_oriented_minus1G_trinuc_percentages.tsv",
  "R3", "damseq/NHF1-4NQO-30m-r3-results/damage_sites/NHF1-4NQO-30m-r3.fragment_oriented_minus1G_trinuc_percentages.tsv"
)

atl_files <- tribble(
  ~timepoint, ~replicate, ~file,

  "0.5h", "R1",
  "/work/users/a/r/arslank/NHF1-4NQO-30m-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-30m-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",
  "0.5h", "R2",
  "/work/users/a/r/arslank/NHF1-4NQO-30m-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-30m-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "2h", "R1",
  "/work/users/a/r/arslank/NHF1-4NQO-2h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-2h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",
  "2h", "R2",
  "/work/users/a/r/arslank/NHF1-4NQO-2h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-2h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "4h", "R1",
  "/work/users/a/r/arslank/NHF1-4NQO-4h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-4h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",
  "4h", "R2",
  "/work/users/a/r/arslank/NHF1-4NQO-4h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-4h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",

  "8h", "R1",
  "/work/users/a/r/arslank/NHF1-4NQO-8h-r1_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-8h-r1_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv",
  "8h", "R2",
  "/work/users/a/r/arslank/NHF1-4NQO-8h-r2_R1_001_mismatch_pipeline/02_filtered_events/NHF1-4NQO-8h-r2_R1_001_singleMismatch_G_to_T_6to13nt_from3prime_20to30mers_trinucleotide_percentages.csv"
)

CONTEXT_ORDER <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

time_order <- c("0.5h", "2h", "4h", "8h")
pseudo <- 1e-6

revcomp <- function(x) {
  sapply(x, function(seq) {
    paste(rev(strsplit(chartr("ACGT", "TGCA", seq), "")[[1]]), collapse = "")
  }, USE.NAMES = FALSE)
}

cosine_similarity <- function(a, b) {
  denom <- sqrt(sum(a^2, na.rm = TRUE) * sum(b^2, na.rm = TRUE))
  if (!is.finite(denom) || denom == 0) return(NA_real_)
  sum(a * b, na.rm = TRUE) / denom
}

read_atl <- function(timepoint, replicate, file) {
  if (!file.exists(file)) stop("Missing ATL file: ", file)

  df <- read_csv(file, show_col_types = FALSE)

  context_col <- intersect(
    c("Trinucleotide_Context_RC", "Trinucleotide_Context",
      "trinucleotide", "Trinucleotide", "context", "Context"),
    colnames(df)
  )[1]

  pct_col <- intersect(
    c("Percentage", "Percent", "percentage", "percent"),
    colnames(df)
  )[1]

  if (is.na(context_col)) stop("No context column in: ", file)
  if (is.na(pct_col)) stop("No percentage column in: ", file)

  df %>%
    transmute(
      timepoint = timepoint,
      replicate = replicate,
      context = toupper(as.character(.data[[context_col]])),
      ATL_percent = as.numeric(.data[[pct_col]])
    ) %>%
    filter(context %in% CONTEXT_ORDER)
}

read_sbs4 <- function(file) {
  if (!file.exists(file)) stop("Missing SBS4 file: ", file)

  df <- read_csv(file, show_col_types = FALSE)

  context_col <- intersect(
    c("Trinucleotide_Context", "Trinucleotide_Context_RC",
      "trinucleotide", "Trinucleotide", "context", "Context"),
    colnames(df)
  )[1]

  pct_col <- intersect(
    c("Percentage", "Percent", "percentage", "percent"),
    colnames(df)
  )[1]

  if (is.na(context_col)) stop("No context column in SBS4 file: ", file)
  if (is.na(pct_col)) stop("No percent column in SBS4 file: ", file)

  df %>%
    transmute(
      context = toupper(as.character(.data[[context_col]])),
      SBS4_percent = as.numeric(.data[[pct_col]])
    ) %>%
    filter(context %in% CONTEXT_ORDER) %>%
    complete(context = CONTEXT_ORDER, fill = list(SBS4_percent = 0)) %>%
    mutate(
      SBS4_percent = 100 * SBS4_percent / sum(SBS4_percent, na.rm = TRUE),
      context = factor(context, levels = CONTEXT_ORDER)
    )
}

read_damage_context_raw <- function(file) {
  if (!file.exists(file)) stop("Missing Damage-seq file: ", file)

  raw_df <- read_tsv(file, show_col_types = FALSE, col_names = FALSE)

  raw_df %>%
    select(1:3) %>%
    setNames(c("Trinucleotide_raw", "Count", "Percentage")) %>%
    filter(Trinucleotide_raw != "Trinucleotide_raw") %>%
    mutate(
      Trinucleotide_raw = toupper(as.character(Trinucleotide_raw)),
      Count = as.numeric(Count)
    ) %>%
    filter(str_detect(Trinucleotide_raw, "^[ACGT]{3}$")) %>%
    mutate(context = revcomp(Trinucleotide_raw)) %>%
    group_by(context) %>%
    summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
    mutate(Damage_percent = 100 * Count / sum(Count, na.rm = TRUE)) %>%
    complete(context = CONTEXT_ORDER, fill = list(Count = 0, Damage_percent = 0)) %>%
    mutate(context = factor(context, levels = CONTEXT_ORDER))
}

read_damage_one <- function(replicate, file) {
  read_damage_context_raw(file) %>%
    mutate(replicate = replicate)
}

read_damage_for_normalization <- function(file) {
  if (!file.exists(file)) stop("Missing Damage-seq normalization file: ", file)

  df <- read_tsv(file, show_col_types = FALSE)

  context_col <- intersect(
    c("RC", "context", "Context", "Trinucleotide", "Trinucleotide_RC",
      "Trinucleotide_Context", "Trinucleotide_Context_RC"),
    colnames(df)
  )[1]

  pct_col <- intersect(
    c("Percentage", "Percent", "percentage", "percent"),
    colnames(df)
  )[1]

  count_col <- intersect(c("Count", "count"), colnames(df))[1]

  if (is.na(context_col)) context_col <- colnames(df)[1]

  if (!is.na(pct_col)) {
    out <- df %>%
      transmute(
        context = toupper(as.character(.data[[context_col]])),
        Damage_percent = as.numeric(.data[[pct_col]])
      )
  } else if (!is.na(count_col)) {
    out <- df %>%
      transmute(
        context = toupper(as.character(.data[[context_col]])),
        Damage_count = as.numeric(.data[[count_col]])
      ) %>%
      mutate(Damage_percent = 100 * Damage_count / sum(Damage_count, na.rm = TRUE)) %>%
      select(context, Damage_percent)
  } else {
    stop("No percentage or count column in Damage-seq normalization file: ", file)
  }

  out %>%
    filter(context %in% CONTEXT_ORDER) %>%
    group_by(context) %>%
    summarise(Damage_percent = sum(Damage_percent, na.rm = TRUE), .groups = "drop") %>%
    complete(context = CONTEXT_ORDER, fill = list(Damage_percent = 0)) %>%
    mutate(
      Damage_percent = 100 * Damage_percent / sum(Damage_percent, na.rm = TRUE),
      context = factor(context, levels = CONTEXT_ORDER)
    )
}

side_label <- function(label, size = 8) {
  ggdraw() +
    draw_label(
      label,
      x = 0.02,
      y = 0.50,
      hjust = 0,
      vjust = 0.5,
      size = size,
      fontface = "bold"
    )
}

blank <- function() ggdraw()

sbs4 <- read_sbs4(sbs4_file)

atl_rep <- pmap_dfr(atl_files, read_atl) %>%
  complete(
    timepoint,
    replicate,
    context = CONTEXT_ORDER,
    fill = list(ATL_percent = 0)
  ) %>%
  mutate(
    timepoint = factor(timepoint, levels = time_order),
    context = factor(context, levels = CONTEXT_ORDER)
  )

atl_summary <- atl_rep %>%
  group_by(timepoint, context) %>%
  summarise(
    mean_ATL_percent = mean(ATL_percent, na.rm = TRUE),
    sem_ATL_percent = sd(ATL_percent, na.rm = TRUE) / sqrt(sum(is.finite(ATL_percent))),
    .groups = "drop"
  ) %>%
  mutate(
    sem_ATL_percent = if_else(is.na(sem_ATL_percent), 0, sem_ATL_percent),
    context = factor(context, levels = CONTEXT_ORDER)
  ) %>%
  left_join(sbs4, by = "context") %>%
  mutate(x_num = as.numeric(context))

damage_rep <- pmap_dfr(damage_files, read_damage_one)

damage_summary <- damage_rep %>%
  group_by(context) %>%
  summarise(
    mean_Damage_percent = mean(Damage_percent, na.rm = TRUE),
    sem_Damage_percent = sd(Damage_percent, na.rm = TRUE) / sqrt(sum(is.finite(Damage_percent))),
    .groups = "drop"
  ) %>%
  mutate(
    sem_Damage_percent = if_else(is.na(sem_Damage_percent), 0, sem_Damage_percent),
    context = factor(context, levels = CONTEXT_ORDER)
  ) %>%
  left_join(sbs4, by = "context") %>%
  mutate(x_num = as.numeric(context))

damage_norm <- read_damage_for_normalization(damage_percent_file)

normalized_rep <- atl_rep %>%
  left_join(damage_norm, by = "context") %>%
  mutate(
    ATL_fraction = ATL_percent / 100,
    Damage_fraction = Damage_percent / 100,
    repair_over_damage = (ATL_fraction + pseudo) / (Damage_fraction + pseudo),
    log2_repair_over_damage = log2(repair_over_damage)
  )

normalized_summary <- normalized_rep %>%
  group_by(timepoint, context) %>%
  summarise(
    mean_log2 = mean(log2_repair_over_damage, na.rm = TRUE),
    sem_log2 = sd(log2_repair_over_damage, na.rm = TRUE) / sqrt(sum(is.finite(log2_repair_over_damage))),
    .groups = "drop"
  ) %>%
  mutate(
    sem_log2 = if_else(is.na(sem_log2), 0, sem_log2),
    timepoint = factor(timepoint, levels = time_order),
    context = factor(context, levels = CONTEXT_ORDER)
  )

cosine_dat <- atl_summary %>%
  group_by(timepoint) %>%
  summarise(
    cosine_similarity_to_SBS4 = cosine_similarity(mean_ATL_percent, SBS4_percent),
    .groups = "drop"
  )

damage_cosine <- damage_summary %>%
  summarise(
    cosine_similarity_to_SBS4 = cosine_similarity(mean_Damage_percent, SBS4_percent)
  ) %>%
  pull(cosine_similarity_to_SBS4)

percent_ymax <- max(
  sbs4$SBS4_percent,
  atl_summary$mean_ATL_percent + atl_summary$sem_ATL_percent,
  damage_summary$mean_Damage_percent + damage_summary$sem_Damage_percent,
  na.rm = TRUE
)
percent_ymax <- ceiling(percent_ymax * 1.12)

log2_ylim <- max(
  abs(normalized_summary$mean_log2 + normalized_summary$sem_log2),
  abs(normalized_summary$mean_log2 - normalized_summary$sem_log2),
  na.rm = TRUE
)
log2_ylim <- ceiling(log2_ylim * 1.15 * 10) / 10

write_tsv(atl_rep, file.path(outdir, "ATL_replicate_percent.tsv"))
write_tsv(atl_summary, file.path(outdir, "ATL_summary_percent.tsv"))
write_tsv(damage_summary, file.path(outdir, "DamageSeq_4NQO_summary_percent.tsv"))
write_tsv(damage_norm, file.path(outdir, "DamageSeq_used_for_repair_over_damage.tsv"))
write_tsv(normalized_rep, file.path(outdir, "ATL_replicate_log2_repair_over_damage.tsv"))
write_tsv(normalized_summary, file.path(outdir, "ATL_summary_log2_repair_over_damage.tsv"))
write_tsv(cosine_dat, file.path(outdir, "ATL_cosine_similarity_to_SBS4.tsv"))

theme_percent <- theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5.7,
      color = "black",
      margin = margin(t = 1)
    ),
    axis.text.y = element_text(size = 5.7, color = "black"),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.32),
    axis.ticks = element_line(color = "black", linewidth = 0.32),
    axis.ticks.length = unit(1.3, "pt"),
    panel.grid = element_blank(),
    plot.title = element_text(
      hjust = 0.5,
      size = 8,
      face = "bold",
      margin = margin(b = 1)
    ),
    plot.margin = margin(2, 3, 4, 3)
  )

theme_log2 <- theme_classic(base_size = 8) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5,
      size = 5.7,
      color = "black",
      margin = margin(t = 1)
    ),
    axis.text.y = element_text(size = 5.7, color = "black"),
    axis.title = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.32),
    axis.ticks = element_line(color = "black", linewidth = 0.32),
    axis.ticks.length = unit(1.3, "pt"),
    panel.grid = element_blank(),
    plot.margin = margin(2, 3, 4, 3)
  )

make_sbs4_panel <- function(show_x = TRUE) {
  ggplot(sbs4, aes(x = context, y = SBS4_percent)) +
    geom_col(
      width = 0.48,
      fill = "skyblue",
      color = "black",
      linewidth = 0.18
    ) +
    scale_x_discrete(drop = FALSE, limits = CONTEXT_ORDER) +
    scale_y_continuous(
      limits = c(0, percent_ymax),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(title = "SBS4") +
    theme_percent +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 5.7, color = "black", margin = margin(t = 1)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line(color = "black", linewidth = 0.32) else element_blank()
    )
}

make_damage_panel <- function(show_x = TRUE) {
  cosine_txt <- paste0("cosine = ", sprintf("%.2f", damage_cosine))

  ggplot(damage_summary, aes(x = context, y = mean_Damage_percent)) +
    geom_col(
      width = 0.48,
      fill = "red3",
      color = "black",
      linewidth = 0.18
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_Damage_percent - sem_Damage_percent, 0),
        ymax = mean_Damage_percent + sem_Damage_percent
      ),
      width = 0.10,
      linewidth = 0.21,
      color = "black"
    ) +
    geom_segment(
      aes(x = x_num, xend = x_num, y = 0, yend = SBS4_percent),
      inherit.aes = FALSE,
      color = "skyblue",
      linewidth = 0.55,
      lineend = "butt"
    ) +
    annotate(
      "text",
      x = 2,
      y = percent_ymax * 0.86,
      label = cosine_txt,
      hjust = 0,
      vjust = 1,
      size = 2.2
    ) +
    scale_x_discrete(drop = FALSE, limits = CONTEXT_ORDER) +
    scale_y_continuous(
      limits = c(0, percent_ymax),
      expand = expansion(mult = c(0, 0.03))
    ) +
    labs(title = "4NQO Damage-seq") +
    theme_percent +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 5.7, color = "black", margin = margin(t = 1)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line(color = "black", linewidth = 0.32) else element_blank()
    )
}

make_atl_percent_panel <- function(tp, show_x = TRUE) {
  df <- atl_summary %>% filter(timepoint == tp)

  cs <- cosine_dat %>%
    filter(timepoint == tp) %>%
    pull(cosine_similarity_to_SBS4)

  cosine_txt <- paste0("cosine = ", sprintf("%.2f", cs))

  ggplot(df, aes(x = context, y = mean_ATL_percent)) +
    geom_col(
      width = 0.48,
      fill = "#4F759B",
      color = "black",
      linewidth = 0.18
    ) +
    geom_errorbar(
      aes(
        ymin = pmax(mean_ATL_percent - sem_ATL_percent, 0),
        ymax = mean_ATL_percent + sem_ATL_percent
      ),
      width = 0.10,
      linewidth = 0.21,
      color = "black"
    ) +
    geom_segment(
      aes(x = x_num, xend = x_num, y = 0, yend = SBS4_percent),
      inherit.aes = FALSE,
      color = "skyblue",
      linewidth = 0.55,
      lineend = "butt"
    ) +
    annotate(
      "text",
      x = 2,
      y = percent_ymax * 0.90,
      label = cosine_txt,
      hjust = 0,
      vjust = 1,
      size = 2.2
    ) +
    scale_x_discrete(drop = FALSE, limits = CONTEXT_ORDER) +
    scale_y_continuous(
      limits = c(0, percent_ymax),
      expand = expansion(mult = c(0, 0.03))
    ) +
    theme_percent +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 5.7, color = "black", margin = margin(t = 1)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line(color = "black", linewidth = 0.32) else element_blank()
    )
}

make_log2_panel <- function(tp, show_x = TRUE) {
  df <- normalized_summary %>% filter(timepoint == tp)

  ggplot(df, aes(x = context, y = mean_log2)) +
    geom_hline(
      yintercept = 0,
      linetype = "dashed",
      linewidth = 0.32,
      color = "black"
    ) +
    geom_errorbar(
      aes(
        ymin = mean_log2 - sem_log2,
        ymax = mean_log2 + sem_log2
      ),
      width = 0.08,
      linewidth = 0.20,
      color = "black"
    ) +
    geom_point(
      size = 1.15,
      color = "#6A3D9A"
    ) +
    scale_x_discrete(drop = FALSE, limits = CONTEXT_ORDER) +
    scale_y_continuous(
      limits = c(-log2_ylim, log2_ylim),
      expand = expansion(mult = c(0.04, 0.04))
    ) +
    theme_log2 +
    theme(
      axis.text.x = if (show_x) element_text(
        angle = 90, hjust = 1, vjust = 0.5,
        size = 5.7, color = "black", margin = margin(t = 1)
      ) else element_blank(),
      axis.ticks.x = if (show_x) element_line(color = "black", linewidth = 0.32) else element_blank()
    )
}

legend_atl <- ggplot() +
  annotate(
    "rect",
    xmin = 0.01,
    xmax = 0.04,
    ymin = 0.50,
    ymax = 0.70,
    fill = "#4F759B",
    color = "black",
    linewidth = 0.1
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 0.80,
    label = "Pinpointed ATL-tXR-Seq signal",
    hjust = 0,
    vjust = 0.5,
    size = 1.6
  ) +
  annotate(
    "segment",
    x = 0.01,
    xend = 0.04,
    y = 0.25,
    yend = 0.25,
    color = "skyblue",
    linewidth = 0.9
  ) +
  annotate(
    "text",
    x = 0.05,
    y = 0.2,
    label = "SBS4",
    hjust = 0,
    vjust = 0.5,
    size = 1.6
  ) +
  coord_cartesian(xlim = c(0, 1), ylim = c(0, 1), clip = "off") +
  theme_void()
  
top_row <- (
  make_sbs4_panel(show_x = TRUE) |
    blank() |
    make_damage_panel(show_x = TRUE) |
    blank()
) +
  plot_layout(widths = c(1, 0.16, 1, 0.055))

atl_legend_row <- legend_atl

make_pair_row <- function(tp, show_x = TRUE) {
  (
    make_atl_percent_panel(tp, show_x = show_x) |
      side_label(tp) |
      make_log2_panel(tp, show_x = show_x) |
      side_label(tp)
  ) +
    plot_layout(widths = c(1, 0.16, 1, 0.055))
}

combined_core <- top_row /
  atl_legend_row /
  make_pair_row("0.5h", show_x = TRUE) /
  make_pair_row("2h", show_x = TRUE) /
  make_pair_row("4h", show_x = TRUE) /
  make_pair_row("8h", show_x = TRUE) +
  plot_layout(heights = c(1, 0.10, 1, 1, 1, 1))

final_plot <- ggdraw() +
  draw_plot(
    combined_core,
    x = 0.085,
    y = 0.055,
    width = 0.88,
    height = 0.91
  ) +
  draw_label(
    "% contribution",
    x = 0.080,
    y = 0.850,
    angle = 90,
    size = 8
  ) +
  draw_label(
    "% contribution",
    x = 0.535,
    y = 0.850,
    angle = 90,
    size = 8
  ) +
  draw_label(
    "% contribution",
    x = 0.075,
    y = 0.435,
    angle = 90,
    size = 8
  ) +
  draw_label(
    expression(log[2]("Repair / Damage")),
    x = 0.545,
    y = 0.435,
    angle = 90,
    size = 8
  ) +
    draw_label(
    "Trinucleotide context",
    x = 0.30,
    y = 0.055,
    size = 8
  ) +
  draw_label(
    "Trinucleotide context",
    x = 0.75,
    y = 0.055,
    size = 8
  )
  
ggsave(
  filename = file.path(outdir, "Fig2_CDEF.pdf"),
  plot = final_plot,
  width = 5.6,
  height = 7.6,
  units = "in",
  device = cairo_pdf
)

ggsave(
  filename = file.path(outdir, "Fig2_CDEF.png"),
  plot = final_plot,
  width = 5.6,
  height = 7.6,
  units = "in",
  dpi = 600
)
message("Done: ", outdir)