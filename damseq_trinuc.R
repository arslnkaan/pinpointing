#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

outdir <- "DamageSeq_CPD_CT_TC_CC_trinuc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

ct_tc_file <- "/work/users/a/r/arslank/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CT_TC_Ccenter_3nt.sequence.tsv"
cc_file    <- "/work/users/a/r/arslank/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CC_Ccenter_3nt.sequence.tsv"

context_order <- c(
  "ACA", "ACC", "ACG", "ACT",
  "CCA", "CCC", "CCG", "CCT",
  "GCA", "GCC", "GCG", "GCT",
  "TCA", "TCC", "TCG", "TCT"
)

read_damage_file <- function(file, fallback_class) {
  if (!file.exists(file)) stop("Missing file: ", file)

  read_tsv(
    file,
    col_names = c("site_id", "context"),
    show_col_types = FALSE
  ) %>%
    mutate(
      context = toupper(context),
      lesion_class = case_when(
        str_detect(site_id, "_CT_Ccenter") ~ "CT",
        str_detect(site_id, "_TC_Ccenter") ~ "TC",
        str_detect(site_id, "_CC_Ccenter") ~ "CC",
        TRUE ~ fallback_class
      )
    ) %>%
    filter(
      lesion_class %in% c("CT", "TC", "CC"),
      context %in% context_order,
      substr(context, 2, 2) == "C"
    )
}

damage <- bind_rows(
  read_damage_file(ct_tc_file, "CT_TC"),
  read_damage_file(cc_file, "CC")
)

counts_by_class <- damage %>%
  count(lesion_class, context, name = "Count") %>%
  complete(
    lesion_class = c("CT", "TC", "CC"),
    context = context_order,
    fill = list(Count = 0)
  ) %>%
  group_by(lesion_class) %>%
  mutate(
    Total = sum(Count),
    Percentage = if_else(Total > 0, 100 * Count / Total, 0)
  ) %>%
  ungroup()

counts_combined <- damage %>%
  count(context, name = "Count") %>%
  complete(context = context_order, fill = list(Count = 0)) %>%
  mutate(
    lesion_class = "CT_TC_CC_combined",
    Total = sum(Count),
    Percentage = if_else(Total > 0, 100 * Count / Total, 0)
  ) %>%
  select(lesion_class, context, Count, Total, Percentage)

write_tsv(counts_by_class, file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_by_class.tsv"))
write_tsv(counts_combined, file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_combined.tsv"))

write_tsv(
  counts_combined %>% transmute(RC = context, Count, Percentage),
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_combined_percentages.tsv")
)

p_combined <- ggplot(
  counts_combined,
  aes(x = factor(context, levels = context_order), y = Percentage)
) +
  geom_col(fill = "#4F759B", color = "black", linewidth = 0.25, width = 0.75) +
  labs(
    title = "Damage-seq CPD trinucleotide context",
    subtitle = "CT + TC + CC combined; TT excluded",
    x = NULL,
    y = "Percent contribution"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_combined.pdf"),
       p_combined, width = 5, height = 3.2, device = cairo_pdf)

ggsave(file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_combined.png"),
       p_combined, width = 5, height = 3.2, dpi = 600)

p_by_class <- ggplot(
  counts_by_class,
  aes(x = factor(context, levels = context_order), y = Percentage)
) +
  geom_col(fill = "#4F759B", color = "black", linewidth = 0.25, width = 0.75) +
  facet_wrap(~lesion_class, ncol = 1) +
  labs(
    title = "Damage-seq CPD trinucleotide context",
    subtitle = "CT, TC, and CC separately; TT excluded",
    x = NULL,
    y = "Percent contribution"
  ) +
  theme_classic(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black"),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_by_class.pdf"),
       p_by_class, width = 5, height = 7, device = cairo_pdf)

ggsave(file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_by_class.png"),
       p_by_class, width = 5, height = 7, dpi = 600)

message("Done: ", outdir)
message("Normalization file: ", file.path(outdir, "DamageSeq_CPD_CT_TC_CC_trinuc_combined_percentages.tsv"))