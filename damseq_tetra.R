#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(tidyverse)
})

outdir <- "DamageSeq_CPD_CT_TC_CC_tetranuc"
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

genome_fasta <- "/proj/seq/data/hg38_UCSC/Sequence/WholeGenomeFasta/genome.fa"

ct_tc_file <- "/work/users/a/r/arslank/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CT_TC_Ccenter_3nt.sequence.tsv"
cc_file    <- "/work/users/a/r/arslank/damseq/NHF1_CPD_0h_r1_results/damage_sites/NHF1_CPD_0h_CC_Ccenter_3nt.sequence.tsv"

bases <- c("A", "C", "G", "T")

context_order <- c(
  expand_grid(n1 = bases, n2 = bases) %>%
    transmute(context = paste0(n1, "CT", n2)) %>%
    pull(context),

  expand_grid(n1 = bases, n2 = bases) %>%
    transmute(context = paste0(n1, "TC", n2)) %>%
    pull(context),

  expand_grid(n1 = bases, n2 = bases) %>%
    transmute(context = paste0(n1, "CC", n2)) %>%
    pull(context)
)

read_damage_file <- function(file, fallback_class) {
  if (!file.exists(file)) stop("Missing file: ", file)

  read_tsv(
    file,
    col_names = c("site_id", "trinuc_context"),
    show_col_types = FALSE
  ) %>%
    mutate(
      trinuc_context = toupper(trinuc_context),
      lesion_class = case_when(
        str_detect(site_id, "_CT_Ccenter") ~ "CT",
        str_detect(site_id, "_TC_Ccenter") ~ "TC",
        str_detect(site_id, "_CC_Ccenter") ~ "CC",
        TRUE ~ fallback_class
      ),
      chrom = str_extract(site_id, "chr[^:]+"),
      start = as.integer(str_match(site_id, ":(\\d+)-(\\d+)\\(")[, 2]),
      end   = as.integer(str_match(site_id, ":(\\d+)-(\\d+)\\(")[, 3]),
      strand = str_match(site_id, "\\(([+-])\\)$")[, 2]
    ) %>%
    filter(
      lesion_class %in% c("CT", "TC", "CC"),
      !is.na(chrom),
      !is.na(start),
      !is.na(end),
      !is.na(strand),
      end - start == 3,
      start >= 1
    )
}

damage_sites <- bind_rows(
  read_damage_file(ct_tc_file, "CT_TC"),
  read_damage_file(cc_file, "CC")
)

write_tsv(
  damage_sites,
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_3nt_sites_parsed.tsv")
)

# Extract 5 nt around each original 3-nt interval, strand-aware.
# The original 3-mer is positions 2:4 of the extracted 5-mer.
bed <- damage_sites %>%
  mutate(
    bed_start = start - 1L,
    bed_end = end + 1L,
    bed_name = paste0("site_", row_number(), "|", lesion_class)
  ) %>%
  select(chrom, bed_start, bed_end, bed_name, score = lesion_class, strand)

tmp_bed <- tempfile(fileext = ".bed")
tmp_fa <- tempfile(fileext = ".fa")

write_tsv(bed, tmp_bed, col_names = FALSE)

cmd <- paste(
  "bedtools getfasta",
  "-fi", shQuote(genome_fasta),
  "-bed", shQuote(tmp_bed),
  "-s",
  "-name",
  "-fo", shQuote(tmp_fa)
)

status <- system(cmd)

if (status != 0) {
  stop("bedtools getfasta failed. Make sure bedtools is loaded.")
}

fa_lines <- readLines(tmp_fa)

headers <- fa_lines[grepl("^>", fa_lines)] %>%
  str_remove("^>")

seqs <- fa_lines[!grepl("^>", fa_lines)] %>%
  toupper()

if (length(headers) != length(seqs)) {
  stop("FASTA header/sequence mismatch.")
}

tetranuc_sites <- tibble(
  header = headers,
  seq5 = seqs
) %>%
  separate(header, into = c("site_index", "lesion_class"), sep = "\\|", remove = FALSE) %>%
  mutate(
    lesion_class = str_extract(lesion_class, "CT|TC|CC"),
    tetranuc_context = case_when(
      lesion_class == "CT" ~ substr(seq5, 2, 5), # NCTN
      lesion_class == "TC" ~ substr(seq5, 1, 4), # NTCN
      lesion_class == "CC" ~ case_when(
        str_detect(substr(seq5, 2, 5), "^[ACGT]CC[ACGT]$") ~ substr(seq5, 2, 5), # NCCN
        str_detect(substr(seq5, 1, 4), "^[ACGT]CC[ACGT]$") ~ substr(seq5, 1, 4), # NCCN
        TRUE ~ NA_character_
      ),
      TRUE ~ NA_character_
    )
  ) %>%
  filter(
    !is.na(tetranuc_context),
    tetranuc_context %in% context_order,
    substr(tetranuc_context, 2, 3) %in% c("CT", "TC", "CC")
  )

write_tsv(
  tetranuc_sites,
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_sites.tsv")
)

counts_by_class <- tetranuc_sites %>%
  count(lesion_class, tetranuc_context, name = "Count") %>%
  complete(
    lesion_class = c("CT", "TC", "CC"),
    tetranuc_context = context_order,
    fill = list(Count = 0)
  ) %>%
  group_by(lesion_class) %>%
  mutate(
    Total = sum(Count),
    Percentage = if_else(Total > 0, 100 * Count / Total, 0)
  ) %>%
  ungroup()

counts_combined <- tetranuc_sites %>%
  count(tetranuc_context, name = "Count") %>%
  complete(tetranuc_context = context_order, fill = list(Count = 0)) %>%
  mutate(
    lesion_class = "CT_TC_CC_combined",
    Total = sum(Count),
    Percentage = if_else(Total > 0, 100 * Count / Total, 0)
  ) %>%
  select(lesion_class, tetranuc_context, Count, Total, Percentage)

write_tsv(
  counts_by_class,
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_by_class.tsv")
)

write_tsv(
  counts_combined,
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_combined.tsv")
)

write_tsv(
  counts_combined %>%
    transmute(Context = tetranuc_context, Count, Percentage),
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_combined_percentages.tsv")
)

p_combined <- ggplot(
  counts_combined,
  aes(x = factor(tetranuc_context, levels = context_order), y = Percentage)
) +
  geom_col(fill = "#4F759B", color = "black", linewidth = 0.25, width = 0.75) +
  labs(
    title = "Damage-seq CPD tetranucleotide context",
    subtitle = "CT + TC + CC; TT excluded",
    x = NULL,
    y = "Percent contribution"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 6),
    axis.text.y = element_text(color = "black"),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_combined.pdf"),
  p_combined,
  width = 7,
  height = 3.4,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_combined.png"),
  p_combined,
  width = 7,
  height = 3.4,
  dpi = 600
)

p_by_class <- ggplot(
  counts_by_class,
  aes(x = factor(tetranuc_context, levels = context_order), y = Percentage)
) +
  geom_col(fill = "#4F759B", color = "black", linewidth = 0.25, width = 0.75) +
  facet_wrap(~lesion_class, ncol = 1) +
  labs(
    title = "Damage-seq CPD tetranucleotide context",
    subtitle = "CT, TC, and CC separately; TT excluded",
    x = NULL,
    y = "Percent contribution"
  ) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, color = "black", size = 6),
    axis.text.y = element_text(color = "black"),
    strip.background = element_rect(fill = "white", color = NA),
    strip.text = element_text(size = 11),
    plot.title = element_text(hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5)
  )

ggsave(
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_by_class.pdf"),
  p_by_class,
  width = 7,
  height = 8,
  device = cairo_pdf
)

ggsave(
  file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_by_class.png"),
  p_by_class,
  width = 7,
  height = 8,
  dpi = 600
)

message("Done: ", outdir)
message("Normalization file: ", file.path(outdir, "DamageSeq_CPD_CT_TC_CC_tetranuc_combined_percentages.tsv"))