library(tidyverse)

infile <- "results/damage_sites/NHF1-4NQO-30m-r1.readstart_minus3_to_0.sequence.tsv"

out_tsv <- "NHF1-4NQO-30m-r1_readstart_minus3_to_0_monomer_percentages.tsv"

df <- read_tsv(
  infile,
  col_names = c("ReadID", "Sequence"),
  show_col_types = FALSE
) %>%
  transmute(Sequence = toupper(Sequence)) %>%
  filter(nchar(Sequence) == 4)

mono_df <- tibble(
  Position = rep(c("-3", "-2", "-1", "Read start"), each = nrow(df)),
  Base = c(
    substr(df$Sequence, 1, 1),
    substr(df$Sequence, 2, 2),
    substr(df$Sequence, 3, 3),
    substr(df$Sequence, 4, 4)
  )
) %>%
  filter(Base %in% c("A", "C", "G", "T")) %>%
  mutate(
    Position = factor(Position, levels = c("-3", "-2", "-1", "Read start")),
    Base = factor(Base, levels = c("C", "T", "A", "G"))
  ) %>%
  count(Position, Base, name = "Count", .drop = FALSE) %>%
  group_by(Position) %>%
  mutate(Percentage = 100 * Count / sum(Count)) %>%
  ungroup()

write_tsv(mono_df, out_tsv)

p <- ggplot(mono_df, aes(x = Position, y = Percentage, fill = Base)) +
  geom_col(width = 0.75, color = "black", linewidth = 0.15) +
  scale_fill_manual(values = c(
    "C" = "dodgerblue4",
    "T" = "orange",
    "A" = "green4",
    "G" = "darkorchid4"
  )) +
  labs(
    title = "Damage-seq sequence context at read start",
    subtitle = "4-NQO, positions -3 to read start",
    x = NULL,
    y = "Base percentage",
    fill = "Base"
  ) +
  theme_classic(base_size = 14)

ggsave("NHF1-4NQO-30m-r1_readstart_minus3_to_0_monomer.png", p, width = 5.5, height = 4.5, dpi = 300)
ggsave("NHF1-4NQO-30m-r1_readstart_minus3_to_0_monomer.pdf", p, width = 5.5, height = 4.5, useDingbats = FALSE)