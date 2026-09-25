# ============================================================
# 04_relative_abundance.R
# ============================================================

library(phyloseq)
library(ggplot2)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

# Transform to relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Phylum-level
ps.phylum <- tax_glom(
  ps.rel,
  taxrank = "Phylum"
)

# Convert to dataframe
df.phylum <- psmelt(ps.phylum)

# Plot
p <- ggplot(
  df.phylum,
  aes(
    x = Sample,
    y = Abundance,
    fill = Phylum
  )
) +
  geom_bar(
    stat = "identity"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_blank()
  )

ggsave(
  "results/figures/phylum_relative_abundance.png",
  p,
  width = 12,
  height = 7
)