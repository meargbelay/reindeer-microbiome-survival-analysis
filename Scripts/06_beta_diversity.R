# ============================================================
# 06_beta_diversity.R
# ============================================================

library(phyloseq)
library(vegan)
library(ggplot2)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

# Relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Bray-Curtis
bray <- phyloseq::distance(
  ps.rel,
  method = "bray"
)

# PCoA
pcoa <- ordinate(
  ps.rel,
  method = "PCoA",
  distance = bray
)

# NMDS
nmds <- ordinate(
  ps.rel,
  method = "NMDS",
  distance = bray
)

# PCoA plot
p <- plot_ordination(
  ps.rel,
  pcoa,
  color = "Survival"
) +
  geom_point(
    size = 3
  ) +
  theme_bw()

ggsave(
  "results/figures/PCoA_survival.png",
  p,
  width = 7,
  height = 5
)

# NMDS
p2 <- plot_ordination(
  ps.rel,
  nmds,
  color = "Survival"
) +
  geom_point(
    size = 3
  ) +
  theme_bw()

ggsave(
  "results/figures/NMDS_survival.png",
  p2,
  width = 7,
  height = 5
)

saveRDS(
  bray,
  "data/processed/bray_curtis.rds"
)

saveRDS(
  pcoa,
  "data/processed/pcoa.rds"
)

saveRDS(
  nmds,
  "data/processed/nmds.rds"
)