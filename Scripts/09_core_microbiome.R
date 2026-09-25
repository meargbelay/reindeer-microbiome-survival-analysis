# ============================================================
# 09_core_microbiome.R
# ============================================================

library(phyloseq)
library(microbiome)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

# Relative abundance
ps.rel <- transform_sample_counts(
  ps,
  function(x) x / sum(x)
)

# Example: mouth survivors
mouth_survivors <- subset_samples(
  ps.rel,
  Site == "Mouth" &
    Survival == 1
)

core_survivors <- core(
  mouth_survivors,
  detection = 0.0001,
  prevalence = 0.70
)

core_survivors