# ============================================================
# 12_differential_abundance_ancombc.R
# ============================================================

library(phyloseq)
library(ANCOMBC)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

mouth <- subset_samples(
  ps,
  Site == "Mouth"
)

mouth_genus <- tax_glom(
  mouth,
  taxrank = "Genus"
)

result <- ancombc(
  data = mouth_genus,
  formula = "Survival",
  p_adj_method = "BH",
  prv_cut = 0.10,
  lib_cut = 0,
  group = "Survival",
  struc_zero = TRUE,
  neg_lb = TRUE
)

writeRDS(
  result,
  "results/tables/ANCOMBC_mouth.rds"
)