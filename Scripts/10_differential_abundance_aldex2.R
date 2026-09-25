# ============================================================
# 10_differential_abundance_aldex2.R
# ============================================================

library(phyloseq)
library(ALDEx2)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

# Select mouth samples
mouth <- subset_samples(
  ps,
  Site == "Mouth"
)

# Aggregate to genus
mouth_genus <- tax_glom(
  mouth,
  taxrank = "Genus"
)

# Extract count matrix
counts <- as(otu_table(mouth_genus), "matrix")

if (!taxa_are_rows(mouth_genus)) {
  counts <- t(counts)
}

metadata <- data.frame(
  sample_data(mouth_genus)
)

condition <- metadata$Survival

# ALDEx2
aldex <- aldex(
  counts,
  condition,
  test = "wilcox",
  effect = TRUE,
  denom = "all",
  verbose = TRUE
)

write.csv(
  aldex,
  "results/tables/ALDEx2_mouth.csv"
)