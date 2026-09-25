# ============================================================
# 11_differential_abundance_maaslin2.R
# ============================================================

library(Maaslin2)
library(phyloseq)

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

# Convert to abundance matrix
abundance <- as.data.frame(
  t(as(otu_table(mouth_genus), "matrix"))
)

metadata <- data.frame(
  sample_data(mouth_genus)
)

# Ensure matching sample IDs
metadata <- metadata[
  rownames(abundance),
  ,
  drop = FALSE
]

fit_data <- cbind(
  metadata,
  abundance
)

fit <- Maaslin2(
  input_data = abundance,
  input_metadata = metadata,
  output = "results/Maaslin2_mouth",
  fixed_effects = c("Survival")
)