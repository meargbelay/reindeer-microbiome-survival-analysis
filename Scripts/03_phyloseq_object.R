# ============================================================
# 03_phyloseq_object.R
# ============================================================

library(phyloseq)

seqtab <- readRDS(
  "data/processed/seqtab_nochim.rds"
)

tax <- readRDS(
  "data/processed/taxonomy_silva.rds"
)

# Read metadata
metadata <- read.csv(
  "metadata/sample_metadata.csv",
  row.names = 1,
  check.names = FALSE
)

# Create objects
OTU <- otu_table(
  seqtab,
  taxa_are_rows = FALSE
)

TAX <- tax_table(
  as.matrix(tax)
)

SAM <- sample_data(metadata)

# Combine
ps <- phyloseq(
  OTU,
  TAX,
  SAM
)

# Remove zero-depth samples
ps <- prune_samples(
  sample_sums(ps) > 0,
  ps
)

# Save
saveRDS(
  ps,
  "data/processed/phyloseq_object.rds"
)

ps