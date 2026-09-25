# ============================================================
# 08_survival_glm.R
# ============================================================

library(phyloseq)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

pcoa <- readRDS(
  "data/processed/pcoa.rds"
)

# Extract PCoA coordinates
scores <- as.data.frame(
  pcoa$vectors
)

scores$SampleID <- rownames(scores)

# Metadata
metadata <- data.frame(
  sample_data(ps)
)

metadata$SampleID <- rownames(metadata)

# Combine
df <- merge(
  scores,
  metadata,
  by = "SampleID"
)

# Rename coordinates
df$PCoA1 <- df$Axis.1
df$PCoA2 <- df$Axis.2
df$PCoA3 <- df$Axis.3

# Logistic regression
model <- glm(
  Survival ~ PCoA1 + PCoA2 + PCoA3,
  data = df,
  family = binomial
)

summary(model)

# Odds ratios
OR <- exp(coef(model))

results <- data.frame(
  coefficient = coef(model),
  odds_ratio = OR,
  p_value = summary(model)$coefficients[, 4]
)

write.csv(
  results,
  "results/tables/survival_glm.csv"
)