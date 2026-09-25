# ============================================================
# 05_alpha_diversity.R
# ============================================================

library(phyloseq)
library(ggplot2)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

alpha <- estimate_richness(
  ps,
  measures = c(
    "Observed",
    "Shannon"
  )
)

alpha$SampleID <- rownames(alpha)

metadata <- data.frame(
  sample_data(ps)
)

metadata$SampleID <- rownames(metadata)

alpha <- merge(
  alpha,
  metadata,
  by = "SampleID"
)

write.csv(
  alpha,
  "results/tables/alpha_diversity.csv",
  row.names = FALSE
)

# Shannon by survival
p <- ggplot(
  alpha,
  aes(
    x = factor(Survival),
    y = Shannon
  )
) +
  geom_boxplot() +
  geom_jitter(width = 0.1) +
  theme_bw() +
  labs(
    x = "Survival",
    y = "Shannon diversity"
  )

ggsave(
  "results/figures/shannon_survival.png",
  p,
  width = 7,
  height = 5
)