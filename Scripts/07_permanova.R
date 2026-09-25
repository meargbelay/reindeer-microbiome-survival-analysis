# ============================================================
# 07_permanova.R
# ============================================================

library(vegan)
library(phyloseq)

ps <- readRDS(
  "data/processed/phyloseq_object.rds"
)

bray <- readRDS(
  "data/processed/bray_curtis.rds"
)

metadata <- data.frame(
  sample_data(ps)
)

# PERMANOVA
permanova <- adonis2(
  bray ~ Site + Sex + Survival,
  data = metadata,
  permutations = 999,
  by = "margin"
)

print(permanova)

write.csv(
  as.data.frame(permanova),
  "results/tables/permanova.csv"
)

# Dispersion
disp_site <- betadisper(
  bray,
  metadata$Site
)

disp_sex <- betadisper(
  bray,
  metadata$Sex
)

disp_survival <- betadisper(
  bray,
  metadata$Survival
)

anova(disp_site)
anova(disp_sex)
anova(disp_survival)