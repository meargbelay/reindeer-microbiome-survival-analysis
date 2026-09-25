# ============================================================
# 13_method_comparison.R
# ============================================================

aldex_taxa <- c(
  "Capnocytophaga",
  "Pseudomonas",
  "Sphingomonas"
)

maaslin_taxa <- c(
  "Capnocytophaga",
  "Xylophilus",
  "Olsenella",
  "Monoglobus",
  "Roseisolibacter",
  "Prevotellaceae UCG-001",
  "Stakelama",
  "Parafilimonas",
  "[Ruminococcus] torques group",
  "Desulfovibrio"
)

ancom_taxa <- c(
  "Monoglobus",
  "Roseisolibacter",
  "Parafilimonas",
  "Phascolarctobacterium",
  "Chitinophaga",
  "Flexivirga",
  "Pedococcus-Phycicoccus",
  "Nocardia",
  "Amnipila",
  "Colidextribacter",
  "Alysiella",
  "Candidatus Berkiella"
)

all_taxa <- unique(
  c(
    aldex_taxa,
    maaslin_taxa,
    ancom_taxa
  )
)

comparison <- data.frame(
  Genus = all_taxa,
  ALDEx2 = all_taxa %in% aldex_taxa,
  MaAsLin2 = all_taxa %in% maaslin_taxa,
  ANCOMBC = all_taxa %in% ancom_taxa
)

write.csv(
  comparison,
  "results/tables/differential_abundance_method_comparison.csv",
  row.names = FALSE
)