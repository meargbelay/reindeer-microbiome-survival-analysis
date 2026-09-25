# ============================================================
# 01_quality_control_dada2.R
# Reindeer microbiome survival project
# ============================================================

library(dada2)

# ------------------------------------------------------------
# 1. Input/output directories
# ------------------------------------------------------------

path <- "data/raw"

fwd <- sort(list.files(
  path,
  pattern = "_R1_001.fastq.gz",
  full.names = TRUE
))

rev <- sort(list.files(
  path,
  pattern = "_R2_001.fastq.gz",
  full.names = TRUE
))

sample_names <- sub(
  "_R1_001.fastq.gz",
  "",
  basename(fwd)
)

# ------------------------------------------------------------
# 2. Quality profiles
# ------------------------------------------------------------

plotQualityProfile(fwd[1:2])
plotQualityProfile(rev[1:2])

# ------------------------------------------------------------
# 3. Filter and trim
# ------------------------------------------------------------

filt_path <- "data/processed/filtered"

dir.create(
  filt_path,
  recursive = TRUE,
  showWarnings = FALSE
)

filt_fwd <- file.path(
  filt_path,
  paste0(sample_names, "_F_filt.fastq.gz")
)

filt_rev <- file.path(
  filt_path,
  paste0(sample_names, "_R_filt.fastq.gz")
)

out <- filterAndTrim(
  fwd,
  filt_fwd,
  rev,
  filt_rev,
  maxN = 0,
  maxEE = c(2, 2),
  truncQ = 2,
  truncLen = c(240, 200),
  rm.phix = TRUE,
  compress = TRUE,
  multithread = TRUE
)

write.csv(
  out,
  "results/tables/filtering_summary.csv"
)

# ------------------------------------------------------------
# 4. Learn error rates
# ------------------------------------------------------------

errF <- learnErrors(
  filt_fwd,
  multithread = TRUE
)

errR <- learnErrors(
  filt_rev,
  multithread = TRUE
)

plotErrors(errF, nominalQ = TRUE)
plotErrors(errR, nominalQ = TRUE)

# ------------------------------------------------------------
# 5. Dereplication
# ------------------------------------------------------------

derepF <- derepFastq(filt_fwd)
derepR <- derepFastq(filt_rev)

names(derepF) <- sample_names
names(derepR) <- sample_names

# ------------------------------------------------------------
# 6. Denoising
# ------------------------------------------------------------

dadaF <- dada(
  derepF,
  err = errF,
  multithread = TRUE
)

dadaR <- dada(
  derepR,
  err = errR,
  multithread = TRUE
)

# ------------------------------------------------------------
# 7. Merge paired reads
# ------------------------------------------------------------

mergers <- mergePairs(
  dadaF,
  derepF,
  dadaR,
  derepR,
  verbose = TRUE
)

# ------------------------------------------------------------
# 8. Construct ASV table
# ------------------------------------------------------------

seqtab <- makeSequenceTable(mergers)

# ------------------------------------------------------------
# 9. Remove chimeras
# ------------------------------------------------------------

seqtab.nochim <- removeBimeraDenovo(
  seqtab,
  method = "consensus",
  multithread = TRUE,
  verbose = TRUE
)

# ------------------------------------------------------------
# 10. Save ASV table
# ------------------------------------------------------------

saveRDS(
  seqtab.nochim,
  "data/processed/seqtab_nochim.rds"
)

cat(
  "Number of samples:",
  nrow(seqtab.nochim),
  "\n"
)

cat(
  "Number of ASVs:",
  ncol(seqtab.nochim),
  "\n"
)