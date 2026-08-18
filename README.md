# reindeer-microbiome-survival-analysis
Bioinformatics analysis of the gut microbiome of reindeer calves in relation to survival.
# Reindeer Microbiome Survival Analysis

## Overview

This project presents a bioinformatics analysis of the gut microbiome of reindeer calves in relation to survival.

The analysis investigates whether gut microbial community composition differs between surviving and non-surviving individuals and identifies microbial taxa associated with survival.

## Research Questions

- Does gut microbiome composition differ between surviving and non-surviving reindeer calves?
- Are alpha-diversity measures associated with survival?
- Does beta-diversity differ between survival groups?
- Which microbial taxa are associated with survival?
- Do different differential abundance methods identify overlapping microbial taxa?

## Bioinformatics Workflow

The analysis includes:

1. Sequence data preprocessing
2. ASV inference using DADA2
3. Chimera removal
4. Taxonomic classification
5. Phyloseq object construction
6. Alpha-diversity analysis
7. Beta-diversity analysis
8. PCoA
9. PERMANOVA
10. Differential abundance analysis

## Tools and Packages

- R
- DADA2
- Phyloseq
- SILVA
- vegan
- ALDEx2
- MaAsLin2
- ANCOM-BC

## Statistical Analysis

### Alpha Diversity

Alpha diversity was evaluated using measures such as:

- Observed ASVs
- Shannon diversity

### Beta Diversity

Community composition was evaluated using:

- Bray-Curtis dissimilarity
- Principal Coordinates Analysis (PCoA)
- PERMANOVA

### Differential Abundance

Differential abundance was evaluated using three complementary approaches:

- ALDEx2
- MaAsLin2
- ANCOM-BC

The comparison of these approaches was used to identify taxa repeatedly detected across statistical methods.

## Project Structure

```text
reindeer-microbiome-survival-analysis/
│
├── README.md
│
├── R/
│   ├── 01_data_preprocessing.R
│   ├── 02_phyloseq_analysis.R
│   ├── 03_alpha_diversity.R
│   ├── 04_beta_diversity.R
│   ├── 05_PERMANOVA.R
│   └── 06_differential_abundance.R
│
├── figures/
│
├── results/
│
└── docs/
```

## Reproducibility

The scripts are organized according to the main stages of the bioinformatics and statistical workflow.

Raw sequencing data and confidential metadata are not included in this repository.

## Research Context

This project was conducted as part of Master's studies in Bioinformatics and focuses on the relationship between the gut microbiome and survival in reindeer calves.

## Author

**Mearg Belay Brhane**

MSc Bioinformatics  
Gothenburg, Sweden
