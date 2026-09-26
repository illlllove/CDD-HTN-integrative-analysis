# CDD-HTN integrative analysis
## Human Genomics submission release v1.1.0

This directory contains the frozen derived data, final figures, and
submission-specific reconstruction scripts accompanying the Human Genomics
manuscript:

**Prior-sensitive regional colocalization between cervical disc disorders
and hypertension highlights two candidate regions**

## Scope

This release documents the computational evidence used in the revised
Human Genomics submission.

The main genetic findings are organized as two candidate genomic regions:

- chromosome 10, containing the competing transcript hypotheses NT5C2,
  CNNM2 and MARCKSL1P1;
- chromosome 18, containing RMC1 as a locus-specific follow-up hypothesis
  within the broader NPC1-RMC1 region.

The release does not establish a common causal variant across cis-eQTL,
hypertension and cervical disc disorders, and it does not assign a resolved
causal effector transcript.

## Directory structure

### frozen_data/Figure2

Frozen derived tables used for regional association architecture and
conditional H4 localization.

### frozen_data/Figure3

Frozen derived tables used for shared-variant prior-sensitivity analyses.

### frozen_data/Figure4

Frozen derived tables used for external genetic evaluation with MVP SBP
and UK Biobank M50.

### figures

Contains the final manuscript PDFs.

### figures/reconstructed

Contains PDFs produced by the submission-specific reconstruction scripts
from the frozen derived tables.

These reconstructed figures are provided as reproducibility checks. They
are not intended to be pixel-identical replacements for the frozen final
manuscript figures.

### scripts

Contains submission-specific R reconstruction scripts generated from the
frozen derived tables:

- Figure2_regional_architecture_RECONSTRUCTION.R
- Figure3_prior_sensitivity_RECONSTRUCTION.R
- Figure4_external_genetic_evaluation_RECONSTRUCTION.R

These scripts are not presented as the original historical interactive
analysis scripts. They provide a documented reconstruction of the final
figure-level results from the frozen derived tables.

## Software

The reconstruction scripts were tested with:

- R 4.5.0
- data.table
- ggplot2
- patchwork

## Reproducibility boundary

The archived files permit reconstruction of the retained candidate set,
regional result tables, prior-sensitivity summaries and external genetic
evaluation used in the Human Genomics submission.

The archive does not imply that every analytical threshold, execution step
or analysis timestamp was prospectively preregistered. Where historical
analysis chronology could not be independently established, the manuscript
and repository describe the evidence according to its documented
inferential role rather than asserting a fully prospective sequence.

## Data restrictions

Raw third-party data are not redistributed.

Relevant source datasets include:

- FinnGen Release 13 M13_CERVICDISC;
- hypertension GWAS GCST90044351;
- eQTLGen Phase I blood cis-eQTL data;
- Million Veteran Program systolic blood pressure summary statistics;
- UK Biobank M50 summary statistics;
- GEO bulk-expression datasets GSE153761, GSE28360 and GSE34095;
- GEO single-cell dataset GSE230809.

Controlled-access MVP summary statistics are not included in this release.

## Version

Human Genomics submission release: v1.1.0