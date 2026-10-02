# CDD-HTN integrative analysis: reproducibility code

Repository: <https://github.com/illlllove/CDD-HTN-integrative-analysis>

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22923192.svg)](https://doi.org/10.5281/zenodo.22923192)

This repository contains the custom code and frozen derived results retained for the study integrating cervical disc disorder (CDD), hypertension, cis-eQTL, bulk-transcriptomic and single-cell evidence.

## Current release

Version `v1.3.1` adds the submission-final update under [`submission/v1.3.1/`](submission/v1.3.1/). It contains formal HTN-MVP-SBP colocalization results, the eQTL quantitative-trait sdY sensitivity analysis, updated Figures 1 and 4, Supplementary Tables S1-S2, analysis scripts, a file manifest and SHA-256 checksums.

The earlier frozen Human Genomics submission package remains available under [`submission/v1.1.0/`](submission/v1.1.0/). Existing tags and archived releases are retained unchanged.

## Scope

The repository provides code and documented outputs for:

- environment and dataset auditing;
- reconstruction of the archived shared-SMR screen and prior-sensitivity candidate counts;
- allele harmonization and pairwise colocalization helper functions;
- documented post hoc regional reconstruction;
- generation of Figure 1 and Supplementary Figure S3;
- reconstruction of the frozen Figure 2 regional architecture and Figure 3 prior-sensitivity results;
- formal HTN-MVP-SBP colocalization and the updated Figure 4 external genetic evaluation;
- post hoc eQTL quantitative-trait sdY sensitivity analysis.

The repository does not claim to reconstruct undocumented historical timestamps or prospective threshold selection. Scripts are labelled according to whether they reconstruct archived results, perform a documented post hoc audit, or generate figures.

## Data access

Raw third-party data are not redistributed.

- FinnGen R13 CDD endpoint: `M13_CERVICDISC`.
- Hypertension GWAS: GWAS Catalog accession `GCST90044351`.
- cis-eQTL: eQTLGen Consortium Phase 1 blood cis-eQTL summary data.
- UK Biobank M50 sensitivity data: publicly released Neale Lab UK Biobank Round 2 summary statistics.
- MVP SBP: controlled-access summary statistics associated with dbGaP accession `phs001672.v11.p1`; these data are not included.
- MVP IVDD: summary statistics associated with dbGaP accession `phs002453.v1.p1`; restricted source data are not redistributed.
- Bulk expression: GEO accessions `GSE153761`, `GSE28360`, and `GSE34095`.
- Single-cell RNA-seq: GEO accession `GSE230809`.

See `data/README.md`, `submission/v1.1.0/README_submission_v1.1.0.md`, and `submission/v1.3.1/README_submission_v1.3.1.md` for the expected inputs and frozen submission packages.

## Running the code

R 4.5 or later is recommended. Set the project root before running the original retained scripts:

```powershell
$env:CDD_HTN_ROOT = "D:\\path\\to\\CDD_HTN"
Rscript scripts/01_dataset_audit.R
```

Run the v1.3.1 submission-specific analyses from the repository root after providing the locally authorized inputs described by the package documentation:

```powershell
Rscript submission/v1.3.1/scripts/01_HTN_MVP_SBP_formal_coloc.R
Rscript submission/v1.3.1/scripts/02_sdY_sensitivity.R
Rscript submission/v1.3.1/scripts/03_Figure4_reconstruction.R
```

## Dependencies

Core R packages used by the retained scripts include `data.table`, `coloc`, `ggplot2`, `patchwork`, `dplyr`, `readr`, `tidyr`, `scales`, `grid`, and `ragg`. The Python reconstruction script uses only the Python 3 standard library.

## Reproducibility limits

Several historical analyses were originally performed interactively and no complete original script was retained. The deposited scripts reproduce or audit only the components explicitly described in their package documentation. Third-party data-use restrictions and the absence of some historical execution records prevent a one-command reconstruction from raw data.

## Citation

Please cite the associated article and the archived Zenodo release:

> Liang W, Yao H, Tu W, Zhai Y. Reproducibility code for the CDD-HTN integrative analysis. Version v1.3.1. Zenodo (2026). <https://doi.org/10.5281/zenodo.22923192>

The DOI above represents all versions and resolves to the latest Zenodo release.

## License

Custom code in this repository is released under the MIT License. Third-party data remain subject to their original terms and licenses.
