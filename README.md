# CDD-HTN integrative analysis: reproducibility code

Repository: https://github.com/illlllove/CDD-HTN-integrative-analysis

This repository contains the custom code retained for the study integrating cervical disc disorder (CDD), hypertension, cis-eQTL, bulk-transcriptomic and single-cell evidence.

## Scope

The repository provides code for:

- environment and dataset auditing;
- reconstruction of the archived shared-SMR screen and prior-sensitivity candidate counts;
- allele harmonization and pairwise colocalization helper functions;
- the documented post hoc regional reconstruction used for Figure 3;
- generation of Figure 1;
- generation of Supplementary Figure S3.

The repository does not claim to reconstruct undocumented historical timestamps or prospective threshold selection. Scripts are labelled according to whether they reconstruct archived results, perform a documented post hoc audit, or generate figures.

## Data access

Raw third-party data are not redistributed.

- FinnGen R13 CDD endpoint: `M13_CERVICDISC`.
- Hypertension GWAS: GWAS Catalog accession `GCST90044351`.
- cis-eQTL: eQTLGen Consortium Phase 1 blood cis-eQTL summary data.
- UK Biobank M50 sensitivity data: publicly released Neale Lab UK Biobank Round 2 summary statistics.
- MVP SBP: controlled-access summary statistics associated with dbGaP accession `phs001672.v11.p1`; these data are not included.
- Bulk expression: GEO accessions `GSE153761`, `GSE28360`, and `GSE34095`.
- Single-cell RNA-seq: GEO accession `GSE230809`.

See `data/README.md` for the expected input layout.

## Running the code

R 4.5 or later is recommended. Set the project root before running scripts:

```powershell
$env:CDD_HTN_ROOT = "D:\path\to\CDD_HTN"
Rscript scripts/01_dataset_audit.R
```

An optional output directory can be set for plotting scripts:

```powershell
$env:CDD_HTN_OUTPUT_DIR = "D:\path\to\outputs\figures"
Rscript scripts/05_plot_figure1.R
Rscript scripts/06_plot_supplementary_figure_s3.R
```

The shared-SMR reconstruction scripts operate on the lightweight audit-package layout described in `data/README.md`:

```powershell
Rscript scripts/02_shared_SMR_coloc_reconstruction.R D:\path\to\audit_package
python scripts/02_shared_SMR_coloc_reconstruction.py --root D:\path\to\audit_package
```

## Dependencies

Core R packages used by the retained scripts include `data.table`, `coloc`, `ggplot2`, `patchwork`, `dplyr`, `readr`, `tidyr`, `scales`, `grid`, and `ragg`. The Python reconstruction script uses only the Python 3 standard library.

## Reproducibility limits

Several historical analyses were originally performed interactively and no complete original script was retained. They are therefore not presented here as fully scripted workflows. The deposited scripts reproduce or audit only the components explicitly described above. Third-party data-use restrictions and the absence of some historical execution records prevent a one-command reconstruction from raw data.

## Citation

Please cite the associated article and the archived Zenodo release. The Zenodo DOI will be added after the first public release.

## License

Custom code in this repository is released under the MIT License. Third-party data remain subject to their original terms and licenses.
