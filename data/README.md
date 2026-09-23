# Input data layout

No raw or controlled-access data are included in this repository.

For the full local project, set `CDD_HTN_ROOT` to a directory containing the study data in the layout below:

```text
01_genetics/
  GWAS/HTN/GCST90044351_buildGRCh37.tsv.gz
03_colocalization/
04_replication/
  HTN_BP_MVP/
    raw/MVP_R4.1000G_AGR.Systolic_Mean_INT.EUR.GIA.dbGaP.txt.gz
    MVP_replication_plan_frozen_2026-09-13.tsv
    MVP_SBP_replication_harmonized_2026-09-13.tsv
    MVP_SBP_replication_FINAL_summary_2026-09-13.tsv
05_expression/
  results/GSE34095_frozen_candidate_secondary_validation_2026-09-13.tsv
07_integration/
  S3_bulk_plotData_extraction_FINAL_2026-09-19_v3/
```

The shared-SMR reconstruction scripts expect a separate lightweight audit package with:

```text
01_primary_inputs/
  CDD_eQTLGen_all.smr
  HTN_eQTLGen_all.smr
03_historical_selection/
  batch_coloc_412_all_priors.tsv
  candidate_freeze_2026-09-13.tsv
05_results_recalculated/
```

Users must obtain all third-party datasets from their original providers and comply with the applicable terms. Do not commit MVP controlled-access files, raw GWAS archives, credentials, access tokens, participant-level data, or local download logs.
