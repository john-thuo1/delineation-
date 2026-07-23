# MedalCare-XL ECG Delineation Project

Pipeline to generate, quality-control, curate, and validate P/QRS/T fiducial labels for
MedalCare-XL simulated 12-lead ECGs, and to prepare a model-ready training set for a global
(12-lead-in → one-delineation-out) ECG-delineation model.

## Layout

```
config/               paths.yaml + pipeline_settings.yaml (single source of truth)
ECG_TOOL/ECGdeli/     the ECGdeli MATLAB delineation toolbox
WP2_largeDataset_Noise/            raw 12x5000 waveforms (source data)
WP2_largeDataset_ParameterFiles/   simulation parameter files (source data)

ecgdeli_labelling/    STEP 1 - delineate with ECGdeli + QC + NeuroKit cross-check
  scripts/            run_ecgdeli_medalcare.m, qc_review_list.py, crosscheck_neurokit(_peak).py, ecgdeli_config_sweep.m
  data/input|primary|qc|neurokit_crosscheck/   manifest, master fiducials, QC, cross-check outputs
  logs/

dataset_curation/     STEP 2 - assemble the training set + QC columns + review worklist
  scripts/            build_master, add_qc_status, add_crosscheck_qc, build_manual_worklist, build_reconciled_global, extract_signal
  data/assembled|global|review|audit_samples/

manual_labelling/     refinement of STEP 2 - manually correct the tier-1 critical pseudo-labels
  MANUAL_LABELLING_PROTOCOL.md · tool/ (browser corrector) · scripts/ (worklist + merge-back) · data/ (plan, batches, corrections)

statistics/           STEP 3 - population validation vs Table 6 and the preprint
  scripts/ data/ figures/

reports/              the labelling report (+ figures)
ml_modelling/         STEP 4 - the delineation model (configs, scripts, data, checkpoints, results, figures)
archive/              superseded artifacts
Dissertation/         the thesis
```

## Pipeline order

1. **ecgdeli_labelling** — `run_ecgdeli_medalcare.m` produces `data/primary/medalcare_fiducials_ecgdeli.csv`
   (one row per record×lead×beat). `qc_review_list.py` tiers each beat critical/minor/clean.
   `crosscheck_neurokit(_peak).py` compare against an independent delineator.
   *This step is NOT rerun routinely — the master fiducials CSV is the fixed source for everything below.*
2. **dataset_curation** — `build_master.py` → `master_labels.csv`, `add_qc_status.py` + `add_crosscheck_qc.py`
   add the QC columns, `build_manual_worklist.py` → the review lists, `build_reconciled_global.py` → the global table.
3. **statistics** — `build_per_signal_stats.py` etc. reproduce the population statistics and compare to Table 6 / the preprint.
4. **ml_modelling** — trains the delineation model on `dataset_curation/data/assembled/master_labels.csv`.

All input/output locations are declared in `config/paths.yaml`, pipeline parameters (QC thresholds,
aggregation, split, reconciliation) in `config/pipeline_settings.yaml`.

## Note

The ECGdeli labelling (the master fiducials CSV) is treated as fixed input. Everything downstream
(curation, statistics) is recomputed from it and is fully reproducible via the scripts above.
