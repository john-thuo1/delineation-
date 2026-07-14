# ml_modelling

The ECG-delineation model. Global design: 12 leads in → one delineation per beat (sequence
labelling / 1-D semantic segmentation), trained on the curated pseudo-labels.

## Inputs
- Labels: `../dataset_curation/data/assembled/master_labels.csv` (per record×lead×beat fiducials +
  `split` ∈ {train, val} + QC columns `qc_status`, `xmethod_flag`, `needs_review`).
- Waveforms: `../dataset_curation/data/review/signals_index.csv` maps each `record_id` to its
  raw 12×5000 signal under `../WP2_largeDataset_Noise/`.
- Testing is on the separate internal (MonoAlg3D) dataset, so MedalCare-XL is all training/val.

## Folders
```
configs/       model + training configs (yaml)
scripts/       data loaders, model, train/eval
data/          cached tensors / splits (derived from master_labels.csv)
checkpoints/   saved model weights
results/       metrics (point error + biomarker error)
figures/       plots
```

## Notes
- Keep the val split by record (already encoded in `master_labels.csv`), never by beat, to avoid
  near-duplicate-beat leakage.
- Consider filtering/weighting by `qc_status` / `needs_review`, and reserving the manually-labelled
  `require_manual_label_priority.csv` units for a clean evaluation slice.
