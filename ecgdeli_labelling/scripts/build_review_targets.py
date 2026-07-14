"""
build_review_targets.py  -  filter the NeuroKit2 disagreements down to the cross-method review set.

The cross-check (crosscheck_neurokit.py) writes consensus_disagreements.csv: every beat where the
two tools differ by more than the per-fiducial tolerance. This script keeps only the beats that
differ by more than a stricter threshold (default 60 ms) on an actual PEAK (P/R/T) or the T-OFFSET
— the landmarks that matter for delineation — and writes them as consensus_review_targets.csv, the
compact table that add_crosscheck_qc.py merges into master_labels.csv as the xmethod_* columns.

Run (after crosscheck_neurokit.py):
    python3 build_review_targets.py            # 60 ms threshold, DWT disagreements
Outputs:
    consensus_review_targets.csv   record_id, disease_class, lead, r_peak, issues
        issues e.g. "p_peak:+68ms;t_offset:-66ms"  (delta = NeuroKit - ECGdeli, in ms)
"""
import os, csv, sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT, "config", "paths.yaml")):
    ROOT = os.path.dirname(ROOT)
XINT = os.path.join(ROOT, "ecgdeli_labelling", "data", "neurokit_crosscheck", "intermediates")
SRC  = os.path.join(XINT, "consensus_disagreements.csv")        # from crosscheck_neurokit.py (DWT)
OUT  = os.path.join(XINT, "consensus_review_targets.csv")

THRESH_MS = 60                                                   # stricter than the per-fiducial tolerance
LANDMARKS = ["p_peak", "r_peak", "t_peak", "t_offset"]           # peaks + the T-offset

if not os.path.isfile(SRC):
    sys.exit(f"missing {SRC} — run crosscheck_neurokit.py first")

n_in = n_out = 0
with open(SRC, newline="") as fin, open(OUT, "w", newline="") as fout:
    r = csv.DictReader(fin)
    w = csv.writer(fout); w.writerow(["record_id", "disease_class", "lead", "r_peak", "issues"])
    for row in r:
        n_in += 1
        issues = []
        for k in LANDMARKS:
            dv = row.get(k + "_delta_ms", "")
            if dv in ("", "None"):
                continue
            d = float(dv)
            if abs(d) > THRESH_MS:
                issues.append(f"{k}:{d:+.0f}ms")
        if issues:
            w.writerow([row["record_id"], row.get("disease_class", ""), row["lead"],
                        row.get("r_peak", ""), ";".join(issues)])
            n_out += 1
print(f"read {n_in:,} disagreement beats -> wrote {n_out:,} review targets (> {THRESH_MS} ms on {'/'.join(LANDMARKS)})")
print(f"  {OUT}")
