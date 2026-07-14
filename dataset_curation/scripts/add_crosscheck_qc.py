#!/usr/bin/env python3
"""
add_crosscheck_qc.py  -  Add the cross-method (ECGdeli vs NeuroKit) agreement signal to
master_labels.csv, so quality is judged by two independent delineators rather than one.

NeuroKit's independent cross-check disagreements — the beats where NeuroKit differs from ECGdeli
by more than 60 ms on a peak (P/R/T) or the T-offset (consensus_review_targets.csv) — are merged
in, matched on (record_id, lead, r_peak). Adds three columns:
  xmethod_flag    1 if an independent tool disagrees on a P/Q/R/S/T landmark for this beat, else 0
  xmethod_issues  the specific disagreements (e.g. "t_peak:-208ms;t_offset:-66ms"), "" if none
  needs_review    1 if qc_status=='critical' OR xmethod_flag==1  (the comprehensive review pool)

This is QC tier 2: tier 1 (qc_status) is ECGdeli's own structural/physiological plausibility;
tier 2 catches "plausible but wrong" beats where an independent method disagrees. Narrow the
manual-labelling effort from there (e.g. beats that are BOTH critical and xmethod-flagged, or
the genuine per-(record,lead) units).

Run:  python3 add_crosscheck_qc.py     (rewrites master_labels.csv in place via a temp file)
Prereq: consensus_review_targets.csv from the NeuroKit cross-check must exist.
"""
import csv, os
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT,"config","paths.yaml")):
    ROOT = os.path.dirname(ROOT)
SRC  = os.path.join(ROOT, "dataset_curation","data","assembled","master_labels.csv"); TMP = SRC + ".tmp"
RT   = os.path.join(ROOT, "ecgdeli_labelling","data","neurokit_crosscheck","intermediates","consensus_review_targets.csv")

# 1) cross-method landmark disagreements keyed by (record_id, lead, r_peak)
iss = {}
if os.path.isfile(RT):
    for r in csv.DictReader(open(RT, newline="")):
        # strip any stray CR/whitespace so it can never leak into the LF-terminated master file
        iss[(r["record_id"], r["lead"], r["r_peak"])] = (r["issues"] or "").replace("\r", "").strip()
    print(f"cross-method targets loaded: {len(iss):,}")
else:
    print(f"WARNING: {RT} not found — run crosscheck_neurokit_peak.py first. "
          "xmethod columns will be empty.")

# 2) stream master_labels and append the columns
with open(SRC, newline="") as fin, open(TMP, "w", newline="") as fout:
    r = csv.DictReader(fin)
    w = csv.DictWriter(fout, fieldnames=r.fieldnames + ["xmethod_flag", "xmethod_issues", "needs_review"])
    w.writeheader()
    nflag = nreview = 0
    for row in r:
        is_ = iss.get((row["record_id"], row["lead"], row["r_peak_sample"]), "")
        f = 1 if is_ else 0
        needs = 1 if (row.get("qc_status") == "critical" or f == 1) else 0
        row["xmethod_flag"] = f; row["xmethod_issues"] = is_; row["needs_review"] = needs
        nflag += f; nreview += needs
        w.writerow(row)
os.replace(TMP, SRC)
print(f"xmethod_flag=1: {nflag:,} | needs_review=1: {nreview:,}")
