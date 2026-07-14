#!/usr/bin/env python3
"""
build_manual_worklist.py  -  Build require_manual_label.csv: the (record, lead) units that need
manual fiducial labelling, tiered by confidence, from master_labels.csv.

A unit (record, lead) is included when a MAJORITY (>=50%) of its beats are flagged — either by the
rule-based QC (qc_status=='critical') or the cross-method check (xmethod_flag==1). Majority filtering
separates genuine problems from lone transient glitches (each ECG is one beat repeated ~13x).

priority_tier:
  1_critical  - majority of beats are rule-based CRITICAL (definite delineation errors) -> label first
  2_both      - flagged by BOTH the rule-based and cross-method signals
  3_xmethod   - majority flagged only by the independent cross-method check (weaker; often
                low-amplitude P/T ambiguity in aVR/aVL) -> lower priority

Load the signal via path_raw (or Master/extract_signal.py), review the flagged lead, correct the
fiducials from master_labels.csv. Run:  python3 build_manual_worklist.py
"""
import csv, os
from collections import defaultdict, Counter
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT,"config","paths.yaml")):
    ROOT = os.path.dirname(ROOT)
SRC = os.path.join(ROOT, "dataset_curation","data","assembled","master_labels.csv"); IDX = os.path.join(ROOT, "dataset_curation","data","review","signals_index.csv")
OUT = os.path.join(ROOT, "dataset_curation","data","review","require_manual_label.csv")

u = defaultdict(lambda: {"tot":0,"crit":0,"xm":0,"need":0,"dc":"","exb":"","exf":"","exi":""})
for r in csv.DictReader(open(SRC)):
    d = u[(r["record_id"], r["lead"])]; d["tot"] += 1; d["dc"] = r["disease_class"]
    if r["qc_status"] == "critical": d["crit"] += 1
    if r["xmethod_flag"] == "1":     d["xm"]  += 1
    if r["needs_review"] == "1":
        d["need"] += 1
        if not d["exb"]: d["exb"], d["exf"], d["exi"] = r["beat_id"], r["qc_flags"], r["xmethod_issues"]

path = {r["record_id"]: r["path_raw"] for r in csv.DictReader(open(IDX))}
rows = []
for (rid, ld), d in u.items():
    fc, fx, fn = d["crit"]/d["tot"], d["xm"]/d["tot"], d["need"]/d["tot"]
    if fc < 0.5 and fx < 0.5: continue
    tier = "1_critical" if fc >= 0.5 else ("2_both" if (d["crit"] > 0 and d["xm"] > 0) else "3_xmethod")
    rows.append({"record_id":rid,"disease_class":d["dc"],"lead":ld,"path_raw":path.get(rid,""),
        "n_beats":d["tot"],"n_critical":d["crit"],"n_xmethod":d["xm"],"n_needs":d["need"],
        "frac_critical":round(fc,2),"frac_xmethod":round(fx,2),"frac_needs":round(fn,2),
        "priority_tier":tier,"example_beat":d["exb"],"example_qc_flags":d["exf"],
        "example_xmethod_issues":d["exi"]})
rows.sort(key=lambda d: (d["priority_tier"], -d["frac_critical"], -d["frac_needs"]))
with open(OUT, "w", newline="") as f:
    w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
print(f"{len(rows)} units, {len({r['record_id'] for r in rows})} signals")
print(dict(Counter(r["priority_tier"] for r in rows)))
