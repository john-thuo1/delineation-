#!/usr/bin/env python3
"""Step 2/2: per-class amplitude panels for Figure 6 — healthy (sinus) vs disease, in the
disease-relevant leads. amp = signal at peak minus pre-QRS baseline; per-signal median then
cohort mean over sampled recordings. Writes fig6_amplitudes_by_class.csv."""
import pandas as pd, numpy as np, os, csv
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT,"config","paths.yaml")):
    ROOT = os.path.dirname(ROOT)
DATA = os.path.join(HERE, "..", "data")
FIG  = os.path.join(HERE, "..", "figures")
POS  = os.path.join(DATA,"fig6_amp_positions.csv")
MAN  = os.path.join(ROOT,"ecgdeli_labelling","data","input","medalcare_manifest.csv")
LROW = {"II":1,"V2":7,"aVL":4,"V6":11}; N = 100
pos = pd.read_csv(POS, na_values=["","None"])
path = {r["record_id"]: r["path_raw"] for r in csv.DictReader(open(MAN))}
CACHE = {}
def load(rid):
    if rid in CACHE: return CACHE[rid]
    try:
        M = pd.read_csv(os.path.join(ROOT, path[rid]), header=None).to_numpy()
        if M.shape[0] != 12: M = M.T
    except Exception: M = None
    CACHE[rid] = M; return M

def amp_for(cls, lead, peakcol):
    sub = pos[(pos.disease_class==cls)&(pos.lead==lead)]
    recs = [r for r in sub.record_id.unique() if r in path][:N]
    per_rec=[]
    for rid,g in sub[sub.record_id.isin(recs)].groupby("record_id"):
        M = load(rid)
        if M is None: continue
        sig = M[LROW[lead]]; vals=[]
        for _,row in g.iterrows():
            qon=row["qrs_onset_sample"]; pk=row[peakcol]
            if pd.isna(qon) or pd.isna(pk): continue
            qon=int(qon); base=np.median(sig[max(0,qon-11):qon+1])
            if 0<=int(pk)<len(sig): vals.append(sig[int(pk)]-base)
        if vals: per_rec.append(np.median(vals))
    return float(np.mean(per_rec)) if per_rec else np.nan, len(per_rec)

# (feature, lead, peak column, disease class)  -- the amplitude panels Figure 6 shows
PANELS = [("Qamp","II","q_peak_sample","mi"),
          ("Ramp","V2","r_peak_sample","mi"),
          ("Ramp","II","r_peak_sample","lae"),
          ("Pamp","aVL","p_peak_sample","fam"),
          ("Pamp","V6","p_peak_sample","fam")]
rows=[]
print(f"{'feat':5s} {'lead':4s} {'class':5s} {'healthy':>8s} {'disease':>8s} {'shift':>7s}")
for feat,lead,pk,cls in PANELS:
    h,_ = amp_for("sinus",lead,pk); d,nd = amp_for(cls,lead,pk)
    rows.append({"feature":feat,"lead":lead,"class":cls,
                 "healthy_mV":round(h,3),"disease_mV":round(d,3),"shift_mV":round(d-h,3),"n":nd})
    print(f"{feat:5s} {lead:4s} {cls:5s} {h:8.3f} {d:8.3f} {d-h:+7.3f}")
pd.DataFrame(rows).to_csv(os.path.join(DATA,"fig6_amplitudes_by_class.csv"), index=False)
print("\nsaved fig6_amplitudes_by_class.csv")
