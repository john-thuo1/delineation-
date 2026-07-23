#!/usr/bin/env python3
# NOTE: retained ONLY to regenerate the ECGdeli-vs-NeuroKit2 comparison table in the
# dissertation Appendix. NeuroKit2 is not part of the QC pipeline; QC uses the rule-based
# screen and the cross-lead consistency check (see build_crosslead_priority.py).
"""
crosscheck_neurokit.py  -  Independent second-tool verification of the ECGdeli labels.

Runs NeuroKit2 (an independent DWT delineator) and compares its P/QRS/T fiducials
against the ECGdeli labels, beat by beat, matched on the R-peak. Agreement between the
two tools provides convergent supporting evidence, while disagreement beyond tolerance
identifies cases requiring further inspection, neither tool is treated as ground truth.
Disagreeing beats are written to a disagreement list for manual review.

FULL-RUN configuration
    N_PER_CLASS = None   -> every record (all 16,848), not a sample
    ALL_LEADS   = True   -> all 12 leads (else lead II only)
Expect a long run (all records x 12 leads is ~200k delineations). Progress prints
every 200 records so you can see it is alive. It is safe to reduce ALL_LEADS to
False (lead II only) for a much faster but still complete-record verification.

Run
    pip install neurokit2 numpy
    python3 crosscheck_neurokit.py

Outputs (this folder)
    consensus_agreement_summary.txt   per-fiducial agreement + median/mean offset
    consensus_disagreements.csv       every beat where the tools differ > tolerance (with lead)
"""
import os, csv, sys, time
import numpy as np
try:
    import neurokit2 as nk
except ImportError:
    sys.exit("Install NeuroKit2 first:  pip install neurokit2")

# ---------------- CONFIG ----------------
HERE      = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT,"config","paths.yaml")):
    ROOT = os.path.dirname(ROOT)
DATA_ROOT = ROOT                                               # repo root (for signal path_raw)
MANIFEST  = os.path.join(ROOT,"ecgdeli_labelling","data","input","medalcare_manifest.csv")
FIDUCIALS = os.path.join(ROOT,"ecgdeli_labelling","data","primary","medalcare_fiducials_ecgdeli.csv")
XSUM      = os.path.join(ROOT,"ecgdeli_labelling","data","neurokit_crosscheck","summaries")
XINT      = os.path.join(ROOT,"ecgdeli_labelling","data","neurokit_crosscheck","intermediates")
os.makedirs(XSUM, exist_ok=True); os.makedirs(XINT, exist_ok=True)   # never lose a long run at write time
N_PER_CLASS = None            # None = ALL records, or an int for a sample per class
ALL_LEADS   = True            # True = all 12 leads, False = lead II only
SIGNAL_VER  = "raw"
METHOD      = "dwt"           # wavelet (DWT) delineator. The non-wavelet check lives in crosscheck_neurokit_peak.py
FS          = 500
LEADS = ["I","II","III","aVR","aVL","aVF","V1","V2","V3","V4","V5","V6"]
LEAD_ROW = {l:i for i,l in enumerate(LEADS)}

#Heuristic Agreement Thresholds
TOL_MS = {"p_onset":40,"p_peak":25,"p_offset":40,"qrs_onset":25,"q_peak":20,"r_peak":15,
          "s_peak":20,"qrs_offset":25,"t_onset":50,"t_peak":30,"t_offset":50}
leads_run = LEADS if ALL_LEADS else ["II"]
leads_set = set(leads_run)

NK_KEY = {"p_onset":"ECG_P_Onsets","p_peak":"ECG_P_Peaks","p_offset":"ECG_P_Offsets",
          "qrs_onset":"ECG_R_Onsets","q_peak":"ECG_Q_Peaks","r_peak":"ECG_R_Peaks","s_peak":"ECG_S_Peaks",
          "qrs_offset":"ECG_R_Offsets","t_onset":"ECG_T_Onsets","t_peak":"ECG_T_Peaks","t_offset":"ECG_T_Offsets"}
FIDS = list(TOL_MS.keys())

def gi(v): return int(v) if v not in ("","None",None) else None

# --- 1. records to process ---
from collections import defaultdict
byc = defaultdict(list)
for r in csv.DictReader(open(MANIFEST)):
    byc[r["disease_class"]].append(r)
sample = {}
for c, lst in byc.items():
    for r in (lst if N_PER_CLASS is None else lst[:N_PER_CLASS]):
        sample[r["record_id"]] = (os.path.join(DATA_ROOT, r["path_"+SIGNAL_VER]), c)
ids = set(sample)
print(f"records to process: {len(ids)}  |  leads: {len(leads_run)}  |  {'FULL' if N_PER_CLASS is None else N_PER_CLASS} per class")

# --- 2. ECGdeli fiducials for those records/leads ---
edeli = defaultdict(list)   # (rid,lead) -> list of beat dicts
for r in csv.DictReader(open(FIDUCIALS)):
    if r["lead"] in leads_set and r["record_id"] in ids:
        edeli[(r["record_id"], r["lead"])].append({k: gi(r.get(k+"_sample","")) for k in FIDS})
print(f"loaded ECGdeli fiducials for {len(edeli)} (record,lead) pairs")

# --- 3+4. run NeuroKit2, match beats by R-peak, compare ---
deltas = {k: [] for k in FIDS}
agree  = {k: [0, 0] for k in FIDS}
disagreements = []
done = 0; t0 = time.time()
for rid, (path, cls) in sample.items():
    try:
        M = np.loadtxt(path, delimiter=",")
    except Exception:
        continue
    for ld in leads_run:
        if (rid, ld) not in edeli:
            continue
        try:
            ecg = nk.ecg_clean(M[LEAD_ROW[ld], :], sampling_rate=FS)
            _, rp = nk.ecg_peaks(ecg, sampling_rate=FS)
            rpk = rp["ECG_R_Peaks"]
            _, waves = nk.ecg_delineate(ecg, rpeaks=rpk, sampling_rate=FS, method=METHOD)
        except Exception:
            continue
        nkR = list(rpk)
        nk_beats = []
        for i in range(len(nkR)):
            beat = {"r_peak": int(nkR[i]) if nkR[i] == nkR[i] else None}
            for k in FIDS:
                if k == "r_peak": continue
                arr = waves.get(NK_KEY[k], [])
                v = arr[i] if i < len(arr) else np.nan
                beat[k] = int(v) if (v == v and v is not None) else None
            nk_beats.append(beat)
        used = set()                          # enforce one-to-one each NeuroKit beat matches once
        for eb in edeli[(rid, ld)]:
            er = eb["r_peak"]
            if er is None or not nk_beats: continue
            cand = [(abs(nk_beats[i]["r_peak"] - er), i) for i in range(len(nk_beats))
                    if i not in used and nk_beats[i]["r_peak"] is not None]
            if not cand: continue
            dist, bi = min(cand)
            if dist > 30:  # 60 ms
                continue
            used.add(bi); nb = nk_beats[bi]
            rowdis = {"record_id": rid, "disease_class": cls, "lead": ld, "r_peak": er}
            flagged = False
            for k in FIDS:
                if eb[k] is None or nb[k] is None: continue
                d = nb[k] - eb[k]; deltas[k].append(d); agree[k][1] += 1
                if abs(d) * (1000/FS) <= TOL_MS[k]:
                    agree[k][0] += 1
                else:
                    rowdis[k+"_ecgdeli"] = eb[k]; rowdis[k+"_neurokit"] = nb[k]; rowdis[k+"_delta_ms"] = round(d*1000/FS,1)
                    flagged = True
            if flagged:
                disagreements.append(rowdis)
    done += 1
    if done % 200 == 0:
        el = time.time()-t0
        print(f"  {done}/{len(ids)} records  ({el:.0f}s, ~{el/done*len(ids)/60:.0f} min total)")

# --- 5. write outputs ---
# CSE Working Party two-standard-deviation acceptance limits (Eur. Heart J. 1985,6:815-825, Table 2
# as used by Martinez et al. 2004, IEEE TBME 51(4)). Defined for wave BOUNDARIES only, peaks and the
# T-onset have no CSE localization standard. Used here purely as a citable reference bound for the
# inter-tool difference SD -- NOT to flag the review set (that stays on the coarse screening tolerance).
CSE_2SD = {"p_onset": 10.2, "p_offset": 12.7, "qrs_onset": 6.5, "qrs_offset": 11.6, "t_offset": 30.6}
with open(os.path.join(XSUM, "consensus_agreement_summary.txt"), "w") as f:
    def line(s=""): print(s); f.write(s+"\n")
    scope = f"{len(leads_run)} lead(s), {done} records"
    line(f"NeuroKit2 ({METHOD}) vs ECGdeli  -  {scope}")
    line("per-fiducial timing difference (NeuroKit2 minus ECGdeli). 2s_CSE = CSE Working Party 1985")
    line("acceptance limit (boundaries only); SD<=CSE marks whether the inter-tool SD is within it.\n")
    line(f"{'fiducial':11s} {'n':>8s} {'meanD':>7s} {'SD':>7s} {'medianD':>8s} {'2s_CSE':>7s} {'SD<=CSE':>8s}  (ms)")
    for k in FIDS:
        w, n = agree[k]
        if n == 0: line(f"{k:11s} {0:>8d}"); continue
        d = np.array(deltas[k]) * (1000/FS)
        cse = CSE_2SD.get(k)
        cstr = f"{cse:7.1f}" if cse else f"{'--':>7s}"
        wstr = ("yes" if (cse and d.std() <= cse) else "no") if cse else "--"
        line(f"{k:11s} {n:>8d} {d.mean():7.1f} {d.std():7.1f} {np.median(d):8.1f} {cstr} {wstr:>8s}")
    line(f"\nbeats with >=1 disagreement (coarse screening tolerance, for the review set only): "
         f"{len(disagreements)}  -> consensus_disagreements.csv")

if disagreements:
    keys = sorted({k for d in disagreements for k in d})
    with open(os.path.join(XINT, "consensus_disagreements.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(disagreements)
print("done. See consensus_agreement_summary.txt")
