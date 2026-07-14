"""
crosscheck_neurokit_peak.py  -  NON-WAVELET third-opinion verification of the ECGdeli labels.

This is the sibling of crosscheck_neurokit.py. That script uses NeuroKit2's DWT
(wavelet) delineator; ECGdeli is also wavelet-based, so the two agreeing cannot
rule out a bias the wavelet family shares. This script instead runs NeuroKit2's
"peak" delineator - a DERIVATIVE / local-extrema method from a different algorithm
family (no wavelets) - and compares it to the ECGdeli labels the same way. Three
independent method families (ECGdeli-wavelet, NeuroKit-wavelet, NeuroKit-derivative)
agreeing on a fiducial is far stronger evidence than two wavelet tools agreeing.

The peak method returns fewer fiducials than DWT: the P/Q/S/T PEAKS, the P-onset
and the T-offset. It does NOT return the QRS onset/offset or the T-onset (those
boundaries are wavelet-only in NeuroKit2). That is fine here - the T-offset (the
one biomarker where we differ from the paper) IS returned, so this check
triangulates exactly the fiducial that matters most.

Configuration (kept light on purpose - this is a method-family check, not a full
per-beat consensus, so lead II across every record is sufficient and fast):
    N_PER_CLASS = None   -> every record (all 16,848)
    ALL_LEADS   = False  -> lead II only (set True for all 12 leads, much slower)

Run:
    pip install neurokit2 numpy
    python3 crosscheck_neurokit_peak.py

Outputs (this folder, named distinctly so they never collide with the DWT run):
    consensus_peak_agreement_summary.txt   per-fiducial agreement + median/mean offset
    consensus_peak_disagreements.csv       every beat where the tools differ > tolerance (with lead)
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
N_PER_CLASS = None            # None = ALL records; or an int for a sample per class
ALL_LEADS   = False           # False = lead II only (recommended); True = all 12 leads
SIGNAL_VER  = "raw"
METHOD      = "peak"          # NON-wavelet (derivative / local-extrema) delineator
FS          = 500
LEADS = ["I","II","III","aVR","aVL","aVF","V1","V2","V3","V4","V5","V6"]
LEAD_ROW = {l:i for i,l in enumerate(LEADS)}

# Only the fiducials the peak method actually returns (peaks + P-onset + T-offset)
TOL_MS = {"p_onset":40,"p_peak":25,"q_peak":20,"r_peak":15,"s_peak":20,"t_peak":30,"t_offset":50}
NK_KEY = {"p_onset":"ECG_P_Onsets","p_peak":"ECG_P_Peaks","q_peak":"ECG_Q_Peaks",
          "r_peak":"ECG_R_Peaks","s_peak":"ECG_S_Peaks","t_peak":"ECG_T_Peaks","t_offset":"ECG_T_Offsets"}
FIDS = list(TOL_MS.keys())

leads_run = LEADS if ALL_LEADS else ["II"]
leads_set = set(leads_run)

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
print(f"[peak/non-wavelet] records: {len(ids)}  |  leads: {len(leads_run)}  |  {'FULL' if N_PER_CLASS is None else N_PER_CLASS} per class")

# --- 2. ECGdeli fiducials for those records/leads ---
edeli = defaultdict(list)   # (rid,lead) -> list of beat dicts
for r in csv.DictReader(open(FIDUCIALS)):
    if r["lead"] in leads_set and r["record_id"] in ids:
        edeli[(r["record_id"], r["lead"])].append({k: gi(r.get(k+"_sample","")) for k in (FIDS+["r_peak"])})
print(f"loaded ECGdeli fiducials for {len(edeli)} (record,lead) pairs")

# --- 3+4. run NeuroKit2 (peak method), match beats by R-peak, compare ---
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
        used = set()                          # enforce one-to-one: each NeuroKit beat matches once
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
                if k == "r_peak" or eb.get(k) is None or nb.get(k) is None: continue
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

# --- 5. write outputs (distinct peak-tagged names) ---
with open(os.path.join(XSUM, "consensus_peak_agreement_summary.txt"), "w") as f:
    def line(s=""): print(s); f.write(s+"\n")
    scope = f"{len(leads_run)} lead(s), {done} records"
    line(f"NeuroKit2 (peak / NON-wavelet) vs ECGdeli  -  {scope}\n")
    line("peak method returns peaks + P-onset + T-offset only (no QRS on/off, no T-onset)\n")
    line(f"{'fiducial':11s} {'agree%':>7s} {'n':>8s} {'medianD':>8s} {'meanD':>7s}  (ms)")
    for k in FIDS:
        w, n = agree[k]
        if n == 0: line(f"{k:11s} {'--':>7s} {0:>8d}"); continue
        d = np.array(deltas[k]) * (1000/FS)
        line(f"{k:11s} {100*w/n:6.1f}% {n:>8d} {np.median(d):7.1f} {d.mean():6.1f}")
    line(f"\nbeats with >=1 disagreement: {len(disagreements)}  -> consensus_peak_disagreements.csv")

if disagreements:
    keys = sorted({k for d in disagreements for k in d})
    with open(os.path.join(XINT, "consensus_peak_disagreements.csv"), "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(disagreements)
print("done. See consensus_peak_agreement_summary.txt")
