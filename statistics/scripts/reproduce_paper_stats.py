#!/usr/bin/env python3
"""
reproduce_paper_stats.py  -  Reproduce the MedalCare-XL paper's population statistics
(Table 6, Figures 5-6) from our ECGdeli labels, as a consistency / reproducibility check.

This is the Python counterpart of reproduce_paper_stats.m. It
  * aggregates each feature to ONE value per record (median across the record's repeated
    beats) to match the paper's per-ECG extraction
  * compares the healthy timing means and standard deviations, and the amplitude features,
    with the Table 6 "sim" column
  * checks the per-class shifts against Figure 6
  * draws the Figure-5-style density overlays used in the report
        fig_repro_timing_leadII.png / fig_repro_amp_leadII.png  (per record, vs paper)
        fig_perbeat_timing_leadII.png / fig_perbeat_amp_leadII.png (all beats, no overlay)

NOTE consistency check, not accuracy -- the paper's features were also ECGdeli-derived.

Run pip install pandas numpy matplotlib, python3 reproduce_paper_stats.py
"""
import pandas as pd, numpy as np, os, csv
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = HERE
while ROOT != os.path.dirname(ROOT) and not os.path.isfile(os.path.join(ROOT,"config","paths.yaml")):
    ROOT = os.path.dirname(ROOT)
DATA = os.path.join(HERE, "..", "data")
FIG  = os.path.join(HERE, "..", "figures")
FID  = os.path.join(ROOT,"ecgdeli_labelling","data","primary","medalcare_fiducials_ecgdeli.csv")
MAN  = os.path.join(ROOT,"ecgdeli_labelling","data","input","medalcare_manifest.csv")
MS   = 2.0                                    # 500 Hz -> 2 ms/sample
LEADS = ["I","II","III","aVR","aVL","aVF","V1","V2","V3","V4","V5","V6"]; LROW={l:i for i,l in enumerate(LEADS)}
N_AMP = 200                                   # records sampled for amplitude stats

# Table 6 "sim" per lead (mu,sd) for [Pdur,QRSdur,Tdur,PQint,QTint,RRint] then amplitudes
T6T = {"I":[(124.06,18.37),(131.31,17.56),(178.12,31.08),(128.07,28.06),(310.54,33.23),(758.15,54.97)],
 "II":[(128.09,14.00),(126.10,13.73),(182.33,25.94),(127.18,22.97),(317.08,23.61),(758.02,54.41)],
 "III":[(164.52,24.12),(126.80,14.32),(183.16,28.07),(171.88,30.36),(306.94,26.02),(757.99,54.41)],
 "aVR":[(127.42,16.05),(128.81,15.48),(179.38,24.98),(126.19,23.97),(318.73,22.35),(757.97,54.83)],
 "aVL":[(154.50,25.17),(128.62,16.29),(182.76,28.06),(169.37,34.88),(299.19,24.85),(758.08,54.46)],
 "aVF":[(141.00,19.09),(125.02,12.77),(184.05,26.18),(142.60,24.74),(310.43,24.82),(758.06,54.33)],
 "V1":[(140.78,29.86),(129.05,16.66),(180.72,29.58),(160.57,38.39),(303.89,28.46),(758.06,54.63)],
 "V2":[(155.93,32.59),(136.12,14.12),(176.15,28.01),(181.65,47.68),(287.71,21.43),(757.90,54.24)],
 "V3":[(154.23,27.31),(132.80,14.15),(179.01,26.69),(174.20,39.24),(285.50,21.16),(758.01,54.27)],
 "V4":[(140.60,20.53),(127.32,13.68),(180.05,26.44),(148.55,25.18),(290.09,23.79),(758.03,54.05)],
 "V5":[(128.69,15.28),(123.89,13.28),(177.39,23.83),(126.73,23.52),(310.74,23.04),(758.12,54.32)],
 "V6":[(123.44,11.93),(126.55,13.99),(174.63,20.55),(118.63,21.76),(320.51,20.90),(758.06,54.20)]}
T6A = {"I":[0.09,0.02,0.37,-0.05,-0.03],"II":[0.09,0.06,0.59,-0.20,0.49],"III":[0.03,0.08,0.11,-0.10,0.51],
 "aVR":[-0.09,-0.03,-0.49,0.13,-0.24],"aVL":[0.03,-0.03,0.14,0.02,-0.27],"aVF":[0.05,0.08,0.34,-0.15,0.49],
 "V1":[-0.05,0.13,-0.48,-0.05,0.65],"V2":[-0.02,0.21,-1.55,-0.02,1.47],"V3":[0.06,0.21,-0.95,-0.05,1.09],
 "V4":[0.06,0.21,-0.36,-0.11,0.79],"V5":[0.09,0.08,0.54,-0.23,0.46],"V6":[0.11,0.03,0.72,-0.17,0.30]}
TIMING = ["Pdur","QRSdur","Tdur","PQint","QTint","RRint"]
AMP = ["Pamp","Qamp","Ramp","Samp","Tamp"]; PK={"Pamp":"p_peak_sample","Qamp":"q_peak_sample",
       "Ramp":"r_peak_sample","Samp":"s_peak_sample","Tamp":"t_peak_sample"}

def kde(data, xs):
    d = np.asarray(data, float); d = d[np.isfinite(d)]; sd = d.std()
    h = 1.06*sd*len(d)**(-0.2) if sd>0 else 1.0
    u = (xs[:,None]-d[None,:])/h
    return np.exp(-0.5*u*u).sum(1)/(len(d)*h*np.sqrt(2*np.pi))

# ---- load labels + derive timing features ----
c = ["disease_class","lead","record_id","beat_id","p_onset_sample","p_offset_sample","qrs_onset_sample",
     "qrs_offset_sample","t_onset_sample","t_offset_sample","q_peak_sample","r_peak_sample",
     "p_peak_sample","s_peak_sample","t_peak_sample"]
df = pd.read_csv(FID, usecols=c, na_values=["","None"]).sort_values(["record_id","lead","beat_id"])
for a,b,nm in [("p_offset_sample","p_onset_sample","Pdur"),("qrs_offset_sample","qrs_onset_sample","QRSdur"),
               ("t_offset_sample","t_onset_sample","Tdur"),("qrs_onset_sample","p_onset_sample","PQint"),
               ("t_offset_sample","qrs_onset_sample","QTint")]:
    df[nm] = (df[a]-df[b])*MS
df["RRint"] = df.groupby(["record_id","lead"])["r_peak_sample"].diff()*MS
for k in TIMING: df.loc[df[k] <= 0, k] = np.nan

# per-record aggregation (one median value per record-lead)
rec = df.groupby(["disease_class","record_id","lead"])[TIMING].median().reset_index()

# ---- timing healthy sinus vs Table 6 ----
print("=== TIMING (sinus, per record) vs Table 6 (sim), averaged over 12 leads ===")
sin = rec[rec.disease_class=="sinus"]
for i,f in enumerate(TIMING):
    dmu=[]; sd=[]
    for ld in LEADS:
        v = sin[sin.lead==ld][f]
        dmu.append(v.mean()-T6T[ld][i][0]); sd.append(v.std())
    print(f"  {f:7s} mean|Δμ|={np.mean(np.abs(dmu)):5.1f}  ourσ={np.mean(sd):5.1f}  paperσ={np.mean([T6T[ld][i][1] for ld in LEADS]):5.1f}")

# ---- per-class shift vs sinus ----
print("\n=== PER-CLASS shift vs sinus (ms, per record) ===")
base = {f: rec[rec.disease_class=='sinus'].groupby('lead')[f].mean().mean() for f in TIMING}
for cls in ["avblock","lbbb","rbbb","iab","lae","fam","mi"]:
    sub = rec[rec.disease_class==cls]
    row = {f: sub.groupby('lead')[f].mean().mean()-base[f] for f in ["Pdur","QRSdur","PQint","QTint"]}
    print(f"  {cls:8s} ΔPdur {row['Pdur']:+5.1f}  ΔQRS {row['QRSdur']:+5.1f}  ΔPQ {row['PQint']:+5.1f}  ΔQT {row['QTint']:+5.1f}")

# ---- amplitudes (sampled) ----
man = {r["record_id"]: r["path_raw"] for r in csv.DictReader(open(MAN)) if r["disease_class"]=="sinus"}
d2 = df[df.lead=="II"]; recs = [r for r in d2.record_id.unique() if r in man][:N_AMP]
amp_beat = {f:[] for f in AMP}; amp_rec = {f:[] for f in AMP}
for rid, g in d2[d2.record_id.isin(recs)].groupby("record_id"):
    try: M = pd.read_csv(os.path.join(ROOT, man[rid]), header=None).to_numpy()
    except Exception: continue
    if M.shape[0] != 12: M = M.T
    sig = M[LROW["II"]]; tmp = {f:[] for f in AMP}
    for _, row in g.iterrows():
        qon = row["qrs_onset_sample"]
        if pd.isna(qon): continue
        qon = int(qon); base_mV = np.median(sig[max(0,qon-11):qon+1])
        for f in AMP:
            v = row[PK[f]]
            if pd.notna(v) and 0 <= int(v) < len(sig):
                amp_beat[f].append(sig[int(v)]-base_mV); tmp[f].append(sig[int(v)]-base_mV)
    for f in AMP:
        if tmp[f]: amp_rec[f].append(np.median(tmp[f]))

# ---- figures ----
def fig(data_by_feat, T6, feats, titles, unit, fname, suptitle, overlay=True):
    n = len(feats); fig, ax = plt.subplots(2, 3, figsize=(14,7))
    for a,(f,ti) in zip(ax.ravel(), zip(feats,titles)):
        v = np.asarray(data_by_feat[f], float); v = v[np.isfinite(v)]
        if overlay:
            mu,sd = T6[f]; lo,hi = min(v.min(),mu-4*sd), max(np.quantile(v,0.99),mu+4*sd)
        else:
            lo,hi = np.quantile(v,0.005), np.quantile(v,0.995)
        xs = np.linspace(lo,hi,300)
        a.fill_between(xs, kde(v,xs), color="#4C78A8", alpha=0.5, label="Our labels")
        a.axvline(np.nanmean(v), color="#4C78A8", ls="--", lw=1)
        if overlay:
            a.plot(xs, np.exp(-0.5*((xs-mu)/sd)**2)/(sd*np.sqrt(2*np.pi)), color="#D62728", lw=2, label="Paper Table 6")
            a.axvline(mu, color="#D62728", ls="--", lw=1); a.legend(fontsize=7)
        a.set_title(ti, fontsize=10); a.set_xlabel(unit); a.set_yticks([])
    if n < 6: ax.ravel()[5].axis("off")
    fig.suptitle(suptitle, fontsize=12); plt.tight_layout()
    plt.savefig(os.path.join(HERE,fname), dpi=130, bbox_inches="tight"); plt.close()

recII = {f: sin[sin.lead=="II"][f].dropna().values for f in TIMING}
beatII = {f: d2[f].dropna().values for f in TIMING}
T6II  = {f:(T6T["II"][i][0],T6T["II"][i][1]) for i,f in enumerate(TIMING)}
T6AII = {f:(T6A["II"][i], [0.05,0.17,0.51,0.21,0.22][i]) for i,f in enumerate(AMP)}
tt = ["P duration","QRS duration","T duration","PQ interval","QT interval","RR interval"]
aa = ["P amplitude","Q amplitude","R amplitude","S amplitude","T amplitude"]
fig(recII,  T6II,  TIMING, [t+" (lead II)" for t in tt], "ms", "fig_repro_timing_leadII.png",
    "Healthy sinus timing (one value per record): our labels vs Table 6", overlay=True)
fig(amp_rec,T6AII, AMP,    [t+" (lead II)" for t in aa], "mV", "fig_repro_amp_leadII.png",
    "Healthy sinus amplitude (one value per record): our labels vs Table 6", overlay=True)
fig(beatII, T6II,  TIMING, [t+" (lead II)" for t in tt], "ms", "fig_perbeat_timing_leadII.png",
    "Healthy sinus timing - all beats pooled (no paper overlay)", overlay=False)
fig(amp_beat,T6AII,AMP,    [t+" (lead II)" for t in aa], "mV", "fig_perbeat_amp_leadII.png",
    "Healthy sinus amplitude - all beats pooled (no paper overlay)", overlay=False)
print("\nsaved 4 figures (fig_repro_*, fig_perbeat_*).")
