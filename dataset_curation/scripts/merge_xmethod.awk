# Merge cross-method (NeuroKit non-wavelet peak) landmark-disagreement flags into master_labels.
# File 1 = consensus_review_targets.csv (record_id,disease_class,lead,r_peak,issues).
# File 2 = master_labels data rows. Hardcoded cols: record_id=$1, lead=$9, r_peak_sample=$18, qc_status=$65.
NR==FNR { if(FNR>1) iss[$1 SUBSEP $3 SUBSEP $4]=$5; next }
{
  k=$1 SUBSEP $9 SUBSEP $18
  if(k in iss){ f=1; is=iss[k] } else { f=0; is="" }
  needs=( $65=="critical" || f==1 )?1:0
  print $0","f","is","needs
}
