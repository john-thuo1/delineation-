# Per-beat QC classifier (critical/minor/clean) — reproduces qc_review_list.py rules.
# Hardcoded column indices for the 64-col master schema fs_hz=6, fiducials p_onset..t_offset
# = cols 13..23 (in order), qrs_present=42. Reads DATA rows only, appends ",qc_status,qc_flags".
{
 fs=$6+0; if(fs<=0)fs=500; msf=1000.0/fs;
 struct=0;gross=0;mild=0;prevv="";prevp=0;
 for(i=1;i<=11;i++){ c=12+i; val=$c;
   if(val==""||val=="None")continue; val=val+0;
   if(prevv!=""){ if(prevv>val){ m=prevv-val;
     if(prevp>=4&&prevp<=8&&i>=4&&i<=8&&m>3)struct=1;
     else if(m>20){if(m>gross)gross=m}
     else if(m>3){if(m>mild)mild=m} } }
   prevv=val;prevp=i }
 reasons="";
 if(struct)reasons="QRS_structural";
 if(gross)reasons=reasons (reasons==""?"":";") "gross_boundary:" gross "smp";
 qon=$16;qoff=$20;pon=$13;toff=$23;rpk=$18;qp=$42;
 if(qon!=""&&qon!="None"&&qoff!=""&&qoff!="None"){d=(qoff-qon)*msf;if(d<40||d>200)reasons=reasons (reasons==""?"":";") sprintf("QRSdur=%.0fms",d)}
 if(qon!=""&&qon!="None"&&toff!=""&&toff!="None"){d=(toff-qon)*msf;if(d<250||d>700)reasons=reasons (reasons==""?"":";") sprintf("QT=%.0fms",d)}
 if(pon!=""&&pon!="None"&&qon!=""&&qon!="None"){d=(qon-pon)*msf;if(d<80||d>400)reasons=reasons (reasons==""?"":";") sprintf("PR=%.0fms",d)}
 if(qp!="1"||rpk==""||rpk=="None")reasons=reasons (reasons==""?"":";") "R_missing";
 if(reasons!=""){status="critical";flags=reasons}else if(mild>0){status="minor";flags=""}else{status="clean";flags=""}
 print $0","status","flags
}
