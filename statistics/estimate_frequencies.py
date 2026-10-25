"""Reproduce selection and estimate probe-set evidence-pattern frequencies."""
from pathlib import Path
import csv,hashlib,json,math
ROOT=Path(__file__).resolve().parents[1]
def read(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def wilson(k,n,z=1.959963984540054):
    p=k/n;d=1+z*z/n;c=(p+z*z/(2*n))/d
    h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/d
    return c-h,c+h
def main():
    frame=read(ROOT/'legacy_annotation/eligibility_frame.csv')
    assert len(frame)==2209 and len({r['probe_id'] for r in frame})==2209
    for r in frame:
        assert r['sampling_hash']==hashlib.sha256(('ProbeAudit-R6-20261005|'+r['probe_id']).encode()).hexdigest()
    chosen=sorted(frame,key=lambda r:(r['sampling_hash'],r['probe_id']))[:100]
    rows=[r for r in read(ROOT/'evidence/results/probe_metrics.csv') if r['sample_group']=='extension']
    assert len(rows)==100 and {r['probe_id'] for r in chosen}=={r['probe_id'] for r in rows}
    a={r['probe_id'] for r in rows if int(r['alternative_max_ensembl'])>0}
    b={r['probe_id'] for r in rows if int(r['alternative_max_refseq'])>0}
    full={r['probe_id'] for r in rows if max(int(r['alternative_max_ensembl']),int(r['alternative_max_refseq']))==11}
    sets={'either_reference':a|b,'ensembl':a,'refseq':b,'both_references':a&b,'full_alternative':full,'partial_alternative':(a|b)-full}
    assert [len(v) for v in sets.values()]==[15,13,11,9,5,10]
    out=[]
    for name,ids in sets.items():
        lo,hi=wilson(len(ids),100)
        out.append(dict(pattern=name,count=len(ids),denominator=100,estimate=len(ids)/100,lower95=lo,upper95=hi,probe_ids='|'.join(sorted(ids))))
    path=ROOT/'statistics/frequency_estimates.csv'
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
    (ROOT/'statistics/verification.json').write_text(json.dumps(dict(eligible=2209,selected=100,all_hashes_verified=True,selected_ids_reproduced=True,method='95% Wilson, binomial approximation, no finite-population correction',assumption='Hash ordering approximately unrelated to sequence outcomes',unit='probe set'),indent=2),encoding='utf-8')
    print(json.dumps(out,indent=2))
if __name__=='__main__':main()
