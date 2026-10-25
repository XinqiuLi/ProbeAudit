"""ProbeAudit 0.1.0: reference sequence audit and paired internal benchmarks.
Research software; no clinical use or external validation claim.
"""
import argparse,csv,gzip,hashlib,json,platform,re,sys
from pathlib import Path
import numpy as np,pandas as pd
VERSION='0.1.0'
def rc(s):return s.translate(str.maketrans('ACGTN','TGCAN'))[::-1]
def fasta(path):
 op=gzip.open if str(path).endswith('.gz') else open
 with op(path,'rt') as f:
  head=None;seq=[]
  for line in f:
   if line.startswith('>'):
    if head is not None:yield head,''.join(seq).upper()
    head=line[1:].strip();seq=[]
   else:seq.append(line.strip())
  if head is not None:yield head,''.join(seq).upper()
def fingerprint(path):
 h=hashlib.sha256()
 with open(path,'rb') as f:
  while b:=f.read(1048576):h.update(b)
 return {'path':str(path),'sha256':h.hexdigest(),'bytes':Path(path).stat().st_size}
def provenance(out,args,files):
 info={'tool':'ProbeAudit','version':VERSION,'python':platform.python_version(),'arguments':vars(args),'inputs':[fingerprint(f) for f in files]}
 (out/'provenance.json').write_text(json.dumps(info,indent=2),encoding='utf-8')
def make_matcher(sequences,mismatches):
 import ahocorasick
 patterns={}
 for s in set(sequences):
  for strand,q in [('+',s),('-',rc(s))]:
   seeds=[(0,q)] if mismatches==0 else [(0,q[:len(q)//2]),(len(q)//2,q[len(q)//2:])]
   for offset,seed in seeds:patterns.setdefault(seed,[]).append((s,strand,offset,q))
 ac=ahocorasick.Automaton()
 for seed,items in patterns.items():ac.add_word(seed,(len(seed),items))
 ac.make_automaton();return ac
def matches(ac,sequence,limit):
 seen=set()
 for end,(length,items) in ac.iter(sequence):
  for original,strand,offset,query in items:
   start=end-length+1-offset
   if start<0 or start+len(query)>len(sequence):continue
   key=(original,strand,start)
   if key in seen:continue
   seen.add(key);target=sequence[start:start+len(query)]
   if set(target)-set('ACGT'):continue
   distance=sum(a!=b for a,b in zip(query,target))
   if distance<=limit:yield original,strand,start,distance
def sequence_audit(a):
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 probes=pd.read_csv(a.probes,dtype=str).fillna('')
 required={'probe_set','oligo_id','sequence','expected_symbol'}
 if not required<=set(probes):raise ValueError('Probe CSV requires '+str(sorted(required)))
 if probes.oligo_id.duplicated().any():raise ValueError('Duplicate oligo_id')
 probes.sequence=probes.sequence.str.upper()
 if any(not s or set(s)-set('ACGT') for s in probes.sequence):raise ValueError('Only nonempty ACGT oligos supported')
 ac=make_matcher(probes.sequence,a.mismatches);seq_to_rows={s:g.to_dict('records') for s,g in probes.groupby('sequence')}
 fields=['probe_set','oligo_id','sequence','expected_symbol','transcript','gene_id','gene_symbol','biotype','reference_location','orientation','transcript_start_1based','mismatches']
 genes={};isoforms={};alltrans={};ntrans=0
 with open(out/'sequence_hits.csv','w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
  for file in a.fasta:
   for header,seq in fasta(file):
    meta=dict(re.findall(r'(gene|gene_symbol|gene_biotype|transcript_biotype):([^ ]+)',header));tid=header.split()[0];gid=meta.get('gene','');symbol=meta.get('gene_symbol',gid);bio=meta.get('transcript_biotype',meta.get('gene_biotype',''))
    if not gid:raise ValueError('FASTA requires Ensembl gene: metadata')
    loc=next((x for x in header.split() if x.startswith(('chromosome:','scaffold:','supercontig:'))),'')
    ntrans+=1;alltrans.setdefault(symbol,set()).add(tid)
    for orig,strand,start,dist in matches(ac,seq,a.mismatches):
     for row in seq_to_rows[orig]:
      w.writerow({**{k:row[k] for k in ['probe_set','oligo_id','sequence','expected_symbol']},'transcript':tid,'gene_id':gid,'gene_symbol':symbol,'biotype':bio,'reference_location':loc,'orientation':strand,'transcript_start_1based':start+1,'mismatches':dist})
      for cutoff in range(dist,a.mismatches+1):
       key=(row['oligo_id'],cutoff,strand);genes.setdefault(key,set()).add((gid,symbol));isoforms.setdefault(key,set()).add((tid,symbol))
    if ntrans%100000==0:print('Scanned transcripts',ntrans,flush=True)
 records=[]
 for row in probes.to_dict('records'):
  intended=set(row['expected_symbol'].split('|'))
  for cutoff in range(a.mismatches+1):
   g=genes.get((row['oligo_id'],cutoff,'+'),set());t=isoforms.get((row['oligo_id'],cutoff,'+'),set());anti=genes.get((row['oligo_id'],cutoff,'-'),set())
   records.append({**row,'max_mismatches':cutoff,'sense_gene_count':len(set(x[0] for x in g)),'sense_symbols':'|'.join(sorted(set(x[1] for x in g))),'intended_sense_hit':any(sym in intended for _,sym in g),'offtarget_sense_hit':any(sym not in intended for _,sym in g),'sense_transcript_count':len(t),'intended_transcripts_hit':'|'.join(sorted(tid for tid,sym in t if sym in intended)),'antisense_gene_count':len(set(x[0] for x in anti))})
 detail=pd.DataFrame(records);detail.to_csv(out/'oligo_summary.csv',index=False)
 summaries=[]
 for (ps,cut),d in detail.groupby(['probe_set','max_mismatches']):
  intended=set(d.expected_symbol.iloc[0].split('|'));t=set('|'.join(d.intended_transcripts_hit).split('|'))-{''};den=set().union(*(alltrans.get(s,set()) for s in intended))
  summaries.append({'probe_set':ps,'expected_symbol':d.expected_symbol.iloc[0],'max_mismatches':cut,'n_oligos':len(d),'n_sense_matched':int((d.sense_gene_count>0).sum()),'n_intended_sense':int(d.intended_sense_hit.sum()),'n_offtarget_sense':int(d.offtarget_sense_hit.sum()),'n_unique_intended_sense':int(((d.sense_gene_count==1)&d.intended_sense_hit).sum()),'n_antisense_matched':int((d.antisense_gene_count>0).sum()),'intended_isoforms_hit':len(t),'intended_isoforms_in_reference':len(den),'all_sense_symbols':'|'.join(sorted(set('|'.join(d.sense_symbols).split('|'))-{''}))})
 pd.DataFrame(summaries).to_csv(out/'probe_set_summary.csv',index=False)
 (out/'reference_summary.json').write_text(json.dumps({'transcripts_scanned':ntrans,'reference_gene_symbols':len(alltrans),'interpretation':'Sense is exact/specified mismatch similarity to the transcript orientation, not experimentally established hybridization. Antisense matches reported separately. Unique is within supplied reference only; no whole-genome uniqueness, indels, variants, thermodynamics, or intensity deconvolution.'},indent=2))
 provenance(out,a,[a.probes]+a.fasta)
def genome_audit(a):
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True);probes=pd.read_csv(a.probes,dtype=str);ac=make_matcher(probes.sequence,0);seqrows={s:g.to_dict('records') for s,g in probes.groupby('sequence')};counts={};total=0
 with open(out/'genome_exact_hits.csv','w',newline='',encoding='utf-8') as f:
  w=csv.DictWriter(f,fieldnames=['oligo_id','probe_set','contig','start_1based','end_1based','strand']);w.writeheader()
  for header,seq in fasta(a.fasta):
   contig=header.split()[0];total+=len(seq)
   for orig,strand,start,dist in matches(ac,seq,0):
    for row in seqrows[orig]:
     w.writerow(dict(oligo_id=row['oligo_id'],probe_set=row['probe_set'],contig=contig,start_1based=start+1,end_1based=start+len(orig),strand=strand));counts.setdefault(row['oligo_id'],set()).add((contig,start,start+len(orig)))
   print('Genome scanned',contig,len(seq),flush=True)
 probes['exact_genomic_sites']=probes.oligo_id.map({k:len(v) for k,v in counts.items()}).fillna(0).astype(int);probes.to_csv(out/'genome_oligo_summary.csv',index=False)
 summary=probes.groupby('probe_set').exact_genomic_sites.agg(n_oligos='count',n_unmatched=lambda x:(x==0).sum(),n_unique_site=lambda x:(x==1).sum(),n_multiple_sites=lambda x:(x>1).sum());summary.to_csv(out/'genome_probe_summary.csv')
 (out/'scope.json').write_text(json.dumps({'bases_scanned':total,'scope':'Exact ungapped matches on both strands of supplied primary assembly. Does not include alternate assemblies, variants or mismatches; splice-junction oligos may lack a contiguous genomic hit.'},indent=2));provenance(out,a,[a.probes,a.fasta])
def fold_design(expr,members,candidates,train):
 # Scaling uses training observations only, including gene-to-program summaries.
 scores=[]
 for _,group in members.groupby('program',sort=True):
  x=expr.loc[group.probe].T.to_numpy(float);mean=x[train].mean(0);sd=x[train].std(0,ddof=1);sd=np.where(sd>0,sd,1)
  scores.append(((x-mean)/sd).mean(1))
 programs=np.column_stack(scores);probes=expr.loc[candidates.probe].T.to_numpy(float)
 return programs,probes
def benchmark(a):
 from sklearn.model_selection import StratifiedKFold
 from sklearn.preprocessing import StandardScaler
 from sklearn.linear_model import LogisticRegression
 from sklearn.metrics import roc_auc_score,brier_score_loss
 from scipy.stats import t
 out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
 expr=pd.read_csv(a.expression,index_col=0);lab=pd.read_csv(a.labels,dtype={'sample':str});cand=pd.read_csv(a.candidates,dtype=str).fillna('');mem=pd.read_csv(a.programs,dtype=str)
 if not expr.index.is_unique or not expr.columns.is_unique or lab['sample'].duplicated().any() or cand.probe.duplicated().any() or mem.duplicated(['program','probe']).any():raise ValueError('Duplicate IDs or program membership')
 if not {'sample','outcome'}<=set(lab) or not {'probe','symbol','annotation'}<=set(cand) or not {'program','probe'}<=set(mem):raise ValueError('Invalid input schema')
 if set(lab['sample'])!=set(expr.columns):raise ValueError('Expression/label sample mismatch')
 y=lab.set_index('sample').loc[expr.columns,'outcome'].to_numpy();
 if set(y)!={0,1} or min(np.bincount(y.astype(int)))<a.folds:raise ValueError('Binary outcomes and enough members per fold required')
 if not np.isfinite(expr.to_numpy(float)).all():raise ValueError('Non-finite expression')
 if a.repeats<2 or a.folds<2 or a.C<=0:raise ValueError('Invalid benchmark parameters')
 if not set(cand.probe)|set(mem.probe)<=set(expr.index):raise ValueError('Requested probe missing')
 single=(~cand.annotation.str.contains('///',regex=False)&~cand.annotation.isin(['','---'])).to_numpy()
 if not single.any():raise ValueError('No single-symbol candidates')
 pred=[];foldscores=[];coefs=[]
 for repeat in range(a.repeats):
  seed=a.seed+repeat
  for fold,(train,test) in enumerate(StratifiedKFold(a.folds,shuffle=True,random_state=seed).split(np.zeros(len(y)),y)):
   pg,pr=fold_design(expr,mem,cand,train)
   designs={'programs':pg,'probes_single':pr[:,single],'probes_all':pr,'programs_single':np.column_stack([pg,pr[:,single]]),'programs_all':np.column_stack([pg,pr])}
   names={'programs':sorted(mem.program.unique()),'probes_single':cand.symbol[single].tolist(),'probes_all':cand.symbol.tolist()};names['programs_single']=names['programs']+names['probes_single'];names['programs_all']=names['programs']+names['probes_all']
   for model,x in designs.items():
    scaler=StandardScaler().fit(x[train]);fit=LogisticRegression(C=a.C,solver='lbfgs',max_iter=2000).fit(scaler.transform(x[train]),y[train]);prob=fit.predict_proba(scaler.transform(x[test]))[:,1]
    if fit.n_iter_[0]>=2000:raise RuntimeError('Logistic model did not converge')
    foldscores.append(dict(seed=seed,fold=fold,model=model,n_train=len(train),n_test=len(test),auc=roc_auc_score(y[test],prob),brier=brier_score_loss(y[test],prob)))
    for idx,prb in zip(test,prob):pred.append(dict(seed=seed,fold=fold,model=model,sample=expr.columns[idx],outcome=int(y[idx]),probability=prb))
    for feature,coef in zip(names[model],fit.coef_[0]):coefs.append(dict(seed=seed,fold=fold,model=model,feature=feature,coefficient=coef))
 p=pd.DataFrame(pred);f=pd.DataFrame(foldscores);p.to_csv(out/'held_out_predictions.csv',index=False);f.to_csv(out/'fold_scores.csv',index=False);pd.DataFrame(coefs).to_csv(out/'coefficients.csv',index=False)
 rep=[]
 for (seed,model),d in p.groupby(['seed','model']):rep.append(dict(seed=seed,model=model,auc=roc_auc_score(d.outcome,d.probability),brier=brier_score_loss(d.outcome,d.probability)))
 r=pd.DataFrame(rep);r.to_csv(out/'repetition_scores.csv',index=False);r.groupby('model')[['auc','brier']].agg(['median','min','max']).to_csv(out/'performance_summary.csv')
 pair=[];interval=[]
 for model in ['programs_single','programs_all']:
  a1=f[f.model==model].merge(f[f.model=='programs'],on=['seed','fold'],suffixes=('_added','_base'),validate='one_to_one')
  for metric in ['auc','brier']:
   v=a1[metric+'_added']-a1[metric+'_base'];j=len(v);ratio=(a1.n_test_added/a1.n_train_added).mean();se=np.sqrt((1/j+ratio)*v.var(ddof=1));delta=v.mean();half=t.ppf(.975,j-1)*se
   interval.append(dict(model=model,metric=metric,mean_fold_difference=delta,ci_low=delta-half,ci_high=delta+half,se=se,folds=j,mean_test_train_ratio=ratio,interpretation='Approximate corrected resampling; not selection-adjusted or external validation'))
  b=r[r.model==model].merge(r[r.model=='programs'],on='seed',suffixes=('_added','_base'))
  for row in b.to_dict('records'):pair.append(dict(seed=row['seed'],model=model,delta_auc=row['auc_added']-row['auc_base'],delta_brier=row['brier_added']-row['brier_base']))
 pd.DataFrame(pair).to_csv(out/'paired_repetition_changes.csv',index=False);pd.DataFrame(interval).to_csv(out/'corrected_intervals.csv',index=False)
 cand.assign(single_symbol=single).to_csv(out/'annotation_risk.csv',index=False)
 provenance(out,a,[a.expression,a.labels,a.candidates,a.programs])
 print('Completed',len(p),'predictions',flush=True)
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',action='version',version=VERSION);sub=p.add_subparsers(dest='command',required=True)
 s=sub.add_parser('sequence');s.add_argument('--probes',required=True);s.add_argument('--fasta',nargs='+',required=True);s.add_argument('--mismatches',type=int,choices=[0,1],default=1);s.add_argument('--out',required=True)
 g=sub.add_parser('genome');g.add_argument('--probes',required=True);g.add_argument('--fasta',required=True);g.add_argument('--out',required=True)
 b=sub.add_parser('benchmark')
 for name in ['expression','labels','candidates','programs','out']:b.add_argument('--'+name,required=True)
 b.add_argument('--repeats',type=int,default=20);b.add_argument('--folds',type=int,default=5);b.add_argument('--seed',type=int,default=20260907);b.add_argument('--C',type=float,default=1.)
 a=p.parse_args();{'sequence':sequence_audit,'genome':genome_audit,'benchmark':benchmark}[a.command](a)
if __name__=='__main__':main()
