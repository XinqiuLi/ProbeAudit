"""Run with Python 3 in the staged repository folder; uses only the standard library."""
from pathlib import Path
import csv,json
P=Path(__file__).resolve().parent
def read(name):
 return list(csv.DictReader((P/name).open(encoding='utf-8-sig')))
rows=read('Supplementary_Annotation_Record.csv')
raw={n:{(r['probe_set'],int(r['max_mismatches'])):r for r in read(n+'_probe_summary.csv')} for n in ['historical_14','historical_ref','extension_100','extension_ref']}
assert len(rows)==114 and len({r['probe_id'] for r in rows})==114
out={}
for ref in ['ensembl','refseq']:
 for r in rows:
  k=r['sample_group'] if ref=='ensembl' else ('historical_ref' if r['sample_group']=='historical_14' else 'extension_ref')
  x=raw[k][r['probe_id'],0];y=raw[k][r['probe_id'],1]
  cond={
   'platform_multiple_symbols':int('///' in r['previous_annotation']),
   'different_gene_exact':int(int(x['n_offtarget_sense'])>0),
   'incomplete_intended_exact':int(int(x['n_intended_sense'])<11),
   'different_gene_at_most_one_substitution':int(int(y['n_offtarget_sense'])>0),
   'not_all_genomic_sites_unique':int(int(r['genome_unique_oligos'])<11)}
  for name,val in cond.items():assert int(r[name+'_'+ref])==val,(r['probe_id'],ref,name)
  assert r['intended_exact_'+ref]==x['n_intended_sense']+'/11'
  assert r['other_gene_exact_'+ref]==x['n_offtarget_sense']+'/11'
  category='A' if any(cond[n] for n in ['platform_multiple_symbols','different_gene_exact','incomplete_intended_exact']) else 'B' if any(cond[n] for n in ['different_gene_at_most_one_substitution','not_all_genomic_sites_unique']) else 'C'
  assert r['category_'+ref].startswith(category),(r['probe_id'],ref,'category')
  represented=ref=='ensembl' or r['target_represented_refseq']=='True'
  assert int(r['coverage_only_independent_'+ref])==int(represented and cond['incomplete_intended_exact'] and sum(cond.values())==1)
  assert int(r['any_condition_excluding_incomplete_'+ref])==int(represented and any(v for k,v in cond.items() if k!='incomplete_intended_exact'))
 extension=[r for r in rows if r['sample_group']=='extension_100' and int(r['reference_assessable_'+ref])]
 names=list(cond)+['coverage_only_independent','coverage_only_archived_priority_flags','any_condition_excluding_incomplete']
 out[ref]={'denominator':len(extension),**{n:sum(int(r[n+'_'+ref]) for r in extension) for n in names},'any_review_flag':sum(not r['category_'+ref].startswith('C') for r in extension)}
assert out['ensembl']['any_review_flag']==61 and out['refseq']['any_review_flag']==62
assert out['ensembl']['coverage_only_independent']==3 and out['refseq']['coverage_only_independent']==7
(P/'recomputed_rule_counts.json').write_text(json.dumps(out,indent=2),encoding='utf-8')
print(json.dumps(out,indent=2));print('All 114 rows and independent rule conditions verified.')
