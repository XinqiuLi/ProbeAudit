"""ProbeAudit evidence explorer 1.1.0. Python standard library only.
Summarizes retained sequence hits; does not infer binding or annotation error.
"""
import argparse,csv,json,collections,html,math
from pathlib import Path

def read(path):
    with Path(path).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
def write(path,rows,fields=None):
    with Path(path).open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields or list(rows[0]));w.writeheader();w.writerows(rows)
def canon(s):return 'CYTOR' if s=='LINC00152' else s
def intended(symbol,label):return canon(symbol) in {canon(s) for s in label.split('|')}

def aggregate(inputs):
    inputs=Path(inputs);annotation=read(inputs/'annotation.csv'); byps={r['probe_id']:r for r in annotation}
    genomes={}; universe={};groups={}
    for group in ('historical','extension'):
        for r in read(inputs/f'{group}_genome.csv'):
            ps=r['probe_set'];universe.setdefault(ps,{})[r['oligo_id']]=r
            genomes.setdefault(ps,{})[r['oligo_id']]=int(r['exact_genomic_sites']);groups[ps]=group
    if set(byps)!=set(universe):raise ValueError('Annotation and oligo inventories differ')
    support=collections.defaultdict(set);meta={};tx=collections.defaultdict(set);loc=collections.defaultdict(set)
    for group in ('historical','extension'):
        for ref in ('Ensembl','RefSeq'):
            for r in read(inputs/f'{group}_{ref}_hits.csv'):
                ps,oid=r['probe_set'],r['oligo_id']
                if ps not in universe or oid not in universe[ps]:raise ValueError('Hit outside inventory')
                if r['sequence']!=universe[ps][oid]['sequence']:raise ValueError('Sequence mismatch')
                mismatch=int(r['mismatches'])
                if mismatch not in (0,1):raise ValueError('Expected zero or one substitution')
                for cutoff in range(mismatch,2):
                    k=(ps,ref,r['orientation'],cutoff,r['gene_id'])
                    if k in meta and meta[k]!=r['gene_symbol']:raise ValueError('Inconsistent gene symbol')
                    support[k].add(oid);meta[k]=r['gene_symbol'];tx[k].add(r['transcript']);loc[k].add(r['reference_location'])
    gene_rows=[]
    for k,oligos in sorted(support.items()):
        ps,ref,ori,cut,gid=k;label=byps[ps]['historical_candidate_label']
        gene_rows.append(dict(probe_id=ps,sample_group=groups[ps],candidate_label=label,reference=ref,orientation=ori,max_mismatches=cut,gene_id=gid,gene_symbol=meta[k],role='intended' if intended(meta[k],label) else 'alternative',n_oligos=len(oligos),n_total=len(universe[ps]),oligo_ids='|'.join(sorted(oligos)),transcript_ids='|'.join(sorted(tx[k])),reference_locations='|'.join(sorted(loc[k]-{''}))))
    metrics=[]
    for ps,inventory in sorted(universe.items()):
        ar=byps[ps];label=ar['historical_candidate_label'];n=len(inventory)
        if n!=11:raise ValueError('This study inventory requires 11 oligos per set')
        m=dict(probe_id=ps,sample_group=groups[ps],candidate_label=label,n_oligos=n,platform_annotation=ar['previous_annotation'])
        alts={}
        for ref in ('Ensembl','RefSeq'):
            def collect(ori,cut,role):
                return [(k,v) for k,v in support.items() if k[0:4]==(ps,ref,ori,cut) and intended(meta[k],label)==(role=='intended')]
            ins=collect('+',0,'intended');alt=collect('+',0,'alternative');rev=collect('-',0,'intended');near=collect('+',1,'alternative')
            union=lambda entries:set().union(*(v for k,v in entries))
            i,a,rv=union(ins),union(alt),union(rev)
            sym=collections.defaultdict(set)
            for k,v in alt:
                # Only explicit symbols are cross-reference comparable; unnamed IDs remain reference-specific.
                if not meta[k].startswith(('ENSG','GeneID_')):sym[canon(meta[k])]|=v
            alts[ref]=sym
            f=ref.lower();m.update({
                'intended_'+f:len(i),'alternative_union_'+f:len(a),
                'alternative_max_'+f:max((len(v) for k,v in alt),default=0),
                'alternative_genes_'+f:len(alt),'shared_max_'+f:max((len(i&v) for k,v in alt),default=0),
                'reverse_intended_'+f:len(rv),'reverse_only_intended_'+f:len(rv-i),
                'one_substitution_only_'+f:len(union(near)-a),
                'target_represented_'+f:int(ar.get('reference_assessable_'+f,'1')),
                'category_'+f:ar['category_'+f].split(':')[0],
                'top_alternatives_'+f:'; '.join(f'{meta[k]} [{k[4]}]:{len(v)}' for k,v in sorted(alt,key=lambda kv:(-len(kv[1]),kv[0][4])))})
            # Check against the previously verified aggregate record.
            if len(i)!=int(ar['intended_exact_'+f].split('/')[0]) or len(a)!=int(ar['other_gene_exact_'+f].split('/')[0]):raise ValueError('Retained count mismatch: '+ps+' '+ref)
        common=set(alts['Ensembl'])&set(alts['RefSeq'])
        m['consensus_alternative']=max((len(alts['Ensembl'][s]&alts['RefSeq'][s]) for s in common),default=0)
        m['consensus_symbols']='; '.join(f'{s}:{len(alts["Ensembl"][s]&alts["RefSeq"][s])}' for s in sorted(common))
        m['genome_multiple']=sum(x>1 for x in genomes[ps].values());m['genome_missing']=sum(x==0 for x in genomes[ps].values())
        metrics.append(m)
    return gene_rows,metrics,universe

MODES=('alternative','deficit','ambiguity','consensus','orientation','genome-multiple','genome-missing','combined')
def components(r,reference):
    refs=('ensembl','refseq') if reference=='both' else (reference,)
    return dict(alternative=max(r['alternative_max_'+f] for f in refs),
        deficit=max((r['n_oligos']-r['intended_'+f] for f in refs if r['target_represented_'+f]),default=0),
        ambiguity=max(r['shared_max_'+f] for f in refs),consensus=r['consensus_alternative'],
        orientation=max(r['reverse_only_intended_'+f] for f in refs),
        **{'genome-multiple':r['genome_multiple'],'genome-missing':r['genome_missing']})
def ranked(metrics,mode,reference='both',weights=None):
    weights=weights or {'alternative':1,'deficit':1,'ambiguity':1}
    rows=[]
    for m in metrics:
        c=components(m,reference)
        score=sum(c[k]*w for k,w in weights.items()) if mode=='combined' else c[mode]
        rows.append(dict(m,score=score,ranking_mode=mode,ranking_reference=reference))
    return sorted(rows,key=lambda r:(-r['score'],r['probe_id']))

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--inputs',type=Path,default=Path(__file__).parent/'inputs');p.add_argument('--out',type=Path,default=Path('results'))
    p.add_argument('--mode',choices=MODES,default='alternative');p.add_argument('--reference',choices=['both','ensembl','refseq'],default='both')
    p.add_argument('--group',choices=['all','historical','extension'],default='all');p.add_argument('--probe',default='')
    p.add_argument('--min-alternative',type=int,default=0);p.add_argument('--max-intended',type=int,default=11)
    p.add_argument('--min-consensus',type=int,default=0);p.add_argument('--weights',default='alternative=1,deficit=1,ambiguity=1')
    a=p.parse_args();a.out.mkdir(parents=True,exist_ok=True)
    weights={k:float(v) for k,v in (x.split('=') for x in a.weights.split(','))}
    if set(weights)-set(MODES[:-1]) or any(v<0 or not math.isfinite(v) for v in weights.values()):p.error('Weights require known finite nonnegative components')
    genes,metrics,inventory=aggregate(a.inputs)
    write(a.out/'gene_support.csv',genes);write(a.out/'probe_metrics.csv',metrics)
    hist=[]
    for group in ('historical','extension'):
        for ref in ('Ensembl','RefSeq'):
            g=[r for r in genes if r['sample_group']==group and r['reference']==ref and r['role']=='alternative' and r['orientation']=='+' and r['max_mismatches']==0]
            for n in range(11,0,-1):
                hist.append(dict(sample_group=group,reference=ref,n_supporting_oligos=n,probe_gene_pairs=sum(r['n_oligos']==n for r in g),distinct_probe_sets=len({r['probe_id'] for r in g if r['n_oligos']==n})))
    write(a.out/'support_distribution.csv',hist)
    rows=ranked(metrics,a.mode,a.reference,weights)
    refs=('ensembl','refseq') if a.reference=='both' else (a.reference,)
    rows=[r for r in rows if (a.group=='all' or r['sample_group']==a.group) and (not a.probe or a.probe in r['probe_id']) and max(r['alternative_max_'+f] for f in refs)>=a.min_alternative and min(r['intended_'+f] for f in refs)<=a.max_intended and r['consensus_alternative']>=a.min_consensus]
    write(a.out/'ranking.csv',rows,fields=list(ranked(metrics,a.mode,a.reference,weights)[0]))
    examples=[]
    for mode in MODES:
        for rank,r in enumerate(ranked(metrics,mode),1):examples.append(dict(mode=mode,rank=rank,probe_id=r['probe_id'],candidate_label=r['candidate_label'],score=r['score']))
    write(a.out/'all_rankings.csv',examples)
    # Self-contained evidence browser. All values are escaped before HTML insertion.
    body=['<!doctype html><meta charset="utf-8"><title>ProbeAudit evidence explorer</title><style>body{font:15px Arial;margin:30px;color:#182b3a}table{border-collapse:collapse;width:100%;font-size:12px}td,th{border:1px solid #ccd5df;padding:7px;text-align:left}th{background:#e8f0f5}summary{cursor:pointer;font-weight:bold}details{margin:12px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere}</style><h1>ProbeAudit evidence explorer</h1><p>Scores are transparent sorting counts, not probabilities of error or measures of biological impact. Gene support counts distinct oligo IDs, deduplicating transcripts and positions. A single oligo may support several genes. Reverse complements and one-substitution matches remain separate.</p>']
    body.append('<p>Selected ranking: '+html.escape(a.mode)+'; reference: '+a.reference+'. Change ranking and filters by running prioritize.py --help.</p>')
    for r in rows:
        ps=r['probe_id'];body.append('<details><summary>'+html.escape(f'{ps} | {r["candidate_label"]} | score {r["score"]}')+'</summary><pre>'+html.escape(json.dumps(r,indent=2))+'</pre>')
        body.append('<table><tr><th>Reference / orientation / substitutions</th><th>Gene / role</th><th>Oligos</th><th>Supporting evidence</th></tr>')
        for g in genes:
            if g['probe_id']!=ps:continue
            body.append('<tr><td>'+html.escape(f'{g["reference"]} / {g["orientation"]} / ≤{g["max_mismatches"]}')+'</td><td>'+html.escape(f'{g["gene_symbol"]} [{g["gene_id"]}] / {g["role"]}')+'</td><td>'+str(g['n_oligos'])+'/11</td><td><details><summary>Oligos, transcripts and locations</summary><pre>'+html.escape(g['oligo_ids']+'\n'+g['transcript_ids']+'\n'+g['reference_locations'])+'</pre></details></td></tr>')
        body.append('</table><details><summary>Oligo sequences and genomic site counts</summary><pre>'+html.escape(json.dumps(list(inventory[ps].values()),indent=2))+'</pre></details></details>')
    (a.out/'report.html').write_text('\n'.join(body),encoding='utf-8')
    (a.out/'run.json').write_text(json.dumps(dict(version='1.1.0',mode=a.mode,reference=a.reference,weights=weights,group=a.group,probe=a.probe,min_alternative=a.min_alternative,max_intended=a.max_intended,min_consensus=a.min_consensus,total_sets=len(metrics),shown_sets=len(rows)),indent=2),encoding='utf-8')
    print(f'Generated {len(genes)} gene-support rows and {len(metrics)} probe metrics; displayed {len(rows)} sets.')
if __name__=='__main__':main()
