import unittest,tempfile,csv
from pathlib import Path
from prioritize import aggregate,ranked,components,write

class EvidenceTests(unittest.TestCase):
    def test_synthetic_deduplication_and_separation(self):
        with tempfile.TemporaryDirectory() as td:
            p=Path(td);ids=[f'p_{i}' for i in range(11)]
            genomic=[dict(probe_set='p',oligo_id=x,sequence='A'*25,expected_symbol='TARGET',exact_genomic_sites=1) for x in ids]
            write(p/'historical_genome.csv',genomic);write(p/'extension_genome.csv',[],list(genomic[0]))
            a=dict(probe_id='p',historical_candidate_label='TARGET',previous_annotation='TARGET',category_ensembl='A',category_refseq='A',reference_assessable_ensembl=1,reference_assessable_refseq=1,intended_exact_ensembl='2/11',intended_exact_refseq='1/11',other_gene_exact_ensembl='3/11',other_gene_exact_refseq='1/11')
            write(p/'annotation.csv',[a])
            def hit(i,g,trans='t',ori='+',mm=0):return dict(probe_set='p',oligo_id=ids[i],sequence='A'*25,expected_symbol='TARGET',transcript=trans,gene_id=g,gene_symbol=g,reference_location='',orientation=ori,mismatches=mm)
            ens=[hit(0,'TARGET'),hit(1,'TARGET'),hit(0,'ALT','t1'),hit(0,'ALT','t2'),hit(1,'ALT'),hit(2,'OTHER'),hit(3,'ALT',mm=1),hit(4,'TARGET',ori='-')]
            ref=[hit(0,'TARGET'),hit(1,'ALT')]
            for name,rr in [('historical_Ensembl',ens),('historical_RefSeq',ref),('extension_Ensembl',[]),('extension_RefSeq',[])]:write(p/(name+'_hits.csv'),rr,list(ens[0]))
            g,m,_=aggregate(p);r=m[0]
            self.assertEqual(r['alternative_union_ensembl'],3)
            self.assertEqual(r['alternative_max_ensembl'],2)
            self.assertEqual(r['shared_max_ensembl'],2)
            self.assertEqual(r['consensus_alternative'],1)
            self.assertEqual(r['reverse_only_intended_ensembl'],1)
            self.assertEqual(r['one_substitution_only_ensembl'],1)
            exact=next(x for x in g if x['reference']=='Ensembl' and x['gene_id']=='ALT' and x['max_mismatches']==0 and x['orientation']=='+')
            self.assertEqual(exact['n_oligos'],2)
    def test_archived_counts_and_priority_examples(self):
        g,m,_=aggregate(Path(__file__).parent/'inputs');d={r['probe_id']:r for r in m}
        self.assertEqual(len(m),114)
        for ref,n in [('Ensembl',88),('RefSeq',42)]:
            self.assertEqual(sum(r['sample_group']=='extension' and r['reference']==ref and r['orientation']=='+' and r['max_mismatches']==0 and r['role']=='alternative' for r in g),n)
        self.assertEqual(d['207719_x_at']['consensus_alternative'],11)
        self.assertEqual(d['225799_at']['shared_max_refseq'],11)
        self.assertEqual(d['216249_at']['reverse_only_intended_refseq'],7)
        self.assertEqual(d['216111_x_at']['alternative_union_ensembl'],10)
        self.assertEqual(d['216111_x_at']['alternative_max_ensembl'],5)
        self.assertEqual([r['probe_id'] for r in ranked(m,'consensus')[:3]],['207719_x_at','212746_s_at','225799_at'])
        self.assertEqual(ranked(m,'orientation')[1]['probe_id'],'216249_at')
        self.assertEqual(components(d['1561050_a_at'],'refseq')['deficit'],0)
    def test_ties_and_user_weights(self):
        _,m,_=aggregate(Path(__file__).parent/'inputs')
        a=ranked(m,'alternative');b=ranked(m,'combined',weights={'alternative':2})
        self.assertEqual([r['probe_id'] for r in a],[r['probe_id'] for r in b])
        self.assertEqual(b[0]['score'],2*a[0]['score'])
if __name__=='__main__':unittest.main()
