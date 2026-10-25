import unittest,random
import numpy as np,pandas as pd
from probe_audit_tool import make_matcher,matches,rc,fold_design
class AuditTests(unittest.TestCase):
 def test_matcher_matches_exhaustive_both_strands(self):
  rng=random.Random(9);qs=[''.join(rng.choices('ACGT',k=25)) for _ in range(10)]
  seq='NN'+''.join(q+'N'+q[:8]+('A' if q[8]!='A' else 'C')+q[9:]+'N'+rc(q)+'N' for q in qs)
  for cutoff in [0,1]:
   actual=set(matches(make_matcher(qs,cutoff),seq,cutoff));expected=set()
   for q in qs:
    for strand,t in [('+',q),('-',rc(q))]:
     for i in range(len(seq)-24):
      s=seq[i:i+25];dist=sum(x!=y for x,y in zip(t,s))
      if 'N' not in s and dist<=cutoff:expected.add((q,strand,i,dist))
   self.assertEqual(actual,expected)
 def test_training_scores_ignore_heldout_extremes(self):
  e=pd.DataFrame([[1,2,3,4],[3,5,2,8]],index=['p1','p2']);m=pd.DataFrame({'program':['P','P'],'probe':['p1','p2']});c=pd.DataFrame({'probe':['p1']});tr=np.array([0,1,2])
  a,_=fold_design(e,m,c,tr);e.iloc[:,3]=1e9;b,_=fold_design(e,m,c,tr)
  np.testing.assert_array_equal(a[tr],b[tr])
 def test_reverse_complement_and_palindrome(self):
  self.assertEqual(rc('ACGTN'),'NACGT');q='ACGT'*6+'A';self.assertEqual(rc(rc(q)),q)
if __name__=='__main__':unittest.main()
