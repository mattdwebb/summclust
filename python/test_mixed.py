"""Independent deterministic checks for the explicit CV31 option."""
import csv
import importlib.util
import math
from pathlib import Path
import unittest
from unittest.mock import patch
import numpy as np
import summclust as sc
from test_summclust import fixture

ROOT = Path(__file__).resolve().parents[1]


class MixedTests(unittest.TestCase):
    def fit(self, d=None, cluster=('g', 'h'), **kw):
        return sc.summclust(fixture() if d is None else d, 'y', 'x', ['z', 'w'], cluster, **kw)

    def test_defaults_and_components(self):
        oldspec = importlib.util.spec_from_file_location('baseline', ROOT/'tests/baseline/python/summclust.py')
        old = importlib.util.module_from_spec(oldspec)
        oldspec.loader.exec_module(old)
        for cluster in ['g', ['g','h']]:
            a = old.summclust(fixture(), 'y','x',['z','w'],cluster)
            b = self.fit(cluster=cluster)
            c = self.fit(cluster=cluster,mixed=False)
            for key in a:
                if isinstance(a[key], np.ndarray):
                    np.testing.assert_array_equal(a[key], b[key])
                    np.testing.assert_array_equal(a[key], c[key])
        a = self.fit()
        b = self.fit(mixed=True, table=True, rankcheck=True)
        for key in ['b','cv1','cv1_intersection','cv3_g','cv3_h']:
            np.testing.assert_array_equal(a[key], b[key])
        expected = a['cv3_g'] + a['cv3_h'] - a['cv1_intersection']
        np.testing.assert_array_equal(expected, b['cv31_conventional'])
        self.assertTrue(np.isnan(b['cv3']).all())
        self.assertTrue(np.isnan(b['betanog_i']).all())
        self.assertTrue(np.isfinite(b['rankfail_i']).all())
        self.assertEqual(b['intersection_delete_one_solves'], 0)

    def test_direct_scores_and_deletion_fits(self):
        d=fixture(); x=np.column_stack([d['x'],np.ones(len(d['x'])),d['z'],d['w']]); y=d['y']
        b=np.linalg.lstsq(x,y,rcond=None)[0]; u=y-x@b; bread=np.linalg.inv(x.T@x)
        vg=[]
        for dim in ['g','h']:
            deleted=np.array([np.linalg.lstsq(x[d[dim]!=g],y[d[dim]!=g],rcond=None)[0] for g in np.unique(d[dim])])
            delta=deleted-b
            vg.append((len(deleted)-1)/len(deleted)*delta.T@delta)
        ids=list(zip(d['g'],d['h'])); labels=sorted(set(ids)); meat=np.zeros((x.shape[1],)*2)
        for label in labels:
            use=np.array([v==label for v in ids]); score=x[use].T@u[use]; meat+=np.outer(score,score)
        n,k=x.shape; I=len(labels); v1=I*(n-1)/((I-1)*(n-k))*bread@meat@bread
        out=self.fit(mixed=True)
        np.testing.assert_allclose(out['cv1_intersection'],v1,atol=1e-13,rtol=1e-10)
        np.testing.assert_allclose(out['cv31'],vg[0]+vg[1]-v1,atol=1e-13,rtol=1e-9)

    def test_no_intersection_solver_calls(self):
        original=sc._component; counts={}
        def observed(*args, **kwargs):
            with patch.object(sc, '_solve', wraps=sc._solve) as solve:
                result=original(*args, **kwargs)
                counts[args[-1]]=solve.call_count
                return result
        with patch.object(sc,'_component',side_effect=observed):
            out=self.fit(mixed=True,rankcheck=True,table=True)
        self.assertEqual(counts['I'],0)
        self.assertEqual(counts['G'],out['G']); self.assertEqual(counts['H'],out['H'])

    def test_inference_and_identities(self):
        out=self.fit(mixed=True)
        var=max(out['cv31'][0,0],out['cv3_g'][0,0],out['cv3_h'][0,0])
        self.assertEqual(out['cv31_max_var'],var)
        inf=sc._inference(out['beta'],math.sqrt(var),min(out['G'],out['H'])-1)
        for key in ['se','t','p','lci','uci']:
            self.assertEqual(out['cv31_max_'+key],inf[key])
        self.assertEqual(out['inference_label'],'CV31 (mixed, max)')
        same=self.fit(cluster=['g','g'],mixed=True)
        np.testing.assert_allclose(same['cv31'],2*same['cv3_g']-same['cv1_intersection'])
        d=fixture(); d['coarse']=d['g']//3
        nested=self.fit(d,cluster=['g','coarse'],mixed=True)
        np.testing.assert_allclose(nested['cv31'],nested['cv3_g']+nested['cv3_h']-nested['cv1_g'])
        swap=self.fit(cluster=['h','g'],mixed=True)
        np.testing.assert_allclose(out['cv31'],swap['cv31'],atol=1e-13)

    def test_variations(self):
        d=fixture(); d['a']=d['g']*100+d['h']; d['zero']=np.zeros(len(d['y']))
        for opts in [dict(fevar=['g','h']),dict(absorb='a'),dict(sample=np.arange(len(d['y']))%7!=0),dict(cluster=['gs','h'])]:
            pure=self.fit(d,**opts); mixed=self.fit(d,mixed=True,**opts)
            np.testing.assert_allclose(mixed['cv31'],pure['cv3_g']+pure['cv3_h']-pure['cv1_intersection'])
        perm=np.arange(len(d['y']))[::-1]
        out=self.fit(d,mixed=True); rev=self.fit({k:v[perm] for k,v in d.items()},mixed=True)
        np.testing.assert_allclose(out['cv31'],rev['cv31'],atol=1e-12)
        for dim in 'ghi':
            a=dict(zip(out['cluster_ids'][dim],out['components'][dim]['partlev']))
            b=dict(zip(rev['cluster_ids'][dim],rev['components'][dim]['partlev']))
            for key in a: self.assertAlmostEqual(a[key],b[key],places=12)

    def test_invalid(self):
        for value in [None,0,1,'true',[],np.array([True])]:
            with self.assertRaises(sc.SummclustError): self.fit(mixed=value)
        with self.assertRaises(sc.SummclustError): self.fit(cluster='g',mixed=True)
        with self.assertRaises(sc.SummclustError): self.fit(mixed=True,jackknife=True)

    def test_negative_conventional_and_focal(self):
        d=sc.read_csv(ROOT/'validation/fixtures/negative.csv')
        out=sc.summclust(d,'y','x',['z'],['g','h'],mixed=True)
        self.assertLess(out['cv31_conventional_var'],0)
        self.assertTrue(math.isnan(out['cv31_conventional_se']))
        self.assertTrue(math.isnan(out['inference_conventional']['cv31']['p']))
        self.assertTrue(math.isfinite(out['cv31_max_p']))
        d=fixture()
        pure=sc.summclust(d,'y','z',['x','w'],['g','h'])
        mixed=sc.summclust(d,'y','z',['x','w'],['g','h'],mixed=True)
        np.testing.assert_allclose(mixed['cv31'],pure['cv3_g']+pure['cv3_h']-pure['cv1_intersection'])


if __name__=='__main__': unittest.main()
