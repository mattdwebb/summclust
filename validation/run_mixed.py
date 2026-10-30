"""Export shared deterministic mixed fixtures, including cluster identities."""
import csv
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
(ROOT/'output').mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'python'))
from summclust import summclust,read_csv,export_rows

fields=['specification','statistic','row','column','value','row_label','column_label']
def label(v): return format(v,'.17g') if isinstance(v,(float,np.floating)) else str(v)
with (ROOT/'output/python_mixed.csv').open('w',newline='') as f, (ROOT/'output/python_mixed_keys.csv').open('w',newline='') as k:
    writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader()
    keys=csv.DictWriter(k,fieldnames=['specification','dimension','index','identifier1','identifier2']);keys.writeheader()
    for case in csv.DictReader((ROOT/'validation/fixtures/cases.csv').open()):
        if not case['specification'].startswith('two_'):continue
        d=read_csv(ROOT/'validation/fixtures'/case['dataset'],string_columns=['gs','hs'])
        opts={key:case[key]=='1' for key in ['rankcheck','table','addmeans','regtable']}
        opts.update(mixed=True,rho=float(case['rho']))
        if case['fevar']:opts['fevar']=case['fevar'].split('|')
        if case['absorb']:opts['absorb']=case['absorb']
        if case['sample']:opts['sample']=d[case['sample']]==1
        out=summclust(d,case['outcome'],case['focal'],case['controls'].split('|'),case['cluster'].split('|'),**opts)
        writer.writerows(export_rows(out,case['specification']))
        for dim,comp in out['components'].items():
            for idx,identity in enumerate(comp['cluster_ids'],1):
                a,b=identity if dim=='i' else (identity,'')
                keys.writerow(dict(specification=case['specification'],dimension=dim,index=idx,identifier1=label(a),identifier2=label(b)))
        assert out['intersection_delete_one_solves']==0
        print(case['specification'],'CV31',out['cv31_conventional_var'],out['cv31_max_se'])
