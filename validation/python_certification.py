"""Save Python unit and exact untouched-baseline fixture parity evidence."""
import csv
import importlib.util
import json
import subprocess
import sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
(ROOT/'output').mkdir(exist_ok=True)
sys.path.insert(0,str(ROOT/'python'))
import summclust as new
spec=importlib.util.spec_from_file_location('old',ROOT/'tests/baseline/python/summclust.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
def equal(a,b):
    if isinstance(a,dict):
        for k in a:equal(a[k],b[k])
    elif isinstance(a,np.ndarray):np.testing.assert_array_equal(a,b)
    elif isinstance(a,float) and np.isnan(a):assert np.isnan(b)
    else:assert a==b,(a,b)
statuses=[]
for case in csv.DictReader((ROOT/'validation/fixtures/cases.csv').open()):
    d=new.read_csv(ROOT/'validation/fixtures'/case['dataset'],string_columns=['gs','hs'])
    opts={k:case[k]=='1' for k in ['jackknife','rankcheck','table','addmeans','regtable']}
    opts['rho']=float(case['rho'])
    if case['fevar']:opts['fevar']=case['fevar'].split('|')
    if case['absorb']:opts['absorb']=case['absorb']
    if case['sample']:opts['sample']=d[case['sample']]==1
    args=(d,case['outcome'],case['focal'],case['controls'].split('|'),case['cluster'].split('|'))
    results=[]
    for module,kw in [(old,opts),(new,opts),(new,dict(opts,mixed=False))]:
        try:results.append(module.summclust(*args,**kw))
        except module.SummclustError as e:results.append(e.code)
    equal(results[0],results[1]);equal(results[0],results[2])
    statuses.append(dict(specification=case['specification'],passed=True))
(ROOT/'output/python_default_parity.json').write_text(json.dumps(statuses,indent=2))
p=subprocess.run([sys.executable,'-m','unittest','discover','-s','python','-p','test*.py','-v'],cwd=ROOT,text=True,capture_output=True)
log=f'Baseline fixture parity: {len(statuses)} cases PASS\n'+p.stdout+p.stderr
(ROOT/'output/python_certification.log').write_text(log)
print(log);raise SystemExit(p.returncode)
