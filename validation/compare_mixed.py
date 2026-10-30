"""Compare mixed results with coefficient and actual cluster ID alignment."""
import csv
import json
import math
from pathlib import Path
from compare_helpers import read_results,COVARIANCES
ROOT=Path(__file__).resolve().parents[1]
COVARIANCES.update(['cv31','cv31_conventional'])
loaded={lang:read_results(ROOT/f'output/{lang}_mixed.csv',ROOT/f'output/{lang}_mixed_keys.csv') for lang in ['stata','r','python']}
required={'b','cv1','cv31_conventional','cv3_g','cv3_h','cv1_intersection','cv31_conventional_var','cv31_conventional_se','cv31_max_var','cv31_max_se','cv31_max_t','cv31_max_p','cv31_max_lci','cv31_max_uci','df_r','intersection_delete_one_solves'}
required.update(f'{s}_{d}' for s in ['ng','lever','partlev','betanog'] for d in ['g','h','intersection'])
rows=[];errors=[]
specs=sorted({key[0] for key in loaded['stata']})
for spec in specs:
    for stat in sorted(required):
        keys={k for k in loaded['stata'] if k[:2]==(spec,stat)}
        if not keys:errors.append(f'Stata missing {spec}:{stat}')
        for lang in ['r','python']:
            worst=0.;count=0
            for key in sorted(keys):
                if key not in loaded[lang]:errors.append(f'{lang} missing {key}');continue
                a,b=loaded['stata'][key],loaded[lang][key]
                if math.isnan(a) and math.isnan(b):continue
                count+=1
                if not(math.isfinite(a) and math.isfinite(b)):
                    errors.append(f'missing mismatch {lang}:{key}');continue
                delta=abs(a-b);worst=max(worst,delta)
                if delta>2e-10+2e-8*max(abs(a),abs(b)):
                    errors.append(f'{lang}:{key} delta {delta:.17g}')
            rows.append(dict(specification=spec,statistic=stat,comparison='stata-'+lang,checked_cells=count,max_abs_difference=worst))
with (ROOT/'output/mixed_cross_language_comparison.csv').open('w',newline='') as f:
    w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
summary=dict(cases=len(specs),checked_cells=sum(r['checked_cells'] for r in rows),max_abs_difference=max(r['max_abs_difference'] for r in rows),errors=errors)
(ROOT/'output/mixed_cross_language_summary.json').write_text(json.dumps(summary,indent=2))
print(json.dumps(summary,indent=2));raise SystemExit(bool(errors))
