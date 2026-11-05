"""Worked example for the Python guide, using fictional school-program data."""
from pathlib import Path
import sys

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from summclust import read_csv, summclust

data = read_csv(ROOT / "examples/school_program.csv", string_columns=["school"])
fit = summclust(
    data, outcome="test_score", focal="program", controls=["baseline_score"],
    cluster="school", fevar="year",
)

print(f"{fit['N']} observations in {fit['G']} schools")
for method in ("cv1", "cv3"):
    r = fit["inference"][method]
    print(f"{method.upper()}: estimate={r['beta']:.3f}, SE={r['se']:.3f}, "
          f"p={r['p']:.3f}, 95% CI=[{r['lci']:.3f}, {r['uci']:.3f}]")

print()
print("School       N   Partial leverage   Estimate without school")
for j in np.argsort(-fit["partlev"])[:5]:
    print(f"{fit['cluster_ids'][j]:<10} {int(fit['ng'][j]):3d} "
          f"{fit['partlev'][j]:18.3f} {fit['betanog'][j]:25.3f}")
