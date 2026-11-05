"""Shared string fixed-effect fits and comparison with Stata and R."""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "python"))
from summclust import summclust, read_csv, export_rows
from compare_helpers import read_results, COVARIANCES

CASES = [f"{mode}_{design}" for mode in ("one", "two", "mixed")
         for design in ("fe", "absorb", "both")]


def export():
    (ROOT / "output").mkdir(exist_ok=True)
    fields = ["specification", "statistic", "row", "column", "value", "row_label", "column_label"]
    with (ROOT / "output/python_string_fe.csv").open("w", newline="") as out, \
            (ROOT / "output/python_string_fe_keys.csv").open("w", newline="") as ids:
        writer = csv.DictWriter(out, fieldnames=fields)
        writer.writeheader()
        keys = csv.DictWriter(ids, fieldnames=["specification", "dimension", "index", "identifier1", "identifier2"])
        keys.writeheader()
        for case in CASES:
            mode, design = case.split("_")
            data = read_csv(ROOT / "validation/fixtures/synthetic.csv", string_columns=["gs", "hs"])
            for name, source, prefix, width, modulus in [
                    ("fstr", "fe", "Effect / ", 2, 19),
                    ("hstr", "h", "Period ", 2, 23),
                    ("astr", "a", "Cell: ", 4, 29)]:
                data[name] = np.array([f"{prefix}{int(v):0{width}d}" for v in data[source]], dtype=object)
                data[name][data["row"] % modulus == 0] = ""
            data["fstr"][data["keep"] == 0] = "Excluded level"
            fit = summclust(data, "y", "x", ["z", "w"], ["gs"] if mode == "one" else ["gs", "hs"],
                            fevar=["fstr", "hstr"] if design == "fe" else ["fstr"] if design == "both" else None,
                            absorb="astr" if design != "fe" else None, sample=data["keep"] == 1,
                            rho=.5, rankcheck=True, table=True, jackknife=mode == "one", mixed=mode == "mixed")
            writer.writerows(export_rows(fit, case))
            for dim, comp in fit["components"].items():
                for index, identity in enumerate(comp["cluster_ids"], 1):
                    first, second = identity if dim == "i" else (identity, "")
                    keys.writerow(dict(specification=case, dimension=dim, index=index,
                                       identifier1=first, identifier2=second))
            print("PASS Python string FE export", case)


def compare():
    COVARIANCES.update(["cv31", "cv31_conventional"])
    results = {lang: read_results(ROOT / f"output/{lang}_string_fe.csv", ROOT / f"output/{lang}_string_fe_keys.csv")
               for lang in ("stata", "r", "python")}
    errors, cells, maximum = [], 0, 0.0
    for case in CASES:
        required = {"b", "cv1"}
        if case.startswith("one_"):
            required.update(["cv3", "cv3j", "ng", "lever", "partlev", "betanog", "rankfail"])
        else:
            required.add("cv31_conventional" if case.startswith("mixed_") else "cv3_conventional")
            required.update(f"{stat}_{dim}" for stat in ("ng", "lever", "partlev", "betanog", "rankfail")
                            for dim in ("g", "h", "intersection"))
            required.update(["N", "df_r", "cv3_g", "cv3_h", "cv1_intersection"])
        for stat in required:
            expected = {key for key in results["stata"] if key[:2] == (case, stat)}
            if not expected:
                errors.append(f"Stata missing {case}:{stat}")
            for lang in ("r", "python"):
                actual = {key for key in results[lang] if key[:2] == (case, stat)}
                if actual != expected:
                    errors.append(f"{lang} result keys differ for {case}:{stat}")
                for key in expected & actual:
                    a, b = results["stata"][key], results[lang][key]
                    if math.isnan(a) and math.isnan(b):
                        continue
                    cells += 1
                    if not (math.isfinite(a) and math.isfinite(b)):
                        errors.append(f"Missing mismatch {lang}:{key}")
                        continue
                    delta = abs(a - b)
                    maximum = max(maximum, delta)
                    if delta > 2e-10 + 2e-8 * max(abs(a), abs(b)):
                        errors.append(f"{lang}:{key} difference={delta}")
    summary = dict(cases=len(CASES), checked_cells=cells, max_abs_difference=maximum, errors=errors)
    (ROOT / "output/string_fe_cross_language_summary.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return bool(errors)


if __name__ == "__main__":
    if "--compare" in sys.argv:
        raise SystemExit(compare())
    export()
