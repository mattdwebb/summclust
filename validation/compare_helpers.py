"""Numeric comparison helpers with coefficient and actual-cluster-ID alignment."""
import csv
import json
import math
import re

COVARIANCES = {"cv1", "cv3", "cv3j", "cv3drop", "cv3jdrop",
               "cv1_conventional", "cv3_conventional"}
for estimator in ("cv1", "cv3"):
    COVARIANCES.update(estimator + "_" + dimension for dimension in ("g", "h", "intersection"))
    COVARIANCES.update("V_" + estimator + "_" + dimension for dimension in ("g", "h", "i"))

def number(s):
    return math.nan if s.strip() in ("", ".", "NA", "nan", "NaN") else float(s)

def identifier(value):
    value = str(value).strip()
    try:
        x = float(value)
        return str(int(x)) if x.is_integer() else format(x, ".17g")
    except ValueError:
        return value

def read_results(path, keypath):
    if not path.is_file() or not keypath.is_file():
        raise ValueError(f"Missing results or cluster keys: {path.name}, {keypath.name}")
    clusters = {}
    with keypath.open(newline="", encoding="utf-8-sig") as stream:
        for row in csv.DictReader(stream):
            dimension = "i" if row["dimension"] == "intersection" else row["dimension"]
            key = (row["specification"], dimension, int(row["index"]))
            ids = (identifier(row["identifier1"]), identifier(row["identifier2"]))
            if key in clusters and clusters[key] != ids:
                raise ValueError(f"Conflicting cluster identity: {key}")
            clusters[key] = ids
    with path.open(newline="", encoding="utf-8-sig") as stream:
        entries = list(csv.DictReader(stream))
    coefficients = {}
    for row in entries:
        label = row.get("row_label", "")
        if row["statistic"] == "cv1" and label and not re.fullmatch(r"r\d+", label):
            coefficients[(row["specification"], row["row"])] = label
    values = {}
    for row in entries:
        spec, statistic = row["specification"], row["statistic"]
        if statistic == "clustsum":
            statistic = "summary"
        i, j = row["row"], row["column"]
        if statistic in COVARIANCES:
            target = (coefficients.get((spec, i), row.get("row_label") or i),
                      coefficients.get((spec, j), row.get("column_label") or j))
        elif statistic == "b":
            target = ("", row.get("column_label") or j)
        elif statistic.startswith(("ng", "lever", "partlev", "betanog", "rankfail", "legacyflag")) or statistic in ("scall", "table_g", "table_h"):
            dim = "h" if statistic.endswith("_h") else "i" if statistic.endswith(("_i", "_intersection")) else "g"
            identity_key = (spec, dim, int(i))
            if identity_key not in clusters:
                raise ValueError(f"No actual cluster identifier for {identity_key}")
            target = (json.dumps(clusters[identity_key], separators=(",", ":")), j)
        elif statistic in ("dimsummary", "summary", "bonus", "regsum", "regresstab"):
            target = (i, j)
        else:
            target = ("", "") if i in ("", "1") and j in ("", "1") else (i, j)
        key = (spec, statistic, *target)
        value = number(row["value"])
        if key in values and not ((math.isnan(value) and math.isnan(values[key])) or value == values[key]):
            raise ValueError(f"Conflicting duplicate numeric result: {key}")
        values[key] = value
    return values
