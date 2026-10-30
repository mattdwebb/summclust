"""Independent NumPy implementation of the local Stata summclust 5.1.1 API.

Copyright (c) 2026 Matthew Webb.  See the bundle LICENSE for attribution.
Only Python's standard library and NumPy are required. No bootstrap is run.
"""
from __future__ import annotations

import csv
import math
from collections.abc import Mapping, Sequence
from pathlib import Path
import numpy as np

__version__ = "0.1.0"
_EPS = np.finfo(float).eps


class SummclustError(ValueError):
    """Invalid specification or unidentified native two-way focal coefficient."""
    def __init__(self, message, code=198):
        super().__init__(message)
        self.code = code


def _betacf(a, b, x):
    """Modified Lentz evaluation of the incomplete-beta continued fraction."""
    tiny = 1e-300
    qab, qap, qam = a + b, a + 1.0, a - 1.0
    c, d = 1.0, 1.0 - qab * x / qap
    if abs(d) < tiny:
        d = tiny
    d = 1.0 / d
    h = d
    for m in range(1, 10001):
        m2 = 2 * m
        aa = m * (b - m) * x / ((qam + m2) * (a + m2))
        d, c = 1.0 + aa * d, 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        h *= d * c
        aa = -(a + m) * (qab + m) * x / ((a + m2) * (qap + m2))
        d, c = 1.0 + aa * d, 1.0 + aa / c
        if abs(d) < tiny:
            d = tiny
        if abs(c) < tiny:
            c = tiny
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) <= 4 * _EPS:
            return h
    raise ArithmeticError("Incomplete beta continued fraction did not converge")


def _ibeta(a, b, x):
    if x <= 0:
        return 0.0
    if x >= 1:
        return 1.0
    logbt = math.lgamma(a + b) - math.lgamma(a) - math.lgamma(b)
    logbt += a * math.log(x) + b * math.log1p(-x)
    bt = math.exp(logbt)
    if x < (a + 1) / (a + b + 2):
        return min(1.0, max(0.0, bt * _betacf(a, b, x) / a))
    return min(1.0, max(0.0, 1.0 - bt * _betacf(b, a, 1 - x) / b))


def student_t_sf(t: float, df: float) -> float:
    """Student-t survival probability P(T > t); no normal approximation.

    Uses the regularized incomplete beta and handles both tails, including
    infinite t. A nonpositive or nonfinite degrees of freedom is an error.
    """
    if not math.isfinite(df) or df <= 0:
        raise ValueError("df must be finite and positive")
    if math.isnan(t):
        return math.nan
    if math.isinf(t):
        return 0.0 if t > 0 else 1.0
    if t == 0:
        return 0.5
    if abs(t) <= 0.5:
        z = t * t / (df + t * t)
        central = 0.5 * _ibeta(0.5, df / 2.0, z)
        return 0.5 - central if t > 0 else 0.5 + central
    z = df / (df + t * t)
    half_tail = 0.5 * _ibeta(df / 2.0, 0.5, z)
    return half_tail if t > 0 else 1.0 - half_tail


def student_t_isf(p: float, df: float) -> float:
    """Inverse survival probability, including p=0 or 1 endpoints."""
    if not math.isfinite(df) or df <= 0:
        raise ValueError("df must be finite and positive")
    if not 0 <= p <= 1:
        raise ValueError("p must be between zero and one")
    if p == 0:
        return math.inf
    if p == 1:
        return -math.inf
    if p == 0.5:
        return 0.0
    if p > 0.5:
        return -student_t_isf(1.0 - p, df)
    lo, hi = 0.0, 1.0
    while student_t_sf(hi, df) > p:
        hi *= 2.0
        if math.isinf(hi):
            raise ArithmeticError("Student-t quantile exceeds floating-point range")
    for _ in range(160):
        mid = lo + (hi - lo) / 2
        if student_t_sf(mid, df) > p:
            lo = mid
        else:
            hi = mid
        if hi - lo <= 8 * _EPS * max(abs(mid), 1e-300):
            return lo + (hi - lo) / 2
    return lo + (hi - lo) / 2


def _rank(a):
    s = np.linalg.svd(a, compute_uv=False)
    return int(np.sum(s > (len(s) * s[0] * _EPS))) if len(s) else 0


def _column_lattice(values):
    """Binary quantum and magnitude bound, for identifying exact dot products."""
    maximum = float(np.max(np.abs(values), initial=0))
    if maximum == 0:
        return 0, 0
    bound = math.ceil(math.log2(maximum))
    if np.all(values == np.rint(values)):
        return 0, bound
    bits = np.ascontiguousarray(values).view(np.uint64) & np.uint64(0x7fffffffffffffff)
    bits = bits[bits != 0]
    exponent = (bits >> np.uint64(52)).astype(np.int64)
    mantissa = bits & np.uint64((1 << 52) - 1)
    mantissa |= np.where(exponent != 0, np.uint64(1 << 52), np.uint64(0))
    lowbit = mantissa & (~mantissa + np.uint64(1))
    zeros = np.log2(lowbit).astype(np.int64)
    quantum = int(np.min(np.where(exponent == 0, -1074, exponent - 1075) + zeros))
    return quantum, bound


def _accurate_dot(a, b):
    # Mata quadcross accumulates binary64-rounded products at higher precision.
    # Retaining the exact multiplication error would change its convention.
    return math.fsum(a * b)


def _cross(a, b):
    """Accurate crossproducts, retaining BLAS where every addition is exact."""
    same = a is b
    vector = b.ndim == 1
    b = b[:, None] if vector else b
    result = a.T @ b
    astats = [_column_lattice(a[:, j]) for j in range(a.shape[1])]
    bstats = astats if same else [_column_lattice(b[:, j]) for j in range(b.shape[1])]
    nbits = math.ceil(math.log2(max(1, len(a))))
    for j, (qa, ba) in enumerate(astats):
        for k, (qb, bb) in enumerate(bstats):
            if same and k < j:
                continue
            product_bits = ba + bb - qa - qb
            if nbits + product_bits > 53:
                result[j, k] = _accurate_dot(a[:, j], b[:, k])
                if same:
                    result[k, j] = result[j, k]
    return result[:, 0] if vector else result


def _invsym(a):
    """Symmetric, zero-omission generalized inverse rather than pinv.

    A pivoted symmetric sweep chooses an independent principal submatrix.
    Omitted coordinates are exactly zero, matching the invsym representation.
    See README for the numerical convention and its validation limits.
    """
    a = np.asarray(a, dtype=float)
    n = len(a)
    if n == 0:
        return a.copy()
    work = (np.triu(a) + np.triu(a, 1).T).copy()
    original_diag = np.diag(work).copy()
    remaining = list(range(n))
    selected = []
    while remaining:
        j = max(remaining, key=lambda i: work[i, i])
        pivot = work[j, j]
        remaining.remove(j)
        if pivot <= 1e-9 * abs(original_diag[j]):
            continue
        selected.append(j)
        if remaining:
            w = work[remaining, j].copy()
            # Divide the pivot column first, then form its outer product.
            # Algebraically equivalent orders can pick a different omitted
            # nuisance column when residual pivots tie to rounding precision.
            work[np.ix_(remaining, remaining)] -= np.outer(w / pivot, w)
            # invsym treats a symmetric matrix as its upper triangle. The two
            # computed triangles can differ by one ulp after division-first
            # updates; use the same triangle for subsequent pivot columns.
            work = np.triu(work) + np.triu(work, 1).T
    result = np.zeros_like(a)
    if selected:
        block = a[np.ix_(selected, selected)]
        result[np.ix_(selected, selected)] = np.linalg.inv(block)
    return (result + result.T) * 0.5


def _solve(a, b, force=False, check_condition=True):
    if force or (check_condition and np.linalg.cond(a) > 1e8):
        return _invsym(a) @ b
    try:
        # The authoritative Mata solver uses a right-looking lower Cholesky,
        # with reciprocal multiplication before each symmetric rank update.
        # LAPACK can make a different decision in singular nuisance directions.
        work = (np.tril(a) + np.tril(a, -1).T).copy()
        lower = np.zeros_like(work)
        for j in range(len(a)):
            if work[j, j] <= 0:
                raise np.linalg.LinAlgError("Nonpositive Cholesky pivot")
            lower[j, j] = math.sqrt(work[j, j])
            v = work[j + 1:, j] * (1.0 / lower[j, j])
            lower[j + 1:, j] = v
            work[j + 1:, j + 1:] -= np.outer(v, v)
        tolerance = 1e-13 * float(np.mean(np.abs(np.diag(lower))))
        if np.any(np.diag(lower) <= tolerance):
            raise np.linalg.LinAlgError("Cholesky solve tolerance")
        result = np.asarray(b, dtype=float).copy()
        for j in range(len(a)):
            result[j] = (result[j] - _mata_dot(lower[j, :j], result[:j])) / lower[j, j]
        for j in range(len(a) - 1, -1, -1):
            result[j] = (result[j] - _mata_dot(lower[j + 1:, j], result[j + 1:])) / lower[j, j]
        return result
    except np.linalg.LinAlgError:
        return _invsym(a) @ b


def _mata_dot(a, b):
    """Mata matrix-product accumulation for its triangular row solves."""
    odd = 0.0
    even = 0.0
    for j in range(0, len(a), 2):
        odd = odd + a[j] * b[j]
    for j in range(1, len(a), 2):
        even = even + a[j] * b[j]
    return odd + even


def _identified(a):
    target = np.zeros(len(a))
    target[0] = 1
    return np.max(np.abs(_invsym(a)[0] @ a - target)) < 1e-7


def _names(value):
    if value is None:
        return []
    return [value] if isinstance(value, str) else list(value)


def _missing(value):
    if value is None or (isinstance(value, str) and value == ""):
        return True
    try:
        return not math.isfinite(float(value))
    except (ValueError, TypeError):
        return False


def _encode(values):
    # Cluster types must be homogeneous, as each Stata variable has one type.
    vals = list(values)
    try:
        levels = sorted(set(vals))
    except TypeError as exc:
        raise SummclustError("Each identifier column must have a consistent type") from exc
    index = {label: i for i, label in enumerate(levels)}
    return np.asarray([index[v] for v in vals], dtype=int), levels


def _nested(fine, coarse):
    return all(len(set(coarse[fine == j])) == 1 for j in np.unique(fine))


def _gstar(gamma, twoway=False):
    m = float(np.mean(gamma))
    if not math.isfinite(m) or m == 0 or (twoway and abs(m) < 1e-14):
        return math.nan
    return len(gamma) / (1 + float(np.mean(((gamma - m) / m) ** 2)))


def _cv(values, denominator_n=None, twoway=False):
    z = np.asarray(values, dtype=float)
    m = float(np.mean(z))
    n = len(z) if denominator_n is None else denominator_n
    if n < 2 or m == 0 or not math.isfinite(m) or (twoway and abs(m) < 1e-14):
        return math.nan
    return math.sqrt(float(np.sum((z - m) ** 2)) / ((n - 1) * m * m))


def _inference(beta, se, df):
    if df <= 0 or not np.isfinite(se) or se <= 0:
        t, p, low, high = (math.nan,) * 4
    else:
        t = beta / se
        p = 2 * student_t_sf(abs(t), df)
        critical = student_t_isf(0.025, df)
        low, high = beta - critical * se, beta + critical * se
    return {"beta": float(beta), "se": float(se), "t": float(t), "p": float(p),
            "lci": float(low), "uci": float(high), "df": float(df)}


def _safe_sqrt(value):
    return math.sqrt(value) if value >= 0 else math.nan


def _stata_default_number(value):
    """Numeric value of Mata strofreal(value), whose default is %9.0g.

    Its precision is a format convention, not floating-point accuracy. Keep
    this restricted to one-way legacy aliases, never covariance calculations.
    """
    value = float(value)
    if not math.isfinite(value) or value == 0:
        return value
    absolute = abs(value)
    if absolute < 1e-5 or absolute >= 1e7:
        return float(format(value, ".2e"))
    decimals = 7 if absolute < 1 else max(0, 6 - math.floor(math.log10(absolute)))
    rounded = float(format(value, f".{decimals}f"))
    if abs(rounded) >= 1e7:
        return float(format(value, ".2e"))
    return rounded


def _legacy_oneway(result):
    result["beta_full"] = result["beta"]
    result["beta"] = _stata_default_number(result["beta"])
    result["legacy_inference"] = {}
    for estimator in ("cv1", "cv3", "cv3j", "cv3drop", "cv3jdrop"):
        if estimator not in result["inference"]:
            continue
        full = result["inference"][estimator]
        se = _stata_default_number(full["se"])
        inf = _inference(result["beta"], se, full["df"])
        if math.isfinite(inf["t"]) and full["df"] > 0:
            tail = student_t_sf(inf["t"], full["df"])
            inf["p"] = 2 * min(tail, 1 - tail)
        result["legacy_inference"][estimator] = inf
        if "drop" not in estimator:
            prefix = "cv3J" if estimator == "cv3j" else estimator
            for suffix in ("se", "t", "p", "lci", "uci"):
                result[prefix + suffix + "_full"] = full[suffix]
                result[prefix + suffix] = inf[suffix]
        elif estimator == ("cv3jdrop" if "cv3jdrop" in result["inference"] else "cv3drop"):
            for suffix in ("se", "t", "p", "lci", "uci"):
                result["cv3" + suffix + "drop"] = inf[suffix]
    for key in ("gstarzero", "gstarone", "gstarrho"):
        if key in result:
            result[key + "_full"] = result[key]
            result[key] = _stata_default_number(result[key])


def _summaries(vectors, original_g, addmeans, stata_float=True):
    # Stata's summarize, detail percentiles use type-2 empirical quantiles.
    columns, bonus = [], []
    for j, z in enumerate(vectors):
        z = np.asarray(z, dtype=float)
        if not len(z):
            columns.append(np.full(7, math.nan))
            bonus.append(np.full(6, math.nan))
            continue
        display_z = z.astype(np.float32).astype(float) if stata_float else z
        q1, med, q3 = np.quantile(display_z, [0.25, 0.5, 0.75], method="averaged_inverted_cdf")
        m = float(np.mean(display_z))
        cv = math.sqrt(float(np.sum((z - m) ** 2)) / ((original_g - 1) * m * m)) if m and original_g > 1 else math.nan
        columns.append([np.min(display_z), q1, med, m, q3, np.max(display_z), cv])
        if addmeans:
            with np.errstate(divide="ignore", invalid="ignore"):
                invz = np.where(z != 0, 1 / z, np.nan)
                logz = np.where(z > 0, np.log(z), np.nan)
                hm = 1 / np.nanmean(invz) if np.any(np.isfinite(invz)) else math.nan
                gm = np.exp(np.nanmean(logz)) if np.any(np.isfinite(logz)) else math.nan
                qm = math.sqrt(float(np.mean(z * z)))
                b = [hm, hm / m if m else math.nan, gm, gm / m if m else math.nan,
                     qm, qm / m if m else math.nan]
            if j == 3:  # The authoritative display suppresses beta's H/G means.
                b[:4] = [math.nan] * 4
            bonus.append(b)
    return np.asarray(columns).T, np.asarray(bonus).T if addmeans else None


def _component(x, y, u, code, labels, xx, xy, inverse, beta, xt, xw,
               coeffnum, rank, rankcheck, force, twoway, dimension, dojack=True):
    n, k = x.shape
    groups = len(labels)
    ng = np.zeros(groups)
    lever = np.zeros(groups)
    partlev = np.zeros(groups)
    omitted = np.full((k, groups), np.nan)
    delete_one_solves = 0
    rankfail = np.full(groups, np.nan)
    legacy = np.full(groups, np.nan)
    gam0, gam1 = np.zeros(groups), np.zeros(groups)
    meat = np.zeros((k, k))
    denom = float(xt @ xt)
    contiguous = bool(np.all(code[1:] >= code[:-1]))
    boundaries = np.r_[0, np.flatnonzero(np.diff(code)) + 1, n] if contiguous else None
    for j in range(groups):
        mask = slice(boundaries[j], boundaries[j + 1]) if contiguous else code == j
        xj, yj, uj = x[mask], y[mask], u[mask]
        xjxj = _cross(xj, xj)
        score = xj.T @ uj
        meat += np.outer(score, score)
        a = xx - xjxj
        # CV31 needs the intersection scores, not intersection deletion fits.
        if dojack:
            if twoway and not _identified(a):
                raise SummclustError(f"Focal coefficient is not identified after deleting {dimension} cluster {labels[j]!r}; CV3 is undefined", 498)
            omitted[:, j] = _solve(a, xy - _cross(xj, yj), force, k <= 50)
            delete_one_solves += 1
        ng[j] = boundaries[j + 1] - boundaries[j] if contiguous else int(np.sum(mask))
        lever[j] = np.trace(xjxj @ inverse)
        partlev[j] = float(xt[mask] @ xt[mask]) / denom if denom > 0 else math.nan
        gam0[j] = float(xw[mask] @ xw[mask])
        gam1[j] = float(np.sum(xw[mask])) ** 2
        if dojack:
            legacy[j] = np.any(np.abs(omitted[:, j]) <= 1e-8)
        if rankcheck:
            rankfail[j] = _rank(a) < rank
    adjustment = groups * (n - 1) / ((groups - 1) * (n - coeffnum))
    v1 = adjustment * (inverse @ meat @ inverse)
    if dojack:
        diff = beta[:, None] - omitted
        v3 = (groups - 1) / groups * (diff @ diff.T)
    else:
        v3 = np.full((k, k), np.nan)
    return {"cluster_ids": list(labels), "ng": ng, "lever": lever, "partlev": partlev,
            "betanog": omitted[0], "betas_deleted": omitted, "rankfail": rankfail,
            "legacyflag": legacy, "cv1": v1, "cv3": v3, "gamma0": gam0, "gamma1": gam1,
            "delete_one_solves": delete_one_solves}


def summclust(data: Mapping, outcome: str, focal: str, controls: Sequence[str],
              cluster, *, fevar=None, absorb=None, sample=None, gstar=False,
              rho=None, jackknife=False, rankcheck=False, table=False,
              addmeans=False, regtable=False, mixed=False):
    """Fit OLS and return cluster diagnostics and finite-sample inference.

    Parameters
    ----------
    data : mapping of column names to equal-length one-dimensional sequences
    outcome, focal : names of numeric outcome and coefficient of interest
    controls : ordered names of one or more numeric nuisance regressors
    cluster : one name or a sequence of one/two numeric or string identifier names
    fevar : optional identifier name or ordered sequence of names
        Explicit indicators: first effect retains all levels without absorption;
        subsequent effects omit their last sorted level.
    absorb : optional identifier name
        One within-transformed effect, nested in every clustering dimension.
    sample : optional boolean mask over incoming rows
        Intersected with complete cases of every model/identifier column.
    gstar, rho : optional effective-cluster diagnostics; rho in [0.00001, 1]
        Omitted rho uses Stata's legacy .1234567 when gstar is requested.
    jackknife : include centered CV3J (one-way only)
    rankcheck : compute true deleted-design rank-loss flags
    table, addmeans, regtable : include numerical presentation tables
    mixed : bool, default False
        Two-way CV31 = CV3_G + CV3_H - CV1_I. Skips intersection deletion
        fits; their diagnostics and pure CV3 covariance are NaN. rankcheck
        still computes ranks without deletion coefficient solves. Reported
        CV31 inference uses max(CV31, CV3_G, CV3_H) with min(G,H)-1 df.

    Returns a dict. Standard numerical fields follow Stata r() names; metadata
    includes coefficient_names, cluster_ids, sample_mask and components. One-way
    b and N/rank metadata are useful additions. Two-way cv1/cv3 are conventional
    inclusion-exclusion covariances; their max-SE inference is stored separately.
    """
    if not isinstance(data, Mapping):
        raise SummclustError("data must be a mapping of named columns")
    controls, clusters, fes = _names(controls), _names(cluster), _names(fevar)
    if not controls:
        raise SummclustError("At least one control is required by the Stata varlist interface")
    if len(clusters) not in (1, 2):
        raise SummclustError("Specify one or two clustering columns")
    two = len(clusters) == 2
    if not isinstance(mixed, (bool, np.bool_)):
        raise SummclustError("mixed must be a boolean")
    if mixed and not two:
        raise SummclustError("mixed is available only with two-way clustering")
    if two and jackknife:
        raise SummclustError("Centered two-way CV3J is not supported; use CV3")
    if rho is not None and (not np.isfinite(rho) or not 0.00001 <= rho <= 1):
        raise SummclustError("rho must be in [0.00001, 1]")
    dogstar = gstar or rho is not None
    rho_value = 0.1234567 if rho is None else float(rho)
    columns = list(dict.fromkeys([outcome, focal] + controls + clusters + fes + ([absorb] if absorb else [])))
    absent = [c for c in columns if c not in data]
    if absent:
        raise SummclustError(f"Missing data columns: {absent}")
    nraw = len(data[outcome])
    arrays = {c: np.asarray(data[c]) for c in columns}
    if any(a.ndim != 1 or len(a) != nraw for a in arrays.values()):
        raise SummclustError("Model columns must be equal-length one-dimensional sequences")
    mask = np.ones(nraw, dtype=bool)
    for c in columns:
        if arrays[c].dtype.kind in "biuf":
            mask &= np.isfinite(arrays[c])
        else:
            mask &= np.asarray([not _missing(v) for v in arrays[c]])
    if sample is not None:
        sample = np.asarray(sample)
        if sample.ndim != 1 or len(sample) != nraw or sample.dtype.kind != "b":
            raise SummclustError("sample must be a boolean mask with one entry per input row")
        mask &= sample
    if not np.any(mask):
        raise SummclustError("No complete observations in the requested sample", 2000)
    numeric = {}
    for c in [outcome, focal] + controls:
        try:
            numeric[c] = np.asarray(arrays[c][mask], dtype=float)
        except (ValueError, TypeError) as exc:
            raise SummclustError(f"Regression column {c!r} must be numeric") from exc
    ids, levels = {}, {}
    for c in clusters + fes + ([absorb] if absorb else []):
        ids[c], levels[c] = _encode(arrays[c][mask])
    for c in clusters:
        if len(levels[c]) < 2:
            raise SummclustError("At least two nonempty clusters are required in each dimension", 498)
    if absorb and any(not _nested(ids[absorb], ids[c]) for c in clusters):
        raise SummclustError("absorb groups must be nested within every clustering dimension; use fevar")
    y = numeric[outcome].copy()
    names, xs = [focal], [numeric[focal]]
    if not fes and not absorb:
        names.append("_cons")
        xs.append(np.ones(len(y)))
    names += controls
    xs += [numeric[c] for c in controls]
    for j, f in enumerate(fes):
        keep = len(levels[f]) - int(bool(absorb) or j > 0)
        for ell in range(keep):
            xs.append((ids[f] == ell).astype(float))
            value = ell + 1 if two else levels[f][ell]
            if isinstance(value, (float, np.floating)) and value.is_integer():
                value = int(value)
            names.append(f"{f}_{value}")
    x = np.column_stack(xs)
    del xs, numeric
    if absorb:
        for j in range(len(levels[absorb])):
            use = ids[absorb] == j
            x[use] -= np.mean(x[use], axis=0)
            y[use] -= np.mean(y[use])
    # Stable grouping also matches one-way accumulation/order conventions.
    order_keys = [ids[clusters[0]]]
    if absorb:
        order_keys.append(ids[absorb])
    order = np.lexsort(tuple(reversed(order_keys))) if not two else np.arange(len(y))
    if not two:
        x, y = x[order], y[order]
        ids = {c: z[order] for c, z in ids.items()}
    n, k = x.shape
    rank = _rank(x)
    coeffnum = rank + (len(levels[absorb]) if absorb else 0)
    if n <= coeffnum:
        raise SummclustError("Residual degrees of freedom are not positive after fixed effects", 2001)
    xx, xy = _cross(x, x), _cross(x, y)
    force = np.linalg.cond(xx) > 1e7
    if force and not two:
        xx, xy = np.zeros((k, k)), np.zeros(k)
        for j in range(len(levels[clusters[0]])):
            use = ids[clusters[0]] == j
            xg = x[use]
            xx += _cross(xg, xg)
            xy += _cross(xg, y[use])
    if two and not _identified(xx):
        raise SummclustError("Focal coefficient is not numerically identified in the full model", 498)
    inverse = _invsym(xx)
    beta = inverse @ xy
    residual = y - x @ beta
    other = x[:, 1:]
    gamma = _solve(_cross(other, other), _cross(other, x[:, 0]), check_condition=False)
    xt = x[:, 0] - other @ gamma
    xw = x @ inverse[:, 0]
    dimensions = [("g", ids[clusters[0]], levels[clusters[0]])]
    if two:
        intersection, ilabels = _encode(list(zip(ids[clusters[0]], ids[clusters[1]])))
        actual = [(levels[clusters[0]][a], levels[clusters[1]][b]) for a, b in ilabels]
        dimensions += [("h", ids[clusters[1]], levels[clusters[1]]), ("i", intersection, actual)]
    components = {}
    for dim, code, labels in dimensions:
        components[dim] = _component(x, y, residual, code, labels, xx, xy, inverse,
                                     beta, xt, xw, coeffnum, rank, rankcheck, force, two, dim.upper(),
                                     dojack=not (mixed and dim == "i"))
    result = {"b": beta, "beta": float(beta[0]), "N": n, "k": k, "rank": rank,
              "coefficient_names": names, "sample_mask": mask, "sample_indices": np.flatnonzero(mask),
              "coeffnum": coeffnum, "components": components, "rho": rho_value,
              "cluster_names": clusters, "warnings": [], "inference": {},
              "mixed": bool(mixed), "selected_estimator": "CV31" if mixed else ("CV3J" if jackknife else "CV3")}
    if not two:
        g = components["g"]
        groups = len(g["ng"])
        result.update({key: g[key] for key in ("ng", "lever", "partlev", "betanog", "cv1", "cv3")})
        result.update(G=groups, df_r=groups - 1, cluster_ids=g["cluster_ids"])
        if rankcheck:
            result.update(rankfail=g["rankfail"], legacyflag=g["legacyflag"])
        if jackknife:
            centered = g["betas_deleted"] - np.mean(g["betas_deleted"], axis=1)[:, None]
            result["cv3j"] = (groups - 1) / groups * centered @ centered.T
        good = g["legacyflag"] == 0
        gdrop = int(np.sum(good))
        result["gdrop"] = gdrop
        result["singclustfull"] = groups - gdrop
        if gdrop < groups:
            result["warnings"].append("Legacy near-zero coefficient flags affect dropped-subsample summaries; they are not rank tests")
            deleted_good = g["betas_deleted"][:, good]
            if gdrop:
                diff = beta[:, None] - deleted_good
                drop = (gdrop - 1) / gdrop * diff @ diff.T
                center = deleted_good - np.mean(deleted_good, axis=1)[:, None]
                jdrop = (gdrop - 1) / gdrop * center @ center.T
            else:
                drop = jdrop = np.full((k, k), math.nan)
            result["drop_cluster_ids"] = [label for label, keep in zip(g["cluster_ids"], good) if keep]
            result["inference"]["cv3drop"] = _inference(beta[0], _safe_sqrt(drop[0, 0]), gdrop - 1)
            if jackknife:
                result["inference"]["cv3jdrop"] = _inference(beta[0], _safe_sqrt(jdrop[0, 0]), gdrop - 1)
            if regtable:
                result["cv3jdrop" if jackknife else "cv3drop"] = jdrop if jackknife else drop
            # Stata's legacy macros are overwritten by CV3J-drop with jackknife.
            dinf = result["inference"]["cv3jdrop" if jackknife else "cv3drop"]
            for suffix in ("se", "t", "p", "lci", "uci"):
                result["cv3" + suffix + "drop"] = dinf[suffix]
        vec = [g["ng"], g["lever"], g["partlev"], g["betanog"]]
        if 0 < gdrop < groups:
            vec.append(g["betanog"][good])
        result["clustsum"], bonus = _summaries(vec, groups, addmeans)
        result["clustsum_full"], _ = _summaries(vec, groups, False, stata_float=False)
        if addmeans:
            result["bonus"] = bonus
        if table and groups <= 52:
            result["scall"] = np.column_stack(vec[:4])
        if dogstar:
            result["gstarzero"] = _gstar(g["gamma0"])
            if not absorb and not any(_nested(ids[f], ids[clusters[0]]) for f in fes):
                result["gstarone"] = _gstar(g["gamma1"])
                result["gstarrho"] = _gstar(rho_value * g["gamma1"] + (1 - rho_value) * g["gamma0"])
        estimators = ["cv1", "cv3"] + (["cv3j"] if jackknife else [])
        for estimator in estimators:
            inf = _inference(beta[0], _safe_sqrt(result[estimator][0, 0]), groups - 1)
            result["inference"][estimator] = inf
            prefix = "cv3J" if estimator == "cv3j" else estimator
            for suffix in ("se", "t", "p", "lci", "uci"):
                result[prefix + suffix] = inf[suffix]
    else:
        g, h, i = (components[d] for d in "ghi")
        G, H, I = (len(z["ng"]) for z in (g, h, i))
        result.update(G=G, H=H, I=I, N_clust_g=G, N_clust_h=H, N_clust_intersection=I,
                      df_r=min(G, H) - 1, max_tie_tolerance=1e-10,
                      cluster_ids={d: c["cluster_ids"] for d, c in components.items()}, vce="two-way max-se")
        dimsummary = []
        for dim, code, _ in dimensions:
            c = components[dim]
            for stat in ("ng", "lever", "partlev", "betanog", "rankfail"):
                result[f"{stat}_{dim}"] = c[stat]
                if dim == "i":
                    result[f"{stat}_intersection"] = c[stat]
            gs0 = _gstar(c["gamma0"], True)
            result[f"gstarzero_{dim}"] = gs0
            if dim == "i":
                result["gstarzero_intersection"] = gs0
            if dogstar and not absorb and not any(_nested(ids[f], code) for f in fes):
                result[f"gstarone_{dim}"] = _gstar(c["gamma1"], True)
                result[f"gstarrho_{dim}"] = _gstar(rho_value * c["gamma1"] + (1 - rho_value) * c["gamma0"], True)
                if dim == "i":
                    result["gstarone_intersection"] = result["gstarone_i"]
                    result["gstarrho_intersection"] = result["gstarrho_i"]
            dimsummary.append([_cv(c[s], twoway=True) for s in ("ng", "lever", "partlev", "betanog")] + [len(c["ng"]), gs0])
            if table and dim != "i":
                result[f"table_{dim}"] = np.column_stack([c[s] for s in ("ng", "lever", "partlev", "betanog")])
        result["dimsummary"] = np.asarray(dimsummary)
        result["intersection_delete_one_solves"] = i["delete_one_solves"]
        result["intersection_delete_solves"] = i["delete_one_solves"]
        result["estimator_label"] = "CV31 (mixed)" if mixed else "CV3"
        result["inference_label"] = "CV31 (mixed, max)" if mixed else "CV3 (max)"
        result["inference_conventional"] = {}
        # Preserve the meanings of ordinary CV3 fields when it is uncomputed.
        if mixed:
            result["cv3"] = result["cv3_conventional"] = np.full((k, k), np.nan)
            for dim, c in components.items():
                result[f"V_cv3_{dim}"] = c["cv3"]
                result[f"cv3_{'intersection' if dim == 'i' else dim}"] = c["cv3"]
        ses = []
        for estimator in ("cv1", "cv31" if mixed else "cv3"):
            primary = "cv3" if estimator == "cv31" else estimator
            intersection_v = i["cv1" if estimator == "cv31" else estimator]
            conventional = g[primary] + h[primary] - intersection_v
            if not all(np.all(np.isfinite(v)) for v in (g[primary], h[primary], intersection_v)):
                raise SummclustError("A two-way covariance component is missing", 498)
            result[estimator] = result[estimator + "_conventional"] = conventional
            if estimator != "cv31":
                for dim, c in components.items():
                    result[f"V_{estimator}_{dim}"] = c[estimator]
                    result[f"{estimator}_{'intersection' if dim == 'i' else dim}"] = c[estimator]
            candidates = np.array([np.diag(conventional), np.diag(g[primary]), np.diag(h[primary])])
            maximum = np.max(candidates, axis=0)
            sese = np.asarray([_safe_sqrt(v) for v in maximum])
            ses.append(sese)
            variance = float(maximum[0])
            which = int(np.argmax(candidates[:, 0]))
            ties = sum(bit for bit, value in zip((1, 2, 4), candidates[:, 0]) if abs(value - variance) <= 1e-10 * max(abs(variance), 1e-300))
            result[estimator + "max_component"] = ("two-way", "G", "H")[which]
            result[estimator + "max_ties"] = ties
            result[estimator + "_conventional_var"] = float(conventional[0, 0])
            result[estimator + "_conventional_se"] = _safe_sqrt(conventional[0, 0])
            result["inference_conventional"][estimator] = _inference(beta[0], _safe_sqrt(conventional[0, 0]), result["df_r"])
            for alias in ("_max_var", "var_max"):
                result[estimator + alias] = variance
            for alias in ("_max_se", "se_max"):
                result[estimator + alias] = float(sese[0])
            inf = _inference(beta[0], sese[0], result["df_r"])
            result["inference"][estimator] = inf
            for suffix in ("se", "t", "p", "lci", "uci"):
                result[estimator + suffix] = inf[suffix]
                if estimator == "cv31":
                    result[estimator + "_max_" + suffix] = inf[suffix]
            if conventional[0, 0] < 0:
                result["warnings"].append(f"{estimator.upper()} conventional focal variance is negative; max inference uses a separate component")
        result["regsum"] = np.column_stack([beta] + ses)
        result["regsum_columns"] = ["beta", "CV1 (max)", result["inference_label"]]
    if regtable and not two:
        estimator = "cv3j" if jackknife else "cv3"
        result["regresstab"] = np.array([[inf[s] for s in ("beta", "se", "t", "p", "lci", "uci")]
                                        for inf in (_inference(b, _safe_sqrt(v), result["df_r"])
                                                    for b, v in zip(beta, np.diag(result[estimator])))])
    if not two:
        _legacy_oneway(result)
    return result


def read_csv(path, *, string_columns=()):
    """Read one shared CSV using stdlib csv; empty/dot numeric cells become NaN.

    Declare string identifier columns explicitly; all other columns are parsed
    as float64. Numeric precision is retained by Python's float conversion.
    """
    string_columns = set(string_columns)
    with Path(path).open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not reader.fieldnames:
            raise SummclustError("CSV has no column header")
        data = {name: [] for name in reader.fieldnames}
        for row in reader:
            for name, value in row.items():
                if name in string_columns:
                    data[name].append(value)
                else:
                    data[name].append(math.nan if value in (None, "", ".") else float(value))
    return {name: np.asarray(values, dtype=object if name in string_columns else float)
            for name, values in data.items()}


def flatten_summclust(result):
    """Yield numeric results as statistic,row,column,value dictionaries.

    Coefficient matrices use coefficient names; cluster vectors use actual
    identifiers (intersection labels are repr tuples). Supplemental working
    arrays and nested objects are excluded; scalars are not display-rounded.
    """
    skip = {"sample_mask", "sample_indices", "coefficient_names", "cluster_ids", "cluster_names", "warnings", "components", "inference", "drop_cluster_ids"}
    coeffs = result["coefficient_names"]
    for stat, value in result.items():
        if stat in skip or isinstance(value, (str, Mapping)):
            continue
        arr = np.asarray(value)
        if arr.dtype.kind not in "biuf":
            continue
        if arr.ndim == 0:
            yield {"statistic": stat, "row": "", "column": "", "value": float(arr)}
        elif arr.ndim == 1:
            if stat == "b":
                labels = coeffs
            elif stat in ("ng", "lever", "partlev", "betanog", "legacyflag", "rankfail"):
                labels = result["cluster_ids"]
            else:
                dim = "i" if stat.endswith("_intersection") else stat.rsplit("_", 1)[-1]
                labels = result["cluster_ids"].get(dim, list(range(1, len(arr) + 1))) if isinstance(result["cluster_ids"], Mapping) else list(range(1, len(arr) + 1))
            for label, number in zip(labels, arr):
                yield {"statistic": stat, "row": str(label), "column": "", "value": float(number)}
        elif arr.ndim == 2:
            covariance = stat.startswith("cv") or stat.startswith("V_cv")
            rows = coeffs if covariance or stat in ("regsum", "regresstab") else list(range(1, arr.shape[0] + 1))
            cols = coeffs if covariance else list(range(1, arr.shape[1] + 1))
            for a, row in enumerate(rows):
                for b, col in enumerate(cols):
                    yield {"statistic": stat, "row": str(row), "column": str(col), "value": float(arr[a, b])}


def export_rows(result, specification):
    """Yield the shared comparison CSV schema with one-based matrix indexes.

    Fields are specification,statistic,row,column,value,row_label,column_label.
    Coefficient labels accompany b and covariance arrays. Cluster identity is
    available in result['components'][dimension]['cluster_ids']; indexes alone
    must not be used to align unrelated cluster sort conventions.
    """
    skip = {"sample_mask", "sample_indices", "coefficient_names", "cluster_ids", "cluster_names", "warnings", "components", "inference", "drop_cluster_ids"}
    result = dict(result)
    if "H" not in result:
        result["summary"] = result["clustsum"]
        result["cv1se_full"] = result["inference"]["cv1"]["se"]
        result["cv3se_full"] = result["inference"]["cv3"]["se"]
    for stat, value in result.items():
        if stat in skip or isinstance(value, (str, Mapping)):
            continue
        arr = np.asarray(value)
        if arr.dtype.kind not in "biuf":
            continue
        if arr.ndim == 0:
            arr = arr.reshape(1, 1)
        elif arr.ndim == 1:
            arr = arr.reshape(1, -1) if stat == "b" else arr[:, None]
        if arr.ndim != 2:
            continue
        for (i, j), number in np.ndenumerate(arr):
            covariance = stat.startswith(("cv", "V_cv")) and arr.shape == (result["k"], result["k"])
            rowname = result["coefficient_names"][i] if covariance else ""
            colname = result["coefficient_names"][j] if covariance or stat == "b" else ""
            yield dict(specification=specification, statistic=stat, row=i + 1, column=j + 1,
                       value=float(number), row_label=rowname, column_label=colname)
