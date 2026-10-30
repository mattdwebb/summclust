"""Deterministic, dependency-light numerical checks. Run: python -m unittest -v."""
import math
import unittest
import numpy as np
from summclust import (summclust, SummclustError, student_t_sf, student_t_isf,
                       _invsym, _rank, _stata_default_number, _cross,
                       _mata_dot, _solve)


def fixture():
    j = np.arange(180)
    g = j // 15
    h = j % 9
    x = np.sin(.23 * j) + .1 * g
    z = np.cos(.17 * j) + .03 * h
    w = np.sin(.07 * j + .2)
    y = 2 + .7 * x - .4 * z + .12 * w + .4 * np.sin(.73 * j) + .03 * g
    return dict(y=y, x=x, z=z, w=w, g=g, h=h,
                a=j // 5, f=j % 4, gs=np.array([f"group{n}" for n in g]))


class StudentTTests(unittest.TestCase):
    def test_known_distributions(self):
        for t in [-100, -5, -1, -.0001, 0, .0001, 1, 5, 100]:
            self.assertAlmostEqual(student_t_sf(t, 1), .5 - math.atan(t) / math.pi, places=14)
            self.assertAlmostEqual(student_t_sf(t, 2), .5 - t / (2 * math.sqrt(t * t + 2)), places=14)

    def test_tail_roundtrip(self):
        for df in [1, 2, 5, 19, 42, 99, 499, 1000]:
            for p in [1e-12, 1e-8, .0001, .001, .025, .25, .5, .975]:
                quantile = student_t_isf(p, df)
                self.assertLessEqual(abs(student_t_sf(quantile, df) - p), 2e-12 * p + 1e-15)

    def test_known_quantiles(self):
        for df, critical in [(1, 12.706204736432095), (2, 4.302652729696142),
                             (5, 2.570581835636305), (42, 2.018081702818444)]:
            self.assertAlmostEqual(student_t_isf(.025, df), critical, places=8)

    def test_small_tail_avoids_cancellation(self):
        self.assertGreater(student_t_sf(10, 499), 0)
        self.assertAlmostEqual(student_t_sf(10, 499) / 6.985962852635183e-22, 1, places=11)

    def test_stata_legacy_precision(self):
        self.assertEqual(_stata_default_number(.65955623145), .6595562)
        self.assertEqual(_stata_default_number(-.0679898123), -.0679898)
        self.assertEqual(_stata_default_number(12345678), 12300000)
        self.assertEqual(_stata_default_number(.0000012345678), .00000123)

    def test_invalid(self):
        for df in [0, -1, math.inf]:
            with self.assertRaises(ValueError):
                student_t_sf(0, df)


class SummclustTests(unittest.TestCase):
    def setUp(self):
        self.d = fixture()

    def test_oneway_against_direct_ols(self):
        d = self.d
        r = summclust(d, "y", "x", ["z", "w"], "g", jackknife=True, gstar=True, rankcheck=True)
        x = np.column_stack((d["x"], np.ones(180), d["z"], d["w"]))
        b = np.linalg.lstsq(x, d["y"], rcond=None)[0]
        u = d["y"] - x @ b
        inverse = np.linalg.inv(x.T @ x)
        omit, meat = [], np.zeros((4, 4))
        for g in sorted(set(d["g"])):
            m = d["g"] == g
            omit.append(np.linalg.lstsq(x[~m], d["y"][~m], rcond=None)[0])
            score = x[m].T @ u[m]
            meat += np.outer(score, score)
        omit = np.array(omit).T
        v1 = (12 * 179 / (11 * 176)) * inverse @ meat @ inverse
        diff = b[:, None] - omit
        center = omit - omit.mean(axis=1)[:, None]
        np.testing.assert_allclose(r["b"], b, atol=2e-13, rtol=1e-12)
        np.testing.assert_allclose(r["cv1"], v1, atol=1e-14, rtol=1e-11)
        np.testing.assert_allclose(r["cv3"], 11 / 12 * diff @ diff.T, atol=1e-14, rtol=1e-10)
        np.testing.assert_allclose(r["cv3j"], 11 / 12 * center @ center.T, atol=1e-14, rtol=1e-10)
        np.testing.assert_allclose(r["betanog"], omit[0], atol=2e-13)
        self.assertAlmostEqual(sum(r["partlev"]), 1, places=13)
        self.assertAlmostEqual(sum(r["lever"]), 4, places=13)
        xt = x[:, 0] - x[:, 1:] @ np.linalg.lstsq(x[:, 1:], x[:, 0], rcond=None)[0]
        expected = np.array([sum(xt[d["g"] == j] ** 2) / sum(xt ** 2) for j in range(12)])
        np.testing.assert_allclose(r["partlev"], expected, atol=2e-14)
        self.assertAlmostEqual(r["gstarzero_full"], 1 / np.sum(expected ** 2), places=12)
        weights = x @ inverse[:, 0]
        gam1 = np.array([np.sum(weights[d["g"] == j]) ** 2 for j in range(12)])
        self.assertAlmostEqual(r["gstarone_full"], np.sum(gam1) ** 2 / np.sum(gam1 ** 2), places=12)

    def test_focal_changes_do_not_change_model(self):
        a = summclust(self.d, "y", "x", ["z", "w"], "g")
        b = summclust(self.d, "y", "z", ["x", "w"], "g")
        ix = [2, 1, 0, 3]
        np.testing.assert_allclose(a["b"][ix], b["b"], atol=1e-13)
        np.testing.assert_allclose(a["cv3"][np.ix_(ix, ix)], b["cv3"], atol=1e-13)

    def test_fixed_effects_absorption(self):
        a = summclust(self.d, "y", "x", ["z", "w"], "g", absorb="a", rankcheck=True, gstar=True)
        b = summclust(self.d, "y", "x", ["z", "w"], "g", fevar="a", rankcheck=True, gstar=True)
        np.testing.assert_allclose(a["b"], b["b"][:3], atol=3e-13)
        np.testing.assert_allclose(a["cv1"], b["cv1"][:3, :3], atol=3e-13)
        np.testing.assert_allclose(a["cv3"], b["cv3"][:3, :3], atol=1e-12)
        np.testing.assert_allclose(a["partlev"], b["partlev"], atol=2e-13)
        self.assertNotIn("gstarone", a)
        self.assertTrue(np.all(b["rankfail"] == 1))
        self.assertTrue(np.all(a["rankfail"] == 0))

    def test_missing_string_and_sample(self):
        d = {key: value.copy() for key, value in self.d.items()}
        d["y"][2] = np.nan
        d["z"][4] = np.nan
        d["gs"][7] = ""
        select = np.arange(180) % 3 != 0
        a = summclust(d, "y", "x", ["z", "w"], "gs", sample=select)
        complete = select & np.isfinite(d["y"]) & np.isfinite(d["z"]) & (d["gs"] != "")
        sub = {key: value[complete] for key, value in d.items()}
        b = summclust(sub, "y", "x", ["z", "w"], "gs")
        np.testing.assert_array_equal(a["sample_mask"], complete)
        np.testing.assert_allclose(a["cv3"], b["cv3"], atol=1e-14)
        self.assertEqual(a["cluster_ids"], sorted(set(sub["gs"])))

    def test_observation_order_and_labels(self):
        d = self.d
        perm = np.random.default_rng(892).permutation(len(d["y"]))
        reordered = {key: value[perm] for key, value in d.items()}
        a = summclust(d, "y", "x", ["z", "w"], "gs", absorb="a", jackknife=True)
        b = summclust(reordered, "y", "x", ["z", "w"], "gs", absorb="a", jackknife=True)
        self.assertEqual(a["cluster_ids"], b["cluster_ids"])
        for key in ("b", "cv1", "cv3", "cv3j", "partlev", "betanog"):
            np.testing.assert_allclose(a[key], b[key], atol=2e-13, rtol=1e-10)

    def test_twoway_components_and_ties(self):
        d = self.d
        d["i"] = d["g"] * 100 + d["h"]
        r = summclust(d, "y", "x", ["z", "w"], ["g", "h"], rankcheck=True, table=True, gstar=True)
        comp = [summclust(d, "y", "x", ["z", "w"], c) for c in ("g", "h", "i")]
        for v in ("cv1", "cv3"):
            expected = comp[0][v] + comp[1][v] - comp[2][v]
            np.testing.assert_allclose(r[v], expected, atol=1e-13)
            self.assertAlmostEqual(r[v + "_max_var"], max(expected[0, 0], comp[0][v][0, 0], comp[1][v][0, 0]), places=13)
        self.assertEqual(r["df_r"], 8)
        same = summclust(d, "y", "x", ["z", "w"], ["g", "g"])
        self.assertEqual(same["cv1max_ties"], 7)
        self.assertEqual(same["cv3max_ties"], 7)
        self.assertEqual(same["cv3max_component"], "two-way")

    def test_singular_nuisance_generalized_inverse(self):
        d = self.d
        d["zdup"] = 2 * d["z"]
        r = summclust(d, "y", "x", ["z", "zdup", "w"], ["g", "h"], rankcheck=True)
        full = summclust(d, "y", "x", ["z", "w"], ["g", "h"])
        self.assertEqual(r["rank"], 4)
        self.assertEqual(r["b"][2], 0)
        self.assertAlmostEqual(r["beta"], full["beta"], places=13)
        self.assertAlmostEqual(r["cv3se"], full["cv3se"], places=13)
        np.testing.assert_array_equal(_invsym(np.array([[1., 2], [2, 4]])), [[0, 0], [0, .25]])
        np.testing.assert_allclose(_invsym(np.diag([1., 1e-17])), np.diag([1., 1e17]))

    def test_guardrails(self):
        with self.assertRaises(SummclustError):
            summclust(self.d, "y", "x", ["z"], ["g", "h"], jackknife=True)
        with self.assertRaises(SummclustError):
            summclust(self.d, "y", "x", ["z"], "g", absorb="h")
        with self.assertRaises(SummclustError):
            summclust(self.d, "y", "x", ["z"], "g", rho=0)
        d = self.d
        d["onlyg"] = (d["g"] == 0).astype(float)
        with self.assertRaises(SummclustError):
            summclust(d, "y", "onlyg", ["z"], ["g", "h"])

    def test_legacy_zero_is_not_rank_failure(self):
        d = self.d
        d["y"] = 1 + d["x"]
        r = summclust(d, "y", "x", ["z"], "g", jackknife=True, rankcheck=True, regtable=True)
        self.assertTrue(np.all(r["legacyflag"] == 1))
        self.assertTrue(np.all(r["rankfail"] == 0))
        self.assertEqual(r["gdrop"], 0)


class NumericalConventionTests(unittest.TestCase):
    def test_compensated_sum_and_rounded_products(self):
        a = np.array([1e16, 1, -1e16])[:, None]
        self.assertEqual(_cross(a, np.ones(3))[0], 1)
        # The products, as in Mata quadcross, are rounded to binary64 before
        # accurate accumulation; exact-product correction would give -2^-54.
        a = np.array([1 + 2 ** -27, 1.0])[:, None]
        b = np.array([1 - 2 ** -27, -1.0])
        self.assertEqual(_cross(a, b)[0], 0)

    def test_triangular_product_accumulation_and_solution(self):
        self.assertEqual(_mata_dot(np.ones(4), np.array([1e16, 1, -1e16, 1])), 2)
        a = np.array([[4., 2, 2], [2, 5, 1], [2, 1, 3]])
        answer = _solve(a, np.array([6., -5., 9.]), check_condition=False)
        np.testing.assert_allclose(answer, [1, -2, 3], atol=3e-15)

    def test_zero_omission_generalized_inverse_identities(self):
        a = np.array([[19., -8, -2, -9, 7], [-8, 8, 0, 8, -2],
                      [-2, 0, 2, -2, 1], [-9, 8, -2, 11, -5],
                      [7, -2, 1, -5, 5]])
        inverse = _invsym(a)
        np.testing.assert_allclose(a @ inverse @ a, a, atol=2e-14)
        np.testing.assert_allclose(inverse @ a @ inverse, inverse, atol=2e-15)
        self.assertEqual(np.count_nonzero(np.diag(inverse)), _rank(a))
        np.testing.assert_array_equal(np.flatnonzero(np.diag(inverse)), [0, 1, 3])


if __name__ == "__main__":
    unittest.main()
