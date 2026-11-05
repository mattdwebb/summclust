"""String identifiers must describe the same regression as numeric group codes."""
import math
import pickle
import unittest

import numpy as np

from summclust import SummclustError, summclust


def fixture():
    j = np.arange(288)
    g, h = j // 36, j % 6
    x = np.sin(.23 * j) + .1 * g
    z = np.cos(.17 * j) + .03 * h
    w = np.sin(.07 * j + .2)
    y = 2 + .7 * x - .4 * z + .12 * w + .4 * np.sin(.73 * j) + .03 * g
    return dict(y=y, x=x, z=z, w=w, g=g, h=h,
                f=(j // 6) % 3 + 1, q=j % 5 + 1, a=g * 6 + h + 1)


def string_copy(data, identifiers):
    result = {key: value.copy() for key, value in data.items()}
    for name in identifiers:
        result[name] = np.asarray([f"group {int(v):03d} / café" for v in data[name]], dtype=object)
    return result


class StringFixedEffectsTests(unittest.TestCase):
    def fit(self, data, cluster="g", **options):
        return summclust(data, "y", "x", ["z", "w"], cluster,
                         rankcheck=True, gstar=True, rho=.4,
                         table=True, regtable=True, **options)

    def assert_result_equal(self, expected, actual, ignored=()):
        """Compare every returned numerical field, including nested diagnostics."""
        if isinstance(expected, dict):
            self.assertEqual(expected.keys(), actual.keys())
            for key in expected:
                if key not in ignored:
                    self.assert_result_equal(expected[key], actual[key], ignored)
        elif isinstance(expected, np.ndarray):
            np.testing.assert_array_equal(expected, actual)
        elif isinstance(expected, (float, np.floating)) and math.isnan(expected):
            self.assertTrue(math.isnan(actual))
        else:
            self.assertEqual(expected, actual)

    def test_equivalent_groups_in_all_modes(self):
        numeric = fixture()
        strings = string_copy(numeric, ("f", "q", "a"))
        before = pickle.dumps(strings)
        designs = [dict(fevar="f"), dict(fevar=["f", "q"]),
                   dict(absorb="a"), dict(fevar=["f", "q"], absorb="a")]
        for cluster, mode in [("g", dict(jackknife=True)),
                              (["g", "h"], {}),
                              (["g", "h"], dict(mixed=True))]:
            for design in designs:
                with self.subTest(cluster=cluster, mode=mode, design=design):
                    expected = self.fit(numeric, cluster, **mode, **design)
                    actual = self.fit(strings, cluster, **mode, **design)
                    self.assert_result_equal(expected, actual)
        self.assertEqual(before, pickle.dumps(strings))

    def test_sample_and_missing_values_precede_level_encoding(self):
        numeric = fixture()
        strings = string_copy(numeric, ("f", "q", "a"))
        for name in ("f", "q", "a"):
            numeric[name] = numeric[name].astype(float)
        for row, missing in [(1, ""), (2, None), (3, np.nan)]:
            strings["f"][row] = missing
            numeric["f"][row] = np.nan
        strings["q"][4] = ""
        numeric["q"][4] = np.nan
        strings["a"][5] = None
        numeric["a"][5] = np.nan
        numeric["y"][6] = strings["y"][6] = np.nan
        sample = np.arange(288) % 13 != 0
        # This level sorts first, but exists only outside the requested sample.
        strings["f"][0] = "000 excluded level"
        numeric["f"][0] = -999
        expected_mask = sample.copy()
        expected_mask[1:7] = False
        before = pickle.dumps(strings)
        for cluster, mode in [("g", {}), (["g", "h"], dict(mixed=True))]:
            with self.subTest(cluster=cluster):
                options = dict(fevar=["f", "q"], absorb="a", sample=sample, **mode)
                expected = self.fit(numeric, cluster, **options)
                actual = self.fit(strings, cluster, **options)
                self.assert_result_equal(expected, actual)
                np.testing.assert_array_equal(actual["sample_mask"], expected_mask)
                self.assertIn("f_1", actual["coefficient_names"])
                self.assertNotIn("f_3", actual["coefficient_names"])
        self.assertEqual(before, pickle.dumps(strings))

    def test_nonempty_numeric_looking_strings_are_categories(self):
        numeric = fixture()
        strings = string_copy(numeric, ("a",))
        strings["f"] = np.asarray([("1e999", "Inf", "NaN")[v - 1] for v in numeric["f"]], dtype=object)
        # Exercise the same missing-value handling for clusters and absorption.
        for name in ("g", "a"):
            labels = [f"label {j:03d}" for j in np.unique(numeric[name])]
            labels[:3] = ["1e999", "Inf", "NaN"]
            strings[name] = np.asarray([labels[int(v - numeric[name].min())] for v in numeric[name]], dtype=object)
            ordered = sorted(set(strings[name]))
            numeric[name] = np.asarray([ordered.index(v) for v in strings[name]])
        for cluster, mode in [("g", {}), (["g", "h"], dict(mixed=True))]:
            for design in (dict(fevar="f"), dict(absorb="a")):
                with self.subTest(cluster=cluster, design=design):
                    actual = self.fit(strings, cluster, **mode, **design)
                    expected = self.fit(numeric, cluster, **mode, **design)
                    self.assertEqual(actual["N"], 288)
                    self.assert_result_equal(expected, actual, ignored=("cluster_ids", "drop_cluster_ids"))

    def test_numeric_regression_missing_values_still_drop(self):
        data = fixture()
        data["f"] = np.asarray(["NaN"] * 288, dtype=object)
        data["y"] = data["y"].astype(object)
        data["y"][0:4] = [None, np.nan, math.inf, "NaN"]
        actual = self.fit(data, fevar="f")
        expected = self.fit({key: value[4:] for key, value in data.items()}, fevar="f")
        self.assertEqual(actual["N"], 284)
        for key in ("b", "cv1", "cv3", "partlev", "betanog"):
            np.testing.assert_array_equal(actual[key], expected[key])

    def test_plain_identifier_lists_preserve_real_missing_values(self):
        for with_none in (False, True):
            data = fixture()
            data["f"] = ["nan", "ordinary label", "z"] * 96
            data["f"][0] = np.nan
            data["f"][1] = ""
            if with_none:
                data["f"][2] = None
            complete = np.ones(288, dtype=bool)
            complete[:3 if with_none else 2] = False
            reference = {key: np.asarray(value, dtype=object if key == "f" else None)[complete]
                         for key, value in data.items()}
            before = pickle.dumps(data)
            for cluster, mode in [("g", {}), (["g", "h"], dict(mixed=True))]:
                with self.subTest(with_none=with_none, cluster=cluster):
                    actual = self.fit(data, cluster, fevar="f", **mode)
                    expected = self.fit(reference, cluster, fevar="f", **mode)
                    np.testing.assert_array_equal(actual["sample_mask"], complete)
                    self.assert_result_equal(expected, actual,
                                             ignored=("sample_mask", "sample_indices"))
                    self.assertEqual(actual["N"], int(complete.sum()))
                    self.assertEqual(actual["coefficient_names"],
                                     ["x", "z", "w", "f_1", "f_2", "f_3"])
            self.assertEqual(before, pickle.dumps(data))

    def test_long_string_effect_names_are_distinct(self):
        data = fixture()
        first = "category_name_with_shared_prefxA"
        second = "category_name_with_shared_prefxB"
        self.assertEqual(len(first), 32)
        data[first] = string_copy(data, ("f",))["f"]
        data[second] = string_copy(data, ("q",))["q"]
        for cluster, mode in [("g", {}), (["g", "h"], dict(mixed=True))]:
            with self.subTest(cluster=cluster):
                out = self.fit(data, cluster, fevar=[first, second], **mode)
                names = out["coefficient_names"]
                self.assertTrue(all(len(name) <= 32 for name in names))
                self.assertEqual(len(names), len(set(names)))
                self.assertTrue(any(name.endswith("_f1_1") for name in names))
                self.assertTrue(any(name.endswith("_f2_1") for name in names))

    def test_numeric_effect_labels_are_preserved(self):
        data = fixture()
        data["f"] = np.asarray([-4.5, 2.5, 9.5])[data["f"] - 1]
        out = self.fit(data, fevar="f")
        self.assertEqual(out["coefficient_names"], ["x", "z", "w", "f_-4.5", "f_2.5", "f_9.5"])

    def test_invalid_nesting_is_still_rejected(self):
        data = fixture()
        data["a"] = string_copy(data, ("f",))["f"]
        for cluster, mode in [("g", {}), (["g", "h"], {}), (["g", "h"], dict(mixed=True))]:
            with self.subTest(cluster=cluster, mode=mode):
                with self.assertRaisesRegex(SummclustError, "nested"):
                    self.fit(data, cluster, absorb="a", **mode)


if __name__ == "__main__":
    unittest.main()
