import unittest
import numpy as np
from scipy.stats import fisher_exact
from kinex.table2x2 import Table2x2


class TestTable2x2(unittest.TestCase):
    def test_initialization_valid(self):
        self.valid_table = np.array([[1, 2], [3, 4]])
        t = Table2x2(self.valid_table, shift_zeros=False)
        self.assertTrue(np.array_equal(t.table, self.valid_table))

    def test_initialization_invalid(self):
        with self.assertRaises(ValueError):
            self.invalid_table_1 = np.array([[1, 2, 3], [4, 5, 6]])  # Wrong dimensions
            Table2x2(self.invalid_table_1, shift_zeros=False)
        with self.assertRaises(ValueError):
            self.invalid_table_2 = np.array([[1], [4, 5]])  # Wrong dimensions
            Table2x2(self.invalid_table_2, shift_zeros=False)
        with self.assertRaises(ValueError):
            self.invalid_table_2 = np.array([[1, 2], [5]])  # Wrong dimensions
            Table2x2(self.invalid_table_2, shift_zeros=False)

    def test_initialization_with_shift(self):
        self.valid_table_with_zeros = np.array([[1, 0], [3, 0]])  # Table with zeros
        with self.assertRaises(ValueError):
            Table2x2(self.valid_table_with_zeros, shift_zeros=False)
        Table2x2(self.valid_table_with_zeros, shift_zeros=True)

    def test_odds_ratio(self):
        self.valid_table = np.array([[1, 2], [3, 4]])
        table = Table2x2(self.valid_table, shift_zeros=False)
        expected_odds_ratio = (1 * 4) / (2 * 3)
        self.assertAlmostEqual(table.odds_ratio(), expected_odds_ratio)

    def test_p_val(self):
        self.valid_table = np.array([[1, 2], [3, 4]])
        table = Table2x2(self.valid_table, shift_zeros=False)
        p_value = fisher_exact(self.valid_table, alternative="two-sided").pvalue
        self.assertAlmostEqual(table.p_val("two-sided"), p_value, places=5)


if __name__ == "__main__":
    unittest.main()
