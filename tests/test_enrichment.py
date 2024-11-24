import unittest

from kinex.kinex import Kinex
import pandas as pd
from tests.data import get_test_input_sites


class TestEnrichment(unittest.TestCase):
    def setUp(self):
        self.kinex = Kinex()
        self.enrich = self.kinex.get_enrichment(
            get_test_input_sites(),
            fc_threshold=1.5,
            phospho_priming=False,
            favorability=True,
            method="max",
        )

    def test_total_upregulated(self):
        self.assertEqual(self.enrich.ser_thr.total_upregulated, 1)

    def test_total_downregulated(self):
        self.assertEqual(self.enrich.ser_thr.total_downregulated, 0)

    def test_total_unregulated(self):
        self.assertEqual(self.enrich.ser_thr.total_unregulated, 4)

    def test_enrichment_table_structure(self):
        expected_columns = [
            "upregulated",
            "downregulated",
            "unregulated",
            "upregulated_enrichment_value",
            "upregulated_enrichment_value_log2",
            "upregulated_p_value",
            "upregulated_p_value_log10_abs",
            "upregulated_adjusted_p_value",
            "upregulated_adjusted_p_value_log10_abs",
            "downregulated_enrichment_value",
            "downregulated_enrichment_value_log2",
            "downregulated_p_value",
            "downregulated_p_value_log10_abs",
            "downregulated_adjusted_p_value",
            "downregulated_adjusted_p_value_log10_abs",
            "dominant_direction",
            "dominant_enrichment_value_log2",
            "dominant_p_value_log10_abs",
            "dominant_adjusted_p_value_log10_abs",
        ]
        self.assertListEqual(
            list(self.enrich.ser_thr.enrichment_table.columns), expected_columns
        )

    def test_adjust_background_sites(self):
        # Test when total_unregulated is zero
        self.enrich.ser_thr.total_upregulated = 4
        self.enrich.ser_thr.total_downregulated = 2
        self.enrich.ser_thr.total_unregulated = 0
        self.enrich.ser_thr.adjust_background_sites()
        self.assertEqual(self.enrich.ser_thr.total_unregulated, 1)

        # Test when total_unregulated is not zero
        self.enrich.ser_thr.total_unregulated = 3
        self.enrich.ser_thr.adjust_background_sites()
        self.assertEqual(self.enrich.ser_thr.total_unregulated, 3)

    def test_calculate_enrichment_for_row(self):
        # Test case 1: No unregulated hits
        self.enrich.ser_thr.enrichment_table = pd.DataFrame(
            {
                "kinase": ["kinase1"],
                "upregulated": [1],
                "downregulated": [0],
                "unregulated": [0],  # 1
            }
        )

        self.enrich.ser_thr._calculate_enrichment_for_row(0)
        self.assertEqual(
            self.enrich.ser_thr.enrichment_table.loc[0, "upregulated_enrichment_value"],
            7,
        )

        # Test case 2: No upregulated hits
        self.enrich.ser_thr.enrichment_table = pd.DataFrame(
            {
                "kinase": ["kinase2"],
                "upregulated": [0],
                "downregulated": [1],
                "unregulated": [0],
            }
        )

        self.enrich.ser_thr._calculate_enrichment_for_row(0)
        self.assertEqual(
            self.enrich.ser_thr.enrichment_table.loc[0, "upregulated_enrichment_value"],
            0,
        )
        self.assertEqual(
            self.enrich.ser_thr.enrichment_table.loc[0, "upregulated_p_value"], 1
        )

        # Test case 3: No unregulated and no downregulatedhits
        self.enrich.ser_thr.enrichment_table = pd.DataFrame(
            {
                "kinase": ["kinase5"],
                "upregulated": [0],
                "downregulated": [0],
                "unregulated": [0],
            }
        )

        self.enrich.ser_thr._calculate_enrichment_for_row(0)
        self.assertEqual(
            self.enrich.ser_thr.enrichment_table.loc[
                0, "downregulated_enrichment_value"
            ],
            0,
        )

    def test_determine_dominant_direction(self):
        self.enrich.ser_thr.enrichment_table = pd.DataFrame(
            {
                "kinase": ["kinase1"],
                "upregulated_enrichment_value": [2],
                "downregulated_enrichment_value": [1],
                "upregulated_enrichment_value_log2": [1],
                "downregulated_enrichment_value_log2": [0.5],
                "upregulated_p_value_log10_abs": [1],
                "downregulated_p_value_log10_abs": [0.5],
                "dominant_direction": [""],
            }
        )
        self.enrich.ser_thr._determine_dominant_direction(0)
        self.assertEqual(
            self.enrich.ser_thr.enrichment_table.loc[0, "dominant_direction"],
            "upregulated set",
        )

    def test_reindex_missing_kinases(self):
        self.enrich.ser_thr.enrichment_table = pd.DataFrame(
            {"kinase": ["kinase1"]}
        ).set_index("kinase")
        self.enrich.ser_thr.all_kinases = {"kinase1", "kinase2"}
        self.enrich.ser_thr._reindex_missing_kinases()
        self.assertIn("kinase2", self.enrich.ser_thr.enrichment_table.index)


if __name__ == "__main__":
    unittest.main()
