import unittest
import pandas as pd
from kinex.score import Score
import numpy as np

from kinex.sequence import get_sequence_object


class TestScore(unittest.TestCase):

    def setUp(self):

        self.sequence_object = get_sequence_object("QSTPQ")
        # ['-1S', '-2Q', '1P', '2Q']
        # self.scores_results = self.sequence_object.get_sequence_scores(get_test_pssm()) # With two kinase AAK1 ACVR2A
        # {'AAK1': 2.463825 = 0.9257 * 2.5792 * 0.7872 * 1.3109 , 'ACVR2A': 0.281407 = 0.7687 * 1.0702 *  0.3404 * 1.0049}

        data1 = {
            "kinase": ["AAK1", "ACVR2A"],
            "score": [2.463825, 0.281407],
            "log_score": [1.300900, -1.829272],
            "percentile_score": [92.818561, 50.004229],
        }

        data2 = {
            "kinase": ["BRAF", "CDK2"],
            "score": [1.234567, 3.456789],
            "log_score": [0.567890, 1.234567],
            "percentile_score": [75.123456, 85.678901],
        }

        self.scores_results = [
            pd.DataFrame(data1).set_index("kinase"),
            pd.DataFrame(data2).set_index("kinase"),
        ]
        self.score_object = Score(self.sequence_object, self.scores_results)

    def test_median_percentile(self):
        expected_median1 = np.array([92.818561, 50.004229]).mean()
        expected_median2 = np.array([75.123456, 85.678901]).mean()
        self.assertAlmostEqual(
            self.score_object.median_percentile[0], expected_median1, places=6
        )
        self.assertAlmostEqual(
            self.score_object.median_percentile[1], expected_median2, places=6
        )

    def test_promiscuity_index(self):
        self.assertEqual(self.score_object.promiscuity_index(limit=80), [1, 1])
        self.assertEqual(self.score_object.promiscuity_index(limit=49), [2, 2])
        self.assertEqual(self.score_object.promiscuity_index(limit=93), [0, 0])

    def test_top(self):
        result = self.score_object.top(number=1)
        self.assertEqual(len(result[0]), 1)
        self.assertEqual(result[0].index[0], "AAK1")
        self.assertEqual(len(result[1]), 1)
        self.assertEqual(result[1].index[0], "BRAF")

    def test_repr(self):
        self.assertEqual(repr(self.score_object), "Scoring results for QSTPQ")


if __name__ == "__main__":
    unittest.main()
