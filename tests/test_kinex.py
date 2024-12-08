import unittest

from kinex import Kinex

from kinex.resources import get_pssm_ser_thr, get_pssm_tyr

from tests.data import get_test_input_sites


class TestKinex(unittest.TestCase):
    def setUp(self):
        self.kinex = Kinex()

    def test_sequence_initialization(self):
        with self.assertRaises(ValueError):
            self.kinex.get_score("GRNSLPVQA")

        with self.assertRaises(ValueError):
            self.kinex.get_score("GRNSLPVQAI")

        with self.assertRaises(ValueError):
            self.kinex.get_score("GRNSL*PVQAI")

        with self.assertRaises(ValueError):
            self.kinex.get_score("GRNSL(ph)PVQAI")

        with self.assertRaises(ValueError):
            self.kinex.get_score("SS")

        # with self.assertRaises(ValueError):
        #     self.kinex.get_score("APQST*PAB")
        # B is not allowed
        # Should be ValueError and invalid seq, but it is KeyError

        self.kinex.get_score("GRNS*PVQA")
        self.kinex.get_score("GRNS(ph)PVQA", phospho_priming=False)

    def test_scoring(self):
        result = self.kinex.get_score(
            "GRNSLs*PVQA", phospho_priming=False, favorability=False, method="all"
        )
        self.assertEqual(len(result.ranking[0]), len(get_pssm_ser_thr()))
        self.assertEqual(len(result.ranking[0].columns), 3)

    def test_enrichment(self):
        result = self.kinex.get_enrichment(get_test_input_sites())
        self.assertEqual(len(result.ser_thr.enrichment_table), len(get_pssm_ser_thr()))
        self.assertEqual(len(result.ser_thr.enrichment_table.columns), 19)
        self.assertEqual(len(result.tyr.enrichment_table), len(get_pssm_tyr()))
        self.assertEqual(len(result.failed_sites), 1)


if __name__ == "__main__":
    unittest.main()
