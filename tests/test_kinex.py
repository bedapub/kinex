import unittest
import pandas as pd

from kinex.kinex import Kinex

from kinex.tests.data import get_test_input_sites, get_test_pssm, get_test_scoring_matrix

class TestKinex(unittest.TestCase):

    def test_scoring(self):
        kinex = Kinex(scoring_matrix=get_test_scoring_matrix(), pssm=get_test_pssm())
        result = kinex.get_score('PSVEPPLs*QETFSDL')

        self.assertEqual(len(result.ranking), 2)
        self.assertEqual(len(result.ranking.columns), 3)
    
    def test_enrichment(self):
        kinex = Kinex(scoring_matrix=get_test_scoring_matrix(), pssm=get_test_pssm())
        input_sites = get_test_input_sites()
        result = kinex.get_enrichment(input_sites=input_sites)

        self.assertEqual(len(result.enrichment_table), 2)
        self.assertEqual(len(result.enrichment_table.columns), 19)