import unittest
import pandas as pd

from src.kinex import Kinex

pssm_table = pd.read_csv('tests/data/test_pssm_table.csv', index_col=0)
scoring_matrix = pd.read_csv('tests/data/test_scoring_matrix.csv', index_col=0)

class TestKinex(unittest.TestCase):

    def test_scoring(self):
        kinex = Kinex(pssm=pssm_table, scoring_matrix=scoring_matrix)
        result = kinex.get_score('PSVEPPLs*QETFSDL')

        self.assertEqual(len(result.ranking), 2)
        self.assertEqual(len(result.ranking.columns), 3)
    
    def test_enrichment(self):
        kinex = Kinex(pssm=pssm_table, scoring_matrix=scoring_matrix)
        input_sites = pd.read_csv('test/data/test_input_sites.csv')
        result = kinex.get_enrichment(input_sites=input_sites)

        self.assertEqual(len(result.enrichment_table), 2)
        self.assertEqual(len(result.enrichment_table.columns), 19)