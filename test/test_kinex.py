import unittest
import pandas as pd

from src.kinex import Kinex

pssm_table = pd.read_csv('data/test_pssm_table.csv', index_col=0)
scoring_matrix = pd.read_scv('data/test_scoring_matrix.csv', index_col=0)

class TestKinex(unittest.TestCase):

    def test_scoring(self):
        kinex = Kinex(pssm=pssm_table, scoring_matrix=scoring_matrix)
        result = kinex.score('PSVEPPLs*QETFSDL')

        self.assertEqual(len(result.scores), 2)
        self.assertEqual(len(result.scores.index), 2)