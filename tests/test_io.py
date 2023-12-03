import unittest

from kinex.functions import check_sequence, get_sequence_format 

class TestIO(unittest.TestCase):

    def test_get_sequence_format(self):
        self.assertEqual(get_sequence_format('PSVEPPLs*QETFSDL'), '*')
        self.assertEqual(get_sequence_format('PSVEXPLs*QXTF___'), '*')
        self.assertEqual(get_sequence_format('PSVEPPLsQETFSDL'), 'central')
        self.assertEqual(get_sequence_format('PSVEPPLs(ph)QETFSDL'), '(ph)')
        self.assertEqual(get_sequence_format('LQVKIPSKEEEsAD'), 'unsupported')

    def test_check_sequence(self):
        self.assertEqual(check_sequence('PSVEPPLs*QETFSDL', sequence_format='*'), True)
        self.assertEqual(check_sequence('PSVEXPLs*QXTF___', sequence_format='*'), True)
        self.assertEqual(check_sequence('PSVEPPLsQETFSDL', sequence_format='central'), True)

if __name__ == '__main__':
    unittest.main()