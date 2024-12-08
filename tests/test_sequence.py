import unittest
import pandas as pd
from kinex.sequence import (
    CentralSequence,
    SeparatorSequence,
    get_sequence_type,
    get_sequence_object,
    SequenceSeparator,
    SequenceType,
    is_central_sequence_valid,
    is_separator_sequence_valid,
    get_score,
)


class TestSequenceProcessing(unittest.TestCase):
    def setUp(self):
        # Sample PSSM DataFrame for testing
        self.pssm_data = {
            "kinase": ["AAK1", "ACVR2A"],
            "-5P": [1.2242, 0.7057],
            "-5G": [0.4165, 0.8178],
        }
        self.pssm = pd.DataFrame(self.pssm_data).set_index("kinase")

    def test_is_central_sequence_valid(self):
        # Test valid sequences
        valid_sequence = "APQSTPQPA"
        self.assertTrue(is_central_sequence_valid(valid_sequence))

        # Test invalid sequences
        invalid_sequence = "APQPPQPA"
        self.assertFalse(is_central_sequence_valid(invalid_sequence))

        # Test not a string
        invalid_sequence = {"APQSTPQPA"}
        self.assertFalse(is_central_sequence_valid(invalid_sequence))

    def test_is_separator_sequence_valid(self):
        # Test valid separator sequence
        valid_sequence = "APQ*S*TPQ*PA"
        separator = "*"
        self.assertTrue(is_separator_sequence_valid(valid_sequence, separator))

        # Test invalid separator sequence
        invalid_sequence = "APQ**STPQ*PA"
        self.assertFalse(is_separator_sequence_valid(invalid_sequence, separator))

        # Test invalid separator
        invalid_sequence = "APQ#STPQ#PA"
        separator = "#"
        self.assertFalse(is_separator_sequence_valid(invalid_sequence, separator))

        # Test consecutive separator
        invalid_sequence = "APQ**STPQ*PA"
        separator = "*"
        self.assertFalse(is_separator_sequence_valid(invalid_sequence, separator))

        # Test invalid patterns
        invalid_sequence = "APQSTPQPA"
        separator = "*"
        self.assertFalse(is_separator_sequence_valid(invalid_sequence, separator))

    def test_get_sequence_type(self):
        # Test valid sequence types
        self.assertEqual(get_sequence_type("S"), SequenceType.SER_THR)

        self.assertEqual(get_sequence_type("T"), SequenceType.SER_THR)

        self.assertEqual(get_sequence_type("Y"), SequenceType.TYR)

        # Test invalid sequence types
        invalid_sequence = "APQPQPA"
        with self.assertRaises(ValueError):
            get_sequence_type(invalid_sequence)

    def test_get_sequence_object(self):
        # Test CentralSequence object
        central_sequence_T = "APQTPQP"
        sequence = get_sequence_object(central_sequence_T)
        self.assertIsInstance(sequence, CentralSequence)

        central_sequence_S = "APQSPQP"
        sequence = get_sequence_object(central_sequence_S)
        self.assertIsInstance(sequence, CentralSequence)

        central_sequence_invalid_type = "APQZPQP"
        with self.assertRaises(ValueError):
            sequence = get_sequence_object(central_sequence_invalid_type)

        # Test SeparatorSequence object
        separator_sequence_asterisk = "APQST*PA"
        sequence = get_sequence_object(separator_sequence_asterisk)
        self.assertIsInstance(sequence, SeparatorSequence)

        separator_sequence_ph = "APQST(ph)PA"
        sequence = get_sequence_object(separator_sequence_asterisk)
        self.assertIsInstance(sequence, SeparatorSequence)

        separator_sequence_invalid_type = "APQZP*QP"
        with self.assertRaises(ValueError):
            sequence = get_sequence_object(separator_sequence_invalid_type)

        separator_sequence_invalid_type_ph = "APQZP(ph)QP"
        with self.assertRaises(ValueError):
            sequence = get_sequence_object(separator_sequence_invalid_type_ph)

    def test_get_score(self):
        # Sample PSSM DataFrame
        self.pssm_data = {
            "kinase": ["AAK1", "ACVR2A"],
            "-5P": [1.2242, 0.7057],
            "-5G": [0.4165, 0.8178],
        }
        self.pssm = pd.DataFrame(self.pssm_data).set_index("kinase")

        # Test valid columns
        columns = ["-5P", "-5G"]
        expected_data = {
            "kinase": ["AAK1", "ACVR2A"],
            "score": [1.2242 * 0.4165, 0.7057 * 0.8178],
        }
        expected_df = pd.DataFrame(expected_data).set_index("kinase")
        expected_df["score"] = expected_df["score"].astype(float)

        result_df = get_score(columns, self.pssm)

        pd.testing.assert_frame_equal(result_df, expected_df)

        # Test empty columns
        columns = []
        expected_data = {"kinase": ["AAK1", "ACVR2A"], "score": [1, 1]}
        expected_df = pd.DataFrame(expected_data).set_index("kinase")
        expected_df["score"] = expected_df["score"].astype(float)

        result_df = get_score(columns, self.pssm)

        pd.testing.assert_frame_equal(result_df, expected_df)


class TestCentralSequence(unittest.TestCase):
    def setUp(self):
        self.valid_sequence = "APQATPQPA"
        self.invalid_sequence = "APQPPQPA"
        self.sequence_type = SequenceType.SER_THR

    def test_validate_sequence_valid(self):
        central_seq = CentralSequence(self.valid_sequence, self.sequence_type)
        central_seq.validate_sequence()

    def test_validate_sequence_invalid(self):
        central_seq = CentralSequence(self.invalid_sequence, self.sequence_type)
        with self.assertRaises(ValueError):
            central_seq.validate_sequence()

    def test_get_split_sequence(self):
        central_seq = CentralSequence(self.valid_sequence, self.sequence_type)
        split_sequence = central_seq.get_split_sequence()
        self.assertEqual(split_sequence, ["APQA", "PQPA"])

    def test_get_sequence_scores(self):
        pssm_data = {
            "kinase": ["AAK1", "ACVR2A"],
            "-1A": [0.5, 0.6],
            "-2Q": [0.7, 0.8],
            "-3P": [0.9, 1.0],
            "-4A": [1.1, 1.2],
            "1P": [1.3, 1.4],
            "2Q": [1.5, 1.6],
            "3P": [1.7, 1.8],
            "4A": [1.9, 2.0],
            "0T": [1.0, 1.0],
        }
        pssm = pd.DataFrame(pssm_data).set_index("kinase")

        # Test favorability
        central_seq = CentralSequence(self.valid_sequence, self.sequence_type)
        columns = ["-1A", "-2Q", "-3P", "-4A", "1P", "2Q", "3P", "4A"]
        expected_data = {
            "kinase": ["AAK1", "ACVR2A"],
            "score": [
                0.5 * 0.7 * 0.9 * 1.1 * 1.3 * 1.5 * 1.7 * 1.9 * 1.0,
                0.6 * 0.8 * 1.0 * 1.2 * 1.4 * 1.6 * 1.8 * 2.0 * 1.0,
            ],
        }
        expected_df = pd.DataFrame(expected_data).set_index("kinase")
        expected_df["score"] = expected_df["score"].astype(float)

        result_df = central_seq.get_sequence_scores(pssm, favorability=True)

        pd.testing.assert_frame_equal(result_df[0], expected_df)

    def test_get_column_list(self):
        central_seq = CentralSequence(self.valid_sequence, self.sequence_type)
        columnList = central_seq.get_columns_list()
        self.assertEqual(
            columnList, ["-1A", "-2Q", "-3P", "-4A", "1P", "2Q", "3P", "4A"]
        )


class TestSeparatorSequence(unittest.TestCase):
    def setUp(self):
        self.valid_sequence = "SGLAAS*AAQQQ"
        self.separator = SequenceSeparator.ASTERISK
        self.sequence_type = SequenceType.SER_THR
        self.separator_sequence = SeparatorSequence(
            self.valid_sequence, self.separator, self.sequence_type
        )

    def test_preprocess_sequence(self):
        self.separator_sequence.preprocess_sequence()
        self.assertEqual(self.separator_sequence.sequences, ["SGLAAs", "AAQQQ"])

    def test_get_split_sequence(self):
        split_sequence = self.separator_sequence.get_split_sequence()
        self.assertEqual(split_sequence, ["SGLAA", "AAQQQ"])

    def test_get_columnList(self):
        columnList = self.separator_sequence.get_columns_list()
        self.assertEqual(
            columnList, ["-1A", "-2A", "-3L", "-4G", "-5S", "1A", "2A", "3Q", "4Q"]
        )


if __name__ == "__main__":
    unittest.main()
