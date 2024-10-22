from enum import Enum
import pandas as pd

allowed_characters = {
    "P", "G", "A", "C", "S", "T", "V", "I", "L", "M", "F", "Y", "W", "H",
    "K", "R", "Q", "N", "D", "E", "s", "t", "y", "X", "_", "*"
}


def is_central_sequence_valid(sequence_string: str) -> bool:
    if sequence_string[len(sequence_string) // 2] not in ("S", "T", "Y", "s", "t", "y"):
        return False
    for aminoacid in sequence_string:
        if aminoacid not in [x for x in allowed_characters if x != "*"]:
            return False
    return True


def is_separator_sequence_valid(sequence_string: str, separator: str) -> bool:
    if separator * 2 in sequence_string:
        return False

    valid_patterns = get_valid_patterns(separator)
    if not any(pattern in sequence_string for pattern in valid_patterns):
        return False

    for aminoacid in sequence_string:
        if aminoacid not in allowed_characters:
            return False
    return True


def get_valid_patterns(separator):
    return [
        f"S{separator}",
        f"T{separator}",
        f"Y{separator}",
        f"s{separator}",
        f"t{separator}",
        f"y{separator}"
    ]


class SequenceSeparator(Enum):
    ASTERISK = "*"
    PH = "(ph)"


class SequenceType(Enum):
    SER_THR = ["S", "T"]
    TYR = ["Y"]


class Sequence:
    def __init__(self, sequence_string: str, sequence_type: SequenceType) -> None:
        self.sequence_string = sequence_string
        self.sequence_type = sequence_type

    def validate_sequence(self):
        pass

    def preprocess_sequence(self):
        pass

    def get_split_sequence(self) -> list:
        pass

    def get_columns_list(self):
        skip_characters = {"_", "X"}
        columns = []

        columns_range = (-5, 4)
        if self.sequence_type == SequenceType.SER_THR:
            columns_range = (-5, 4)
        elif self.sequence_type == SequenceType.TYR:
            columns_range = (-5, 5)

        part_id = 0
        parts = self.get_split_sequence()
        for part in parts:
            if part_id == 0:
                part = part[::-1]
            for position, aminoacid in enumerate(part):
                if aminoacid in skip_characters:
                    continue
                pos = (position + 1) * (-1) if part_id == 0 else position + 1
                if columns_range[0] <= pos <= columns_range[1]:
                    columns.append(f"{pos}{aminoacid}")
            part_id += 1
        return columns

    def get_sequence_scores(self, pssm_table: pd.DataFrame, favorability: bool = False) -> list:
        pass


class SeparatorSequence(Sequence):
    def __init__(self, sequence_string: str, separator: SequenceSeparator, sequence_type: SequenceType) -> None:
        super().__init__(sequence_string, sequence_type)
        self.separator = separator
        self.sequences = []

    def preprocess_sequence(self):
        sequences = self.sequence_string.split(self.separator.value)
        for index, item in enumerate(sequences):
            if index != len(sequences) - 1:
                sequences[index] = item[:-1] + item[-1].lower()
        for index in range(len(sequences) - 1):
            seq = f"{sequences[:index + 1]}{self.separator.value}{sequences[index + 1:]}"
            if not is_separator_sequence_valid(seq, self.separator.value):
                continue
            sequences.append(seq)
        self.sequences = sequences

    def validate_sequence(self) -> None:
        if len(self.sequences) == 0:
            raise ValueError("Invalid sequence")

    def get_split_sequence(self) -> list:
        parts = self.sequence_string.split(self.separator.value)
        parts[0] = parts[0][:-1]
        return parts

    def get_sequence_scores(self, pssm_table: pd.DataFrame, favorability: bool = False) -> list:
        score_results = []
        for sequence in self.sequences:
            columns_list = self.get_columns_list()
            if favorability:
                seq_upper = sequence.upper()
                if f"S{self.separator.value}" in seq_upper:
                    columns_list.append("0S")
                elif f"T{self.separator.value}" in seq_upper:
                    columns_list.append("0T")
            score_results.append(get_score(columns_list, pssm_table))
        return score_results


class CentralSequence(Sequence):
    def __init__(self, sequence_string: str, sequence_type: SequenceType) -> None:
        super().__init__(sequence_string, sequence_type)

    def preprocess_sequence(self):
        pass

    def validate_sequence(self):
        if not is_central_sequence_valid(self.sequence_string):
            raise ValueError("Invalid sequence")

    def get_split_sequence(self) -> list:
        sequence_length = len(self.sequence_string)
        return [
            self.sequence_string[:sequence_length // 2],
            self.sequence_string[sequence_length // 2 + 1:]
        ]

    def get_sequence_scores(self, pssm_table: pd.DataFrame, favorability: bool = False) -> list:
        columns = self.get_columns_list()
        if favorability:
            seq_upper = self.sequence_string.upper()
            if seq_upper[len(seq_upper) // 2] == "S":
                columns.append("0S")
            elif seq_upper[len(seq_upper) // 2] == "T":
                columns.append("0T")
        return [get_score(columns, pssm_table)]


def get_sequence_type(aminoacid: str) -> SequenceType:
    aminoacid = aminoacid.upper()
    for sequence_type in SequenceType:
        if aminoacid in sequence_type.value:
            return sequence_type
    raise ValueError(f"Unsupported sequence type. Supported types: SER, THR, TYR")


def get_sequence_object(sequence: str) -> Sequence:
    for separator in SequenceSeparator:
        position = sequence.find(separator.value)
        if position > 0:
            sequence_type = get_sequence_type(sequence[position - 1])
            return SeparatorSequence(sequence, separator, sequence_type)
    if len(sequence) % 2 == 1:
        sequence_type = get_sequence_type(sequence[int(len(sequence) / 2)])
        return CentralSequence(sequence, sequence_type)
    raise ValueError(
        f"Unsupported sequence format. Supported formats: *, (ph) and central")


def get_score(columns: list, pssm: pd.DataFrame) -> pd.DataFrame:
    pssm = pssm.reset_index()
    columns.append("kinase")
    score_df = pssm[columns]
    score_df.insert(0, "score", score_df.prod(axis=1, numeric_only=True))
    score_df = score_df[["kinase", "score"]]
    score_df = score_df.set_index("kinase")
    return score_df