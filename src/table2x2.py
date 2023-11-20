import numpy as np
from scipy.stats import fisher_exact


class Table2x2:
    """
    A 2x2 contingency table

    Attributes
    ----------
        table: numpy.ndarray
            2x2 contingency table
        shift_zeros: bool
            Represents whether to shift zero values

    Methods
    -------
        odds_ratio(self) -> float
            Returns odds ratio value
        p_val(self, mode) -> float
            Returns p-value using Fisher exact test
    """

    def __init__(self, table, shift_zeros):
        self.table = table
        self.shift_zeros = shift_zeros
        if self.shift_zeros:
            if 0 in self.table:
                self.table = np.add(self.table, 0.5)

    def __repr__(self):
        return f'{self.table}'

    @property
    def table(self):
        return self._table

    @table.setter
    def table(self, table):
        n_rows = len(table)
        if not n_rows == 2:
            raise ValueError("Wrong table format")
        for row in table:
            if not len(row) == 2:
                raise ValueError("Wrong table format")
            for cell in row:
                if not isinstance(cell, (int, float, np.int64, np.float64)):
                    raise ValueError("Wrong table format")
        self._table = np.array(table)

    @property
    def shift_zeros(self):
        return self._shift_zeros

    @shift_zeros.setter
    def shift_zeros(self, value):
        if isinstance(value, bool):
            self._shift_zeros = value
        else:
            raise ValueError("Wrong shift_zeros format")

    def odds_ratio(self) -> float:

        return (self.table[0][0] * self.table[1][1] / (self.table[0][1] * self.table[1][0]))

    def p_val(self, mode) -> float:
        return fisher_exact(self.table, alternative=mode).pvalue
