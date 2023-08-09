import numpy as np
from scipy.stats import hypergeom

class Table2x2:
    """Summary of class here.

    Longer class information...
    Longer class information...

    Attributes:
        table: 
        shirt_zeros: 
    """

    def __init__(self, table, shift_zeros):
        """Initializes the instance based on the table and shirt_zeros arguments

        Args:
          table: 
        """
        self.table = table
        self.shift_zeros = shift_zeros
        if self.shift_zeros:
            if 0 in self.table:
                self.table = np.add(self.table, 0.5)

    def __repr__ (self):
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

    def odds_ratio(self):
        return (self.table[0][0] * self.table[1][1] / (self.table[0][1] * self.table[1][0]))
    
    def p_val(self, mode):
        x = self.table[0][0]
        M = self.table.sum()
        n = self.table[0].sum()
        N = self.table[:, 0].sum()
        start, end = hypergeom.support(M, n, N)
        if mode == 'greater':
            return hypergeom.sf(x - 1, M, n, N)
        if mode == 'less':
            return hypergeom.cdf(x, M, n, N)
        if mode == 'two-sided':
            pmf = hypergeom.pmf(x, M, n, N)
            res = 0
            for x_val in range(start, end + 1):
                if hypergeom.pmf(x_val, M, n, N) <= pmf:
                    res += hypergeom.pmf(x_val, M, n, N)
            return res
        else:
            raise ValueError("Wrong mode")