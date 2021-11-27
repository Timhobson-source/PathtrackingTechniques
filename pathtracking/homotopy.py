from typing import Callable
import sympy as sp

DEFAULT_GAMMA = complex(0.44, 0.77)


class Homotopy:
    """
    Class representing a homotopy between two polynomials.

    Only works for single variable polynomials.

    TODO - add docs on homotopy.

    """

    def __init__(self, start_poly: Callable, end_poly: Callable, gamma=DEFAULT_GAMMA):
        self.start_poly = start_poly
        self.end_poly = end_poly
        self.gamma = gamma

    def path(self, x, t):
        return self.gamma * t * self.end_poly(x) + (1-t) * self.start_poly(x)

    @staticmethod
    def _setup_symbols():
        return sp.symbols('x t')

    def path_symbolised(self):
        args = self._setup_symbols()
        return self.path(*args)

    def __repr__(self):
        return str(self.path_symbolised)
