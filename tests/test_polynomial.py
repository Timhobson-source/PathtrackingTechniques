import unittest

from pathtracking.polynomial import PolynomialSystem


class TestPolynomialSystem(unittest.TestCase):
    def test_hello_is_none(self):
        polysys = PolynomialSystem()
        self.assertIsNone(polysys._hello())
