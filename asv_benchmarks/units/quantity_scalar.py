import astropy.units as u
from astropy.units import Quantity


class QuantityScalarBench:
    """
    ASV benchmarks for scalar Quantity arithmetic.

    These benchmarks track performance of Quantity * scalar and
    Quantity / scalar, which are safe fast-path cases.
    """

    def setup(self):
        self.q = Quantity(1.234, u.m)
        self.scalar = 2.5

    def time_mul_scalar(self):
        self.q * self.scalar

    def time_truediv_scalar(self):
        self.q / self.scalar

