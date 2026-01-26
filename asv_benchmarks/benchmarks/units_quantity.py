import astropy.units as u


class QuantityScalarArithmetic:
    """
    Benchmarks for scalar Quantity arithmetic.

    These benchmarks specifically target builtin scalar multiplication
    and division, which are safe fast-path cases.
    """

    def setup(self):
        self.q = 1.234 * u.m

    def time_mul_float(self):
        self.q * 2.0

    def time_div_float(self):
        self.q / 2.0

