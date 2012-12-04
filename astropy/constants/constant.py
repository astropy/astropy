# Licensed under a 3-clause BSD style license - see LICENSE.rst

from ..units.quantity import Quantity

__all__ = ['Constant']


class Constant(Quantity):
    """A physical or astronomical constant.

    These objects are quantities that are meant to represent physical constants

    Attributes
    ----------
    name : str
        The name of this constant.
    error : float
        The uncertainty in the value of this constant.
    origin : str
        The source used for the value of this constant.
    units : `astropy.units.UnitBase` instance
        The units of this constant. Can be set either as a string or
        `astropy.units.UnitBase`.
    """

    def __init__(self, value, unit, error, name, origin):
        super(Constant, self).__init__(value, unit)
        self.error = error
        self.name = name
        self.origin = origin

    def __new__(cls, value, error, name, origin, units):
        return super(Constant, cls).__new__(cls, value)
