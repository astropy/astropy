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
    uncertainty : float
        The uncertainty in the value of this constant.
    origin : str
        The source used for the value of this constant.
    units : `astropy.units.UnitBase` instance
        The units of this constant. Can be set either as a string or
        `astropy.units.UnitBase`.
    """

    def __init__(self, value, unit, uncertainty, name, origin):
        super(Constant, self).__init__(value, unit)
        self.uncertainty = uncertainty
        self.name = name
        self.origin = origin

    def __new__(cls, value, uncertainty, name, origin, units):
        return super(Constant, cls).__new__(cls, value)

    def __repr__(self):
        s = "<Constant: "
        s += "name='{0}' ".format(self.name)
        s += "value={0} ".format(self.value)
        s += "error={0} ".format(self.uncertainty)
        s += "units='{0}' ".format(self.unit)
        s += "origin='{0}'".format(self.origin)
        s += ">"
        return s

    def __str__(self):
        s = "  Name   = {0}\n".format(self.name)
        s += "  Value  = {0}\n".format(self.value)
        s += "  Error  = {0}\n".format(self.uncertainty)
        s += "  Units = {0}\n".format(self.unit)
        s += "  Origin = {0}".format(self.origin)
        return s
