# Licensed under a 3-clause BSD style license - see LICENSE.rst


__all__ = ['Constant']


class Constant(float):
    """A physical or astronomical constant.

    Objects of this class behave just like `float` values when used in
    expressions, but they contain additional information

    Attributes
    ----------
    name : str
        The name of this constant.
    error : float
        The uncertainty in the value of this constant.
    origin : str
        The source used for the value of this constant.
    units : `astropy.units.UnitBase` instance or str
        A string describing the units
    """

    def __init__(self, value, error, name, origin, units):
        super(Constant, self).__init__()
        self.error = error
        self.name = name
        self.origin = origin
        self.units = units

    def __new__(cls, value, error, name, origin, units):
        return super(Constant, cls).__new__(cls, value)

    def __repr__(self):
        s = "<{0}: ".format(self.__class__.__name__)
        s += "name='{0}' ".format(self.name)
        s += "value={0} ".format(super(Constant, self).__repr__())
        s += "error={0} ".format(self.error)
        s += "units='{0}' ".format(self.units)
        s += "origin='{0}' ".format(self.origin)
        s += ">"
        return s

    def __str__(self):
        s = "  Name   = {0}\n".format(self.name)
        s += "  Value  = {0}\n".format(super(Constant, self).__repr__())
        s += "  Error  = {0}\n".format(self.error)
        s += "  Units = {0}\n".format(self.units)
        s += "  Origin = {0}".format(self.origin)
        return s

    @property
    def units(self):
        # There is a necessary cyclical dependency between
        # astropy.units and astropy.constants, since some units are
        # defined in terms of constants.  Therefore, we have to create
        # the unit objects in the getter, not the setter.
        from astropy import units as u

        if self._unitobj is None:
            self._unitobj = u.Unit(self._units)
        return self._unitobj

    @units.setter
    def units(self, unit):
        self._units = unit
        self._unitobj = None
