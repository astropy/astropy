# Licensed under a 3-clause BSD style license - see LICENSE.rst


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
    units : str
        A string describing the units

    .. warning ::
        The `units` attribute will be replaced when the Astropy units package
        is implemented. Don't count on it being a simple string in the future.

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
