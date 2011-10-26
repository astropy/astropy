class Constant(float):

    def __init__(self, value, error, name, origin, units):
        self.error = error
        self.name = name
        self.origin = origin
        self.units = units
        return float.__init__(self, value)

    def __new__(self, value, error, name, origin, units):
        return float.__new__(self, value)

    def __repr__(self):
        s = "<Constant: "
        s += "name='{0}' ".format(self.name)
        s += "value={0} ".format(float.__repr__(self))
        s += "error={0} ".format(self.error)
        s += "units='{0}' ".format(self.units)
        s += "origin='{0}' ".format(self.origin)
        s += ">"
        return s

    def __str__(self):
        s = "  Name   = {0}\n".format(self.name)
        s += "  Value  = {0}\n".format(float.__repr__(self))
        s += "  Error  = {0}\n".format(self.error)
        s += "  Units = {0}\n".format(self.units)
        s += "  Origin = {0}".format(self.origin)
        return s
