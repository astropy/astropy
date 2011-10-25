class Constant(float):

    def __init__(self, value, error, name, origin, system):
        self.error = error
        self.name = name
        self.origin = origin
        self.system = system
        return float.__init__(self, value)

    def __new__(self, value, error, name, origin, system):
        return float.__new__(self, value)

    def __repr__(self):
        s = "<Constant: "
        s += "name='{0}' ".format(self.name)
        s += "value={0} ".format(float.__repr__(self))
        s += "error={0} ".format(self.error)
        s += "system='{0}' ".format(self.system)
        s += "origin='{0}' ".format(self.origin)
        s += ">"
        return s

    def __str__(self):
        s = "  Name   = {0}\n".format(self.name)
        s += "  Value  = {0}\n".format(float.__repr__(self))
        s += "  Error  = {0}\n".format(self.error)
        s += "  System = {0}\n".format(self.system)
        s += "  Origin = {0}".format(self.origin)
        return s
