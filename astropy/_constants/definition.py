# Licensed under a 3-clause BSD style license - see LICENSE.rst
class ConstantDefinition(float):
    """
    A float with meta-data related to physical constants
    """
    def __init__(self, value, uncertainty, name, reference, units,
                 system=None):
        super(ConstantDefinition, self).__init__()
        self.value = value
        self.uncertainty = uncertainty
        self.name = name
        self.reference = reference
        self.units = units
        self.system = system

    def __new__(cls, value, uncertainty, name, reference, units, system=None):
        return super(ConstantDefinition, cls).__new__(cls, value)
