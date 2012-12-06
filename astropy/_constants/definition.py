class ConstantDefinition(float):
    """
    A float with meta-data related to physical constants
    """
    def __init__(self, value, uncertainty, name, origin, units):
        super(ConstantDefinition, self).__init__()
        self.value = value
        self.uncertainty = uncertainty
        self.name = name
        self.origin = origin
        self.units = units

    def __new__(cls, value, uncertainty, name, origin, units):
        return super(ConstantDefinition, cls).__new__(cls, value)
