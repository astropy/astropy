class ConstantDefinition(float):
    """
    A float with meta-data related to physical constants
    """
    def __init__(self, value, error, name, origin, units):
        super(ConstantDefinition, self).__init__()
        self.value = value
        self.error = error
        self.name = name
        self.origin = origin
        self.units = units

    def __new__(cls, value, error, name, origin, units):
        return super(ConstantDefinition, cls).__new__(cls, value)
