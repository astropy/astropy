# Collections of different masks for nddata object

import numpy as np

class BoolMask(np.ndarray):

    def __new__(cls, input_array):
        obj = np.asarray(input_array, dtype='bool').view(cls)
        return obj

    def arithmetic_operation(self, operand):
        if not isinstance(operand, BoolMask):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
        return np.logical_or(self, operand)
        
    __add__ = arithmetic_operation
    __sub__ = arithmetic_operation
    __mul__ = arithmetic_operation
    __div__ = arithmetic_operation
    
    __radd__ = arithmetic_operation
    __rsub__ = arithmetic_operation
    __rmul__ = arithmetic_operation
    __rdiv__ = arithmetic_operation
    