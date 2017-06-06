# Collections of different masks for nddata object

import numpy as np

import pyfits

class NDMask(object):
    #Mask superclass
    pass


class BoolMask(NDMask):
    
    
    def __init__(self, mask):
        self.mask = mask
        
    def _check_operand(self, operand):
        if not isinstance(operand, BoolMask):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
            
    def _arithmetic_operation(self, operand):
        self._check_operand(operand)
        return np.logical_and(self.mask, operand.mask)
        
    def mask_add(self, operand):
        return self._arithmetic_operation(operand)
        
    def mask_sub(self, operand):
        return self._arithmetic_operation(operand)
    
    def mask_mul(self, operand):
        return self._arithmetic_operation(operand)
        
    def mask_div(self, operand):
        return self._arithmetic_operation(operand)
