#module for errors in ndata
class NDError(np.ndarray):
    #super class for errors
    pass
    
    
class SDError(NDError):
    
    def __new__(cls, input_array):
        obj = np.asarray(input_array).view(cls)
        return obj
    
    def interpolate(self, old_lookup_table, new_lookup_table):
        return np.interp(new_lookup_table, old_lookup_table)
    
    def _check_operand(self, operand):
        if not isinstance(operand, (SDError, None)):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
    
    def error_add(self, self_data, operand, operand_data, result_data):
        self._check_operand(operand)
        if operand is None:
            return self
        
        return np.sqrt(self**2 + operand**2)
        
    def error_sub(self, self_data, operand, operand_data, result_data):
        self._check_operand(operand)
        if operand is None:
            return self
        
        return np.sqrt(self**2 + operand**2)

    def error_mul(self, self_data, operand, operand_data, result_data):
        self._check_operand(operand)
        if operand is None:
            return self
        
        return np.sqrt((self / self_data)**2 + (operand / operand_data)**2) * result_data
        
    def error_div(self, self_data, operand, operand_data, result_data):
        self._check_operand(operand)
        if operand is None:
            return self
        
        return np.sqrt((self / self_data)**2 + (operand / operand_data)**2) * result_data

