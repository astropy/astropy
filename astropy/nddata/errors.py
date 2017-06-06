#module for errors in ndata
class NDError(object):
    #super class for errors
    def __init__(self, error)
        self.error = error
    
    def interpolate(self, old_wcs, new_wcs):
        raise NotImplemented('This class is only for structure and not meant for usage')

    def convolve(self, kernel):
        raise NotImplemented('This class is only for structure and not meant for usage')

    def error_add(self,
                  self_data,
                  self_mask,
                  self_flag,
                  operand,
                  operand_data,
                  operand_mask,
                  operand_flag,
                  result_data):
        raise NotImplemented('This class is only for structure and not meant for usage')

    def error_sub(self,
                  self_data,
                  self_mask,
                  self_flag,
                  operand,
                  operand_data,
                  operand_mask,
                  operand_flag,
                  result_data):
        raise NotImplemented('This class is only for structure and not meant for usage')

    def error_mul(self,
                  self_data,
                  self_mask,
                  self_flag,
                  operand,
                  operand_data,
                  operand_mask,
                  operand_flag,
                  result_data):
        raise NotImplemented('This class is only for structure and not meant for usage')

    def error_div(self,
                  self_data,
                  self_mask,
                  self_flag,
                  operand,
                  operand_data,
                  operand_mask,
                  operand_flag,
                  result_data):
        raise NotImplemented('This class is only for structure and not meant for usage')


class SDError(NDError):
    
    def __init__(self, error):
        self.error = error

    
    def interpolate(self, old_wcs, new_wcs):
        raise NotImplemented()

    def convolve(self, kernel):
        raise NotImplemented()
    
    def _check_operand(self, operand):
        if not isinstance(operand, (SDError, None)):
            raise ValueError('unsupported operand type(s) for +: %s and %s' %\
                             (type(self), type(operand)))
    
    def error_add(self,
                  self_data,
                  self_mask,
                  self_flag,
                  operand,
                  operand_data,
                  operand_mask,
                  operand_flag,
                  result_data):#should we also do result_mask, result_flag - chicken / egg problem what comes first??
        self._check_operand(operand)

        if operand is None:
            return self
        
        return np.sqrt(self.error**2 + operand**2)
        

