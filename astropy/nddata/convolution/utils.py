

class NormalizationError(Exception):
    """
    Called when normalization of kernels is wrong.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


class KernelSizeError(Exception):
    """
    Called when size of kernels is wrong.
    """
    def __init__(self, message):
        self._message = message

    def __str__(self):
        return self._message


def next_larger_odd(number):
    """
    Round int to next larger odd int.
    """
    number = int(number + 0.5)
    if number % 2 == 0:
        return number + 1
    else:
        return number
