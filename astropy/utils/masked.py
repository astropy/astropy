# Licensed under a 3-clause BSD style license - see LICENSE.rst

__all__ = ['masked_arrays_equal']


def masked_arrays_equal(array1, array2):
    """
    Check if two masked arrays are equal.

    The rules are as follows:

    * If an entry is masked in both arrays, this returns `True`
    * If an entry is masked in only one array, this returns `False`
    * If an entry is masked in neither array, this returns `True` if the
      values are equal

    This differs from Numpy, which returns masked values in the boolean array
    if an entry is masked in either or both arrays, with a fill value of `True`.
    """

    # First use the standard Numpy comparison, using the default fill values
    # on the arrays.
    equal = array1.filled() == array2.filled()

    # We need to now adjust the results - we want entries masked in both arrays
    # to be considered equal (even if the fill value is different) and we want
    # entries masked in one array and not in the other to not match, even if
    # the fill values would make them match. We need to use equality instead of
    # boolean comparisons as these don't work with structured arrays.
    equal[(array1.mask == True) & (array2.mask == True)] = True
    equal[array1.mask != array2.mask] = False

    return equal
