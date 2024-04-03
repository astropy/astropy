# Licensed under a 3-clause BSD style license - see LICENSE.rst


import numpy as np

from astropy.utils.misc import isiterable

__all__ = ["FlagCollection"]


class FlagCollection(dict):
    """
    The purpose of this class is to provide a dictionary for
    containing arrays of flags for the `NDData` class. Flags should be
    stored in Numpy arrays that have the same dimensions as the parent
    data, so the `FlagCollection` class adds shape checking to a
    dictionary.

    The `FlagCollection` should be initialized like a
    dict, but with the addition of a ``shape=``
    keyword argument used to pass the NDData shape.
    """

    def __init__(self, *args, **kwargs):
        if "shape" in kwargs:
            self.shape = kwargs.pop("shape")
            if not isiterable(self.shape):
                raise ValueError("FlagCollection shape should be an iterable object")
        else:
            raise Exception(
                "FlagCollection should be initialized with the shape of the data"
            )

        super().__init__(self, *args, **kwargs)

    def __setitem__(self, item, value, **kwargs):
        if isinstance(value, np.ndarray):
            if value.shape == self.shape:
                super().__setitem__(item, value)
            else:
                raise ValueError(
                    f"flags array shape {value.shape} does not match data shape "
                    f"{self.shape}"
                )
        else:
            raise TypeError("flags should be given as a Numpy array")
