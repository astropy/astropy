from collections import namedtuple
from contextlib import contextmanager
from copy import copy
from typing import Dict

import numpy as np

from astropy.utils import unbroadcast
from astropy.utils.data_info import MixinInfo
from astropy.utils.shapes import ShapedLikeNDArray

__all__ = ["StokesCoord", "custom_stokes_symbol_mapping", "StokesSymbol"]


StokesSymbol = namedtuple("StokesSymbol", ["symbol", "description"], defaults=[""])

# This is table 29 in the FITS 4.0 paper
FITS_STOKES_VALUE_SYMBOL_MAP = {
    1: StokesSymbol("I", "Standard Stokes unpolarized"),
    2: StokesSymbol("Q", "Standard Stokes linear"),
    3: StokesSymbol("U", "Standard Stokes linear"),
    4: StokesSymbol("V", "Standard Stokes circular"),
    -1: StokesSymbol("RR", "Right-right circular"),
    -2: StokesSymbol("LL", "Left-left circular"),
    -3: StokesSymbol("RL", "Right-left cross-circular"),
    -4: StokesSymbol("LR", "Left-right cross-circular"),
    -5: StokesSymbol("XX", "X parallel linear"),
    -6: StokesSymbol("YY", "Y parallel linear"),
    -7: StokesSymbol("XY", "XY cross linear"),
    -8: StokesSymbol("YX", "YX cross linear"),
}

STOKES_VALUE_SYMBOL_MAP = copy(FITS_STOKES_VALUE_SYMBOL_MAP)


@contextmanager
def custom_stokes_symbol_mapping(
    mapping: Dict[int, StokesSymbol], replace: bool = False
):
    """
    Add a custom set of mappings from values to Stokes symbols.

    Parameters
    ----------
    mappings
        A list of dictionaries with custom mappings between values (integers)
        and `.StokesSymbol` classes.
    replace
        Replace all mappings with this one.
    """
    global STOKES_VALUE_SYMBOL_MAP

    original_mapping = copy(STOKES_VALUE_SYMBOL_MAP)
    if not replace:
        STOKES_VALUE_SYMBOL_MAP = {**original_mapping, **mapping}
    else:
        STOKES_VALUE_SYMBOL_MAP = mapping

    yield

    STOKES_VALUE_SYMBOL_MAP = original_mapping


class StokesCoordInfo(MixinInfo):
    _represent_as_dict_attrs = {"_data"}
    _represent_as_dict_primary_data = "_data"

    @property
    def unit(self):
        return None

    @property
    def dtype(self):
        return self._parent._data.dtype


class StokesCoord(ShapedLikeNDArray):
    """
    A representation of stokes coordinates with helpers for converting to profile names.

    Parameters
    ----------
    stokes_values : array-like
        The numeric values representing stokes coordinates.
    """

    info = StokesCoordInfo()

    def __init__(self, stokes, copy=False):
        stokes = np.asanyarray(stokes)
        if stokes.dtype.kind == "U":
            self._data = self._from_symbols(stokes)
        else:
            self._data = stokes

    @property
    def shape(self):
        return self._data.shape

    @property
    def value(self):
        return self._data

    # def __array__(self):
    #     return self._data

    def _apply(self, method, *args, **kwargs):
        cls = type(self)

        if callable(method):
            return cls(method(self._data, *args, **kwargs))
        elif method == "__getitem__":
            return cls(self._data.__getitem__(*args))
        else:
            return getattr(self, method)(*args, **kwargs)

    @staticmethod
    def _from_symbols(symbols):
        """
        Construct a StokesCoord from strings representing the stokes symbols.
        """
        values_array = np.full_like(symbols, np.nan, dtype=float, subok=False)
        for stokes_value, symbol in STOKES_VALUE_SYMBOL_MAP.items():
            values_array[symbols == symbol.symbol] = stokes_value

        if (nan_values := np.isnan(values_array)).any():
            raise ValueError(
                f"Unknown stokes symbols present in the input array: {np.unique(symbols[nan_values])}"
            )

        return values_array

    @property
    def _stokes_values(self):
        """
        A representation of the coordinate as integers.
        """
        # Note we unbroadcast and re-broadcast here to prevent the new array
        # using more memory than the old one.
        return np.broadcast_to(np.round(unbroadcast(self._data)), self.shape)

    @property
    def symbol(self):
        """
        The coordinate represented as strings
        """
        known_symbols = tuple(
            ["?"] + [s.symbol for s in STOKES_VALUE_SYMBOL_MAP.values()]
        )
        max_len = np.max([len(s) for s in known_symbols])

        # Note we unbroadcast and re-broadcast here to prevent the new array
        # using more memory than the old one.
        symbolarr = np.full(unbroadcast(self._data).shape, "?", dtype=f"<U{max_len}")

        for value, symbol in STOKES_VALUE_SYMBOL_MAP.items():
            symbolarr[unbroadcast(self._stokes_values) == value] = symbol.symbol

        return np.broadcast_to(symbolarr, self.shape)

    def __eq__(self, other):
        try:
            other = self.__class__(other)
        except Exception:
            return NotImplemented

        return self._data == other._data

    def __repr__(self):
        arrstr = np.array2string(self.symbol, separator=", ", prefix="  ")
        return f"{type(self).__name__}({arrstr})"
