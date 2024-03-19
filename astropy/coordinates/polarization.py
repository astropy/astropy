from __future__ import annotations

from contextlib import contextmanager
from typing import NamedTuple

import numpy as np

from astropy.utils import unbroadcast
from astropy.utils.compat import COPY_IF_NEEDED
from astropy.utils.data_info import MixinInfo
from astropy.utils.shapes import ShapedLikeNDArray

__all__ = ["StokesCoord", "custom_stokes_symbol_mapping", "StokesSymbol"]


class StokesSymbol(NamedTuple):
    """Symbol for a Stokes coordinate."""

    symbol: str = ""
    description: str = ""


# This is table 29 in the FITS 4.0 paper
FITS_STOKES_VALUE_SYMBOL_MAP = {
    1: StokesSymbol("I", "Standard Stokes unpolarized"),
    2: StokesSymbol("Q", "Standard Stokes linear"),
    3: StokesSymbol("U", "Standard Stokes linear"),
    4: StokesSymbol("V", "Standard Stokes circular"),
    -1: StokesSymbol("RR", "Right-right circular: <RR*>"),
    -2: StokesSymbol("LL", "Left-left circular: <LL*>"),
    -3: StokesSymbol("RL", "Right-left cross-circular: Re(<RL*>))"),
    -4: StokesSymbol("LR", "Left-right cross-circular: Re(<LR*>)=Im(<RL*>)"),
    -5: StokesSymbol("XX", "X parallel linear: <XX*>"),
    -6: StokesSymbol("YY", "Y parallel linear: <YY*>"),
    -7: StokesSymbol("XY", "XY cross linear: Re(<XY*>)"),
    -8: StokesSymbol("YX", "YX cross linear: Im(<XY*>)"),
}

STOKES_VALUE_SYMBOL_MAP = FITS_STOKES_VALUE_SYMBOL_MAP.copy()
UNKNOWN_STOKES_VALUE = -99999


@contextmanager
def custom_stokes_symbol_mapping(
    mapping: dict[int, StokesSymbol], replace: bool = False
) -> None:
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

    original_mapping = STOKES_VALUE_SYMBOL_MAP.copy()
    if not replace:
        STOKES_VALUE_SYMBOL_MAP = {**original_mapping, **mapping}
    else:
        STOKES_VALUE_SYMBOL_MAP = mapping

    yield

    STOKES_VALUE_SYMBOL_MAP = original_mapping


class StokesCoordInfo(MixinInfo):
    # The attributes containing actual information.
    _represent_as_dict_attrs = {"value"}
    # Since there is only one attribute, use a column with the name to represent it
    # (rather than as name.value)
    _represent_as_dict_primary_data = "value"
    # Attributes that should be presented as positional arguments to
    # the class initializer (which takes "stokes" as an argument, not "value").
    _construct_from_dict_args = ("value",)

    @property
    def unit(self):
        return None

    @property
    def dtype(self):
        return self._parent._data.dtype

    @staticmethod
    def default_format(val):
        return f"{val.symbol}"

    def new_like(self, cols, length, metadata_conflicts="warn", name=None):
        """
        Return a new StokesCoord instance which is consistent with the
        input ``cols`` and has ``length`` rows.

        This is intended for creating an empty column object whose elements can
        be set in-place for table operations like join or vstack.

        Parameters
        ----------
        cols : list
            List of input columns
        length : int
            Length of the output column object
        metadata_conflicts : str ('warn'|'error'|'silent')
            How to handle metadata conflicts
        name : str
            Output column name

        Returns
        -------
        col : `~astropy.coordinates.StokesCoord` (or subclass)
            Empty instance of this class consistent with ``cols``

        """
        # Get merged info attributes like shape, dtype, format, description, etc.
        attrs = self.merge_cols_attributes(
            cols, metadata_conflicts, name, ("meta", "format", "description")
        )

        # Make an empty StokesCoord.
        shape = (length,) + attrs.pop("shape")
        data = np.zeros(shape=shape, dtype=attrs.pop("dtype"))
        # Get arguments needed to reconstruct class
        out = self._construct_from_dict({"value": data})

        # Set remaining info attributes
        for attr, value in attrs.items():
            setattr(out.info, attr, value)

        return out

    def get_sortable_arrays(self):
        """
        Return a list of arrays which can be lexically sorted to represent
        the order of the parent column.

        For StokesCoord this is just the underlying values.

        Returns
        -------
        arrays : list of ndarray
        """
        return [self._parent._data]


class StokesCoord(ShapedLikeNDArray):
    """
    A representation of stokes coordinates with helpers for converting to profile names.

    Parameters
    ----------
    stokes : array-like
        The numeric values representing stokes coordinates.
    """

    info = StokesCoordInfo()

    def __init__(self, stokes, copy=False):
        if isinstance(stokes, type(self)):
            data = stokes._data.copy() if copy else stokes._data
            self.info = stokes.info
        else:
            data = np.asanyarray(stokes)

            if data.dtype.kind == "O":
                msg = "StokesCoord objects cannot be initialised with an object array."
                raise ValueError(msg)

            if data.dtype.kind == "U":
                data = self._from_symbols(data)
            else:
                data = data.copy() if copy and data is stokes else data
        self._data = data

    @property
    def shape(self):
        return self._data.shape

    @property
    def value(self):
        return self._data

    @property
    def dtype(self):
        return self._data.dtype

    def __array__(self, dtype=None, copy=COPY_IF_NEEDED):
        return self._data.astype(dtype, copy=copy)

    def _apply(self, method, *args, **kwargs):
        cls = type(self)

        if callable(method):
            new = cls(method(self._data, *args, **kwargs))
        else:
            new = cls(getattr(self._data, method)(*args, **kwargs))

        # Copy other 'info' attr only if it has actually been defined.
        # See PR #3898 for further explanation and justification, along
        # with Quantity.__array_finalize__
        if "info" in self.__dict__:
            new.info = self.info

        return new

    @staticmethod
    def _from_symbols(symbols):
        """
        Convert an array of symbols to an array of values
        """
        values_array = np.full_like(
            symbols, UNKNOWN_STOKES_VALUE, dtype=int, subok=False
        )
        for stokes_value, symbol in STOKES_VALUE_SYMBOL_MAP.items():
            values_array[symbols == symbol.symbol] = stokes_value

        if (unknown_values := np.equal(values_array, UNKNOWN_STOKES_VALUE)).any():
            raise ValueError(
                f"Unknown stokes symbols present in the input array: {np.unique(symbols[unknown_values])}"
            )

        return values_array

    @property
    def symbol(self):
        """The coordinate represented as strings."""
        known_symbols = tuple(
            ["?"] + [s.symbol for s in STOKES_VALUE_SYMBOL_MAP.values()]
        )
        max_len = np.max([len(s) for s in known_symbols])

        # Note we unbroadcast and re-broadcast here to prevent the new array
        # using more memory than the old one.
        unbroadcasted = np.round(unbroadcast(self.value))
        symbolarr = np.full(unbroadcasted.shape, "?", dtype=f"<U{max_len}")
        for value, symbol in STOKES_VALUE_SYMBOL_MAP.items():
            symbolarr[unbroadcasted == value] = symbol.symbol

        return np.broadcast_to(symbolarr, self.shape)

    def __setitem__(self, item, value):
        self._data[item] = type(self)(value)._data

    def __eq__(self, other):
        try:
            other = self.__class__(other)
        except Exception:
            return NotImplemented

        return self._data == other._data

    def __str__(self):
        arrstr = np.array2string(self.symbol, separator=", ", prefix="  ")
        return f"{type(self).__name__}({arrstr})"

    def __repr__(self):
        return self.__str__()
