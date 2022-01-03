# Licensed under a 3-clause BSD style license - see PYFITS.rst

import copy
import operator
import warnings
import weakref

from contextlib import suppress
from functools import reduce

import numpy as np

from numpy import char as chararray

from .column import (ASCIITNULL, FITS2NUMPY, ASCII2NUMPY, ASCII2STR, ColDefs,
                     _AsciiColDefs, _FormatX, _FormatP, _VLF, _get_index,
                     _wrapx, _unwrapx, _makep, Delayed)
from .util import decode_ascii, encode_ascii, _rstrip_inplace
from astropy.utils import lazyproperty


class FITS_record:
    """
    FITS record class.

    `FITS_record` is used to access records of the `FITS_rec` object.
    This will allow us to deal with scaled columns.  It also handles
    conversion/scaling of columns in ASCII tables.  The `FITS_record`
    class expects a `FITS_rec` object as input.
    """

    def __init__(self, input, row=0, start=None, end=None, step=None,
                 base=None, **kwargs):
        """
        Parameters
        ----------
        input : array
            The array to wrap.
        row : int, optional
            The starting logical row of the array.
        start : int, optional
            The starting column in the row associated with this object.
            Used for subsetting the columns of the `FITS_rec` object.
        end : int, optional
            The ending column in the row associated with this object.
            Used for subsetting the columns of the `FITS_rec` object.
        """

        self.array = input
        self.row = row
        if base:
            width = len(base)
        else:
            width = self.array._nfields

        s = slice(start, end, step).indices(width)
        self.start, self.end, self.step = s
        self.base = base

    def __getitem__(self, key):
        if isinstance(key, str):
            indx = _get_index(self.array.names, key)

            if indx < self.start or indx > self.end - 1:
                raise KeyError(f"Key '{key}' does not exist.")
        elif isinstance(key, slice):
            return type(self)(self.array, self.row, key.start, key.stop,
                              key.step, self)
        else:
            indx = self._get_index(key)

            if indx > self.array._nfields - 1:
                raise IndexError('Index out of bounds')

        return self.array.field(indx)[self.row]

    def __setitem__(self, key, value):
        if isinstance(key, str):
            indx = _get_index(self.array.names, key)

            if indx < self.start or indx > self.end - 1:
                raise KeyError(f"Key '{key}' does not exist.")
        elif isinstance(key, slice):
            for indx in range(slice.start, slice.stop, slice.step):
                indx = self._get_indx(indx)
                self.array.field(indx)[self.row] = value
        else:
            indx = self._get_index(key)
            if indx > self.array._nfields - 1:
                raise IndexError('Index out of bounds')

        self.array.field(indx)[self.row] = value

    def __len__(self):
        return len(range(self.start, self.end, self.step))

    def __repr__(self):
        """
        Display a single row.
        """

        outlist = []
        for idx in range(len(self)):
            outlist.append(repr(self[idx]))
        return f"({', '.join(outlist)})"

    def field(self, field):
        """
        Get the field data of the record.
        """

        return self.__getitem__(field)

    def setfield(self, field, value):
        """
        Set the field data of the record.
        """

        self.__setitem__(field, value)

    @lazyproperty
    def _bases(self):
        bases = [weakref.proxy(self)]
        base = self.base
        while base:
            bases.append(base)
            base = base.base
        return bases

    def _get_index(self, indx):
        indices = np.ogrid[:self.array._nfields]
        for base in reversed(self._bases):
            if base.step < 1:
                s = slice(base.start, None, base.step)
            else:
                s = slice(base.start, base.end, base.step)
            indices = indices[s]
        return indices[indx]


class FITS_rec(np.recarray):
    """
    FITS record array class.

    `FITS_rec` is the data part of a table HDU's data part.  This is a layer
    over the `~numpy.recarray`, so we can deal with scaled columns.

    It inherits all of the standard methods from `numpy.ndarray`.
    """

    _record_type = FITS_record
    _character_as_bytes = False

    def __new__(subtype, input):
        """
        Construct a FITS record array from a recarray.
        """

        # input should be a record array
        if input.dtype.subdtype is None:
            self = np.recarray.__new__(subtype, input.shape, input.dtype,
                                       buf=input.data)
        else:
            self = np.recarray.__new__(subtype, input.shape, input.dtype,
                                       buf=input.data, strides=input.strides)

        self._init()
        if self.dtype.fields:
            self._nfields = len(self.dtype.fields)

        return self

    def __setstate__(self, state):
        meta = state[-1]
        column_state = state[-2]
        state = state[:-2]

        super().__setstate__(state)

        self._col_weakrefs = weakref.WeakSet()

        for attr, value in zip(meta, column_state):
            setattr(self, attr, value)

    def __reduce__(self):
        """
        Return a 3-tuple for pickling a FITS_rec. Use the super-class
        functionality but then add in a tuple of FITS_rec-specific
        values that get used in __setstate__.
        """

        reconst_func, reconst_func_args, state = super().__reduce__()

        # Define FITS_rec-specific attrs that get added to state
        column_state = []
        meta = []

        for attrs in ['_converted', '_heapoffset', '_heapsize', '_nfields',
                      '_gap', '_uint', 'parnames', '_coldefs']:

            with suppress(AttributeError):
                # _coldefs can be Delayed, and file objects cannot be
                # picked, it needs to be deepcopied first
                if attrs == '_coldefs':
                    column_state.append(self._coldefs.__deepcopy__(None))
                else:
                    column_state.append(getattr(self, attrs))
                meta.append(attrs)

        state = state + (column_state, meta)

        return reconst_func, reconst_func_args, state

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if isinstance(obj, FITS_rec):
            self._character_as_bytes = obj._character_as_bytes

        if isinstance(obj, FITS_rec) and obj.dtype == self.dtype:
            self._converted = obj._converted
            self._heapoffset = obj._heapoffset
            self._heapsize = obj._heapsize
            self._col_weakrefs = obj._col_weakrefs
            self._coldefs = obj._coldefs
            self._nfields = obj._nfields
            self._gap = obj._gap
            self._uint = obj._uint
        elif self.dtype.fields is not None:
            # This will allow regular ndarrays with fields, rather than
            # just other FITS_rec objects
            self._nfields = len(self.dtype.fields)
            self._converted = {}

            self._heapoffset = getattr(obj, '_heapoffset', 0)
            self._heapsize = getattr(obj, '_heapsize', 0)

            self._gap = getattr(obj, '_gap', 0)
            self._uint = getattr(obj, '_uint', False)
            self._col_weakrefs = weakref.WeakSet()
            self._coldefs = ColDefs(self)

            # Work around chicken-egg problem.  Column.array relies on the
            # _coldefs attribute to set up ref back to parent FITS_rec; however
            # in the above line the self._coldefs has not been assigned yet so
            # this fails.  This patches that up...
            for col in self._coldefs:
                del col.array
                col._parent_fits_rec = weakref.ref(self)
        else:
            self._init()

    def _init(self):
        """Initializes internal attributes specific to FITS-isms."""

        self._nfields = 0
        self._converted = {}
        self._heapoffset = 0
        self._heapsize = 0
        self._col_weakrefs = weakref.WeakSet()
        self._coldefs = None
        self._gap = 0
        self._uint = False

    @classmethod
    def from_columns(cls, columns, nrows=0, fill=False, character_as_bytes=False):
        """
        Given a `ColDefs` object of unknown origin, initialize a new `FITS_rec`
        object.

        .. note::

            This was originally part of the ``new_table`` function in the table
            module but was moved into a class method since most of its
            functionality always had more to do with initializing a `FITS_rec`
            object than anything else, and much of it also overlapped with
            ``FITS_rec._scale_back``.

        Parameters
        ----------
        columns : sequence of `Column` or a `ColDefs`
            The columns from which to create the table data.  If these
            columns have data arrays attached that data may be used in
            initializing the new table.  Otherwise the input columns
            will be used as a template for a new table with the requested
            number of rows.

        nrows : int
            Number of rows in the new table.  If the input columns have data
            associated with them, the size of the largest input column is used.
            Otherwise the default is 0.

        fill : bool
            If `True`, will fill all cells with zeros or blanks.  If
            `False`, copy the data from input, undefined cells will still
            be filled with zeros/blanks.
        """

        if not isinstance(columns, ColDefs):
            columns = ColDefs(columns)

        # read the delayed data
        for column in columns:
            arr = column.array
            if isinstance(arr, Delayed):
                if arr.hdu.data is None:
                    column.array = None
                else:
                    column.array = _get_recarray_field(arr.hdu.data,
                                                       arr.field)
        # Reset columns._arrays (which we may want to just do away with
        # altogether
        del columns._arrays

        # use the largest column shape as the shape of the record
        if nrows == 0:
            for arr in columns._arrays:
                if arr is not None:
                    dim = arr.shape[0]
                else:
                    dim = 0
                if dim > nrows:
                    nrows = dim

        raw_data = np.empty(columns.dtype.itemsize * nrows, dtype=np.uint8)
        raw_data.fill(ord(columns._padding_byte))
        data = np.recarray(nrows, dtype=columns.dtype, buf=raw_data).view(cls)
        data._character_as_bytes = character_as_bytes

        # Previously this assignment was made from hdu.columns, but that's a
        # bug since if a _TableBaseHDU has a FITS_rec in its .data attribute
        # the _TableBaseHDU.columns property is actually returned from
        # .data._coldefs, so this assignment was circular!  Don't make that
        # mistake again.
        # All of this is an artifact of the fragility of the FITS_rec class,
        # and that it can't just be initialized by columns...
        data._coldefs = columns

        # If fill is True we don't copy anything from the column arrays.  We're
        # just using them as a template, and returning a table filled with
        # zeros/blanks
        if fill:
            return data

        # Otherwise we have to fill the recarray with data from the input
        # columns
        for idx, column in enumerate(columns):
            # For each column in the ColDef object, determine the number of
            # rows in that column.  This will be either the number of rows in
            # the ndarray associated with the column, or the number of rows
            # given in the call to this function, which ever is smaller.  If
            # the input FILL argument is true, the number of rows is set to
            # zero so that no data is copied from the original input data.
            arr = column.array

            if arr is None:
                array_size = 0
            else:
                array_size = len(arr)

            n = min(array_size, nrows)

            # TODO: At least *some* of this logic is mostly redundant with the
            # _convert_foo methods in this class; see if we can eliminate some
            # of that duplication.

            if not n:
                # The input column had an empty array, so just use the fill
                # value
                continue

            field = _get_recarray_field(data, idx)
            name = column.name
            fitsformat = column.format
            recformat = fitsformat.recformat

            outarr = field[:n]
            inarr = arr[:n]

            if isinstance(recformat, _FormatX):
                # Data is a bit array
                if inarr.shape[-1] == recformat.repeat:
                    _wrapx(inarr, outarr, recformat.repeat)
                    continue
            elif isinstance(recformat, _FormatP):
                data._cache_field(name, _makep(inarr, field, recformat,
                                               nrows=nrows))
                continue
            # TODO: Find a better way of determining that the column is meant
            # to be FITS L formatted
            elif recformat[-2:] == FITS2NUMPY['L'] and inarr.dtype == bool:
                # column is boolean
                # The raw data field should be filled with either 'T' or 'F'
                # (not 0).  Use 'F' as a default
                field[:] = ord('F')
                # Also save the original boolean array in data._converted so
                # that it doesn't have to be re-converted
                converted = np.zeros(field.shape, dtype=bool)
                converted[:n] = inarr
                data._cache_field(name, converted)
                # TODO: Maybe this step isn't necessary at all if _scale_back
                # will handle it?
                inarr = np.where(inarr == np.False_, ord('F'), ord('T'))
            elif (columns[idx]._physical_values and
                    columns[idx]._pseudo_unsigned_ints):
                # Temporary hack...
                bzero = column.bzero
                converted = np.zeros(field.shape, dtype=inarr.dtype)
                converted[:n] = inarr
                data._cache_field(name, converted)
                if n < nrows:
                    # Pre-scale rows below the input data
                    field[n:] = -bzero

                inarr = inarr - bzero
            elif isinstance(columns, _AsciiColDefs):
                # Regardless whether the format is character or numeric, if the
                # input array contains characters then it's already in the raw
                # format for ASCII tables
                if fitsformat._pseudo_logical:
                    # Hack to support converting from 8-bit T/F characters
                    # Normally the column array is a chararray of 1 character
                    # strings, but we need to view it as a normal ndarray of
                    # 8-bit ints to fill it with ASCII codes for 'T' and 'F'
                    outarr = field.view(np.uint8, np.ndarray)[:n]
                elif arr.dtype.kind not in ('S', 'U'):
                    # Set up views of numeric columns with the appropriate
                    # numeric dtype
                    # Fill with the appropriate blanks for the column format
                    data._cache_field(name, np.zeros(nrows, dtype=arr.dtype))
                    outarr = data._converted[name][:n]

                outarr[:] = inarr
                continue

            if inarr.shape != outarr.shape:
                if (inarr.dtype.kind == outarr.dtype.kind and
                        inarr.dtype.kind in ('U', 'S') and
                        inarr.dtype != outarr.dtype):

                    inarr_rowsize = inarr[0].size
                    inarr = inarr.flatten().view(outarr.dtype)

                # This is a special case to handle input arrays with
                # non-trivial TDIMn.
                # By design each row of the outarray is 1-D, while each row of
                # the input array may be n-D
                if outarr.ndim > 1:
                    # The normal case where the first dimension is the rows
                    inarr_rowsize = inarr[0].size
                    inarr = inarr.reshape(n, inarr_rowsize)
                    outarr[:, :inarr_rowsize] = inarr
                else:
                    # Special case for strings where the out array only has one
                    # dimension (the second dimension is rolled up into the
                    # strings
                    outarr[:n] = inarr.ravel()
            else:
                outarr[:] = inarr

        # Now replace the original column array references with the new
        # fields
        # This is required to prevent the issue reported in
        # https://github.com/spacetelescope/PyFITS/issues/99
        for idx in range(len(columns)):
            columns._arrays[idx] = data.field(idx)

        return data

    def __repr__(self):
        # Force use of the normal ndarray repr (rather than the new
        # one added for recarray in Numpy 1.10) for backwards compat
        return np.ndarray.__repr__(self)

    def __getattribute__(self, attr):
        # First, see if ndarray has this attr, and return it if so. Note that
        # this means a field with the same name as an ndarray attr cannot be
        # accessed by attribute, this is Numpy's default behavior.
        # We avoid using np.recarray.__getattribute__ here because after doing
        # this check it would access the columns without doing the conversions
        # that we need (with .field, see below).
        try:
            return object.__getattribute__(self, attr)
        except AttributeError:
            pass

        # attr might still be a fieldname.  If we have column definitions,
        # we should access this via .field, as the data may have to be scaled.
        if self._coldefs is not None and attr in self.columns.names:
            return self.field(attr)

        # If not, just let the usual np.recarray override deal with it.
        return super().__getattribute__(attr)

    def __getitem__(self, key):
        if self._coldefs is None:
            return super().__getitem__(key)

        if isinstance(key, str):
            return self.field(key)

        # Have to view as a recarray then back as a FITS_rec, otherwise the
        # circular reference fix/hack in FITS_rec.field() won't preserve
        # the slice.
        out = self.view(np.recarray)[key]
        if type(out) is not np.recarray:
            # Oops, we got a single element rather than a view. In that case,
            # return a Record, which has no __getstate__ and is more efficient.
            return self._record_type(self, key)

        # We got a view; change it back to our class, and add stuff
        out = out.view(type(self))
        out._uint = self._uint
        out._coldefs = ColDefs(self._coldefs)
        arrays = []
        out._converted = {}
        for idx, name in enumerate(self._coldefs.names):
            #
            # Store the new arrays for the _coldefs object
            #
            arrays.append(self._coldefs._arrays[idx][key])

            # Ensure that the sliced FITS_rec will view the same scaled
            # columns as the original; this is one of the few cases where
            # it is not necessary to use _cache_field()
            if name in self._converted:
                dummy = self._converted[name]
                field = np.ndarray.__getitem__(dummy, key)
                out._converted[name] = field

        out._coldefs._arrays = arrays
        return out

    def __setitem__(self, key, value):
        if self._coldefs is None:
            return super().__setitem__(key, value)

        if isinstance(key, str):
            self[key][:] = value
            return

        if isinstance(key, slice):
            end = min(len(self), key.stop or len(self))
            end = max(0, end)
            start = max(0, key.start or 0)
            end = min(end, start + len(value))

            for idx in range(start, end):
                self.__setitem__(idx, value[idx - start])
            return

        if isinstance(value, FITS_record):
            for idx in range(self._nfields):
                self.field(self.names[idx])[key] = value.field(self.names[idx])
        elif isinstance(value, (tuple, list, np.void)):
            if self._nfields == len(value):
                for idx in range(self._nfields):
                    self.field(idx)[key] = value[idx]
            else:
                raise ValueError('Input tuple or list required to have {} '
                                 'elements.'.format(self._nfields))
        else:
            raise TypeError('Assignment requires a FITS_record, tuple, or '
                            'list as input.')

    def _ipython_key_completions_(self):
        return self.names

    def copy(self, order='C'):
        """
        The Numpy documentation lies; `numpy.ndarray.copy` is not equivalent to
        `numpy.copy`.  Differences include that it re-views the copied array as
        self's ndarray subclass, as though it were taking a slice; this means
        ``__array_finalize__`` is called and the copy shares all the array
        attributes (including ``._converted``!).  So we need to make a deep
        copy of all those attributes so that the two arrays truly do not share
        any data.
        """

        new = super().copy(order=order)

        new.__dict__ = copy.deepcopy(self.__dict__)
        return new

    @property
    def columns(self):
        """A user-visible accessor for the coldefs."""

        return self._coldefs

    @property
    def _coldefs(self):
        # This used to be a normal internal attribute, but it was changed to a
        # property as a quick and transparent way to work around the reference
        # leak bug fixed in https://github.com/astropy/astropy/pull/4539
        #
        # See the long comment in the Column.array property for more details
        # on this.  But in short, FITS_rec now has a ._col_weakrefs attribute
        # which is a WeakSet of weakrefs to each Column in _coldefs.
        #
        # So whenever ._coldefs is set we also add each Column in the ColDefs
        # to the weakrefs set.  This is an easy way to find out if a Column has
        # any references to it external to the FITS_rec (i.e. a user assigned a
        # column to a variable).  If the column is still in _col_weakrefs then
        # there are other references to it external to this FITS_rec.  We use
        # that information in __del__ to save off copies of the array data
        # for those columns to their Column.array property before our memory
        # is freed.
        return self.__dict__.get('_coldefs')

    @_coldefs.setter
    def _coldefs(self, cols):
        self.__dict__['_coldefs'] = cols
        if isinstance(cols, ColDefs):
            for col in cols.columns:
                self._col_weakrefs.add(col)

    @_coldefs.deleter
    def _coldefs(self):
        try:
            del self.__dict__['_coldefs']
        except KeyError as exc:
            raise AttributeError(exc.args[0])

    def __del__(self):
        try:
            del self._coldefs
            if self.dtype.fields is not None:
                for col in self._col_weakrefs:

                    if col.array is not None:
                        col.array = col.array.copy()

        # See issues #4690 and #4912
        except (AttributeError, TypeError):  # pragma: no cover
            pass

    @property
    def names(self):
        """List of column names."""

        if self.dtype.fields:
            return list(self.dtype.names)
        elif getattr(self, '_coldefs', None) is not None:
            return self._coldefs.names
        else:
            return None

    @property
    def formats(self):
        """List of column FITS formats."""

        if getattr(self, '_coldefs', None) is not None:
            return self._coldefs.formats

        return None

    @property
    def _raw_itemsize(self):
        """
        Returns the size of row items that would be written to the raw FITS
        file, taking into account the possibility of unicode columns being
        compactified.

        Currently for internal use only.
        """

        if _has_unicode_fields(self):
            total_itemsize = 0
            for field in self.dtype.fields.values():
                itemsize = field[0].itemsize
                if field[0].kind == 'U':
                    itemsize = itemsize // 4
                total_itemsize += itemsize
            return total_itemsize
        else:
            # Just return the normal itemsize
            return self.itemsize

    def field(self, key):
        """
        A view of a `Column`'s data as an array.
        """

        # NOTE: The *column* index may not be the same as the field index in
        # the recarray, if the column is a phantom column
        column = self.columns[key]
        name = column.name
        format = column.format

        if format.dtype.itemsize == 0:
            warnings.warn(
                'Field {!r} has a repeat count of 0 in its format code, '
                'indicating an empty field.'.format(key))
            return np.array([], dtype=format.dtype)

        # If field's base is a FITS_rec, we can run into trouble because it
        # contains a reference to the ._coldefs object of the original data;
        # this can lead to a circular reference; see ticket #49
        base = self
        while (isinstance(base, FITS_rec) and
                isinstance(base.base, np.recarray)):
            base = base.base
        # base could still be a FITS_rec in some cases, so take care to
        # use rec.recarray.field to avoid a potential infinite
        # recursion
        field = _get_recarray_field(base, name)

        if name not in self._converted:
            recformat = format.recformat
            # TODO: If we're now passing the column to these subroutines, do we
            # really need to pass them the recformat?
            if isinstance(recformat, _FormatP):
                # for P format
                converted = self._convert_p(column, field, recformat)
            else:
                # Handle all other column data types which are fixed-width
                # fields
                converted = self._convert_other(column, field, recformat)

            # Note: Never assign values directly into the self._converted dict;
            # always go through self._cache_field; this way self._converted is
            # only used to store arrays that are not already direct views of
            # our own data.
            self._cache_field(name, converted)
            return converted

        return self._converted[name]

    def _cache_field(self, name, field):
        """
        Do not store fields in _converted if one of its bases is self,
        or if it has a common base with self.

        This results in a reference cycle that cannot be broken since
        ndarrays do not participate in cyclic garbage collection.
        """

        base = field
        while True:
            self_base = self
            while True:
                if self_base is base:
                    return

                if getattr(self_base, 'base', None) is not None:
                    self_base = self_base.base
                else:
                    break

            if getattr(base, 'base', None) is not None:
                base = base.base
            else:
                break

        self._converted[name] = field

    def _update_column_attribute_changed(self, column, idx, attr, old_value,
                                         new_value):
        """
        Update how the data is formatted depending on changes to column
        attributes initiated by the user through the `Column` interface.

        Dispatches column attribute change notifications to individual methods
        for each attribute ``_update_column_<attr>``
        """

        method_name = f'_update_column_{attr}'
        if hasattr(self, method_name):
            # Right now this is so we can be lazy and not implement updaters
            # for every attribute yet--some we may not need at all, TBD
            getattr(self, method_name)(column, idx, old_value, new_value)

    def _update_column_name(self, column, idx, old_name, name):
        """Update the dtype field names when a column name is changed."""

        dtype = self.dtype
        # Updating the names on the dtype should suffice
        dtype.names = dtype.names[:idx] + (name,) + dtype.names[idx + 1:]

    def _convert_x(self, field, recformat):
        """Convert a raw table column to a bit array as specified by the
        FITS X format.
        """

        dummy = np.zeros(self.shape + (recformat.repeat,), dtype=np.bool_)
        _unwrapx(field, dummy, recformat.repeat)
        return dummy

    def _convert_p(self, column, field, recformat):
        """Convert a raw table column of FITS P or Q format descriptors
        to a VLA column with the array data returned from the heap.
        """

        dummy = _VLF([None] * len(self), dtype=recformat.dtype)
        raw_data = self._get_raw_data()

        if raw_data is None:
            raise OSError(
                "Could not find heap data for the {!r} variable-length "
                "array column.".format(column.name))

        for idx in range(len(self)):
            offset = field[idx, 1] + self._heapoffset
            count = field[idx, 0]

            if recformat.dtype == 'a':
                dt = np.dtype(recformat.dtype + str(1))
                arr_len = count * dt.itemsize
                da = raw_data[offset:offset + arr_len].view(dt)
                da = np.char.array(da.view(dtype=dt), itemsize=count)
                dummy[idx] = decode_ascii(da)
            else:
                dt = np.dtype(recformat.dtype)
                arr_len = count * dt.itemsize
                dummy[idx] = raw_data[offset:offset + arr_len].view(dt)
                dummy[idx].dtype = dummy[idx].dtype.newbyteorder('>')
                # Each array in the field may now require additional
                # scaling depending on the other scaling parameters
                # TODO: The same scaling parameters apply to every
                # array in the column so this is currently very slow; we
                # really only need to check once whether any scaling will
                # be necessary and skip this step if not
                # TODO: Test that this works for X format; I don't think
                # that it does--the recformat variable only applies to the P
                # format not the X format
                dummy[idx] = self._convert_other(column, dummy[idx],
                                                 recformat)

        return dummy

    def _convert_ascii(self, column, field):
        """
        Special handling for ASCII table columns to convert columns containing
        numeric types to actual numeric arrays from the string representation.
        """

        format = column.format
        recformat = getattr(format, 'recformat', ASCII2NUMPY[format[0]])
        # if the string = TNULL, return ASCIITNULL
        nullval = str(column.null).strip().encode('ascii')
        if len(nullval) > format.width:
            nullval = nullval[:format.width]

        # Before using .replace make sure that any trailing bytes in each
        # column are filled with spaces, and *not*, say, nulls; this causes
        # functions like replace to potentially leave gibberish bytes in the
        # array buffer.
        dummy = np.char.ljust(field, format.width)
        dummy = np.char.replace(dummy, encode_ascii('D'), encode_ascii('E'))
        null_fill = encode_ascii(str(ASCIITNULL).rjust(format.width))

        # Convert all fields equal to the TNULL value (nullval) to empty fields.
        # TODO: These fields really should be converted to NaN or something else undefined.
        # Currently they are converted to empty fields, which are then set to zero.
        dummy = np.where(np.char.strip(dummy) == nullval, null_fill, dummy)

        # always replace empty fields, see https://github.com/astropy/astropy/pull/5394
        if nullval != b'':
            dummy = np.where(np.char.strip(dummy) == b'', null_fill, dummy)

        try:
            dummy = np.array(dummy, dtype=recformat)
        except ValueError as exc:
            indx = self.names.index(column.name)
            raise ValueError(
                '{}; the header may be missing the necessary TNULL{} '
                'keyword or the table contains invalid data'.format(
                    exc, indx + 1))

        return dummy

    def _convert_other(self, column, field, recformat):
        """Perform conversions on any other fixed-width column data types.

        This may not perform any conversion at all if it's not necessary, in
        which case the original column array is returned.
        """

        if isinstance(recformat, _FormatX):
            # special handling for the X format
            return self._convert_x(field, recformat)

        (_str, _bool, _number, _scale, _zero, bscale, bzero, dim) = \
            self._get_scale_factors(column)

        indx = self.names.index(column.name)

        # ASCII table, convert strings to numbers
        # TODO:
        # For now, check that these are ASCII columns by checking the coldefs
        # type; in the future all columns (for binary tables, ASCII tables, or
        # otherwise) should "know" what type they are already and how to handle
        # converting their data from FITS format to native format and vice
        # versa...
        if not _str and isinstance(self._coldefs, _AsciiColDefs):
            field = self._convert_ascii(column, field)

        # Test that the dimensions given in dim are sensible; otherwise
        # display a warning and ignore them
        if dim:
            # See if the dimensions already match, if not, make sure the
            # number items will fit in the specified dimensions
            if field.ndim > 1:
                actual_shape = field.shape[1:]
                if _str:
                    actual_shape = actual_shape + (field.itemsize,)
            else:
                actual_shape = field.shape[0]

            if dim == actual_shape:
                # The array already has the correct dimensions, so we
                # ignore dim and don't convert
                dim = None
            else:
                nitems = reduce(operator.mul, dim)
                if _str:
                    actual_nitems = field.itemsize
                elif len(field.shape) == 1:  # No repeat count in TFORMn, equivalent to 1
                    actual_nitems = 1
                else:
                    actual_nitems = field.shape[1]
                if nitems > actual_nitems:
                    warnings.warn(
                        'TDIM{} value {:d} does not fit with the size of '
                        'the array items ({:d}).  TDIM{:d} will be ignored.'
                        .format(indx + 1, self._coldefs[indx].dims,
                                actual_nitems, indx + 1))
                    dim = None

        # further conversion for both ASCII and binary tables
        # For now we've made columns responsible for *knowing* whether their
        # data has been scaled, but we make the FITS_rec class responsible for
        # actually doing the scaling
        # TODO: This also needs to be fixed in the effort to make Columns
        # responsible for scaling their arrays to/from FITS native values
        if not column.ascii and column.format.p_format:
            format_code = column.format.p_format
        else:
            # TODO: Rather than having this if/else it might be nice if the
            # ColumnFormat class had an attribute guaranteed to give the format
            # of actual values in a column regardless of whether the true
            # format is something like P or Q
            format_code = column.format.format

        if (_number and (_scale or _zero) and not column._physical_values):
            # This is to handle pseudo unsigned ints in table columns
            # TODO: For now this only really works correctly for binary tables
            # Should it work for ASCII tables as well?
            if self._uint:
                if bzero == 2**15 and format_code == 'I':
                    field = np.array(field, dtype=np.uint16)
                elif bzero == 2**31 and format_code == 'J':
                    field = np.array(field, dtype=np.uint32)
                elif bzero == 2**63 and format_code == 'K':
                    field = np.array(field, dtype=np.uint64)
                    bzero64 = np.uint64(2 ** 63)
                else:
                    field = np.array(field, dtype=np.float64)
            else:
                field = np.array(field, dtype=np.float64)

            if _scale:
                np.multiply(field, bscale, field)
            if _zero:
                if self._uint and format_code == 'K':
                    # There is a chance of overflow, so be careful
                    test_overflow = field.copy()
                    try:
                        test_overflow += bzero64
                    except OverflowError:
                        warnings.warn(
                            "Overflow detected while applying TZERO{:d}. "
                            "Returning unscaled data.".format(indx + 1))
                    else:
                        field = test_overflow
                else:
                    field += bzero

            # mark the column as scaled
            column._physical_values = True

        elif _bool and field.dtype != bool:
            field = np.equal(field, ord('T'))
        elif _str:
            if not self._character_as_bytes:
                with suppress(UnicodeDecodeError):
                    field = decode_ascii(field)

        if dim:
            # Apply the new field item dimensions
            nitems = reduce(operator.mul, dim)
            if field.ndim > 1:
                field = field[:, :nitems]
            if _str:
                fmt = field.dtype.char
                dtype = (f'|{fmt}{dim[-1]}', dim[:-1])
                field.dtype = dtype
            else:
                field.shape = (field.shape[0],) + dim

        return field

    def _get_heap_data(self):
        """
        Returns a pointer into the table's raw data to its heap (if present).

        This is returned as a numpy byte array.
        """

        if self._heapsize:
            raw_data = self._get_raw_data().view(np.ubyte)
            heap_end = self._heapoffset + self._heapsize
            return raw_data[self._heapoffset:heap_end]
        else:
            return np.array([], dtype=np.ubyte)

    def _get_raw_data(self):
        """
        Returns the base array of self that "raw data array" that is the
        array in the format that it was first read from a file before it was
        sliced or viewed as a different type in any way.

        This is determined by walking through the bases until finding one that
        has at least the same number of bytes as self, plus the heapsize.  This
        may be the immediate .base but is not always.  This is used primarily
        for variable-length array support which needs to be able to find the
        heap (the raw data *may* be larger than nbytes + heapsize if it
        contains a gap or padding).

        May return ``None`` if no array resembling the "raw data" according to
        the stated criteria can be found.
        """

        raw_data_bytes = self.nbytes + self._heapsize
        base = self
        while hasattr(base, 'base') and base.base is not None:
            base = base.base
            if hasattr(base, 'nbytes') and base.nbytes >= raw_data_bytes:
                return base

    def _get_scale_factors(self, column):
        """Get all the scaling flags and factors for one column."""

        # TODO: Maybe this should be a method/property on Column?  Or maybe
        # it's not really needed at all...
        _str = column.format.format == 'A'
        _bool = column.format.format == 'L'

        _number = not (_bool or _str)
        bscale = column.bscale
        bzero = column.bzero

        _scale = bscale not in ('', None, 1)
        _zero = bzero not in ('', None, 0)

        # ensure bscale/bzero are numbers
        if not _scale:
            bscale = 1
        if not _zero:
            bzero = 0

        # column._dims gives a tuple, rather than column.dim which returns the
        # original string format code from the FITS header...
        dim = column._dims

        return (_str, _bool, _number, _scale, _zero, bscale, bzero, dim)

    def _scale_back(self, update_heap_pointers=True):
        """
        Update the parent array, using the (latest) scaled array.

        If ``update_heap_pointers`` is `False`, this will leave all the heap
        pointers in P/Q columns as they are verbatim--it only makes sense to do
        this if there is already data on the heap and it can be guaranteed that
        that data has not been modified, and there is not new data to add to
        the heap.  Currently this is only used as an optimization for
        CompImageHDU that does its own handling of the heap.
        """

        # Running total for the new heap size
        heapsize = 0

        for indx, name in enumerate(self.dtype.names):
            column = self._coldefs[indx]
            recformat = column.format.recformat
            raw_field = _get_recarray_field(self, indx)

            # add the location offset of the heap area for each
            # variable length column
            if isinstance(recformat, _FormatP):
                # Irritatingly, this can return a different dtype than just
                # doing np.dtype(recformat.dtype); but this returns the results
                # that we want.  For example if recformat.dtype is 'a' we want
                # an array of characters.
                dtype = np.array([], dtype=recformat.dtype).dtype

                if update_heap_pointers and name in self._converted:
                    # The VLA has potentially been updated, so we need to
                    # update the array descriptors
                    raw_field[:] = 0  # reset
                    npts = [len(arr) for arr in self._converted[name]]

                    raw_field[:len(npts), 0] = npts
                    raw_field[1:, 1] = (np.add.accumulate(raw_field[:-1, 0]) *
                                        dtype.itemsize)
                    raw_field[:, 1][:] += heapsize

                heapsize += raw_field[:, 0].sum() * dtype.itemsize
                # Even if this VLA has not been read or updated, we need to
                # include the size of its constituent arrays in the heap size
                # total

            if isinstance(recformat, _FormatX) and name in self._converted:
                _wrapx(self._converted[name], raw_field, recformat.repeat)
                continue

            _str, _bool, _number, _scale, _zero, bscale, bzero, _ = \
                self._get_scale_factors(column)

            field = self._converted.get(name, raw_field)

            # conversion for both ASCII and binary tables
            if _number or _str:
                if _number and (_scale or _zero) and column._physical_values:
                    dummy = field.copy()
                    if _zero:
                        dummy -= bzero
                    if _scale:
                        dummy /= bscale
                    # This will set the raw values in the recarray back to
                    # their non-physical storage values, so the column should
                    # be mark is not scaled
                    column._physical_values = False
                elif _str or isinstance(self._coldefs, _AsciiColDefs):
                    dummy = field
                else:
                    continue

                # ASCII table, convert numbers to strings
                if isinstance(self._coldefs, _AsciiColDefs):
                    self._scale_back_ascii(indx, dummy, raw_field)
                # binary table string column
                elif isinstance(raw_field, chararray.chararray):
                    self._scale_back_strings(indx, dummy, raw_field)
                # all other binary table columns
                else:
                    if len(raw_field) and isinstance(raw_field[0],
                                                     np.integer):
                        dummy = np.around(dummy)

                    if raw_field.shape == dummy.shape:
                        raw_field[:] = dummy
                    else:
                        # Reshaping the data is necessary in cases where the
                        # TDIMn keyword was used to shape a column's entries
                        # into arrays
                        raw_field[:] = dummy.ravel().view(raw_field.dtype)

                del dummy

            # ASCII table does not have Boolean type
            elif _bool and name in self._converted:
                choices = (np.array([ord('F')], dtype=np.int8)[0],
                           np.array([ord('T')], dtype=np.int8)[0])
                raw_field[:] = np.choose(field, choices)

        # Store the updated heapsize
        self._heapsize = heapsize

    def _scale_back_strings(self, col_idx, input_field, output_field):
        # There are a few possibilities this has to be able to handle properly
        # The input_field, which comes from the _converted column is of dtype
        # 'Un' so that elements read out of the array are normal str
        # objects (i.e. unicode strings)
        #
        # At the other end the *output_field* may also be of type 'S' or of
        # type 'U'.  It will *usually* be of type 'S' because when reading
        # an existing FITS table the raw data is just ASCII strings, and
        # represented in Numpy as an S array.  However, when a user creates
        # a new table from scratch, they *might* pass in a column containing
        # unicode strings (dtype 'U').  Therefore the output_field of the
        # raw array is actually a unicode array.  But we still want to make
        # sure the data is encodable as ASCII.  Later when we write out the
        # array we use, in the dtype 'U' case, a different write routine
        # that writes row by row and encodes any 'U' columns to ASCII.

        # If the output_field is non-ASCII we will worry about ASCII encoding
        # later when writing; otherwise we can do it right here
        if input_field.dtype.kind == 'U' and output_field.dtype.kind == 'S':
            try:
                _ascii_encode(input_field, out=output_field)
            except _UnicodeArrayEncodeError as exc:
                raise ValueError(
                    "Could not save column '{}': Contains characters that "
                    "cannot be encoded as ASCII as required by FITS, starting "
                    "at the index {!r} of the column, and the index {} of "
                    "the string at that location.".format(
                        self._coldefs[col_idx].name,
                        exc.index[0] if len(exc.index) == 1 else exc.index,
                        exc.start))
        else:
            # Otherwise go ahead and do a direct copy into--if both are type
            # 'U' we'll handle encoding later
            input_field = input_field.flatten().view(output_field.dtype)
            output_field.flat[:] = input_field

        # Ensure that blanks at the end of each string are
        # converted to nulls instead of spaces, see Trac #15
        # and #111
        _rstrip_inplace(output_field)

    def _scale_back_ascii(self, col_idx, input_field, output_field):
        """
        Convert internal array values back to ASCII table representation.

        The ``input_field`` is the internal representation of the values, and
        the ``output_field`` is the character array representing the ASCII
        output that will be written.
        """

        starts = self._coldefs.starts[:]
        spans = self._coldefs.spans
        format = self._coldefs[col_idx].format

        # The the index of the "end" column of the record, beyond
        # which we can't write
        end = super().field(-1).itemsize
        starts.append(end + starts[-1])

        if col_idx > 0:
            lead = starts[col_idx] - starts[col_idx - 1] - spans[col_idx - 1]
        else:
            lead = 0

        if lead < 0:
            warnings.warn('Column {!r} starting point overlaps the previous '
                          'column.'.format(col_idx + 1))

        trail = starts[col_idx + 1] - starts[col_idx] - spans[col_idx]

        if trail < 0:
            warnings.warn('Column {!r} ending point overlaps the next '
                          'column.'.format(col_idx + 1))

        # TODO: It would be nice if these string column formatting
        # details were left to a specialized class, as is the case
        # with FormatX and FormatP
        if 'A' in format:
            _pc = '{:'
        else:
            _pc = '{:>'

        fmt = ''.join([_pc, format[1:], ASCII2STR[format[0]], '}',
                       (' ' * trail)])

        # Even if the format precision is 0, we should output a decimal point
        # as long as there is space to do so--not including a decimal point in
        # a float value is discouraged by the FITS Standard
        trailing_decimal = (format.precision == 0 and
                            format.format in ('F', 'E', 'D'))

        # not using numarray.strings's num2char because the
        # result is not allowed to expand (as C/Python does).
        for jdx, value in enumerate(input_field):
            value = fmt.format(value)
            if len(value) > starts[col_idx + 1] - starts[col_idx]:
                raise ValueError(
                    "Value {!r} does not fit into the output's itemsize of "
                    "{}.".format(value, spans[col_idx]))

            if trailing_decimal and value[0] == ' ':
                # We have some extra space in the field for the trailing
                # decimal point
                value = value[1:] + '.'

            output_field[jdx] = value

        # Replace exponent separator in floating point numbers
        if 'D' in format:
            output_field[:] = output_field.replace(b'E', b'D')

    def tolist(self):
        # Override .tolist to take care of special case of VLF

        column_lists = [self[name].tolist() for name in self.columns.names]

        return [list(row) for row in zip(*column_lists)]


def _get_recarray_field(array, key):
    """
    Compatibility function for using the recarray base class's field method.
    This incorporates the legacy functionality of returning string arrays as
    Numeric-style chararray objects.
    """

    # Numpy >= 1.10.dev recarray no longer returns chararrays for strings
    # This is currently needed for backwards-compatibility and for
    # automatic truncation of trailing whitespace
    field = np.recarray.field(array, key)
    if (field.dtype.char in ('S', 'U') and
            not isinstance(field, chararray.chararray)):
        field = field.view(chararray.chararray)
    return field


class _UnicodeArrayEncodeError(UnicodeEncodeError):
    def __init__(self, encoding, object_, start, end, reason, index):
        super().__init__(encoding, object_, start, end, reason)
        self.index = index


def _ascii_encode(inarray, out=None):
    """
    Takes a unicode array and fills the output string array with the ASCII
    encodings (if possible) of the elements of the input array.  The two arrays
    must be the same size (though not necessarily the same shape).

    This is like an inplace version of `np.char.encode` though simpler since
    it's only limited to ASCII, and hence the size of each character is
    guaranteed to be 1 byte.

    If any strings are non-ASCII an UnicodeArrayEncodeError is raised--this is
    just a `UnicodeEncodeError` with an additional attribute for the index of
    the item that couldn't be encoded.
    """

    out_dtype = np.dtype((f'S{inarray.dtype.itemsize // 4}',
                         inarray.dtype.shape))
    if out is not None:
        out = out.view(out_dtype)

    op_dtypes = [inarray.dtype, out_dtype]
    op_flags = [['readonly'], ['writeonly', 'allocate']]
    it = np.nditer([inarray, out], op_dtypes=op_dtypes,
                   op_flags=op_flags, flags=['zerosize_ok'])

    try:
        for initem, outitem in it:
            outitem[...] = initem.item().encode('ascii')
    except UnicodeEncodeError as exc:
        index = np.unravel_index(it.iterindex, inarray.shape)
        raise _UnicodeArrayEncodeError(*(exc.args + (index,)))

    return it.operands[1]


def _has_unicode_fields(array):
    """
    Returns True if any fields in a structured array have Unicode dtype.
    """

    dtypes = (d[0] for d in array.dtype.fields.values())
    return any(d.kind == 'U' for d in dtypes)
