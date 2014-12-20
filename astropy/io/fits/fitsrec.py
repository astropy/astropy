# Licensed under a 3-clause BSD style license - see PYFITS.rst

import copy
import operator
import warnings
import weakref

from functools import reduce

import numpy as np

from numpy import char as chararray

from .column import (ASCIITNULL, FITS2NUMPY, ASCII2NUMPY, ASCII2STR, ColDefs,
                     _AsciiColDefs, _FormatX, _FormatP, _VLF, _get_index,
                     _wrapx, _unwrapx, _makep, _convert_ascii_format, Delayed)
from .util import decode_ascii, encode_ascii
from ...extern.six import string_types
from ...extern.six.moves import xrange, map
from ...utils import lazyproperty
from ...utils.compat import ignored
from ...utils.exceptions import AstropyDeprecationWarning, AstropyUserWarning


class FITS_record(object):
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

        # For backward compatibility...
        for arg in [('startColumn', 'start'), ('endColumn', 'end')]:
            if arg[0] in kwargs:
                warnings.warn('The %s argument to FITS_record is deprecated; '
                              'use %s instead' % arg, AstropyDeprecationWarning)
                if arg[0] == 'startColumn':
                    start = kwargs[arg[0]]
                elif arg[0] == 'endColumn':
                    end = kwargs[arg[0]]

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
        if isinstance(key, string_types):
            indx = _get_index(self.array.names, key)

            if indx < self.start or indx > self.end - 1:
                raise KeyError("Key '%s' does not exist." % key)
        elif isinstance(key, slice):
            return type(self)(self.array, self.row, key.start, key.stop,
                              key.step, self)
        else:
            indx = self._get_index(key)

            if indx > self.array._nfields - 1:
                raise IndexError('Index out of bounds')

        return self.array.field(indx)[self.row]

    def __setitem__(self, key, value):
        if isinstance(key, string_types):
            indx = _get_index(self.array._coldefs.names, key)

            if indx < self.start or indx > self.end - 1:
                raise KeyError("Key '%s' does not exist." % key)
        elif isinstance(key, slice):
            for indx in xrange(slice.start, slice.stop, slice.step):
                indx = self._get_indx(indx)
                self.array.field(indx)[self.row] = value
        else:
            indx = self._get_index(key)
            if indx > self.array._nfields - 1:
                raise IndexError('Index out of bounds')

        self.array.field(indx)[self.row] = value

    def __getslice__(self, start, end):
        return self[slice(start, end)]

    def __len__(self):
        return len(xrange(self.start, self.end, self.step))

    def __repr__(self):
        """
        Display a single row.
        """

        outlist = []
        for idx in xrange(len(self)):
            outlist.append(repr(self[idx]))
        return '(%s)' % ', '.join(outlist)

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

        self._nfields = len(self.dtype.names)
        self._convert = [None] * len(self.dtype.names)
        self._heapoffset = 0
        self._heapsize = 0
        self._coldefs = None
        self._gap = 0
        self._uint = False
        self.formats = None
        return self

    def __setstate__(self, state):
        meta = state[-1]
        column_state = state[-2]
        state = state[:-2]

        super(FITS_rec, self).__setstate__(state)

        for attr, value in zip(meta, column_state):
            setattr(self, attr, value)

    def __reduce__(self):
        """
        Return a 3-tuple for pickling a FITS_rec. Use the super-class
        functionality but then add in a tuple of FITS_rec-specific
        values that get used in __setstate__.
        """

        reconst_func, reconst_func_args, state = \
            super(FITS_rec, self).__reduce__()

        # Define FITS_rec-specific attrs that get added to state
        column_state = []
        meta = []

        for attrs in ['_convert', '_heapoffset', '_heapsize', '_nfields',
                      '_gap', '_uint', 'formats', 'parnames', '_coldefs']:

            with ignored(AttributeError):
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
            self._convert = obj._convert
            self._heapoffset = obj._heapoffset
            self._heapsize = obj._heapsize
            self._coldefs = obj._coldefs
            self._nfields = obj._nfields
            self._gap = obj._gap
            self._uint = obj._uint
            self.formats = obj.formats
        else:
            # This will allow regular ndarrays with fields, rather than
            # just other FITS_rec objects
            self._nfields = len(obj.dtype.names)
            self._convert = [None] * len(obj.dtype.names)

            self._heapoffset = getattr(obj, '_heapoffset', 0)
            self._heapsize = getattr(obj, '_heapsize', 0)

            self._coldefs = None
            self._gap = getattr(obj, '_gap', 0)
            self._uint = getattr(obj, '_uint', False)

            # Bypass setattr-based assignment to fields; see #86
            self.formats = None

            attrs = ['_convert', '_coldefs', '_gap']
            for attr in attrs:
                if hasattr(obj, attr):
                    value = getattr(obj, attr, None)
                    if value is None:
                        warnings.warn('Setting attribute %s as None' % attr, AstropyUserWarning)
                    setattr(self, attr, value)

            if self._coldefs is None:
                self._coldefs = ColDefs(self)
            self.formats = self._coldefs.formats

    @classmethod
    def from_columns(cls, columns, nrows=0, fill=False):
        """
        Given a `ColDefs` object of unknown origin, initialize a new `FITS_rec`
        object.

        .. note::

            This was originally part of the `new_table` function in the table
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
        for idx in range(len(columns)):
            arr = columns._arrays[idx]
            if isinstance(arr, Delayed):
                if arr.hdu.data is None:
                    columns._arrays[idx] = None
                else:
                    columns._arrays[idx] = np.rec.recarray.field(arr.hdu.data,
                                                                 arr.field)

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

        # Previously this assignment was made from hdu.columns, but that's a
        # bug since if a _TableBaseHDU has a FITS_rec in its .data attribute
        # the _TableBaseHDU.columns property is actually returned from
        # .data._coldefs, so this assignment was circular!  Don't make that
        # mistake again.
        # All of this is an artifact of the fragility of the FITS_rec class,
        # and that it can't just be initialized by columns...
        data._coldefs = columns
        data.formats = columns.formats

        # If fill is True we don't copy anything from the column arrays.  We're
        # just using them as a template, and returning a table filled with
        # zeros/blanks
        if fill:
            return data

        # Otherwise we have to fill the recarray with data from the input
        # columns
        for idx in range(len(columns)):
            # For each column in the ColDef object, determine the number of
            # rows in that column.  This will be either the number of rows in
            # the ndarray associated with the column, or the number of rows
            # given in the call to this function, which ever is smaller.  If
            # the input FILL argument is true, the number of rows is set to
            # zero so that no data is copied from the original input data.
            arr = columns._arrays[idx]

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

            field = np.rec.recarray.field(data, idx)
            fitsformat = columns.formats[idx]
            recformat = columns._recformats[idx]

            outarr = field[:n]
            inarr = arr[:n]

            if isinstance(recformat, _FormatX):
                # Data is a bit array
                if inarr.shape[-1] == recformat.repeat:
                    _wrapx(inarr, outarr, recformat.repeat)
                    continue
            elif isinstance(recformat, _FormatP):
                data._convert[idx] = _makep(inarr, field, recformat,
                                            nrows=nrows)
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
                data._convert[idx] = np.zeros(field.shape, dtype=bool)
                data._convert[idx][:n] = inarr
                # TODO: Maybe this step isn't necessary at all if _scale_back
                # will handle it?
                inarr = np.where(inarr == False, ord('F'), ord('T'))
            elif (columns[idx]._physical_values and
                    columns[idx]._pseudo_unsigned_ints):
                # Temporary hack...
                bzero = columns[idx].bzero
                data._convert[idx] = np.zeros(field.shape, dtype=inarr.dtype)
                data._convert[idx][:n] = inarr
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
                elif not isinstance(arr, chararray.chararray):
                    # Fill with the appropriate blanks for the column format
                    data._convert[idx] = np.zeros(nrows, dtype=arr.dtype)
                    outarr = data._convert[idx][:n]

                outarr[:] = inarr
                continue

            if inarr.shape != outarr.shape:
                if inarr.dtype != outarr.dtype:
                    inarr = inarr.view(outarr.dtype)

                # This is a special case to handle input arrays with
                # non-trivial TDIMn.
                # By design each row of the outarray is 1-D, while each row of
                # the input array may be n-D
                if outarr.ndim > 1:
                    # The normal case where the first dimension is the rows
                    inarr_rowsize = inarr[0].size
                    inarr = inarr.reshape((n, inarr_rowsize))
                    outarr[:, :inarr_rowsize] = inarr
                else:
                    # Special case for strings where the out array only has one
                    # dimension (the second dimension is rolled up into the
                    # strings
                    outarr[:n] = inarr.ravel()
            else:
                outarr[:] = inarr

        return data

    def __repr__(self):
        return np.recarray.__repr__(self)

    def __getattribute__(self, attr):
        # See the comment in __setattr__
        if attr in ('names', 'formats'):
            return object.__getattribute__(self, attr)
        else:
            return super(FITS_rec, self).__getattribute__(attr)

    def __setattr__(self, attr, value):
        # Overrides the silly attribute-based assignment to fields supported by
        # recarrays for our two built-in public attributes: names and formats
        # Otherwise, the default behavior, bad as it is, is preserved.  See
        # ticket #86
        if attr in ('names', 'formats'):
            return object.__setattr__(self, attr, value)
        else:
            return super(FITS_rec, self).__setattr__(attr, value)

    def __getitem__(self, key):
        if isinstance(key, string_types):
            return self.field(key)
        elif isinstance(key, (slice, np.ndarray, tuple, list)):
            # Have to view as a recarray then back as a FITS_rec, otherwise the
            # circular reference fix/hack in FITS_rec.field() won't preserve
            # the slice
            subtype = type(self)
            out = self.view(np.recarray).__getitem__(key).view(subtype)
            out._coldefs = ColDefs(self._coldefs)
            arrays = []
            out._convert = [None] * len(self.dtype.names)
            for idx in range(len(self.dtype.names)):
                #
                # Store the new arrays for the _coldefs object
                #
                arrays.append(self._coldefs._arrays[idx][key])

                # touch all fields to expand the original ._convert list
                # so the sliced FITS_rec will view the same scaled columns as
                # the original
                dummy = self.field(idx)
                if self._convert[idx] is not None:
                    out._convert[idx] = \
                        np.ndarray.__getitem__(self._convert[idx], key)
            del dummy

            out._coldefs._arrays = arrays
            out._coldefs._shape = len(arrays[0])

            return out

        # if not a slice, do this because Record has no __getstate__.
        # also more efficient.
        else:
            if isinstance(key, int) and key >= len(self):
                raise IndexError("Index out of bounds")

            newrecord = self._record_type(self, key)
            return newrecord

    def __setitem__(self, key, value):
        if isinstance(key, string_types):
            self[key][:] = value
            return

        if isinstance(key, slice):
            end = min(len(self), key.stop or len(self))
            end = max(0, end)
            start = max(0, key.start or 0)
            end = min(end, start + len(value))

            for idx in xrange(start, end):
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
                raise ValueError('Input tuple or list required to have %s '
                                 'elements.' % self._nfields)
        else:
            raise TypeError('Assignment requires a FITS_record, tuple, or '
                            'list as input.')

    def __getslice__(self, start, end):
        return self[slice(start, end)]

    def __setslice__(self, start, end, value):
        self[slice(start, end)] = value

    def copy(self, order='C'):
        """
        The Numpy documentation lies; `numpy.ndarray.copy` is not equivalent to
        `numpy.copy`.  Differences include that it re-views the copied array as
        self's ndarray subclass, as though it were taking a slice; this means
        ``__array_finalize__`` is called and the copy shares all the array
        attributes (including ``._convert``!).  So we need to make a deep copy
        of all those attributes so that the two arrays truly do not share any
        data.
        """

        new = super(FITS_rec, self).copy(order=order)

        new.__dict__ = copy.deepcopy(self.__dict__)
        return new

    @property
    def columns(self):
        """
        A user-visible accessor for the coldefs.

        See https://aeon.stsci.edu/ssb/trac/pyfits/ticket/44
        """

        return self._coldefs

    @property
    def names(self):
        """List of column names."""

        if hasattr(self, '_coldefs') and self._coldefs is not None:
            return self._coldefs.names
        else:
            return list(self.dtype.names)


    def field(self, key):
        """
        A view of a `Column`'s data as an array.
        """

        # NOTE: The *column* index may not be the same as the field index in
        # the recarray, if the column is a phantom column
        col_indx = _get_index(self.columns.names, key)
        if self.columns[col_indx]._phantom:
            warnings.warn(
                'Field %r has a repeat count of 0 in its format code, '
                'indicating an empty field.' % key)
            recformat = self.columns._recformats[col_indx].lstrip('0')
            return np.array([], dtype=recformat)
        # Ignore phantom columns in determining the physical field number
        n_phantom = len([c for c in self.columns[:col_indx] if c._phantom])
        field_indx = col_indx - n_phantom

        recformat = self._coldefs._recformats[col_indx]

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
        field = np.recarray.field(base, field_indx)

        if self._convert[field_indx] is None:
            if isinstance(recformat, _FormatP):
                # for P format
                converted = self._convert_p(col_indx, field, recformat)
            else:
                # Handle all other column data types which are fixed-width
                # fields
                converted = self._convert_other(col_indx, field, recformat)

            self._convert[field_indx] = converted
            return converted

        return self._convert[field_indx]

    def _convert_x(self, field, recformat):
        """Convert a raw table column to a bit array as specified by the
        FITS X format.
        """

        dummy = np.zeros(self.shape + (recformat.repeat,), dtype=np.bool_)
        _unwrapx(field, dummy, recformat.repeat)
        return dummy

    def _convert_p(self, indx, field, recformat):
        """Convert a raw table column of FITS P or Q format descriptors
        to a VLA column with the array data returned from the heap.
        """

        dummy = _VLF([None] * len(self), dtype=recformat.dtype)
        raw_data = self._get_raw_data()

        if raw_data is None:
            raise IOError(
                "Could not find heap data for the %r variable-length "
                "array column." % self.columns.names[indx])

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
                dummy[idx] = self._convert_other(indx, dummy[idx], recformat)

        return dummy

    def _convert_ascii(self, indx, field):
        """Special handling for ASCII table columns to convert columns
        containing numeric types to actual numeric arrays from the string
        representation.
        """

        format = self._coldefs.formats[indx]
        recformat = ASCII2NUMPY[format[0]]
        # if the string = TNULL, return ASCIITNULL
        nullval = str(self._coldefs.nulls[indx]).strip().encode('ascii')
        if len(nullval) > format.width:
            nullval = nullval[:format.width]

        # Before using .replace make sure that any trailing bytes in each
        # column are filled with spaces, and *not*, say, nulls; this causes
        # functions like replace to potentially leave gibberish bytes in the
        # array buffer.
        dummy = np.char.ljust(field, format.width)
        dummy = np.char.replace(dummy, encode_ascii('D'), encode_ascii('E'))
        null_fill = encode_ascii(str(ASCIITNULL).rjust(format.width))
        dummy = np.where(np.char.strip(dummy) == nullval, null_fill, dummy)

        try:
            dummy = np.array(dummy, dtype=recformat)
        except ValueError as exc:
            raise ValueError(
                '%s; the header may be missing the necessary TNULL%d '
                'keyword or the table contains invalid data' %
                (exc, indx + 1))

        return dummy

    def _convert_other(self, indx, field, recformat):
        """Perform conversions on any other fixed-width column data types.

        This may not perform any conversion at all if it's not necessary, in
        which case the original column array is returned.
        """

        if isinstance(recformat, _FormatX):
            # special handling for the X format
            return self._convert_x(field, recformat)

        (_str, _bool, _number, _scale, _zero, bscale, bzero, dim) = \
            self._get_scale_factors(indx)

        # ASCII table, convert strings to numbers
        # TODO:
        # For now, check that these are ASCII columns by checking the coldefs
        # type; in the future all columns (for binary tables, ASCII tables, or
        # otherwise) should "know" what type they are already and how to handle
        # converting their data from FITS format to native format and vice
        # versa...
        if not _str and isinstance(self._coldefs, _AsciiColDefs):
            field = self._convert_ascii(indx, field)

        # Test that the dimensions given in dim are sensible; otherwise
        # display a warning and ignore them
        if dim:
            # See if the dimensions already match, if not, make sure the
            # number items will fit in the specified dimensions
            if field.ndim > 1:
                actual_shape = field.shape[1:]
                if _str:
                    actual_shape = (field.itemsize,) + actual_shape
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
                else:
                    actual_nitems = field.shape[1]
                if nitems > actual_nitems:
                    warnings.warn(
                        'TDIM%d value %s does not fit with the size of '
                        'the array items (%d).  TDIM%d will be ignored.'
                        % (indx + 1, self._coldefs.dims[indx],
                           actual_nitems, indx + 1))
                    dim = None

        # further conversion for both ASCII and binary tables
        # For now we've made columns responsible for *knowing* whether their
        # data has been scaled, but we make the FITS_rec class responsible for
        # actually doing the scaling
        # TODO: This also needs to be fixed in the effort to make Columns
        # responsible for scaling their arrays to/from FITS native values
        column = self._coldefs[indx]
        if (_number and (_scale or _zero) and not column._physical_values):
            # This is to handle pseudo unsigned ints in table columns
            # TODO: For now this only really works correctly for binary tables
            # Should it work for ASCII tables as well?
            if self._uint:
                if bzero == 2**15 and 'I' in self._coldefs.formats[indx]:
                    field = np.array(field, dtype=np.uint16)
                elif bzero == 2**31 and 'J' in self._coldefs.formats[indx]:
                    field = np.array(field, dtype=np.uint32)
                elif bzero == 2**63 and 'K' in self._coldefs.formats[indx]:
                    field = np.array(field, dtype=np.uint64)
                    bzero64 = np.uint64(2 ** 63)
                else:
                    field = np.array(field, dtype=np.float64)
            else:
                field = np.array(field, dtype=np.float64)

            if _scale:
                np.multiply(field, bscale, field)
            if _zero:
                if self._uint and 'K' in self._coldefs.formats[indx]:
                    # There is a chance of overflow, so be careful
                    test_overflow = field.copy()
                    try:
                        test_overflow += bzero64
                    except OverflowError:
                        warnings.warn(
                            "Overflow detected while applying TZERO{0:d}. "
                            "Returning unscaled data.".format(indx + 1))
                    else:
                        field = test_overflow
                else:
                    field += bzero
        elif _bool and field.dtype != bool:
            field = np.equal(field, ord('T'))
        elif _str:
            with ignored(UnicodeDecodeError):
                field = decode_ascii(field)

        if dim:
            # Apply the new field item dimensions
            nitems = reduce(operator.mul, dim)
            if field.ndim > 1:
                field = field[:, :nitems]
            if _str:
                fmt = field.dtype.char
                dtype = ('|%s%d' % (fmt, dim[-1]), dim[:-1])
                field.dtype = dtype
            else:
                field.shape = (field.shape[0],) + dim

        return field

    def _clone(self, shape):
        """
        Overload this to make mask array indexing work properly.
        """

        from .hdu.table import new_table

        hdu = new_table(self._coldefs, nrows=shape[0])
        return hdu.data

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

    def _get_scale_factors(self, indx):
        """
        Get the scaling flags and factors for one field.

        `indx` is the index of the field.
        """

        if isinstance(self._coldefs, _AsciiColDefs):
            _str = self._coldefs.formats[indx][0] == 'A'
            _bool = False  # there is no boolean in ASCII table
        else:
            _str = 'a' in self._coldefs._recformats[indx]
            # TODO: Determine a better way to determine if the column is bool
            # formatted
            _bool = self._coldefs._recformats[indx][-2:] == FITS2NUMPY['L']

        _number = not (_bool or _str)
        bscale = self._coldefs.bscales[indx]
        bzero = self._coldefs.bzeros[indx]
        _scale = bscale not in ('', None, 1)
        _zero = bzero not in ('', None, 0)
        # ensure bscale/bzero are numbers
        if not _scale:
            bscale = 1
        if not _zero:
            bzero = 0

        dim = self._coldefs._dims[indx]

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

        for indx in range(len(self.dtype.names)):
            recformat = self._coldefs._recformats[indx]
            field = super(FITS_rec, self).field(indx)

            # add the location offset of the heap area for each
            # variable length column
            if isinstance(recformat, _FormatP):
                # Irritatingly, this can return a different dtype than just
                # doing np.dtype(recformat.dtype); but this returns the results
                # that we want.  For example if recformat.dtype is 'a' we want
                # an array of characters.
                dtype = np.array([], dtype=recformat.dtype).dtype

                if update_heap_pointers and self._convert[indx] is not None:
                    # The VLA has potentially been updated, so we need to
                    # update the array descriptors
                    field[:] = 0  # reset
                    npts = [len(arr) for arr in self._convert[indx]]

                    field[:len(npts), 0] = npts
                    field[1:, 1] = (np.add.accumulate(field[:-1, 0]) *
                                    dtype.itemsize)
                    field[:, 1][:] += heapsize

                heapsize += field[:, 0].sum() * dtype.itemsize
                # Even if this VLA has not been read or updated, we need to
                # include the size of its constituent arrays in the heap size
                # total

            if self._convert[indx] is None:
                continue

            if isinstance(recformat, _FormatX):
                _wrapx(self._convert[indx], field, recformat.repeat)
                continue

            _str, _bool, _number, _scale, _zero, bscale, bzero, _ = \
                self._get_scale_factors(indx)

            # conversion for both ASCII and binary tables
            if _number or _str:
                column = self._coldefs[indx]
                if _number and (_scale or _zero) and column._physical_values:
                    dummy = self._convert[indx].copy()
                    if _zero:
                        dummy -= bzero
                    if _scale:
                        dummy /= bscale
                    # This will set the raw values in the recarray back to
                    # their non-physical storage values, so the column should
                    # be mark is not scaled
                    column._physical_values = False
                elif _str:
                    dummy = self._convert[indx]
                elif isinstance(self._coldefs, _AsciiColDefs):
                    dummy = self._convert[indx]
                else:
                    continue

                # ASCII table, convert numbers to strings
                if isinstance(self._coldefs, _AsciiColDefs):
                    starts = self._coldefs.starts[:]
                    spans = self._coldefs.spans
                    format = self._coldefs.formats[indx].strip()

                    # The the index of the "end" column of the record, beyond
                    # which we can't write
                    end = super(FITS_rec, self).field(-1).itemsize
                    starts.append(end + starts[-1])

                    if indx > 0:
                        lead = (starts[indx] - starts[indx - 1] -
                                spans[indx - 1])
                    else:
                        lead = 0

                    if lead < 0:
                        warnings.warn(
                            'Column %r starting point overlaps the '
                            'previous column.' % (indx + 1))

                    trail = starts[indx + 1] - starts[indx] - spans[indx]

                    if trail < 0:
                        warnings.warn(
                            'Column %r ending point overlaps the next '
                            'column.' % (indx + 1))

                    # TODO: It would be nice if these string column formatting
                    # details were left to a specialized class, as is the case
                    # with FormatX and FormatP
                    if 'A' in format:
                        _pc = '%-'
                    else:
                        _pc = '%'

                    fmt = ''.join([_pc, format[1:], ASCII2STR[format[0]],
                                   (' ' * trail)])

                    # not using numarray.strings's num2char because the
                    # result is not allowed to expand (as C/Python does).
                    for jdx in xrange(len(dummy)):
                        x = fmt % dummy[jdx]
                        if len(x) > starts[indx + 1] - starts[indx]:
                            raise ValueError(
                                "Value %r does not fit into the output's "
                                "itemsize of %s." % (x, spans[indx]))
                        else:
                            field[jdx] = x
                    # Replace exponent separator in floating point numbers
                    if 'D' in format:
                        field.replace(encode_ascii('E'), encode_ascii('D'))
                # binary table
                else:
                    if len(field) and isinstance(field[0], np.integer):
                        dummy = np.around(dummy)
                    elif isinstance(field, np.chararray):
                        # Ensure that blanks at the end of each string are
                        # converted to nulls instead of spaces, see Trac #15
                        # and #111
                        itemsize = dummy.itemsize
                        if dummy.dtype.kind == 'U':
                            pad = self._coldefs._padding_byte
                        else:
                            pad = self._coldefs._padding_byte.encode('ascii')

                        for idx in range(len(dummy)):
                            val = dummy[idx]
                            dummy[idx] = val + (pad * (itemsize - len(val)))

                        # Encode *after* handling the padding byte or else
                        # Numpy will complain about trying to append bytes to
                        # an array
                        if dummy.dtype.kind == 'U':
                            dummy = dummy.encode('ascii')

                    if field.shape == dummy.shape:
                        field[:] = dummy
                    else:
                        # Reshaping the data is necessary in cases where the
                        # TDIMn keyword was used to shape a column's entries
                        # into arrays
                        field[:] = dummy.ravel().view(field.dtype)

                del dummy

            # ASCII table does not have Boolean type
            elif _bool:
                field[:] = np.choose(self._convert[indx],
                                     (np.array([ord('F')], dtype=np.int8)[0],
                                      np.array([ord('T')], dtype=np.int8)[0]))

        # Store the updated heapsize
        self._heapsize = heapsize
