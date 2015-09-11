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
                     _wrapx, _unwrapx, _makep, Delayed)
from .util import decode_ascii, encode_ascii
from ...extern.six import string_types
from ...extern.six.moves import xrange
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
        self._converted = {}
        self._heapoffset = 0
        self._heapsize = 0
        self._coldefs = None
        self._gap = 0
        self._uint = False
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

        for attrs in ['_converted', '_heapoffset', '_heapsize', '_nfields',
                      '_gap', '_uint', 'parnames', '_coldefs']:

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
            self._converted = obj._converted
            self._heapoffset = obj._heapoffset
            self._heapsize = obj._heapsize
            self._coldefs = obj._coldefs
            self._nfields = obj._nfields
            self._gap = obj._gap
            self._uint = obj._uint
        else:
            # This will allow regular ndarrays with fields, rather than
            # just other FITS_rec objects
            self._nfields = len(obj.dtype.names)
            self._converted = {}

            self._heapoffset = getattr(obj, '_heapoffset', 0)
            self._heapsize = getattr(obj, '_heapsize', 0)

            self._coldefs = None
            self._gap = getattr(obj, '_gap', 0)
            self._uint = getattr(obj, '_uint', False)

            attrs = ['_converted', '_coldefs', '_gap']
            for attr in attrs:
                if hasattr(obj, attr):
                    value = getattr(obj, attr, None)
                    if value is None:
                        warnings.warn('Setting attribute %s as None' % attr,
                                      AstropyUserWarning)
                    setattr(self, attr, value)

            if self._coldefs is None:
                self._coldefs = ColDefs(self)

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

        # Make sure the data is a listener for changes to the columns
        columns._add_listener(data)

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
                data._converted[name] = _makep(inarr, field, recformat,
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
                data._converted[name] = np.zeros(field.shape, dtype=bool)
                data._converted[name][:n] = inarr
                # TODO: Maybe this step isn't necessary at all if _scale_back
                # will handle it?
                inarr = np.where(inarr == False, ord('F'), ord('T'))
            elif (columns[idx]._physical_values and
                    columns[idx]._pseudo_unsigned_ints):
                # Temporary hack...
                bzero = column.bzero
                data._converted[name] = np.zeros(field.shape,
                                                 dtype=inarr.dtype)
                data._converted[name][:n] = inarr
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
                    data._converted[name] = np.zeros(nrows, dtype=arr.dtype)
                    outarr = data._converted[name][:n]

                outarr[:] = inarr
                continue

            if inarr.shape != outarr.shape:
                if (inarr.dtype.kind == outarr.dtype.kind and
                        inarr.dtype.kind in ('U', 'S') and
                        inarr.dtype != outarr.dtype):
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
            out._converted = {}
            for idx, name in enumerate(self._coldefs.names):
                #
                # Store the new arrays for the _coldefs object
                #
                arrays.append(self._coldefs._arrays[idx][key])

                # touch all fields to expand the original ._converted dict
                # so the sliced FITS_rec will view the same scaled columns as
                # the original
                dummy = self.field(idx)
                if name in self._converted:
                    out._converted[name] = \
                        np.ndarray.__getitem__(self._converted[name], key)
            del dummy

            out._coldefs._arrays = arrays
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
        attributes (including ``._converted``!).  So we need to make a deep
        copy of all those attributes so that the two arrays truly do not share
        any data.
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

    @property
    def formats(self):
        """List of column FITS foramts."""

        return self._coldefs.formats

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
                'Field %r has a repeat count of 0 in its format code, '
                'indicating an empty field.' % key)
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

            self._converted[name] = converted
            return converted

        return self._converted[name]

    def _update_column_attribute_changed(self, column, idx, attr, old_value,
                                         new_value):
        """
        Update how the data is formatted depending on changes to column
        attributes initiated by the user through the `Column` interface.

        Dispatches column attribute change notifications to individual methods
        for each attribute ``_update_column_<attr>``
        """

        method_name = '_update_column_{0}'.format(attr)
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
            raise IOError(
                "Could not find heap data for the %r variable-length "
                "array column." % column.name)

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
        recformat = ASCII2NUMPY[format[0]]
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
        dummy = np.where(np.char.strip(dummy) == nullval, null_fill, dummy)

        try:
            dummy = np.array(dummy, dtype=recformat)
        except ValueError as exc:
            indx = self._coldefs.names.index(column.name)
            raise ValueError(
                '%s; the header may be missing the necessary TNULL%d '
                'keyword or the table contains invalid data' %
                (exc, indx + 1))

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

        indx = self._coldefs.names.index(column.name)

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
                elif len(field.shape) == 1:  # No repeat count in TFORMn, equivalent to 1
                    actual_nitems = 1
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
            field = _get_recarray_field(self, indx)

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
                    field[:] = 0  # reset
                    npts = [len(arr) for arr in self._converted[name]]

                    field[:len(npts), 0] = npts
                    field[1:, 1] = (np.add.accumulate(field[:-1, 0]) *
                                    dtype.itemsize)
                    field[:, 1][:] += heapsize

                heapsize += field[:, 0].sum() * dtype.itemsize
                # Even if this VLA has not been read or updated, we need to
                # include the size of its constituent arrays in the heap size
                # total

            if name not in self._converted:
                continue

            if isinstance(recformat, _FormatX):
                _wrapx(self._converted[name], field, recformat.repeat)
                continue

            _str, _bool, _number, _scale, _zero, bscale, bzero, _ = \
                self._get_scale_factors(column)

            # conversion for both ASCII and binary tables
            if _number or _str:
                if _number and (_scale or _zero) and column._physical_values:
                    dummy = self._converted[name].copy()
                    if _zero:
                        dummy -= bzero
                    if _scale:
                        dummy /= bscale
                    # This will set the raw values in the recarray back to
                    # their non-physical storage values, so the column should
                    # be mark is not scaled
                    column._physical_values = False
                elif _str:
                    dummy = self._converted[name]
                elif isinstance(self._coldefs, _AsciiColDefs):
                    dummy = self._converted[name]
                else:
                    continue

                # ASCII table, convert numbers to strings
                if isinstance(self._coldefs, _AsciiColDefs):
                    self._scale_back_ascii(indx, dummy, field)
                # binary table
                else:
                    if len(field) and isinstance(field[0], np.integer):
                        dummy = np.around(dummy)
                    elif isinstance(field, np.chararray):
                        # Ensure that blanks at the end of each string are
                        # converted to nulls instead of spaces, see Trac #15
                        # and #111
                        _rstrip_inplace(dummy)

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
                field[:] = np.choose(self._converted[name],
                                     (np.array([ord('F')], dtype=np.int8)[0],
                                      np.array([ord('T')], dtype=np.int8)[0]))

        # Store the updated heapsize
        self._heapsize = heapsize

    def _scale_back_ascii(self, col_idx, input_field, output_field):
        """
        Convert internal array values back to ASCII table representation.

        The ``input_field`` is the internal representation of the values, and
        the ``output_field`` is the character array representing the ASCII
        output that will be written.
        """

        starts = self._coldefs.starts[:]
        spans = self._coldefs.spans
        format = self._coldefs.formats[col_idx]

        # The the index of the "end" column of the record, beyond
        # which we can't write
        end = super(FITS_rec, self).field(-1).itemsize
        starts.append(end + starts[-1])

        if col_idx > 0:
            lead = starts[col_idx] - starts[col_idx - 1] - spans[col_idx - 1]
        else:
            lead = 0

        if lead < 0:
            warnings.warn('Column %r starting point overlaps the previous '
                          'column.' % (col_idx + 1))

        trail = starts[col_idx + 1] - starts[col_idx] - spans[col_idx]

        if trail < 0:
            warnings.warn('Column %r ending point overlaps the next '
                          'column.' % (col_idx + 1))

        # TODO: It would be nice if these string column formatting
        # details were left to a specialized class, as is the case
        # with FormatX and FormatP
        if 'A' in format:
            _pc = '%-'
        else:
            _pc = '%'

        fmt = ''.join([_pc, format[1:], ASCII2STR[format[0]],
                       (' ' * trail)])

        # Even if the format precision is 0, we should output a decimal point
        # as long as there is space to do so--not including a decimal point in
        # a float value is discouraged by the FITS Standard
        trailing_decimal = (format.precision == 0 and
                            format.format in ('F', 'E', 'D'))

        # not using numarray.strings's num2char because the
        # result is not allowed to expand (as C/Python does).
        for jdx, value in enumerate(input_field):
            value = fmt % value
            if len(value) > starts[col_idx + 1] - starts[col_idx]:
                raise ValueError(
                    "Value %r does not fit into the output's itemsize of "
                    "%s." % (value, spans[col_idx]))

            if trailing_decimal and value[0] == ' ':
                # We have some extra space in the field for the trailing
                # decimal point
                value = value[1:] + '.'

            output_field[jdx] = value

        # Replace exponent separator in floating point numbers
        if 'D' in format:
            output_field.replace(encode_ascii('E'), encode_ascii('D'))


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
    if field.dtype.char in ('S', 'U') and not isinstance(field, np.chararray):
        field = field.view(np.chararray)
    return field


def _rstrip_inplace(array, chars=None):
    """
    Performs an in-place rstrip operation on string arrays.
    This is necessary since the built-in `np.char.rstrip` in Numpy does not
    perform an in-place calculation.  This can be removed if ever
    https://github.com/numpy/numpy/issues/6303 is implemented (however, for
    the purposes of this module the only in-place vectorized string function
    we need is rstrip).
    """

    for item in np.nditer(array, flags=['zerosize_ok'],
                                 op_flags=['readwrite']):
        item[...] = item.item().rstrip(chars)
