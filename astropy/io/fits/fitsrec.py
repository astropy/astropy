# Licensed under a 3-clause BSD style license - see PYFITS.rst

import operator
import sys
import warnings
import weakref

import numpy as np

from .column import (ASCIITNULL, FITS2NUMPY, TDIM_RE, Column, ColDefs,
                     _FormatX, _FormatP, _VLF, _get_index, _wrapx, _unwrapx,
                     _convert_format, _convert_ascii_format)
from .util import _array_from_file, decode_ascii

from ...utils import lazyproperty


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
           Used for subsetting the columns of the FITS_rec object.

        end : int, optional
           The ending column in the row associated with this object.
           Used for subsetting the columns of the FITS_rec object.
        """

        # For backward compatibility...
        for arg in [('startColumn', 'start'), ('endColumn', 'end')]:
            if arg[0] in kwargs:
                warnings.warn('The %s argument to FITS_record is deprecated; '
                              'use %s instead' % arg, DeprecationWarning)
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
        if isinstance(key, basestring):
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
        if isinstance(key, basestring):
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

    `FITS_rec` is the data part of a table HDU's data part.  This is a
    layer over the `recarray`, so we can deal with scaled columns.

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
        self._file = None
        self._buffer = None
        self._coldefs = None
        self._gap = 0
        self.names = list(self.dtype.names)
        self.formats = None
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return

        if isinstance(obj, FITS_rec):
            self._convert = obj._convert
            self._heapoffset = obj._heapoffset
            self._heapsize = obj._heapsize
            self._file = obj._file
            self._buffer = obj._buffer
            self._coldefs = obj._coldefs
            self._nfields = obj._nfields
            self._gap = obj._gap
            self.names = obj.names
            self.formats = obj.formats
        else:
            # This will allow regular ndarrays with fields, rather than
            # just other FITS_rec objects
            self._nfields = len(obj.dtype.names)
            self._convert = [None] * len(obj.dtype.names)

            self._heapoffset = getattr(obj, '_heapoffset', 0)
            self._heapsize = getattr(obj, '_heapsize', 0)
            self._file = getattr(obj, '_file', None)
            self._buffer = getattr(obj, '_buffer', None)

            self._coldefs = None
            self._gap = 0

            # Bypass setattr-based assignment to fields; see #86
            self.names = list(obj.dtype.names)
            self.formats = None

            attrs = ['_convert', '_coldefs', '_gap']
            for attr in attrs:
                if hasattr(obj, attr):
                    value = getattr(obj, attr, None)
                    if value is None:
                        warnings.warn('Setting attribute %s as None' % attr)
                    setattr(self, attr, value)

            if self._coldefs is None:
                self._coldefs = ColDefs(self)
            self.formats = self._coldefs.formats

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
        if isinstance(key, basestring):
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

    def __setitem__(self, row, value):
        if isinstance(row, slice):
            end = min(len(self), row.stop or len(self))
            end = max(0, end)
            start = max(0, row.start or 0)
            end = min(end, start + len(value))

            for idx in range(start, end):
                self.__setitem__(idx, value[idx - start])
            return

        if isinstance(value, FITS_record):
            for idx in range(self._nfields):
                self.field(self.names[idx])[row] = value.field(self.names[idx])
        elif isinstance(value, (tuple, list)):
            if self._nfields == len(value):
                for idx in range(self._nfields):
                    self.field(idx)[row] = value[idx]
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

    @property
    def columns(self):
        """
        A user-visible accessor for the coldefs.  See ticket #44.
        """

        return self._coldefs

    def field(self, key):
        """
        A view of a `Column`'s data as an array.
        """

        indx = _get_index(self.names, key)
        recformat = self._coldefs._recformats[indx]

        # If field's base is a FITS_rec, we can run into trouble because it
        # contains a reference to the ._coldefs object of the original data;
        # this can lead to a circular reference; see ticket #49
        base = self
        while isinstance(base, FITS_rec) and \
              isinstance(base.base, np.recarray):
            base = base.base
        # base could still be a FITS_rec in some cases, so take care to
        # use rec.recarray.field to avoid a potential infinite
        # recursion
        field = np.recarray.field(base, indx)

        if self._convert[indx] is None:
            # for X format
            if isinstance(recformat, _FormatX):
                _nx = recformat._nx
                dummy = np.zeros(self.shape + (_nx,), dtype=np.bool_)
                _unwrapx(field, dummy, _nx)
                self._convert[indx] = dummy
                return self._convert[indx]

            (_str, _bool, _number, _scale, _zero, bscale, bzero, dim) = \
                self._get_scale_factors(indx)

            # for P format
            buff = None
            if isinstance(recformat, _FormatP):
                dummy = _VLF([None] * len(self), dtype=recformat.dtype)
                for i in range(len(self)):
                    _offset = field[i, 1] + self._heapoffset

                    if self._file is not None:
                        self._file.seek(_offset)
                        def get_pdata(dtype, count):
                            return _array_from_file(self._file, dtype=dtype,
                                                    count=count, sep='')
                    else:  # There must be a _buffer or something is wrong
                        # Sometimes the buffer is already a Numpy array; in
                        # particular this can occur in compressed HDUs.
                        # Hypothetically other cases as well.
                        if buff is None:
                            buff = self._buffer
                        if not isinstance(buff, np.ndarray):
                            # Go ahead and great a single ndarray from the
                            # buffer if it is not already one; we will then
                            # take slices from it.  This is more efficient than
                            # the previous approach that created separate
                            # arrays for each VLA.
                            buff = np.fromstring(buff, dtype=np.uint8)

                        def get_pdata(dtype, count):
                            dtype = np.dtype(dtype)
                            nbytes = count * dtype.itemsize
                            slc = slice(_offset, _offset + nbytes)
                            return buff[slc].view(dtype=dtype)

                    if recformat.dtype == 'a':
                        count = field[i, 0]
                        dt = recformat.dtype + str(1)
                        da = get_pdata(dt, count)
                        dummy[i] = np.char.array(da, itemsize=count)
                        dummy[i] = decode_ascii(dummy[i])
                    else:
                        count = field[i, 0]
                        dt = recformat.dtype
                        dummy[i] = get_pdata(dt, count)
                        dummy[i].dtype = dummy[i].dtype.newbyteorder('>')

                # scale by TSCAL and TZERO
                if _scale or _zero:
                    for i in range(len(self)):
                        dummy[i][:] = dummy[i] * bscale + bzero

                # Boolean (logical) column
                if recformat.dtype == FITS2NUMPY['L']:
                    for i in range(len(self)):
                        dummy[i] = np.equal(dummy[i], ord('T'))

                self._convert[indx] = dummy
                return self._convert[indx]

            # ASCII table, convert strings to numbers
            if not _str and self._coldefs._tbtype == 'TableHDU':
                _fmap = {'I': np.int32, 'F': np.float32, 'E': np.float32,
                         'D': np.float64}
                _type = _fmap[self._coldefs.formats[indx][0]]

                # if the string = TNULL, return ASCIITNULL
                nullval = self._coldefs.nulls[indx].strip().encode('ascii')
                dummy = field.replace('D'.encode('ascii'),
                                      'E'.encode('ascii'))
                dummy = np.where(dummy.strip() == nullval, str(ASCIITNULL),
                                 dummy)
                dummy = np.array(dummy, dtype=_type)

                self._convert[indx] = dummy
            else:
                dummy = field

            # Test that the dimensions given in dim are sensible; otherwise
            # display a warning and ignore them
            if dim:
                # See if the dimensions already match, if not, make sure the
                # number items will fit in the specified dimensions
                if dummy.ndim > 1:
                    actual_shape = dummy[0].shape
                    if _str:
                        actual_shape = (dummy[0].itemsize,) + actual_shape
                else:
                    actual_shape = len(dummy[0])
                if dim == actual_shape:
                    # The array already has the correct dimensions, so we
                    # ignore dim and don't convert
                    dim = None
                else:
                    nitems = reduce(operator.mul, dim)
                    if _str:
                        actual_nitems = dummy.itemsize
                    else:
                        actual_nitems = dummy.shape[1]
                    if nitems != actual_nitems:
                        warnings.warn(
                        'TDIM%d value %s does not fit with the size of '
                            'the array items (%d).  TDIM%d will be ignored.'
                            % (indx + 1, self._coldefs.dims[indx],
                               actual_nitems, indx + 1))
                        dim = None

            # further conversion for both ASCII and binary tables
            if _number and (_scale or _zero):

                # only do the scaling the first time and store it in _convert
                self._convert[indx] = np.array(dummy, dtype=np.float64)
                if _scale:
                    np.multiply(self._convert[indx], bscale,
                                self._convert[indx])
                if _zero:
                    self._convert[indx] += bzero
            elif _bool:
                self._convert[indx] = np.equal(dummy, ord('T'))
            elif _str:
                try:
                    self._convert[indx] = decode_ascii(dummy)
                except UnicodeDecodeError:
                    pass

            if dim:
                if self._convert[indx] is None:
                    self._convert[indx] = dummy
                if _str:
                    fmt = self._convert[indx].dtype.char
                    dtype = ('|%s%d' % (fmt, dim[-1]), dim[:-1])
                    self._convert[indx].dtype = dtype
                else:
                    self._convert[indx].shape = (dummy.shape[0],) + dim

        if self._convert[indx] is not None:
            return self._convert[indx]
        else:
            return dummy

    def _clone(self, shape):
        """
        Overload this to make mask array indexing work properly.
        """

        from .hdu.table import new_table

        hdu = new_table(self._coldefs, nrows=shape[0])
        return hdu.data

    def _get_scale_factors(self, indx):
        """
        Get the scaling flags and factors for one field.

        `indx` is the index of the field.
        """

        if self._coldefs._tbtype == 'BinTableHDU':
            _str = 'a' in self._coldefs._recformats[indx]
            _bool = self._coldefs._recformats[indx][-2:] == FITS2NUMPY['L']
        else:
            _str = self._coldefs.formats[indx][0] == 'A'
            _bool = False             # there is no boolean in ASCII table
        _number = not(_bool or _str)
        bscale = self._coldefs.bscales[indx]
        bzero = self._coldefs.bzeros[indx]
        _scale = bscale not in ['', None, 1]
        _zero = bzero not in ['', None, 0]
        # ensure bscale/bzero are numbers
        if not _scale:
            bscale = 1
        if not _zero:
            bzero = 0
        dim = self._coldefs.dims[indx]
        m = dim and TDIM_RE.match(dim)
        if m:
            dim = m.group('dims')
            dim = tuple(int(d.strip()) for d in dim.split(','))[::-1]
        else:
            # Ignore any dim values that don't specify a multidimensional
            # column
            dim = ''

        return (_str, _bool, _number, _scale, _zero, bscale, bzero, dim)

    def _scale_back(self):
        """
        Update the parent array, using the (latest) scaled array.
        """

        _fmap = {'A': 's', 'I': 'd', 'J': 'd', 'F': 'f', 'E': 'E', 'D': 'E'}
        # calculate the starting point and width of each field for ASCII table
        # TODO: Ick--fix this _tbtype usage eventually...
        if self._coldefs._tbtype == 'TableHDU':
            loc = self._coldefs.starts
            widths = []

            idx = 0
            for idx in range(len(self.dtype.names)):
                f = _convert_ascii_format(self._coldefs.formats[idx])
                widths.append(f[1])
            loc.append(loc[-1] + super(FITS_rec, self).field(idx).itemsize)

        for indx in range(len(self.dtype.names)):
            recformat = self._coldefs._recformats[indx]
            field = super(FITS_rec, self).field(indx)

            if self._convert[indx] is None:
                continue

            if isinstance(recformat, _FormatX):
                _wrapx(self._convert[indx], field, recformat._nx)
                continue

            (_str, _bool, _number, _scale, _zero, bscale, bzero, dim) = \
                self._get_scale_factors(indx)

            # add the location offset of the heap area for each
            # variable length column
            if isinstance(recformat, _FormatP):
                # Reset the heapsize and recompute it starting from the first P
                # column
                if indx == 0:
                    self._heapsize = 0

                field[:] = 0  # reset
                npts = map(len, self._convert[indx])

                # Irritatingly, this can return a different dtype than just
                # doing np.dtype(recformat.dtype); but this returns the results
                # that we want.  For example if recformat.dtype is 'a' we want
                # an array of characters.
                dtype = np.array([], dtype=recformat.dtype).dtype
                field[:len(npts), 0] = npts
                field[1:, 1] = (np.add.accumulate(field[:-1, 0]) *
                                dtype.itemsize)
                field[:, 1][:] += self._heapsize
                self._heapsize += field[:, 0].sum() * dtype.itemsize

            # conversion for both ASCII and binary tables
            if _number or _str:
                if _number and (_scale or _zero):
                    dummy = self._convert[indx].copy()
                    if _zero:
                        dummy -= bzero
                    if _scale:
                        dummy /= bscale
                elif _str:
                    dummy = self._convert[indx]
                elif self._coldefs._tbtype == 'TableHDU':
                    dummy = self._convert[indx]
                else:
                    continue

                # ASCII table, convert numbers to strings
                if self._coldefs._tbtype == 'TableHDU':
                    format = self._coldefs.formats[indx].strip()
                    lead = self._coldefs.starts[indx] - loc[indx]
                    if lead < 0:
                        raise ValueError(
                            'Column `%s` starting point overlaps to the '
                            'previous column.' % indx + 1)
                    trail = (loc[indx + 1] - widths[indx] -
                             self._coldefs.starts[indx])
                    if trail < 0:
                        raise ValueError(
                            'Column `%s` ending point overlaps to the next '
                            'column.' % indx + 1)
                    if 'A' in format:
                        _pc = '%-'
                    else:
                        _pc = '%'

                    fmt = ''.join([(' ' * lead), _pc, format[1:],
                                   _fmap[format[0]], (' ' * trail)])

                    # not using numarray.strings's num2char because the
                    # result is not allowed to expand (as C/Python does).
                    for jdx in range(len(dummy)):
                        x = fmt % dummy[jdx]
                        if len(x) > (loc[indx + 1] - loc[indx]):
                            raise ValueError(
                                "Number `%s` does not fit into the output's "
                                "itemsize of %s." % (x, widths[indx]))
                        else:
                            field[jdx] = x
                    # Replace exponent separator in floating point numbers
                    if 'D' in format:
                        field.replace('E', 'D')
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

                        for idx in xrange(len(dummy)):
                            val = dummy[idx]
                            dummy[idx] = val + (pad * (itemsize - len(val)))
                    if dummy.dtype != field.dtype:
                        dummy = dummy.astype(field.dtype)
                    field[:] = dummy

                del dummy

            # ASCII table does not have Boolean type
            elif _bool:
                field[:] = np.choose(self._convert[indx],
                                     (np.array([ord('F')], dtype=np.int8)[0],
                                      np.array([ord('T')], dtype=np.int8)[0]))
