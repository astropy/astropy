import re
import sys
import warnings
import weakref

import numpy as np
from numpy import char as chararray

from pyfits.card import Card
from pyfits.util import lazyproperty, pairwise, _is_int, _convert_array, \
                        encode_ascii, deprecated


__all__ = ['Column', 'ColDefs', 'Delayed']


# mapping from TFORM data type to numpy data type (code)
# L: Logical (Boolean)
# B: Unsigned Byte
# I: 16-bit Integer
# J: 32-bit Integer
# K: 64-bit Integer
# E: Single-precision Floating Point
# D: Double-precision Floating Point
# C: Single-precision Complex
# M: Double-precision Complex
# A: Character
FITS2NUMPY = {'L': 'i1', 'B': 'u1', 'I': 'i2', 'J': 'i4', 'K': 'i8', 'E': 'f4',
              'D': 'f8', 'C': 'c8', 'M': 'c16', 'A': 'a'}

# the inverse dictionary of the above
NUMPY2FITS = dict([(val, key) for key, val in FITS2NUMPY.iteritems()])

# This is the order in which values are converted to FITS types
# Note that only double precision floating point/complex are supported
FORMATORDER = ['L', 'B', 'I', 'J', 'K', 'D', 'M', 'A']

# lists of column/field definition common names and keyword names, make
# sure to preserve the one-to-one correspondence when updating the list(s).
# Use lists, instead of dictionaries so the names can be displayed in a
# preferred order.
KEYWORD_NAMES = ['TTYPE', 'TFORM', 'TUNIT', 'TNULL', 'TSCAL', 'TZERO',
                 'TDISP', 'TBCOL', 'TDIM']
KEYWORD_ATTRIBUTES = ['name', 'format', 'unit', 'null', 'bscale', 'bzero',
                      'disp', 'start', 'dim']

# TFORM regular expression
TFORMAT_RE = re.compile(r'(?P<repeat>^[0-9]*)(?P<dtype>[A-Za-z])'
                        r'(?P<option>[!-~]*)')

# table definition keyword regular expression
TDEF_RE = re.compile(r'(?P<label>^T[A-Z]*)(?P<num>[1-9][0-9 ]*$)')

# table dimension keyword regular expression (fairly flexible with whitespace)
TDIM_RE = re.compile(r'\(\s*(?P<dims>(?:\d+,\s*)+\s*\d+)\s*\)\s*')

ASCIITNULL = 0          # value for ASCII table cell with value = TNULL
                        # this can be reset by user.


class Delayed(object):
    """Delayed file-reading data."""

    def __init__(self, hdu=None, field=None):
        self.hdu = weakref.proxy(hdu)
        self.field = field

    def __getitem__(self, key):
        # This forces the data for the HDU to be read, which will replace
        # the corresponding Delayed objects in the Tables Columns to be
        # transformed into ndarrays.  It will also return the value of the
        # requested data element.
        return self.hdu.data[key][self.field]


class _FormatX(str):
    """For X format in binary tables."""

    pass


class _FormatP(str):
    """For P format in variable length table."""

    pass


class Column(object):
    """
    Class which contains the definition of one column, e.g.  `ttype`,
    `tform`, etc. and the array containing values for the column.
    Does not support `theap` yet.
    """

    def __init__(self, name=None, format=None, unit=None, null=None, \
                       bscale=None, bzero=None, disp=None, start=None, \
                       dim=None, array=None):
        """
        Construct a `Column` by specifying attributes.  All attributes
        except `format` can be optional.

        Parameters
        ----------
        name : str, optional
            column name, corresponding to ``TTYPE`` keyword

        format : str, optional
            column format, corresponding to ``TFORM`` keyword

        unit : str, optional
            column unit, corresponding to ``TUNIT`` keyword

        null : str, optional
            null value, corresponding to ``TNULL`` keyword

        bscale : int-like, optional
            bscale value, corresponding to ``TSCAL`` keyword

        bzero : int-like, optional
            bzero value, corresponding to ``TZERO`` keyword

        disp : str, optional
            display format, corresponding to ``TDISP`` keyword

        start : int, optional
            column starting position (ASCII table only), corresponding
            to ``TBCOL`` keyword

        dim : str, optional
            column dimension corresponding to ``TDIM`` keyword
        """

        # any of the input argument (except array) can be a Card or just
        # a number/string
        for attr in KEYWORD_ATTRIBUTES:
            value = locals()[attr]           # get the argument's value

            if isinstance(value, Card):
                setattr(self, attr, value.value)
            else:
                setattr(self, attr, value)

        # if the column data is not ndarray, make it to be one, i.e.
        # input arrays can be just list or tuple, not required to be ndarray
        if format is not None:
            # check format
            try:

                # legit FITS format? convert to record format (e.g. '3J'->'3i4')
                recfmt = _convert_format(format)
            except ValueError:
                try:
                    # legit recarray format?
                    recfmt = format
                    format = _convert_format(recfmt, reverse=True)
                except ValueError:
                    raise ValueError('Illegal format `%s`.' % format)

            self.format = format
            # Zero-length formats are legal in the FITS format, but since they
            # are not supported by numpy we mark columns that use them as
            # "phantom" columns, that are not considered when reading the data
            # as a record array.
            if self.format[0] == '0' or \
               (self.format[-1] == '0' and self.format[-2].isalpha()):
                self._phantom = True
            else:
                self._phantom = False

            # does not include Object array because there is no guarantee
            # the elements in the object array are consistent.
            if not isinstance(array,
                              (np.ndarray, chararray.chararray, Delayed)):
                try: # try to convert to a ndarray first
                    if array is not None:
                        array = np.array(array)
                except:
                    try: # then try to conver it to a strings array
                        array = chararray.array(array,
                                                itemsize=eval(recfmt[1:]))

                    # then try variable length array
                    except:
                        if isinstance(recfmt, _FormatP):
                            try:
                                array = _VLF(array)
                            except:
                                try:
                                    # this handles ['abc'] and [['a','b','c']]
                                    # equally, beautiful!
                                    _func = lambda x: \
                                                chararray.array(x, itemsize=1)
                                    array = _VLF(map(_func, array))
                                except:
                                    raise ValueError('Inconsistent input data '
                                                     'array: %s' % array)
                            array._dtype = recfmt._dtype
                        else:
                            raise ValueError('Data is inconsistent with the '
                                             'format `%s`.' % format)

        else:
            raise ValueError('Must specify format to construct Column.')

        # scale the array back to storage values if there is bscale/bzero
        if isinstance(array, np.ndarray):

            # boolean needs to be scaled too
            if recfmt[-2:] == FITS2NUMPY['L']:
                array = np.where(array==0, ord('F'), ord('T'))

            # make a copy if scaled, so as not to corrupt the original array
            if bzero not in ['', None, 0] or bscale not in ['', None, 1]:
                array = array.copy()
                if bzero not in ['', None, 0]:
                    array += -bzero
                if bscale not in ['', None, 1]:
                    array /= bscale

        array = self._convert_to_valid_data_type(array, self.format)
        self.array = array

    def __repr__(self):
        text = ''
        for attr in KEYWORD_ATTRIBUTES:
            value = getattr(self, attr)
            if value is not None:
                text += attr + ' = ' + repr(value) + '; '
        return text[:-2]

    def copy(self):
        """
        Return a copy of this `Column`.
        """
        tmp = Column(format='I') # just use a throw-away format
        tmp.__dict__ = self.__dict__.copy()
        return tmp

    def _convert_to_valid_data_type(self, array, format):
        # Convert the format to a type we understand
        if isinstance(array, Delayed):
            return array
        elif array is None:
            return array
        else:
            if 'A' in format and 'P' not in format:
                if array.dtype.char in 'SU':
                    fsize = int(_convert_format(format)[1:])
                    return chararray.array(array, itemsize=fsize)
                else:
                    numpy_format = _convert_format(format)
                    return _convert_array(array, np.dtype(numpy_format))
            elif 'X' not in format and 'P' not in format:
                (repeat, fmt, option) = _parse_tformat(format)
                # Preserve byte order of the original array for now; see #77
                numpy_format = array.dtype.byteorder + _convert_format(fmt)
                return _convert_array(array, np.dtype(numpy_format))
            elif 'X' in format:
                return _convert_array(array, np.dtype('uint8'))
            else:
                return array


class ColDefs(object):
    """
    Column definitions class.

    It has attributes corresponding to the `Column` attributes
    (e.g. `ColDefs` has the attribute `~ColDefs.names` while `Column`
    has `~Column.name`). Each attribute in `ColDefs` is a list of
    corresponding attribute values from all `Column` objects.
    """

    _padding_byte = '\x00'

    def __new__(cls, input, tbtype='BinTableHDU'):
        from pyfits.hdu.table import TableHDU

        if tbtype == 'BinTableHDU':
            klass = cls
        elif tbtype == 'TableHDU':
            klass = _ASCIIColDefs
        else:
            raise ValueError('Invalid table type: %s.' % tbtype)

        if isinstance(input, TableHDU):
            klass = _ASCIIColDefs

        return object.__new__(klass)

    def __init__(self, input, tbtype='BinTableHDU'):
        """
        Parameters
        ----------

        input :
            An existing table HDU, an existing ColDefs, or recarray

        **(Deprecated)** tbtype : str, optional
            which table HDU, ``"BinTableHDU"`` (default) or
            ``"TableHDU"`` (text table).
            Now ColDefs for a normal (binary) table by default, but converted
            automatically to ASCII table ColDefs in the appropriate contexts
            (namely, when creating an ASCII table).
        """

        from pyfits.hdu.table import _TableBaseHDU

        self._tbtype = tbtype

        if isinstance(input, ColDefs):
            self.columns = [col.copy() for col in input.columns]

        # if the input is a list of Columns
        elif isinstance(input, (list, tuple)):
            for col in input:
                if not isinstance(col, Column):
                    raise TypeError(
                           'Element %d in the ColDefs input is not a Column.'
                           % input.index(col))
            self.columns = [col.copy() for col in input]

        # Construct columns from the fields of a record array
        elif isinstance(input, np.ndarray) and input.dtype.fields is not None:
            self.columns = []
            for idx in range(len(input.dtype)):
                cname = input.dtype.names[idx]
                ftype = input.dtype.fields[cname][0]
                # String formats should have 'A' first
                if ftype.type == np.string_:
                    format = 'A' + str(ftype.itemsize)
                else:
                    format = _convert_format(ftype, reverse=True)
                # Determine the appropriate dimensions for items in the column
                # (typically just 1D)
                dim = input.dtype[idx].shape[::-1]
                if dim and (len(dim) > 1 or 'A' in format):
                    if 'A' in format:
                        # n x m string arrays must include the max string
                        # length in their dimensions (e.g. l x n x m)
                        dim = (input.dtype[idx].base.itemsize,) + dim
                    dim = repr(dim).replace(' ', '')
                else:
                    dim = None

                c = Column(name=cname, format=format,
                           array=input.view(np.ndarray)[cname], dim=dim)
                self.columns.append(c)

        # Construct columns from fields in an HDU header
        elif isinstance(input, _TableBaseHDU):
            hdr = input._header
            nfields = hdr['TFIELDS']
            self._width = hdr['NAXIS1']
            self._shape = hdr['NAXIS2']

            # go through header keywords to pick out column definition keywords
            # definition dictionaries for each field
            col_attributes = [{} for i in range(nfields)]
            for card in hdr.ascard:
                key = TDEF_RE.match(card.key)
                try:
                    keyword = key.group('label')
                except:
                    continue               # skip if there is no match
                if (keyword in KEYWORD_NAMES):
                    col = int(key.group('num'))
                    if col <= nfields and col > 0:
                        idx = KEYWORD_NAMES.index(keyword)
                        attr = KEYWORD_ATTRIBUTES[idx]
                        col_attributes[col - 1][attr] = card.value

            # data reading will be delayed
            for col in range(nfields):
                col_attributes[col]['array'] = Delayed(input, col)

            # now build the columns
            self.columns = [Column(**attrs) for attrs in col_attributes]
            self._listener = input
        else:
            raise TypeError('Input to ColDefs must be a table HDU or a list '
                            'of Columns.')

        # For ASCII tables, reconstruct string columns and ensure that spaces
        # are used for padding instead of \x00, and do the reverse for binary
        # table columns.
        for col in self.columns:
            array = col.array
            if not isinstance(array, chararray.chararray):
                continue
            for i in range(len(array)):
                al = len(array[i])
                if isinstance(array[i], unicode):
                    pad = self._padding_byte
                else:
                    pad = self._padding_byte.encode('ascii')
                array[i] = array[i] + (pad * (array.itemsize - al))

    def __getattr__(self, name):
        """
        Automatically returns the values for the given keyword attribute for
        all `Column`s in this list.

        Implements for example self.units, self.formats, etc.
        """

        cname = name[:-1]
        if cname in KEYWORD_ATTRIBUTES and name[-1] == 's':
            attr = [''] * len(self)
            for idx in range(len(self)):
                val = getattr(self[idx], cname)
                if val is not None:
                    attr[idx] = val
            self.__dict__[name] = attr
            return self.__dict__[name]
        raise AttributeError(name)

    @property
    @deprecated(message='The %(func)s attribute is deprecated; use the '
                        '%(alternative)s attribute instead.',
                alternative='.columns')
    def data(self):
        """
        What was originally self.columns is now self.data; this provides some
        backwards compatibility.
        """

        return self.columns

    @lazyproperty
    def _arrays(self):
        return [col.array for col in self.columns]

    @lazyproperty
    def _recformats(self):
        return [_convert_format(fmt) for fmt in self.formats]

    def __getitem__(self, key):
        x = self.columns[key]
        if _is_int(key):
            return x
        else:
            return ColDefs(x)

    def __len__(self):
        return len(self.columns)

    def __repr__(self):
        rep = 'ColDefs('
        if self.columns:
            rep += '\n    '
            rep += '\n    '.join([repr(c) for c in self.columns])
            rep += '\n'
        rep += ')'
        return rep

    def __add__(self, other, option='left'):
        if isinstance(other, Column):
            b = [other]
        elif isinstance(other, ColDefs):
            b = list(other.columns)
        else:
            raise TypeError('Wrong type of input.')
        if option == 'left':
            tmp = list(self.columns) + b
        else:
            tmp = b + list(self.columns)
        return ColDefs(tmp)

    def __radd__(self, other):
        return self.__add__(other, 'right')

    def __sub__(self, other):
        if not isinstance(other, (list, tuple)):
            other = [other]
        _other = [_get_index(self.names, key) for key in other]
        indx = range(len(self))
        for x in _other:
            indx.remove(x)
        tmp = [self[i] for i in indx]
        return ColDefs(tmp)

    def _update_listener(self):
        if hasattr(self, '_listener'):
            if self._listener._data_loaded:
                del self._listener.data
            self._listener.columns = self

    def add_col(self, column):
        """
        Append one `Column` to the column definition.

        .. warning::

            *New in pyfits 2.3*: This function appends the new column
            to the `ColDefs` object in place.  Prior to pyfits 2.3,
            this function returned a new `ColDefs` with the new column
            at the end.
        """

        assert isinstance(column, Column)

        for cname in KEYWORD_ATTRIBUTES:
            attr = getattr(self, cname+'s')
            attr.append(getattr(column, cname))

        self._arrays.append(column.array)
        # Obliterate caches of certain things
        del self._recformats

        self.columns.append(column)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()
        return self

    def del_col(self, col_name):
        """
        Delete (the definition of) one `Column`.

        col_name : str or int
            The column's name or index
        """

        indx = _get_index(self.names, col_name)

        for cname in KEYWORD_ATTRIBUTES:
            attr = getattr(self, cname+'s')
            del attr[indx]

        del self._arrays[indx]
        # Obliterate caches of certain things
        del self._recformats

        del self.columns[indx]

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()
        return self

    def change_attrib(self, col_name, attrib, new_value):
        """
        Change an attribute (in the commonName list) of a `Column`.

        col_name : str or int
            The column name or index to change

        attrib : str
            The attribute name

        value : object
            The new value for the attribute
        """

        indx = _get_index(self.names, col_name)
        getattr(self, attrib+'s')[indx] = new_value

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

    def change_name(self, col_name, new_name):
        """
        Change a `Column`'s name.

        col_name : str
            The current name of the column

        new_name : str
            The new name of the column
        """

        if new_name != col_name and new_name in self.names:
            raise ValueError('New name %s already exists.' % new_name)
        else:
            self.change_attrib(col_name, 'name', new_name)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

    def change_unit(self, col_name, new_unit):
        """
        Change a `Column`'s unit.

        col_name : str or int
            The column name or index

        new_unit : str
            The new unit for the column
        """

        self.change_attrib(col_name, 'unit', new_unit)

        # If this ColDefs is being tracked by a Table, inform the
        # table that its data is now invalid.
        self._update_listener()

    def info(self, attrib='all', output=None):
        """
        Get attribute(s) information of the column definition.

        Parameters
        ----------
        attrib : str
            Can be one or more of the attributes listed in
            `KEYWORD_ATTRIBUTES`.  The default is ``"all"`` which will print
            out all attributes.  It forgives plurals and blanks.  If
            there are two or more attribute names, they must be
            separated by comma(s).

        output : file, optional
            File-like object to output to.  Outputs to stdout by default.
            If False, returns the attributes as a dict instead.

        Notes
        -----
        This function doesn't return anything by default; it just prints to
        stdout.
        """

        if output is None:
            output = sys.stdout

        if attrib.strip().lower() in ['all', '']:
            lst = KEYWORD_ATTRIBUTES
        else:
            lst = attrib.split(',')
            for idx in range(len(lst)):
                lst[idx] = lst[idx].strip().lower()
                if lst[idx][-1] == 's':
                    lst[idx]=list[idx][:-1]

        ret = {}

        for attr in lst:
            if output:
                if attr not in KEYWORD_ATTRIBUTES:
                    output.write("'%s' is not an attribute of the column "
                                 "definitions.\n" % attr)
                continue
                output.write("%s:\n" % attr)
                output.write('    %s\n' % getattr(self, attr + 's'))
            else:
                ret[attr] = getattr(self, attr + 's')

        if not output:
            return ret


class _ASCIIColDefs(ColDefs):
    """ColDefs implementation for ASCII tables."""

    _ascii_fmt = {'A':'A1', 'I':'I10', 'J':'I15', 'E':'E15.7', 'F':'F16.7',
                  'D':'D25.17'}

    _padding_byte = ' '

    def __init__(self, input, tbtype='TableHDU'):
        super(_ASCIIColDefs, self).__init__(input, tbtype)

        # if the format of an ASCII column has no width, add one
        if not isinstance(input, _ASCIIColDefs):
            for col in self.columns:
                (type, width) = _convert_ascii_format(col.format)
                if width is None:
                    col.format = self._ascii_fmt[col.format]

    @lazyproperty
    def spans(self):
        # make sure to consider the case that the starting column of
        # a field may not be the column right after the last field
        end = 0
        spans = [0] * len(self)
        for idx in range(len(self)):
            format, width = _convert_ascii_format(self.formats[idx])
            if not self.starts[idx]:
                self.starts[idx] = end + 1
            end = self.starts[idx] + width - 1
            spans[idx] = width
        self._width = end
        return spans

    @lazyproperty
    def _recformats(self):
        if len(self) == 1:
            widths = []
        else:
            widths = [y - x for x, y in pairwise(self.starts)]
        # NOTE: The self._width attribute only exists if this ColDefs was
        # instantiated with a _TableHDU object; make sure that's the only
        # context in which this is used, for now...
        # Touch spans to make sure self.starts is set
        self.spans
        widths.append(self._width - self.starts[-1] + 1)
        return ['a' + str(w) for w in widths]

    def add_col(self, column):
        # Clear existing spans value
        del self.spans
        super(_ASCIIColDefs, self).add_col(column)

    def del_col(self, col_name):
        # Clear existing spans value
        del self.spans
        super(_ASCIIColDefs, self).del_col(col_name)


class _VLF(np.ndarray):
    """Variable length field object."""

    def __new__(cls, args):
        """
        Parameters
        ----------
        args
            a sequence of variable-sized elements.
        """

        a = np.array(args, dtype=np.object)
        self = np.ndarray.__new__(cls, shape=(len(args)), buffer=a,
                                  dtype=np.object)
        self._max = 0
        return self

    def __array_finalize__(self, obj):
        if obj is None:
            return
        self._max = obj._max

    def __setitem__(self, key, value):
        """
        To make sure the new item has consistent data type to avoid
        misalignment.
        """

        if isinstance(value, np.ndarray) and value.dtype == self.dtype:
            pass
        elif isinstance(value, chararray.chararray) and value.itemsize == 1:
            pass
        elif self._dtype == 'a':
            value = chararray.array(value, itemsize=1)
        else:
            value = np.array(value, dtype=self._dtype)
        np.ndarray.__setitem__(self, key, value)
        self._max = max(self._max, len(value))


def _get_index(names, key):
    """
    Get the index of the `key` in the `names` list.

    The `key` can be an integer or string.  If integer, it is the index
    in the list.  If string,

        a. Field (column) names are case sensitive: you can have two
           different columns called 'abc' and 'ABC' respectively.

        b. When you *refer* to a field (presumably with the field
           method), it will try to match the exact name first, so in
           the example in (a), field('abc') will get the first field,
           and field('ABC') will get the second field.

        If there is no exact name matched, it will try to match the
        name with case insensitivity.  So, in the last example,
        field('Abc') will cause an exception since there is no unique
        mapping.  If there is a field named "XYZ" and no other field
        name is a case variant of "XYZ", then field('xyz'),
        field('Xyz'), etc. will get this field.
    """

    if _is_int(key):
        indx = int(key)
    elif isinstance(key, basestring):
        # try to find exact match first
        try:
            indx = names.index(key.rstrip())
        except ValueError:
            # try to match case-insentively,
            _key = key.lower().rstrip()
            names = [n.lower().rstrip() for n in names]
            count = names.count(_key) # occurrence of _key in names
            if count == 1:
                indx = names.index(_key)
            elif count == 0:
                raise KeyError("Key '%s' does not exist." % key)
            else:              # multiple match
                raise KeyError("Ambiguous key name '%s'." % key)
    else:
        raise KeyError("Illegal key '%s'." % repr(key))

    return indx


def _unwrapx(input, output, nx):
    """
    Unwrap the X format column into a Boolean array.

    Parameters
    ----------
    input
        input ``Uint8`` array of shape (`s`, `nbytes`)

    output
        output Boolean array of shape (`s`, `nx`)

    nx
        number of bits
    """

    pow2 = np.array([128, 64, 32, 16, 8, 4, 2, 1], dtype='uint8')
    nbytes = ((nx-1) // 8) + 1
    for i in range(nbytes):
        _min = i*8
        _max = min((i+1)*8, nx)
        for j in range(_min, _max):
            output[...,j] = np.bitwise_and(input[...,i], pow2[j-i*8])


def _wrapx(input, output, nx):
    """
    Wrap the X format column Boolean array into an ``UInt8`` array.

    Parameters
    ----------
    input
        input Boolean array of shape (`s`, `nx`)

    output
        output ``Uint8`` array of shape (`s`, `nbytes`)

    nx
        number of bits
    """

    output[...] = 0 # reset the output
    nbytes = ((nx-1) // 8) + 1
    unused = nbytes*8 - nx
    for i in range(nbytes):
        _min = i*8
        _max = min((i+1)*8, nx)
        for j in range(_min, _max):
            if j != _min:
                np.left_shift(output[...,i], 1, output[...,i])
            np.add(output[...,i], input[...,j], output[...,i])

    # shift the unused bits
    np.left_shift(output[...,i], unused, output[...,i])


def _makep(input, desp_output, dtype, nrows=None):
    """
    Construct the P format column array, both the data descriptors and
    the data.  It returns the output "data" array of data type `dtype`.

    The descriptor location will have a zero offset for all columns
    after this call.  The final offset will be calculated when the file
    is written.

    Parameters
    ----------
    input
        input object array

    desp_output
        output "descriptor" array of data type ``Int32``--must be nrows wide in
        its first dimension

    dtype
        data type of the variable array

    nrows : int, optional
        number of rows to create in the column; defaults to the number of rows
        in the input array
    """

    _offset = 0

    if not nrows:
        nrows = len(input)
    n = min(len(input), nrows)

    data_output = _VLF([None] * nrows)
    data_output._dtype = dtype

    if dtype == 'a':
        _nbytes = 1
    else:
        _nbytes = np.array([], dtype=np.typeDict[dtype]).itemsize

    for idx in range(nrows):
        if idx < len(input):
            rowval = input[idx]
        else:
            if dtype == 'a':
                rowval = ' ' * data_output._max
            else:
                rowval = [0] * data_output._max
        if dtype == 'a':
            data_output[idx] = chararray.array(encode_ascii(rowval),
                                               itemsize=1)
        else:
            data_output[idx] = np.array(rowval, dtype=dtype)

        desp_output[idx,0] = len(data_output[idx])
        desp_output[idx,1] = _offset
        _offset += len(data_output[idx]) * _nbytes

    return data_output


def _parse_tformat(tform):
    """Parse the ``TFORM`` value into `repeat`, `dtype`, and `option`."""

    try:
        (repeat, dtype, option) = TFORMAT_RE.match(tform.strip()).groups()
    except:
        warnings.warn('Format "%s" is not recognized.' % tform)


    if repeat == '':
        repeat = 1
    else:
        repeat = int(repeat)

    return (repeat, dtype, option)


def _scalar_to_format(value):
    """
    Given a scalar value or string, returns the minimum FITS column format
    that can represent that value.  'minimum' is defined by the order given in
    FORMATORDER.
    """

    # TODO: Numpy 1.6 and up has a min_scalar_type() function that can handle
    # this; in the meantime we have to use our own implementation (which for
    # now is pretty naive)

    # First, if value is a string, try to convert to the appropriate scalar
    # value
    for type_ in (int, float, complex):
        try:
            value = type_(value)
            break
        except ValueError:
            continue

    if isinstance(value, int) and value in (0, 1):
        # Could be a boolean
        return 'L'
    elif isinstance(value, int):
        for char in ('B', 'I', 'J', 'K'):
            type_ = np.dtype(FITS2NUMPY[char]).type
            if type_(value) == value:
                return char
    elif isinstance(value, float):
        # For now just assume double precision
        return 'D'
    elif isinstance(value, complex):
        return 'M'
    else:
        return 'A' + str(len(value))

def _cmp_recformats(f1, f2):
    """
    Compares two numpy recformats using the ordering given by FORMATORDER.
    """

    if f1[0] == 'a' and f2[0] == 'a':
        return cmp(int(f1[1:]), int(f2[1:]))
    else:
        f1, f2 = NUMPY2FITS[f1], NUMPY2FITS[f2]
        return cmp(FORMATORDER.index(f1), FORMATORDER.index(f2))

def _convert_fits2record(format):
    """
    Convert FITS format spec to record format spec.
    """

    repeat, dtype, option = _parse_tformat(format)

    if dtype in FITS2NUMPY:                            # FITS format
        if dtype == 'A':
            output_format = FITS2NUMPY[dtype] + str(repeat)
            # to accomodate both the ASCII table and binary table column
            # format spec, i.e. A7 in ASCII table is the same as 7A in
            # binary table, so both will produce 'a7'.
            if format.lstrip()[0] == 'A' and option != '':
                 # make sure option is integer
                output_format = FITS2NUMPY[dtype] + str(int(option))
        else:
            repeat_str = ''
            if repeat != 1:
                repeat_str = str(repeat)
            output_format = repeat_str + FITS2NUMPY[dtype]

    elif dtype == 'X':
        nbytes = ((repeat-1) // 8) + 1
        # use an array, even if it is only ONE u1 (i.e. use tuple always)
        output_format = _FormatX(repr((nbytes,)) + 'u1')
        output_format._nx = repeat

    elif dtype == 'P':
        output_format = _FormatP('2i4')
        output_format._dtype = FITS2NUMPY[option[0]]
    elif dtype == 'F':
        output_format = 'f8'
    else:
        raise ValueError('Illegal format %s.' % format)

    return output_format


def _convert_record2fits(format):
    """
    Convert record format spec to FITS format spec.
    """

    if isinstance(format, np.dtype):
        shape = format.shape
        kind = format.base.kind
        option = str(format.base.itemsize)
        if kind in ('U', 'S'):
            kind = 'a'
        dtype = kind

        ndims = len(shape)
        repeat = 1
        if ndims > 0:
            nel = np.array(shape, dtype='i8').prod()
            if nel > 1:
                repeat = nel
    else:
        repeat, dtype, option = _parse_tformat(format)

    if dtype == 'a':
        # This is a kludge that will place string arrays into a
        # single field, so at least we won't lose data.  Need to
        # use a TDIM keyword to fix this, declaring as (slength,
        # dim1, dim2, ...)  as mwrfits does

        ntot = int(repeat) * int(option)

        output_format = str(ntot) + NUMPY2FITS[dtype]
    elif isinstance(dtype, _FormatX):
        warnings.warn('X format')
    elif dtype + option in NUMPY2FITS: # record format
        if repeat != 1:
            repeat = str(repeat)
        else:
            repeat = ''
        output_format = repeat + NUMPY2FITS[dtype + option]
    else:
        raise ValueError('Illegal format %s.' % format)

    return output_format

def _convert_format(format, reverse=False):
    """
    Convert FITS format spec to record format spec.  Do the opposite if
    reverse=True.
    """

    if reverse:
        return _convert_record2fits(format)
    else:
        return _convert_fits2record(format)


def _convert_ascii_format(input_format):
    """Convert ASCII table format spec to record format spec."""

    ascii2rec = {'A': 'a', 'I': 'i4', 'J': 'i8', 'F': 'f4', 'E': 'f4',
                 'D': 'f8'}
    _re = re.compile(r'(?P<dtype>[AIJFED])(?P<width>[0-9]*)')

    # Parse the TFORM value into data type and width.
    try:
        (dtype, width) = _re.match(input_format.strip()).groups()
        dtype = ascii2rec[dtype]
        if width == '':
            width = None
        else:
            width = int(width)
    except KeyError:
        raise ValueError('Illegal format `%s` for ASCII table.'
                         % input_format)

    return (dtype, width)
