# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division  # confidence high

import csv
import os
import re
import sys
import textwrap
import warnings

import numpy as np
from numpy import char as chararray

from .base import DELAYED, _ValidHDU, ExtensionHDU
from ..column import (FITS2NUMPY, KEYWORD_NAMES, KEYWORD_ATTRIBUTES, TDEF_RE,
                      Column, ColDefs, _AsciiColDefs, _FormatX, _FormatP,
                      _FormatQ, _makep, _VLF, _parse_tformat,
                      _scalar_to_format, _convert_format, _cmp_recformats,
                      _get_index)
from ..fitsrec import FITS_rec
from ..header import Header, _pad_length
from ..util import _is_int, _str_to_num

from ....utils import deprecated, lazyproperty
from ....utils.exceptions import AstropyUserWarning


class FITSTableDumpDialect(csv.excel):
    """
    A CSV dialect for the PyFITS format of ASCII dumps of FITS tables.
    """

    delimiter = ' '
    lineterminator = '\n'
    quotechar = '"'
    quoting = csv.QUOTE_ALL
    skipinitialspace = True


class _TableLikeHDU(_ValidHDU):
    """
    A class for HDUs that have table-like data.  This is used for both
    Binary/ASCII tables as well as Random Access Group HDUs (which are
    otherwise too dissimlary for tables to use _TableBaseHDU directly).
    """

    _data_type = FITS_rec
    _columns_type = ColDefs

    # TODO: Temporary flag representing whether uints are enabled; remove this
    # after restructuring to support uints by default on a per-column basis
    _uint = False

    @classmethod
    def match_header(cls, header):
        """
        This is an abstract HDU type for HDUs that contain table-like data.
        This is even more abstract than _TableBaseHDU which is specifically for
        the standard ASCII and Binary Table types.
        """

        raise NotImplementedError

    @lazyproperty
    def columns(self):
        # The base class doesn't make any assumptions about where the column
        # definitions come from, so just return an empty ColDefs
        return ColDefs([])

    def _get_tbdata(self):
        """Get the table data from an input HDU object."""

        columns = self.columns

        # TODO: Details related to variable length arrays need to be dealt with
        # specifically in the BinTableHDU class, since they're a detail
        # specific to FITS binary tables
        if (any(type(r) in (_FormatP, _FormatQ)
                for r in columns._recformats) and
                self._data_size > self._theap):
            # We have a heap; include it in the raw_data
            raw_data = self._get_raw_data(self._data_size, np.uint8,
                                          self._data_offset)
            data = raw_data[:self._theap].view(dtype=columns.dtype,
                                               type=np.rec.recarray)
        else:
            raw_data = self._get_raw_data(columns._shape, columns.dtype,
                                          self._data_offset)
            data = raw_data.view(np.rec.recarray)

        self._init_tbdata(data)
        return data.view(self._data_type)

    def _init_tbdata(self, data):
        columns = self.columns

        data.dtype = data.dtype.newbyteorder('>')

        # hack to enable pseudo-uint support
        data._uint = self._uint

        # pass datLoc, for P format
        data._heapoffset = self._theap
        data._heapsize = self._header['PCOUNT']
        tbsize = self._header['NAXIS1'] * self._header['NAXIS2']
        data._gap = self._theap - tbsize

        # pass the attributes
        fidx = 0
        for idx in range(len(columns)):
            if not columns[idx]._phantom:
                # get the data for each column object from the rec.recarray
                columns[idx].array = data.field(fidx)
                fidx += 1

        # delete the _arrays attribute so that it is recreated to point to the
        # new data placed in the column object above
        del columns._arrays


class _TableBaseHDU(ExtensionHDU, _TableLikeHDU):
    """
    FITS table extension base HDU class.
    """

    def __init__(self, data=None, header=None, name=None, uint=False):
        """
        Parameters
        ----------
        header : Header instance
            header to be used

        data : array
            data to be used

        name : str
            name to be populated in ``EXTNAME`` keyword

        uint : bool, optional
            set to ``True`` if the table contains unsigned integer columns.
        """

        super(_TableBaseHDU, self).__init__(data=data, header=header,
                                            name=name)

        if header is not None and not isinstance(header, Header):
            raise ValueError('header must be a Header object.')
        self._uint = uint
        if data is DELAYED:
            # this should never happen
            if header is None:
                raise ValueError('No header to setup HDU.')

            # if the file is read the first time, no need to copy, and keep it
            # unchanged
            else:
                self._header = header
        else:
            # construct a list of cards of minimal header
            cards = [
                ('XTENSION',      '', ''),
                ('BITPIX',         8, 'array data type'),
                ('NAXIS',          2, 'number of array dimensions'),
                ('NAXIS1',         0, 'length of dimension 1'),
                ('NAXIS2',         0, 'length of dimension 2'),
                ('PCOUNT',         0, 'number of group parameters'),
                ('GCOUNT',         1, 'number of groups'),
                ('TFIELDS',        0, 'number of table fields')]

            if header is not None:
                # Make a "copy" (not just a view) of the input header, since it
                # may get modified.  the data is still a "view" (for now)
                hcopy = header.copy(strip=True)
                cards.extend(hcopy.cards)

            self._header = Header(cards)

            if isinstance(data, np.ndarray) and data.dtype.fields is not None:
                # self._data_type is FITS_rec.
                if isinstance(data, self._data_type):
                    self.data = data
                else:
                    # Just doing a view on the input data screws up unsigned
                    # columns, so treat those more carefully.
                    # TODO: I need to read this code a little more closely
                    # again, but I think it can be simplified quite a bit with
                    # the use of some appropriate utility functions
                    update_coldefs = {}
                    if 'u' in [data.dtype[k].kind for k in data.dtype.names]:
                        self._uint = True
                        bzeros = {2: np.uint16(2**15), 4: np.uint32(2**31),
                                  8: np.uint64(2**63)}

                        new_dtype = [
                            (k, data.dtype[k].kind.replace('u', 'i') +
                            str(data.dtype[k].itemsize))
                            for k in data.dtype.names]

                        new_data = np.zeros(data.shape, dtype=new_dtype)

                        for k in data.dtype.fields:
                            dtype = data.dtype[k]
                            if dtype.kind == 'u':
                                new_data[k] = data[k] - bzeros[dtype.itemsize]
                                update_coldefs[k] = bzeros[dtype.itemsize]
                            else:
                                new_data[k] = data[k]
                        self.data = new_data.view(self._data_type)
                        # Uck...
                        self.data._uint = True
                    else:
                        self.data = data.view(self._data_type)
                    for k in update_coldefs:
                        indx = _get_index(self.data.names, k)
                        self.data._coldefs[indx].bzero = update_coldefs[k]
                        # This is so bad that we have to update this in
                        # duplicate...
                        self.data._coldefs.bzeros[indx] = update_coldefs[k]
                        # More uck...
                        self.data._coldefs[indx]._physical_values = False
                        self.data._coldefs[indx]._pseudo_unsigned_ints = True

                self._header['NAXIS1'] = self.data.itemsize
                self._header['NAXIS2'] = self.data.shape[0]
                self._header['TFIELDS'] = len(self.data._coldefs)

                self.columns = self.data._coldefs
                self.update()

                try:
                   # Make the ndarrays in the Column objects of the ColDefs
                   # object of the HDU reference the same ndarray as the HDU's
                   # FITS_rec object.
                    for idx in range(len(self.columns)):
                        self.columns[idx].array = self.data.field(idx)

                    # Delete the _arrays attribute so that it is recreated to
                    # point to the new data placed in the column objects above
                    del self.columns._arrays
                except (TypeError, AttributeError) as e:
                    # This shouldn't happen as long as self.columns._arrays
                    # is a lazyproperty
                    pass
            elif data is None:
                pass
            else:
                raise TypeError('Table data has incorrect type.')

        if not (isinstance(self._header[0], basestring) and
                self._header[0].rstrip() == self._extension):
            self._header[0] = (self._extension, self._ext_comment)

        # Ensure that the correct EXTNAME is set on the new header if one was
        # created, or that it overrides the existing EXTNAME if different
        if name:
            self.name = name

    @classmethod
    def match_header(cls, header):
        """
        This is an abstract type that implements the shared functionality of
        the ASCII and Binary Table HDU types, which should be used instead of
        this.
        """

        raise NotImplementedError

    @lazyproperty
    def columns(self):
        if self._has_data and hasattr(self.data, '_coldefs'):
            return self.data._coldefs
        return self._columns_type(self)

    @lazyproperty
    def data(self):
        data = self._get_tbdata()
        data._coldefs = self.columns
        data.formats = self.columns.formats
        # Columns should now just return a reference to the data._coldefs
        del self.columns
        return data

    @data.setter
    def data(self, data):
        if 'data' in self.__dict__:
            if self.__dict__['data'] is data:
                return
            else:
                self._data_replaced = True
        else:
            self._data_replaced = True

        self._modified = True

        if data is None and self.columns:
            # Create a new table with the same columns, but empty rows
            formats = ','.join(self.columns._recformats)
            data = np.rec.array(None, formats=formats,
                                names=self.columns.names,
                                shape=0)

        if isinstance(data, np.ndarray) and data.dtype.fields is not None:
            # Go ahead and always make a view, even if the data is already the
            # correct class (self._data_type) so we can update things like the
            # column defs, if necessary
            data = data.view(self._data_type)

            if not isinstance(data.columns, self._columns_type):
                # This would be the place, if the input data was for an ASCII
                # table and this is binary table, or vice versa, to convert the
                # data to the appropriate format for the table type
                new_columns = self._columns_type(data.columns)
                data = FITS_rec.from_columns(new_columns)

            self.__dict__['data'] = data

            self.columns = self.data.columns
            self.update()

            try:
               # Make the ndarrays in the Column objects of the ColDefs
               # object of the HDU reference the same ndarray as the HDU's
               # FITS_rec object.
                for idx in range(len(self.columns)):
                    self.columns[idx].array = self.data.field(idx)

                # Delete the _arrays attribute so that it is recreated to
                # point to the new data placed in the column objects above
                del self.columns._arrays
            except (TypeError, AttributeError):
                # This shouldn't happen as long as self.columns._arrays
                # is a lazyproperty
                pass
        elif data is None:
            pass
        else:
            raise TypeError('Table data has incorrect type.')

        # returning the data signals to lazyproperty that we've already handled
        # setting self.__dict__['data']
        return data

    @lazyproperty
    def _theap(self):
        size = self._header['NAXIS1'] * self._header['NAXIS2']
        return self._header.get('THEAP', size)

    # TODO: Need to either rename this to update_header, for symmetry with the
    # Image HDUs, or just at some point deprecate it and remove it altogether,
    # since header updates should occur automatically when necessary...
    def update(self):
        """
        Update header keywords to reflect recent changes of columns.
        """

        self._header.set('NAXIS1', self.data.itemsize, after='NAXIS')
        self._header.set('NAXIS2', self.data.shape[0], after='NAXIS1')
        self._header.set('TFIELDS', len(self.columns), after='GCOUNT')

        self._clear_table_keywords()
        self._populate_table_keywords()

    def copy(self):
        """
        Make a copy of the table HDU, both header and data are copied.
        """

        # touch the data, so it's defined (in the case of reading from a
        # FITS file)
        self.data
        return new_table(self.columns, header=self._header,
                         tbtype=self.columns._tbtype)

    def _prewriteto(self, checksum=False, inplace=False):
        if self._has_data:
            self.data._scale_back()
            # check TFIELDS and NAXIS2
            self._header['TFIELDS'] = len(self.data._coldefs)
            self._header['NAXIS2'] = self.data.shape[0]

            # calculate PCOUNT, for variable length tables
            tbsize = self.header['NAXIS1'] * self.header['NAXIS2']
            heapstart = self.header.get('THEAP', tbsize)
            self.data._gap = heapstart - tbsize
            pcount = self.data._heapsize + self.data._gap
            if pcount > 0:
                self.header['PCOUNT'] = pcount

            # update the other T****n keywords
            self._populate_table_keywords()

            # update TFORM for variable length columns
            for idx in range(self.data._nfields):
                format = self.data._coldefs._recformats[idx]
                if isinstance(format, _FormatP):
                    _max = self.data.field(idx).max
                    # May be either _FormatP or _FormatQ
                    format_cls = format.__class__
                    format = format_cls(format.dtype, repeat=format.repeat,
                                        max=_max)
                    self._header['TFORM' + str(idx + 1)] = format.tform
        return super(_TableBaseHDU, self)._prewriteto(checksum, inplace)

    def _verify(self, option='warn'):
        """
        _TableBaseHDU verify method.
        """

        errs = super(_TableBaseHDU, self)._verify(option=option)
        self.req_cards('NAXIS', None, lambda v: (v == 2), 2, option, errs)
        self.req_cards('BITPIX', None, lambda v: (v == 8), 8, option, errs)
        self.req_cards('TFIELDS', 7,
                       lambda v: (_is_int(v) and v >= 0 and v <= 999), 0,
                       option, errs)
        tfields = self._header['TFIELDS']
        for idx in range(tfields):
            self.req_cards('TFORM' + str(idx + 1), None, None, None, option,
                           errs)
        return errs

    def _summary(self):
        """
        Summarize the HDU: name, dimensions, and formats.
        """

        class_name = self.__class__.__name__

        # if data is touched, use data info.
        if self._data_loaded:
            if self.data is None:
                shape, format = (), ''
                nrows = 0
            else:
                nrows = len(self.data)

            ncols = len(self.columns)
            format = self.columns.formats

        # if data is not touched yet, use header info.
        else:
            shape = ()
            nrows = self._header['NAXIS2']
            ncols = self._header['TFIELDS']
            format = ', '.join([self._header['TFORM' + str(j + 1)]
                                for j in range(ncols)])
            format = '[%s]' % format
        dims = "%dR x %dC" % (nrows, ncols)
        ncards = len(self._header)

        return (self.name, class_name, ncards, dims, format)

    def _clear_table_keywords(self):
        """Wipe out any existing table definition keywords from the header."""

        # Go in reverse so as to not confusing indexing while deleting.
        for idx, keyword in reversed(list(enumerate(self._header.keys()))):
            keyword = TDEF_RE.match(keyword)
            try:
                keyword = keyword.group('label')
            except:
                continue                # skip if there is no match
            if keyword in KEYWORD_NAMES:
                del self._header[idx]

    def _populate_table_keywords(self):
        """Populate the new table definition keywords from the header."""

        cols = self.columns

        for idx in range(len(cols)):
            for attr, keyword in zip(KEYWORD_ATTRIBUTES, KEYWORD_NAMES):
                val = getattr(cols, attr + 's')[idx]
                if val:
                    keyword = keyword + str(idx + 1)
                    self._header[keyword] = val


class TableHDU(_TableBaseHDU):
    """
    FITS ASCII table extension HDU class.
    """

    _extension = 'TABLE'
    _ext_comment = 'ASCII table extension'

    _padding_byte = ' '
    _columns_type = _AsciiColDefs

    __format_RE = re.compile(
        r'(?P<code>[ADEFIJ])(?P<width>\d+)(?:\.(?P<prec>\d+))?')

    def __init__(self, data=None, header=None, name=None):
        super(TableHDU, self).__init__(data, header, name=name)

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        return card.keyword == 'XTENSION' and xtension == cls._extension

    def _get_tbdata(self):
        columns = self.columns
        names = [n for idx, n in enumerate(columns.names)
                 if not columns[idx]._phantom]

        # determine if there are duplicate field names and if there
        # are throw an exception
        dup = np.rec.find_duplicate(names)

        if dup:
            raise ValueError("Duplicate field names: %s" % dup)

        # TODO: Determine if this extra logic is necessary--I feel like the
        # _AsciiColDefs class should be responsible for telling the table what
        # its dtype should be...
        itemsize = columns.spans[-1] + columns.starts[-1] - 1
        dtype = {}

        for idx in range(len(columns)):
            data_type = 'S' + str(columns.spans[idx])

            if idx == len(columns) - 1:
                # The last column is padded out to the value of NAXIS1
                if self._header['NAXIS1'] > itemsize:
                    data_type = 'S' + str(columns.spans[idx] +
                                self._header['NAXIS1'] - itemsize)
            dtype[columns.names[idx]] = (data_type, columns.starts[idx] - 1)

        raw_data = self._get_raw_data(columns._shape, dtype, self._data_offset)
        data = raw_data.view(np.rec.recarray)
        self._init_tbdata(data)
        return data.view(self._data_type)

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._has_data:
            # We have the data to be used.
            # We need to pad the data to a block length before calculating
            # the datasum.

            d = np.append(self.data.view(dtype='ubyte'),
                          np.fromstring(_pad_length(self.size) * ' ',
                                        dtype='ubyte'))

            cs = self._compute_checksum(d, blocking=blocking)
            return cs
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a gereric manner.
            return super(TableHDU, self)._calculate_datasum(blocking)

    def _verify(self, option='warn'):
        """
        `TableHDU` verify method.
        """

        errs = super(TableHDU, self)._verify(option=option)
        self.req_cards('PCOUNT', None, lambda v: (v == 0), 0, option, errs)
        tfields = self._header['TFIELDS']
        for idx in range(tfields):
            self.req_cards('TBCOL' + str(idx + 1), None, _is_int, None, option,
                           errs)
        return errs


class BinTableHDU(_TableBaseHDU):
    """
    Binary table HDU class.
    """

    _extension = 'BINTABLE'
    _ext_comment = 'binary table extension'

    @classmethod
    def match_header(cls, header):
        card = header.cards[0]
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        return (card.keyword == 'XTENSION' and
                xtension in (cls._extension, 'A3DTABLE'))

    def _calculate_datasum_from_data(self, data, blocking):
        """
        Calculate the value for the ``DATASUM`` card given the input data
        """

        # Check the byte order of the data.  If it is little endian we
        # must swap it before calculating the datasum.
        for i in range(data._nfields):
            coldata = data.field(i)

            if not isinstance(coldata, chararray.chararray):
                if isinstance(coldata, _VLF):
                    for j, d in enumerate(coldata):
                        if not isinstance(d, chararray.chararray):
                            if d.itemsize > 1:
                                if d.dtype.str[0] != '>':
                                    d[:] = d.byteswap()
                                    d.dtype = d.dtype.newbyteorder('>')
                        field = np.rec.recarray.field(data, i)[j:j + 1]
                        if field.dtype.str[0] != '>':
                            field.byteswap(True)
                else:
                    if coldata.itemsize > 1:
                        if data.field(i).dtype.str[0] != '>':
                            data.field(i)[:] = data.field(i).byteswap()
        data.dtype = data.dtype.newbyteorder('>')

        dout = data.view(dtype='ubyte')

        for i in range(data._nfields):
            if isinstance(data._coldefs._recformats[i], _FormatP):
                for coldata in data.field(i):
                    if len(coldata) > 0:
                        dout = np.append(dout, coldata.view(dtype='ubyte'))

        cs = self._compute_checksum(dout, blocking=blocking)
        return cs

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._has_data:
            # We have the data to be used.
            return self._calculate_datasum_from_data(self.data, blocking)
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
            return super(BinTableHDU, self)._calculate_datasum(blocking)

    def _writedata_internal(self, fileobj):
        size = 0

        if self.data is not None:
            size += self._binary_table_byte_swap(fileobj)
            size += self.data.size * self.data.itemsize

        return size

    def _binary_table_byte_swap(self, fileobj):
        """Prepares data in the native FITS format and writes the raw bytes
        out to the given file object.  This handles byte swapping from native
        to big endian (if necessary).  In addition, however, this also handles
        writing the binary table heap when variable length array columns are
        present.
        """

        to_swap = []
        swapped = []
        nbytes = 0
        if sys.byteorder == 'little':
            swap_types = ('<', '=')
        else:
            swap_types = ('<',)
        try:
            if not fileobj.simulateonly:
                for idx in range(self.data._nfields):
                    field = np.rec.recarray.field(self.data, idx)
                    if isinstance(field, chararray.chararray):
                        continue
                    recformat = self.data.columns._recformats[idx]
                    # only swap unswapped
                    if field.itemsize > 1 and field.dtype.str[0] in swap_types:
                        to_swap.append(field)
                    # deal with var length table
                    if isinstance(recformat, _FormatP):
                        coldata = self.data.field(idx)
                        for c in coldata:
                            if (not isinstance(c, chararray.chararray) and
                                c.itemsize > 1 and
                                    c.dtype.str[0] in swap_types):
                                to_swap.append(c)

                while to_swap:
                    obj = to_swap.pop()
                    obj.byteswap(True)
                    swapped.append(obj)

                fileobj.writearray(self.data)

                # write out the heap of variable length array
                # columns this has to be done after the
                # "regular" data is written (above)
                fileobj.write((self.data._gap * '\0').encode('ascii'))

            nbytes = self.data._gap

            for idx in range(self.data._nfields):
                if isinstance(self.data.columns._recformats[idx], _FormatP):
                    field = self.data.field(idx)
                    for row in field:
                        if len(row) > 0:
                            nbytes += row.nbytes
                            if not fileobj.simulateonly:
                                fileobj.writearray(row)

            self.data._heapsize = nbytes - self.data._gap
        finally:
            for obj in swapped:
                obj.byteswap(True)

        return nbytes

    def __populate_table_keywords(self):
        """Populate the new table definition keywords from the header."""

        cols = self.columns

        for idx in range(len(cols)):
            for attr, keyword in zip(KEYWORD_ATTRIBUTES, KEYWORD_NAMES):
                val = getattr(cols, attr + 's')[idx]
                if val:
                    keyword = keyword + str(idx + 1)
                    if attr == 'format':
                        val = cols._recformats[idx]
                        if isinstance(val, (_FormatX, _FormatP)):
                            val = val.tform
                        else:
                            # There are some cases where the original TFORM and
                            # the one generated from the recformat can have the
                            # same meaning but different string representation;
                            # make sure to use the original representation in
                            # this case
                            orig_val = cols.formats[idx]
                            val = _convert_format(val, reverse=True)
                            if _parse_tformat(orig_val) == _parse_tformat(val):
                                val = orig_val
                    self._header[keyword] = val

    _tdump_file_format = textwrap.dedent("""

        - **datafile:** Each line of the data file represents one row of table
          data.  The data is output one column at a time in column order.  If
          a column contains an array, each element of the column array in the
          current row is output before moving on to the next column.  Each row
          ends with a new line.

          Integer data is output right-justified in a 21-character field
          followed by a blank.  Floating point data is output right justified
          using 'g' format in a 21-character field with 15 digits of
          precision, followed by a blank.  String data that does not contain
          whitespace is output left-justified in a field whose width matches
          the width specified in the ``TFORM`` header parameter for the
          column, followed by a blank.  When the string data contains
          whitespace characters, the string is enclosed in quotation marks
          (``""``).  For the last data element in a row, the trailing blank in
          the field is replaced by a new line character.

          For column data containing variable length arrays ('P' format), the
          array data is preceded by the string ``'VLA_Length= '`` and the
          integer length of the array for that row, left-justified in a
          21-character field, followed by a blank.

          .. note::

              This format does *not* support variable length arrays using the
              ('Q' format) due difficult to overcome ambiguities. What
              this means is that this file format cannot support VLA columns
              in tables stored in files that are over 2 GB in size.

          For column data representing a bit field ('X' format), each bit
          value in the field is output right-justified in a 21-character field
          as 1 (for true) or 0 (for false).

        - **cdfile:** Each line of the column definitions file provides the
          definitions for one column in the table.  The line is broken up into
          8, sixteen-character fields.  The first field provides the column
          name (``TTYPEn``).  The second field provides the column format
          (``TFORMn``).  The third field provides the display format
          (``TDISPn``).  The fourth field provides the physical units
          (``TUNITn``).  The fifth field provides the dimensions for a
          multidimensional array (``TDIMn``).  The sixth field provides the
          value that signifies an undefined value (``TNULLn``).  The seventh
          field provides the scale factor (``TSCALn``).  The eighth field
          provides the offset value (``TZEROn``).  A field value of ``""`` is
          used to represent the case where no value is provided.

        - **hfile:** Each line of the header parameters file provides the
          definition of a single HDU header card as represented by the card
          image.
      """)

    def dump(self, datafile=None, cdfile=None, hfile=None, clobber=False):
        """
        Dump the table HDU to a file in ASCII format.  The table may be dumped
        in three separate files, one containing column definitions, one
        containing header parameters, and one for table data.

        Parameters
        ----------
        datafile : file path, file object or file-like object, optional
            Output data file.  The default is the root name of the
            fits file associated with this HDU appended with the
            extension ``.txt``.

        cdfile : file path, file object or file-like object, optional
            Output column definitions file.  The default is `None`, no
            column definitions output is produced.

        hfile : file path, file object or file-like object, optional
            Output header parameters file.  The default is `None`,
            no header parameters output is produced.

        clobber : bool
            Overwrite the output files if they exist.

        Notes
        -----
        The primary use for the `dump` method is to allow viewing and editing
        the table data and parameters in a standard text editor.
        The `load` method can be used to create a new table from the three
        plain text (ASCII) files.
        """

        # TODO: This is looking pretty long and complicated--might be a few
        # places we can break this up into smaller functions

        # check if the output files already exist
        exist = []
        files = [datafile, cdfile, hfile]

        for f in files:
            if isinstance(f, basestring):
                if os.path.exists(f) and os.path.getsize(f) != 0:
                    if clobber:
                        warnings.warn("Overwriting existing file '%s'." % f, AstropyUserWarning)
                    else:
                        exist.append(f)

        if exist:
            raise IOError('  '.join(["File '%s' already exists." % f
                                     for f in exist]))

        # Process the data
        self._dump_data(datafile)

        # Process the column definitions
        if cdfile:
            self._dump_coldefs(cdfile)

        # Process the header parameters
        if hfile:
            self._header.tofile(hfile, sep='\n', endcard=False, padding=False)

    dump.__doc__ += _tdump_file_format.replace('\n', '\n        ')

    @deprecated('3.1', alternative=':meth:`dump`')
    def tdump(self, datafile=None, cdfile=None, hfile=None, clobber=False):
        self.dump(datafile, cdfile, hfile, clobber)

    def load(cls, datafile, cdfile=None, hfile=None, replace=False,
             header=None):
        """
        Create a table from the input ASCII files.  The input is from up to
        three separate files, one containing column definitions, one containing
        header parameters, and one containing column data.

        The column definition and header parameters files are not required.
        When absent the column definitions and/or header parameters are taken
        from the header object given in the header argument; otherwise sensible
        defaults are inferred (though this mode is not recommended).

        Parameters
        ----------
        datafile : file path, file object or file-like object
            Input data file containing the table data in ASCII format.

        cdfile : file path, file object, file-like object, optional
            Input column definition file containing the names,
            formats, display formats, physical units, multidimensional
            array dimensions, undefined values, scale factors, and
            offsets associated with the columns in the table.  If
            `None`, the column definitions are taken from the current
            values in this object.

        hfile : file path, file object, file-like object, optional
            Input parameter definition file containing the header
            parameter definitions to be associated with the table.  If
            `None`, the header parameter definitions are taken from
            the current values in this objects header.

        replace : bool
            When `True`, indicates that the entire header should be
            replaced with the contents of the ASCII file instead of
            just updating the current header.

        header : Header object
            When the cdfile and hfile are missing, use this Header object in
            the creation of the new table and HDU.  Otherwise this Header
            supercedes the keywords from hfile, which is only used to update
            values not present in this Header, unless replace=True in which
            this Header's values are completely replaced with the values from
            hfile.

        Notes
        -----
        The primary use for the `load` method is to allow the input of ASCII
        data that was edited in a standard text editor of the table data and
        parameters.  The `dump` method can be used to create the initial ASCII
        files.
        """

        # Process the parameter file
        if header is None:
            header = Header()

        if hfile:
            if replace:
                header = Header.fromtextfile(hfile)
            else:
                header.extend(Header.fromtextfile(hfile), update=True,
                              update_first=True)

        coldefs = None
        # Process the column definitions file
        if cdfile:
            coldefs = cls._load_coldefs(cdfile)

        # Process the data file
        data = cls._load_data(datafile, coldefs)
        if coldefs is None:
            coldefs = ColDefs(data)

        # Create a new HDU using the supplied header and data
        hdu = cls(data=data, header=header)
        hdu.columns = coldefs
        return hdu
    load.__doc__ += _tdump_file_format.replace('\n', '\n        ')
    load = classmethod(load)
    # Have to create a classmethod from this here instead of as a decorator;
    # otherwise we can't update __doc__

    @deprecated('3.1', alternative=':meth:`load`')
    @classmethod
    def tcreate(cls, datafile, cdfile=None, hfile=None, replace=False,
                header=None):
        return cls.load(datafile, cdfile, hfile, replace, header)

    def _dump_data(self, fileobj):
        """
        Write the table data in the ASCII format read by BinTableHDU.load()
        to fileobj.
        """

        if not fileobj and self._file:
            root = os.path.splitext(self._file.name)[0]
            fileobj = root + '.txt'

        close_file = False

        if isinstance(fileobj, basestring):
            fileobj = open(fileobj, 'w')
            close_file = True

        linewriter = csv.writer(fileobj, dialect=FITSTableDumpDialect)

        # Process each row of the table and output one row at a time
        def format_value(val, format):
            if format[0] == 'S':
                itemsize = int(format[1:])
                return '%-*s' % (itemsize, val)
            elif format in np.typecodes['AllInteger']:
                # output integer
                return '%21d' % val
            elif format in np.typecodes['Complex']:
                return '%21.15g+%.15gj' % (val.real, val.imag)
            elif format in np.typecodes['Float']:
                # output floating point
                return '%#21.15g' % val

        for row in self.data:
            line = []   # the line for this row of the table

            # Process each column of the row.
            for column in self.columns:
                vla_format = None   # format of data in a variable length array
                                    # where None means it is not a VLA
                format = _convert_format(column.format)

                if isinstance(format, _FormatP):
                    # P format means this is a variable length array so output
                    # the length of the array for this row and set the format
                    # for the VLA data
                    line.append('VLA_Length=')
                    line.append('%-21d' % len(row[column.name]))
                    _, dtype, option = _parse_tformat(column.format)
                    vla_format = FITS2NUMPY[option[0]][0]

                if vla_format:
                    # Output the data for each element in the array
                    for val in row[column.name].flat:
                        line.append(format_value(val, vla_format))
                else:
                    # The column data is a single element
                    dtype = self.data.dtype.fields[column.name][0]
                    array_format = dtype.char
                    if array_format == 'S':
                        array_format += str(dtype.itemsize)
                    line.append(format_value(row[column.name], array_format))
            linewriter.writerow(line)
        if close_file:
            fileobj.close()

    def _dump_coldefs(self, fileobj):
        """
        Write the column definition parameters in the ASCII format read by
        BinTableHDU.load() to fileobj.
        """

        close_file = False

        if isinstance(fileobj, basestring):
            fileobj = open(fileobj, 'w')
            close_file = True

        # Process each column of the table and output the result to the
        # file one at a time
        for column in self.columns:
            line = [column.name, column.format]
            attrs = ['disp', 'unit', 'dim', 'null', 'bscale', 'bzero']
            line += ['%-16s' % (value if value else '""')
                     for value in (getattr(column, attr) for attr in attrs)]
            fileobj.write(' '.join(line))
            fileobj.write('\n')

        if close_file:
            fileobj.close()

    @classmethod
    def _load_data(cls, fileobj, coldefs=None):
        """
        Read the table data from the ASCII file output by BinTableHDU.dump().
        """

        close_file = False

        if isinstance(fileobj, basestring):
            fileobj = open(fileobj, 'r')
            close_file = True

        initialpos = fileobj.tell()  # We'll be returning here later
        linereader = csv.reader(fileobj, dialect=FITSTableDumpDialect)

        # First we need to do some preprocessing on the file to find out how
        # much memory we'll need to reserve for the table.  This is necessary
        # even if we already have the coldefs in order to determine how many
        # rows to reserve memory for
        vla_lengths = []
        recformats = []
        names = []
        nrows = 0
        if coldefs is not None:
            recformats = coldefs._recformats
            names = coldefs.names

        def update_recformats(value, idx):
            fitsformat = _scalar_to_format(value)
            recformat = _convert_format(fitsformat)
            if idx >= len(recformats):
                recformats.append(recformat)
            else:
                if _cmp_recformats(recformats[idx], recformat) < 0:
                    recformats[idx] = recformat

        # TODO: The handling of VLAs could probably be simplified a bit
        for row in linereader:
            nrows += 1
            if coldefs is not None:
                continue
            col = 0
            idx = 0
            while idx < len(row):
                if row[idx] == 'VLA_Length=':
                    if col < len(vla_lengths):
                        vla_length = vla_lengths[col]
                    else:
                        vla_length = int(row[idx + 1])
                        vla_lengths.append(vla_length)
                    idx += 2
                    while vla_length:
                        update_recformats(row[idx], col)
                        vla_length -= 1
                        idx += 1
                    col += 1
                else:
                    if col >= len(vla_lengths):
                        vla_lengths.append(None)
                    update_recformats(row[idx], col)
                    col += 1
                    idx += 1

        # Update the recformats for any VLAs
        for idx, length in enumerate(vla_lengths):
            if length is not None:
                recformats[idx] = str(length) + recformats[idx]

        dtype = np.rec.format_parser(recformats, names, None).dtype

        # TODO: In the future maybe enable loading a bit at a time so that we
        # can convert from this format to an actual FITS file on disk without
        # needing enough physical memory to hold the entire thing at once;
        # new_table() could use a similar feature.
        hdu = new_table(np.recarray(shape=1, dtype=dtype), nrows=nrows,
                        fill=True)
        data = hdu.data
        for idx, length in enumerate(vla_lengths):
            if length is not None:
                arr = data.columns._arrays[idx]
                dt = recformats[idx][len(str(length)):]

                # NOTE: FormatQ not supported here; it's hard to determine
                # whether or not it will be necessary to use a wider descriptor
                # type. The function documentation will have to serve as a
                # warning that this is not supported.
                recformats[idx] = _FormatP(dt, max=length)
                data.columns._recformats[idx] = recformats[idx]
                data._convert[idx] = _makep(arr, arr, recformats[idx])

        # Jump back to the start of the data and create a new line reader
        fileobj.seek(initialpos)
        linereader = csv.reader(fileobj, dialect=FITSTableDumpDialect)
        for row, line in enumerate(linereader):
            col = 0
            idx = 0
            while idx < len(line):
                if line[idx] == 'VLA_Length=':
                    vla_len = vla_lengths[col]
                    idx += 2
                    data[row][col][:] = line[idx:idx + vla_len]
                    idx += vla_len
                else:
                    # TODO: This won't work for complex-valued types; fix this
                    # Kind of silly special handling for bools
                    val = line[idx]
                    if recformats[col] == FITS2NUMPY['L']:
                        val = bool(int(val))
                    elif recformats[col] == FITS2NUMPY['M']:
                        # For some reason, in arrays/fields where numpy expects
                        # a complex it's not happy to take a string
                        # representation (though it's happy to do that in other
                        # contexts), so we have to convert the string
                        # representation for it:
                        val = complex(val)
                    data[row][col] = val
                    idx += 1
                col += 1

        if close_file:
            fileobj.close()

        return data

    @classmethod
    def _load_coldefs(cls, fileobj):
        """
        Read the table column definitions from the ASCII file output by
        BinTableHDU.dump().
        """

        close_file = False

        if isinstance(fileobj, basestring):
            fileobj = open(fileobj, 'r')
            close_file = True

        columns = []

        for line in fileobj:
            words = line[:-1].split()
            kwargs = {}
            for key in ['name', 'format', 'disp', 'unit', 'dim']:
                kwargs[key] = words.pop(0).replace('""', '')

            for key in ['null', 'bscale', 'bzero']:
                word = words.pop(0).replace('""', '')
                if word:
                    word = _str_to_num(word)
                kwargs[key] = word
            columns.append(Column(**kwargs))

        if close_file:
            fileobj.close()

        return ColDefs(columns)


@deprecated('3.2',
            alternative=':meth:`FITS_rec.from_columns` to create a new '
                        ':class:`FITS_rec` data object from the input '
                        'columns to pass into the constructor for '
                        ':class:`BinTableHDU` or :class:`TableHDU`',
            pending=True)
def new_table(input, header=None, nrows=0, fill=False, tbtype=BinTableHDU):
    """
    Create a new table from the input column definitions.

    Warning: Creating a new table using this method creates an in-memory *copy*
    of all the column arrays in the input.  This is because if they are
    separate arrays they must be combined into a single contiguous array.

    If the column data is already in a single contiguous array (such as an
    existing record array) it may be better to create a BinTableHDU instance
    directly.  See the Astropy documentation for more details.

    Parameters
    ----------
    input : sequence of Column or ColDefs objects
        The data to create a table from.

    header : Header instance
        Header to be used to populate the non-required keywords.

    nrows : int
        Number of rows in the new table.

    fill : bool
        If `True`, will fill all cells with zeros or blanks.  If
        `False`, copy the data from input, undefined cells will still
        be filled with zeros/blanks.

    tbtype : str or class
        Table type to be created (BinTableHDU or TableHDU) or the class
        name as a string.  Currently only BinTableHDU and TableHDU (ASCII
        tables) are supported.
    """

    # tbtype defaults to classes now, but in all prior version of PyFITS it was
    # a string, so we still support that use case as well
    if not isinstance(tbtype, basestring):
        cls = tbtype
        tbtype = cls.__name__
    else:
        # Right now the string input must be one of 'TableHDU' or 'BinTableHDU'
        # and nothing else, though we will allow this to be case insensitive
        # This could be done more generically through the HDU registry, but my
        # hope is to deprecate this function anyways so there's not much point
        # in trying to make it more "generic".
        if tbtype.lower() == 'tablehdu':
            tbtype = 'TableHDU'
            cls = TableHDU
        elif tbtype.lower() == 'bintablehdu':
            tbtype = 'BinTableHDU'
            cls = BinTableHDU
        else:
            raise ValueError("tbtype must be one of 'TableHDU' or "
                             "'BinTableHDU'")

    if isinstance(input, FITS_rec):  # input is a FITS_rec
        # Create a new ColDefs object from the input FITS_rec's ColDefs
        # object and assign it to the ColDefs attribute of the new hdu.
        columns = ColDefs(input._coldefs, tbtype)
    else:  # input is a list of Columns or possibly a recarray
        # Create a new ColDefs object from the input list of Columns and
        # assign it to the ColDefs attribute of the new hdu.
        columns = ColDefs(input, tbtype)

    data = FITS_rec.from_columns(columns, nrows=nrows, fill=fill)

    # construct a table HDU of the requested type
    return cls(header=header, data=data)
