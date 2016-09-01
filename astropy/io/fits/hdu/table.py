# Licensed under a 3-clause BSD style license - see PYFITS.rst

from __future__ import division  # confidence high

import contextlib
import csv
import operator
import os
import re
import sys
import textwrap
import warnings

import numpy as np
from numpy import char as chararray

from .base import DELAYED, _ValidHDU, ExtensionHDU
# This module may have many dependencies on astropy.io.fits.column, but
# astropy.io.fits.column has fewer dependencies overall, so it's easier to
# keep table/column-related utilities in astropy.io.fits.column
from .. import _numpy_hacks as nh
from ..column import (FITS2NUMPY, KEYWORD_NAMES, KEYWORD_TO_ATTRIBUTE,
                      ATTRIBUTE_TO_KEYWORD, TDEF_RE, Column, ColDefs,
                      _AsciiColDefs, _FormatP, _FormatQ, _makep,
                      _parse_tformat, _scalar_to_format, _convert_format,
                      _cmp_recformats, _get_index)
from ..fitsrec import FITS_rec, _get_recarray_field, _has_unicode_fields
from ..header import Header, _pad_length
from ..util import _is_int, _str_to_num

from ....extern import six
from ....extern.six import string_types
from ....extern.six.moves import range, zip
from ....utils import lazyproperty
from ....utils.compat import suppress
from ....utils.exceptions import AstropyUserWarning
from ....utils.decorators import deprecated_renamed_argument


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
    otherwise too dissimilar for tables to use _TableBaseHDU directly).
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

    @classmethod
    def from_columns(cls, columns, header=None, nrows=0, fill=False,
                     **kwargs):
        """
        Given either a `ColDefs` object, a sequence of `Column` objects,
        or another table HDU or table data (a `FITS_rec` or multi-field
        `numpy.ndarray` or `numpy.recarray` object, return a new table HDU of
        the class this method was called on using the column definition from
        the input.

        See also `FITS_rec.from_columns`.

        Parameters
        ----------
        columns : sequence of `Column`, `ColDefs`, or other
            The columns from which to create the table data, or an object with
            a column-like structure from which a `ColDefs` can be instantiated.
            This includes an existing `BinTableHDU` or `TableHDU`, or a
            `numpy.recarray` to give some examples.

            If these columns have data arrays attached that data may be used in
            initializing the new table.  Otherwise the input columns will be
            used as a template for a new table with the requested number of
            rows.

        header : `Header`
            An optional `Header` object to instantiate the new HDU yet.  Header
            keywords specifically related to defining the table structure (such
            as the "TXXXn" keywords like TTYPEn) will be overridden by the
            supplied column definitions, but all other informational and data
            model-specific keywords are kept.

        nrows : int
            Number of rows in the new table.  If the input columns have data
            associated with them, the size of the largest input column is used.
            Otherwise the default is 0.

        fill : bool
            If `True`, will fill all cells with zeros or blanks.  If `False`,
            copy the data from input, undefined cells will still be filled with
            zeros/blanks.

        Notes
        -----

        Any additional keyword arguments accepted by the HDU class's
        ``__init__`` may also be passed in as keyword arguments.
        """

        coldefs = cls._columns_type(columns)
        data = FITS_rec.from_columns(coldefs, nrows=nrows, fill=fill)
        hdu = cls(data=data, header=header, **kwargs)
        coldefs._add_listener(hdu)
        return hdu

    @lazyproperty
    def columns(self):
        """
        The :class:`ColDefs` objects describing the columns in this table.
        """

        # The base class doesn't make any assumptions about where the column
        # definitions come from, so just return an empty ColDefs
        return ColDefs([])

    @property
    def _nrows(self):
        """
        Table-like HDUs must provide an attribute that specifies the number of
        rows in the HDU's table.

        For now this is an internal-only attribute.
        """

        raise NotImplementedError

    def _get_tbdata(self):
        """Get the table data from an input HDU object."""

        columns = self.columns

        # TODO: Details related to variable length arrays need to be dealt with
        # specifically in the BinTableHDU class, since they're a detail
        # specific to FITS binary tables
        if (any(type(r) in (_FormatP, _FormatQ)
                for r in columns._recformats) and
                self._data_size is not None and
                self._data_size > self._theap):
            # We have a heap; include it in the raw_data
            raw_data = self._get_raw_data(self._data_size, np.uint8,
                                          self._data_offset)
            data = raw_data[:self._theap].view(dtype=columns.dtype,
                                               type=np.rec.recarray)
        else:
            raw_data = self._get_raw_data(self._nrows, columns.dtype,
                                          self._data_offset)
            if raw_data is None:
                # This can happen when a brand new table HDU is being created
                # and no data has been assigned to the columns, which case just
                # return an empty array
                raw_data = np.array([], dtype=columns.dtype)

            data = raw_data.view(np.rec.recarray)

        self._init_tbdata(data)
        data = data.view(self._data_type)
        columns._add_listener(data)
        return data

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
        for idx, col in enumerate(columns):
            # get the data for each column object from the rec.recarray
            col.array = data.field(idx)

        # delete the _arrays attribute so that it is recreated to point to the
        # new data placed in the column object above
        del columns._arrays

    def _update_column_added(self, columns, column):
        """
        Update the data upon addition of a new column through the `ColDefs`
        interface.
        """

        # TODO: It's not clear that this actually works--it probably does not.
        # This is what the code used to do before introduction of the
        # notifier interface, but I don't believe it actually worked (there are
        # several bug reports related to this...)
        if self._data_loaded:
            del self.data

    def _update_column_removed(self, columns, col_idx):
        """
        Update the data upon removal of a column through the `ColDefs`
        interface.
        """

        # For now this doesn't do anything fancy--it just deletes the data
        # attribute so that it is forced to be recreated again.  It doesn't
        # change anything on the existing data recarray (this is also how this
        # worked before introducing the notifier interface)
        if self._data_loaded:
            del self.data


class _TableBaseHDU(ExtensionHDU, _TableLikeHDU):
    """
    FITS table extension base HDU class.
    """

    _manages_own_heap = False
    """
    This flag implies that when writing VLA tables (P/Q format) the heap
    pointers that go into P/Q table columns should not be reordered or
    rearranged in any way by the default heap management code.

    This is included primarily as an optimization for compressed image HDUs
    which perform their own heap maintenance.
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
            set to `True` if the table contains unsigned integer columns.
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

                # TODO: Too much of the code in this class uses header keywords
                # in making calculations related to the data size.  This is
                # unreliable, however, in cases when users mess with the header
                # unintentionally--code that does this should be cleaned up.
                self._header['NAXIS1'] = self.data._raw_itemsize
                self._header['NAXIS2'] = self.data.shape[0]
                self._header['TFIELDS'] = len(self.data._coldefs)

                self.columns = self.data._coldefs
                self.update()

                with suppress(TypeError, AttributeError):
                    # Make the ndarrays in the Column objects of the ColDefs
                    # object of the HDU reference the same ndarray as the HDU's
                    # FITS_rec object.
                    for idx, col in enumerate(self.columns):
                        col.array = self.data.field(idx)

                    # Delete the _arrays attribute so that it is recreated to
                    # point to the new data placed in the column objects above
                    del self.columns._arrays
            elif data is None:
                pass
            else:
                raise TypeError('Table data has incorrect type.')

        if not (isinstance(self._header[0], string_types) and
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
        """
        The :class:`ColDefs` objects describing the columns in this table.
        """

        if self._has_data and hasattr(self.data, '_coldefs'):
            return self.data._coldefs
        return self._columns_type(self)

    @lazyproperty
    def data(self):
        data = self._get_tbdata()
        data._coldefs = self.columns
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

            with suppress(TypeError, AttributeError):
                # Make the ndarrays in the Column objects of the ColDefs
                # object of the HDU reference the same ndarray as the HDU's
                # FITS_rec object.
                for idx, col in enumerate(self.columns):
                    col.array = self.data.field(idx)

                # Delete the _arrays attribute so that it is recreated to
                # point to the new data placed in the column objects above
                del self.columns._arrays
        elif data is None:
            pass
        else:
            raise TypeError('Table data has incorrect type.')

        # returning the data signals to lazyproperty that we've already handled
        # setting self.__dict__['data']
        return data

    @property
    def _nrows(self):
        if not self._data_loaded:
            return self._header.get('NAXIS2', 0)
        else:
            return len(self.data)

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

        self._header.set('NAXIS1', self.data._raw_itemsize, after='NAXIS')
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
        return self.__class__(data=self.data.copy(),
                              header=self._header.copy())

    def _prewriteto(self, checksum=False, inplace=False):
        if self._has_data:
            self.data._scale_back(
                update_heap_pointers=not self._manages_own_heap)
            # check TFIELDS and NAXIS2
            self._header['TFIELDS'] = len(self.data._coldefs)
            self._header['NAXIS2'] = self.data.shape[0]

            # calculate PCOUNT, for variable length tables
            tbsize = self._header['NAXIS1'] * self._header['NAXIS2']
            heapstart = self._header.get('THEAP', tbsize)
            self.data._gap = heapstart - tbsize
            pcount = self.data._heapsize + self.data._gap
            if pcount > 0:
                self._header['PCOUNT'] = pcount

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
            format = '[{}]'.format(format)
        dims = "{}R x {}C".format(nrows, ncols)
        ncards = len(self._header)

        return (self.name, class_name, ncards, dims, format)

    def _update_column_removed(self, columns, idx):
        super(_TableBaseHDU, self)._update_column_removed(columns, idx)

        # Fix the header to reflect the column removal
        self._clear_table_keywords(index=idx)

    def _update_column_attribute_changed(self, column, col_idx, attr,
                                         old_value, new_value):
        """
        Update the header when one of the column objects is updated.
        """

        # base_keyword is the keyword without the index such as TDIM
        # while keyword is like TDIM1
        base_keyword = ATTRIBUTE_TO_KEYWORD[attr]
        keyword = base_keyword + str(col_idx + 1)

        if keyword in self._header:
            if new_value is None:
                # If the new value is None, i.e. None was assigned to the
                # column attribute, then treat this as equivalent to deleting
                # that attribute
                del self._header[keyword]
            else:
                self._header[keyword] = new_value
        else:
            keyword_idx = KEYWORD_NAMES.index(base_keyword)
            # Determine the appropriate keyword to insert this one before/after
            # if it did not already exist in the header
            for before_keyword in reversed(KEYWORD_NAMES[:keyword_idx]):
                before_keyword += str(col_idx + 1)
                if before_keyword in self._header:
                    self._header.insert(before_keyword, (keyword, new_value),
                                        after=True)
                    break
            else:
                for after_keyword in KEYWORD_NAMES[keyword_idx + 1:]:
                    after_keyword += str(col_idx + 1)
                    if after_keyword in self._header:
                        self._header.insert(after_keyword,
                                            (keyword, new_value))
                        break
                else:
                    # Just append
                    self._header[keyword] = new_value

    def _clear_table_keywords(self, index=None):
        """
        Wipe out any existing table definition keywords from the header.

        If specified, only clear keywords for the given table index (shifting
        up keywords for any other columns).  The index is zero-based.
        Otherwise keywords for all columns.
        """

        # First collect all the table structure related keyword in the header
        # into a single list so we can then sort them by index, which will be
        # useful later for updating the header in a sensible order (since the
        # header *might* not already be written in a reasonable order)
        table_keywords = []

        for idx, keyword in enumerate(self._header.keys()):
            match = TDEF_RE.match(keyword)
            try:
                base_keyword = match.group('label')
            except Exception:
                continue                # skip if there is no match

            if base_keyword in KEYWORD_TO_ATTRIBUTE:
                num = int(match.group('num')) - 1  # convert to zero-base
                table_keywords.append((idx, match.group(0), base_keyword,
                                       num))

        # First delete
        rev_sorted_idx_0 = sorted(table_keywords, key=operator.itemgetter(0),
                                  reverse=True)
        for idx, keyword, _, num in rev_sorted_idx_0:
            if index is None or index == num:
                del self._header[idx]

        # Now shift up remaining column keywords if only one column was cleared
        if index is not None:
            sorted_idx_3 = sorted(table_keywords, key=operator.itemgetter(3))
            for _, keyword, base_keyword, num in sorted_idx_3:
                if num <= index:
                    continue

                old_card = self._header.cards[keyword]
                new_card = (base_keyword + str(num), old_card.value,
                            old_card.comment)
                self._header.insert(keyword, new_card)
                del self._header[keyword]

            # Also decrement TFIELDS
            if 'TFIELDS' in self._header:
                self._header['TFIELDS'] -= 1

    def _populate_table_keywords(self):
        """Populate the new table definition keywords from the header."""

        for idx, column in enumerate(self.columns):
            for keyword, attr in six.iteritems(KEYWORD_TO_ATTRIBUTE):
                val = getattr(column, attr)
                if val is not None:
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
        if isinstance(xtension, string_types):
            xtension = xtension.rstrip()
        return card.keyword == 'XTENSION' and xtension == cls._extension

    def _get_tbdata(self):
        columns = self.columns
        names = [n for idx, n in enumerate(columns.names)]

        # determine if there are duplicate field names and if there
        # are throw an exception
        dup = np.rec.find_duplicate(names)

        if dup:
            raise ValueError("Duplicate field names: {}".format(dup))

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

        raw_data = self._get_raw_data(self._nrows, dtype, self._data_offset)
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
            bytes_array = self.data.view(type=np.ndarray, dtype=np.ubyte)
            padding = np.fromstring(_pad_length(self.size) * b' ',
                                    dtype=np.ubyte)

            d = np.append(bytes_array, padding)

            cs = self._compute_checksum(d, blocking=blocking)
            return cs
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
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
        if isinstance(xtension, string_types):
            xtension = xtension.rstrip()
        return (card.keyword == 'XTENSION' and
                xtension in (cls._extension, 'A3DTABLE'))

    def _calculate_datasum_with_heap(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card given the input data
        """

        with _binary_table_byte_swap(self.data) as data:
            dout = data.view(type=np.ndarray, dtype=np.ubyte)
            csum = self._compute_checksum(dout, blocking=blocking)

            # Now add in the heap data to the checksum (we can skip any gap
            # between the table and the heap since it's all zeros and doesn't
            # contribute to the checksum
            # TODO: The following code may no longer be necessary since it is
            # now possible to get a pointer directly to the heap data as a
            # whole.  That said, it is possible for the heap section to contain
            # data that is not actually pointed to by the table (i.e. garbage;
            # this *shouldn't* happen but it is not disallowed either)--need to
            # double check whether or not the checksum should include such
            # garbage
            for idx in range(data._nfields):
                if isinstance(data.columns._recformats[idx], _FormatP):
                    for coldata in data.field(idx):
                        # coldata should already be byteswapped from the call
                        # to _binary_table_byte_swap
                        if not len(coldata):
                            continue

                        csum = self._compute_checksum(coldata, csum,
                                                      blocking=blocking)

            return csum

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._has_data:
            # This method calculates the datasum while incorporating any
            # heap data, which is obviously not handled from the base
            # _calculate_datasum
            return self._calculate_datasum_with_heap(blocking)
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
            return super(BinTableHDU, self)._calculate_datasum(blocking)

    def _writedata_internal(self, fileobj):
        size = 0

        if self.data is None:
            return size

        with _binary_table_byte_swap(self.data) as data:
            if _has_unicode_fields(data):
                # If the raw data was a user-supplied recarray, we can't write
                # unicode columns directly to the file, so we have to switch
                # to a slower row-by-row write
                self._writedata_by_row(fileobj)
            else:
                fileobj.writearray(data)
                # write out the heap of variable length array columns this has
                # to be done after the "regular" data is written (above)
                fileobj.write((data._gap * '\0').encode('ascii'))

            nbytes = data._gap

            if not self._manages_own_heap:
                # Write the heap data one column at a time, in the order
                # that the data pointers appear in the column (regardless
                # if that data pointer has a different, previous heap
                # offset listed)
                for idx in range(data._nfields):
                    if not isinstance(data.columns._recformats[idx],
                                      _FormatP):
                        continue

                    field = self.data.field(idx)
                    for row in field:
                        if len(row) > 0:
                            nbytes += row.nbytes
                            if not fileobj.simulateonly:
                                fileobj.writearray(row)
            else:
                heap_data = data._get_heap_data()
                if len(heap_data) > 0:
                    nbytes += len(heap_data)
                    if not fileobj.simulateonly:
                        fileobj.writearray(heap_data)

            data._heapsize = nbytes - data._gap
            size += nbytes

        size += self.data.size * self.data._raw_itemsize

        return size

    def _writedata_by_row(self, fileobj):
        fields = [self.data.field(idx)
                  for idx in range(len(self.data.columns))]

        # Creating Record objects is expensive (as in
        # `for row in self.data:` so instead we just iterate over the row
        # indices and get one field at a time:
        for idx in range(len(self.data)):
            for field in fields:
                item = field[idx]
                field_width = None

                if field.dtype.kind == 'U':
                    # Read the field *width* by reading past the field kind.
                    i = field.dtype.str.index(field.dtype.kind)
                    field_width = int(field.dtype.str[i+1:])
                    item = np.char.encode(item, 'ascii')

                fileobj.writearray(item)
                if field_width is not None:
                    j = item.dtype.str.index(item.dtype.kind)
                    item_length = int(item.dtype.str[j+1:])
                    # Fix padding problem (see #5296).
                    padding = '\x00'*(field_width - item_length)
                    fileobj.write(padding.encode('ascii'))

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
              ('Q' format) due to difficult to overcome ambiguities. What this
              means is that this file format cannot support VLA columns in
              tables stored in files that are over 2 GB in size.

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

    @deprecated_renamed_argument('clobber', 'overwrite', '1.3', pending=True)
    def dump(self, datafile=None, cdfile=None, hfile=None, overwrite=False):
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

        overwrite : bool, optional
            If ``True``, overwrite the output file if it exists. Raises an
            ``OSError`` (``IOError`` for Python 2) if ``False`` and the
            output file exists. Default is ``False``.

            .. versionchanged:: 1.3
               ``overwrite`` replaces the deprecated ``clobber`` argument.

        Notes
        -----
        The primary use for the `dump` method is to allow viewing and editing
        the table data and parameters in a standard text editor.
        The `load` method can be used to create a new table from the three
        plain text (ASCII) files.
        """

        # check if the output files already exist
        exist = []
        files = [datafile, cdfile, hfile]

        for f in files:
            if isinstance(f, string_types):
                if os.path.exists(f) and os.path.getsize(f) != 0:
                    if overwrite:
                        warnings.warn(
                            "Overwriting existing file '{}'.".format(f),
                            AstropyUserWarning)
                        os.remove(f)
                    else:
                        exist.append(f)

        if exist:
            raise IOError('  '.join(["File '{}' already exists.".format(f)
                                     for f in exist]))

        # Process the data
        self._dump_data(datafile)

        # Process the column definitions
        if cdfile:
            self._dump_coldefs(cdfile)

        # Process the header parameters
        if hfile:
            self._header.tofile(hfile, sep='\n', endcard=False, padding=False)

    if isinstance(dump.__doc__, string_types):
        dump.__doc__ += _tdump_file_format.replace('\n', '\n        ')

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
            values not present in this Header, unless ``replace=True`` in which
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

    if isinstance(load.__doc__, string_types):
        load.__doc__ += _tdump_file_format.replace('\n', '\n        ')

    load = classmethod(load)
    # Have to create a classmethod from this here instead of as a decorator;
    # otherwise we can't update __doc__

    def _dump_data(self, fileobj):
        """
        Write the table data in the ASCII format read by BinTableHDU.load()
        to fileobj.
        """

        if not fileobj and self._file:
            root = os.path.splitext(self._file.name)[0]
            fileobj = root + '.txt'

        close_file = False

        if isinstance(fileobj, string_types):
            fileobj = open(fileobj, 'w')
            close_file = True

        linewriter = csv.writer(fileobj, dialect=FITSTableDumpDialect)

        # Process each row of the table and output one row at a time
        def format_value(val, format):
            if format[0] == 'S':
                itemsize = int(format[1:])
                return '{:{size}}'.format(val, size=itemsize)
            elif format in np.typecodes['AllInteger']:
                # output integer
                return '{:21d}'.format(val)
            elif format in np.typecodes['Complex']:
                return '{:21.15g}+{:.15g}j'.format(val.real, val.imag)
            elif format in np.typecodes['Float']:
                # output floating point
                # workaround as py2 doesn't support alternate form for format()
                if six.PY2:
                    return '%#21.15g' % val
                else:
                    return '{:#21.15g}'.format(val)

        for row in self.data:
            line = []   # the line for this row of the table

            # Process each column of the row.
            for column in self.columns:
                # format of data in a variable length array
                # where None means it is not a VLA:
                vla_format = None
                format = _convert_format(column.format)

                if isinstance(format, _FormatP):
                    # P format means this is a variable length array so output
                    # the length of the array for this row and set the format
                    # for the VLA data
                    line.append('VLA_Length=')
                    line.append('{:21d}'.format(len(row[column.name])))
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
                    if array_format == 'V':
                        array_format = dtype.base.char
                    if array_format == 'S':
                        array_format += str(dtype.itemsize)

                    if dtype.char == 'V':
                        for value in row[column.name].flat:
                            line.append(format_value(value, array_format))
                    else:
                        line.append(format_value(row[column.name],
                                    array_format))
            linewriter.writerow(line)
        if close_file:
            fileobj.close()

    def _dump_coldefs(self, fileobj):
        """
        Write the column definition parameters in the ASCII format read by
        BinTableHDU.load() to fileobj.
        """

        close_file = False

        if isinstance(fileobj, string_types):
            fileobj = open(fileobj, 'w')
            close_file = True

        # Process each column of the table and output the result to the
        # file one at a time
        for column in self.columns:
            line = [column.name, column.format]
            attrs = ['disp', 'unit', 'dim', 'null', 'bscale', 'bzero']
            line += ['{:16s}'.format(value if value else '""')
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

        if isinstance(fileobj, string_types):
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
        # needing enough physical memory to hold the entire thing at once
        hdu = BinTableHDU.from_columns(np.recarray(shape=1, dtype=dtype),
                                       nrows=nrows, fill=True)

        # TODO: It seems to me a lot of this could/should be handled from
        # within the FITS_rec class rather than here.
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
                name = data.columns.names[idx]
                data._cache_field(name, _makep(arr, arr, recformats[idx]))

        def format_value(col, val):
            # Special formatting for a couple particular data types
            if recformats[col] == FITS2NUMPY['L']:
                return bool(int(val))
            elif recformats[col] == FITS2NUMPY['M']:
                # For some reason, in arrays/fields where numpy expects a
                # complex it's not happy to take a string representation
                # (though it's happy to do that in other contexts), so we have
                # to convert the string representation for it:
                return complex(val)
            else:
                return val

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
                    slice_ = slice(idx, idx + vla_len)
                    data[row][col][:] = line[idx:idx + vla_len]
                    idx += vla_len
                elif dtype[col].shape:
                    # This is an array column
                    array_size = int(np.multiply.reduce(dtype[col].shape))
                    slice_ = slice(idx, idx + array_size)
                    idx += array_size
                else:
                    slice_ = None

                if slice_ is None:
                    # This is a scalar row element
                    data[row][col] = format_value(col, line[idx])
                    idx += 1
                else:
                    data[row][col].flat[:] = [format_value(col, val)
                                              for val in line[slice_]]

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

        if isinstance(fileobj, string_types):
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


@contextlib.contextmanager
def _binary_table_byte_swap(data):
    """
    Ensures that all the data of a binary FITS table (represented as a FITS_rec
    object) is in a big-endian byte order.  Columns are swapped in-place one
    at a time, and then returned to their previous byte order when this context
    manager exits.

    Because a new dtype is needed to represent the byte-swapped columns, the
    new dtype is temporarily applied as well.
    """

    orig_dtype = data.dtype

    names = []
    formats = []
    offsets = []

    to_swap = []

    if sys.byteorder == 'little':
        swap_types = ('<', '=')
    else:
        swap_types = ('<',)

    for idx, name in enumerate(orig_dtype.names):
        field = _get_recarray_field(data, idx)

        field_dtype, field_offset = orig_dtype.fields[name]
        names.append(name)
        formats.append(field_dtype)
        offsets.append(field_offset)

        if isinstance(field, chararray.chararray):
            continue

        # only swap unswapped
        # must use field_dtype.base here since for multi-element dtypes,
        # the .str with be '|V<N>' where <N> is the total bytes per element
        if field.itemsize > 1 and field_dtype.base.str[0] in swap_types:
            to_swap.append(field)
            # Override the dtype for this field in the new record dtype with
            # the byteswapped version
            formats[-1] = field_dtype.newbyteorder()

        # deal with var length table
        recformat = data.columns._recformats[idx]
        if isinstance(recformat, _FormatP):
            coldata = data.field(idx)
            for c in coldata:
                if (not isinstance(c, chararray.chararray) and
                        c.itemsize > 1 and c.dtype.str[0] in swap_types):
                    to_swap.append(c)

    for arr in reversed(to_swap):
        arr.byteswap(True)

    new_dtype = nh.realign_dtype(np.dtype(list(zip(names, formats))),
                                 offsets)

    data.dtype = new_dtype

    yield data

    for arr in to_swap:
        arr.byteswap(True)

    data.dtype = orig_dtype
