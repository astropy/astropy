from __future__ import division # confidence high

import csv
import os
import re
import sys
import textwrap
import warnings

import numpy as np
from numpy import char as chararray

from pyfits.card import Card, CardList
# This module may have many dependencies on pyfits.column, but pyfits.column
# has fewer dependencies overall, so it's easier to keep table/column-related
# utilities in pyfits.column
from pyfits.column import (FITS2NUMPY, KEYWORD_NAMES, KEYWORD_ATTRIBUTES,
                           TDEF_RE, Delayed, Column, ColDefs, _ASCIIColDefs,
                           _FormatX, _FormatP, _wrapx, _makep, _VLF,
                           _parse_tformat, _scalar_to_format, _convert_format,
                           _cmp_recformats)
from pyfits.fitsrec import FITS_rec
from pyfits.hdu.base import DELAYED, _ValidHDU, ExtensionHDU
from pyfits.header import Header
from pyfits.util import (lazyproperty, _is_int, _str_to_num, _pad_length,
                         deprecated)


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

        # TODO: Need to find a way to eliminate the check for phantom columns;
        # this detail really needn't be worried about outside the ColDefs class
        columns = self.columns
        recformats = [f for idx, f in enumerate(columns._recformats)
                      if not columns[idx]._phantom]
        formats = ','.join(recformats)
        names = [n for idx, n in enumerate(columns.names)
                 if not columns[idx]._phantom]
        dtype = np.rec.format_parser(formats, names, None).dtype
        raw_data = self._file.readarray(offset=self._datLoc, dtype=dtype,
                                        shape=columns._shape)
        data = raw_data.view(np.rec.recarray)
        self._init_tbdata(data)
        return data.view(FITS_rec)

    def _init_tbdata(self, data):
        columns = self.columns

        data.dtype = data.dtype.newbyteorder('>')

        # pass datLoc, for P format
        data._heapoffset = self._theap + self._datLoc
        data._file = self._file
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

    def __init__(self, data=None, header=None, name=None):
        """
        Parameters
        ----------
        header : Header instance
            header to be used

        data : array
            data to be used

        name : str
            name to be populated in ``EXTNAME`` keyword
        """

        super(_TableBaseHDU, self).__init__(data=data, header=header,
                                            name=name)

        if header is not None and not isinstance(header, Header):
            raise ValueError('header must be a Header object.')

        if data is DELAYED:
            # this should never happen
            if header is None:
                raise ValueError('No header to setup HDU.')

            # if the file is read the first time, no need to copy, and keep it unchanged
            else:
                self._header = header
        else:
            # construct a list of cards of minimal header
            cards = CardList([
                Card('XTENSION',      '', ''),
                Card('BITPIX',         8, 'array data type'),
                Card('NAXIS',          2, 'number of array dimensions'),
                Card('NAXIS1',         0, 'length of dimension 1'),
                Card('NAXIS2',         0, 'length of dimension 2'),
                Card('PCOUNT',         0, 'number of group parameters'),
                Card('GCOUNT',         1, 'number of groups'),
                Card('TFIELDS',        0, 'number of table fields')
                ])

            if header is not None:
                # Make a "copy" (not just a view) of the input header, since it
                # may get modified.  the data is still a "view" (for now)
                hcopy = header.copy(strip=True)
                cards.extend(hcopy.ascard)

            self._header = Header(cards)

            if isinstance(data, np.ndarray) and data.dtype.fields is not None:
                if isinstance(data, FITS_rec):
                    self.data = data
                else:
                    self.data = data.view(FITS_rec)

                self._header['NAXIS1'] = self.data.itemsize
                self._header['NAXIS2'] = self.data.shape[0]
                self._header['TFIELDS'] = self.data._nfields

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
                except (TypeError, AttributeError), e:
                    # This shouldn't happen as long as self.columns._arrays
                    # is a lazyproperty
                    pass
            elif data is None:
                pass
            else:
                raise TypeError('Table data has incorrect type.')

        if self._header[0].rstrip() != self._extension:
            self._header[0] = self._extension
            self._header.ascard[0].comment = self._ext_comment

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
        if self._data_loaded and hasattr(self.data, '_coldefs'):
            return self.data._coldefs
        return ColDefs(self)

    @lazyproperty
    def data(self):
        data = self._get_tbdata()
        data._coldefs = self.columns
        data.formats = self.columns.formats
        # Columns should now just return a reference to the data._coldefs
        del self.columns
        return data

    @lazyproperty
    def _theap(self):
        size = self._header['NAXIS1'] * self._header['NAXIS2']
        return self._header.get('THEAP', size)

    @deprecated(alternative='the .columns attribute')
    def get_coldefs(self):
        """
        **[Deprecated]** Returns the table's column definitions.
        """

        return self.columns

    def update(self):
        """
        Update header keywords to reflect recent changes of columns.
        """

        update = self._header.update
        update('naxis1', self.data.itemsize, after='naxis')
        update('naxis2', self.data.shape[0], after='naxis1')
        update('tfields', len(self.columns), after='gcount')

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

    def _writeheader(self, fileobj, checksum=False):
        if self._data_loaded and self.data is not None:
            self.data._scale_back()
            # check TFIELDS and NAXIS2
            self._header['TFIELDS'] = self.data._nfields
            self._header['NAXIS2'] = self.data.shape[0]

            # calculate PCOUNT, for variable length tables
            tbsize = self.header['NAXIS1'] * self.header['NAXIS2']
            heapstart = self.header.get('THEAP', tbsize)
            self.data._gap = heapstart - tbsize
            pcount = self.data._heapsize + self.data._gap
            if pcount > 0:
                self.header['PCOUNT'] = pcount

            # update TFORM for variable length columns
            for idx in range(self.data._nfields):
                if isinstance(self.data._coldefs.formats[idx], _FormatP):
                    key = 'TFORM' + str(idx + 1)
                    val = self._header[key]
                    val = val[:val.find('(') + 1] + \
                          repr(self.data.field(idx)._max) + ')'
                    self._header[key] = val
        return super(_TableBaseHDU, self)._writeheader(fileobj, checksum)

    def _writeto(self, fileobj, checksum=False):
        return super(_TableBaseHDU, self)._writeto(fileobj, checksum)

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

            ncols = len(self.columns.formats)
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
        ncards = len(self._header.ascard)

        return (self.name, class_name, ncards, dims, format)

    def _clear_table_keywords(self):
        """Wipe out any existing table definition keywords from the header."""

        # Go in reverse so as to not confusing indexing while deleting.
        for idx, card in enumerate(reversed(self._header.ascard)):
            key = TDEF_RE.match(card.key)
            try:
                keyword = key.group('label')
            except:
                continue                # skip if there is no match
            if (keyword in KEYWORD_NAMES):
                del self._header.ascard[idx]

    def _populate_table_keywords(self):
        """Populate the new table definition keywords from the header."""

        cols = self.columns
        append = self._header.ascard.append

        for idx, col in enumerate(cols):
            for attr, keyword in zip(KEYWORD_ATTRIBUTES, KEYWORD_NAMES):
                val = getattr(cols, attr + 's')[idx]
                if val:
                    keyword = keyword + str(idx + 1)
                    append(Card(keyword, val))


class TableHDU(_TableBaseHDU):
    """
    FITS ASCII table extension HDU class.
    """

    _extension = 'TABLE'
    _ext_comment = 'ASCII table extension'

    _padding_byte = ' '

    __format_RE = re.compile(
        r'(?P<code>[ADEFIJ])(?P<width>\d+)(?:\.(?P<prec>\d+))?')

    def __init__(self, data=None, header=None, name=None):
        super(TableHDU, self).__init__(data, header, name=name)
        if self._data_loaded and self.data is not None and \
           not isinstance(self.data._coldefs, _ASCIIColDefs):
            self.data._coldefs = _ASCIIColDefs(self.data._coldefs)

    @classmethod
    def match_header(cls, header):
        card = header.ascard[0]
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        return card.key == 'XTENSION' and xtension == cls._extension

    def _get_tbdata(self):
        columns = self.columns
        names = [n for idx, n in enumerate(columns.names)
                 if not columns[idx]._phantom]

        # determine if there are duplicate field names and if there
        # are throw an exception
        dup = np.rec.find_duplicate(names)

        if dup:
            raise ValueError("Duplicate field names: %s" % dup)

        itemsize = columns.spans[-1] + columns.starts[-1]-1
        dtype = {}

        for idx in range(len(columns)):
            data_type = 'S' + str(columns.spans[idx])

            if idx == len(columns) - 1:
                # The last column is padded out to the value of NAXIS1
                if self._header['NAXIS1'] > itemsize:
                    data_type = 'S' + str(columns.spans[idx] + \
                                self._header['NAXIS1'] - itemsize)
            dtype[columns.names[idx]] = (data_type, columns.starts[idx] - 1)

        raw_data = self._file.readarray(offset=self._datLoc, dtype=dtype,
                                        shape=columns._shape)
        data = raw_data.view(np.rec.recarray)
        self._init_tbdata(data)
        return data.view(FITS_rec)

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._data_loaded and self.data is not None:
            # We have the data to be used.
            # We need to pad the data to a block length before calculating
            # the datasum.

            if self.size > 0:
                d = np.append(np.fromstring(self.data, dtype='ubyte'),
                              np.fromstring(_pad_length(self.size) * ' ',
                                            dtype='ubyte'))

            cs = self._compute_checksum(np.fromstring(d, dtype='ubyte'),
                                        blocking=blocking)
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
        card = header.ascard[0]
        xtension = card.value
        if isinstance(xtension, basestring):
            xtension = xtension.rstrip()
        return card.key == 'XTENSION' and \
               xtension in (cls._extension, 'A3DTABLE')

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

        dout = np.fromstring(data, dtype='ubyte')

        for i in range(data._nfields):
            if isinstance(data._coldefs._recformats[i], _FormatP):
                for coldata in data.field(i):
                    if len(coldata) > 0:
                        dout = np.append(dout,
                                         np.fromstring(coldata, dtype='ubyte'))

        cs = self._compute_checksum(dout, blocking=blocking)
        return cs

    def _calculate_datasum(self, blocking):
        """
        Calculate the value for the ``DATASUM`` card in the HDU.
        """

        if self._data_loaded and self.data is not None:
            # We have the data to be used.
            return self._calculate_datasum_from_data(self.data, blocking)
        else:
            # This is the case where the data has not been read from the file
            # yet.  We can handle that in a generic manner so we do it in the
            # base class.  The other possibility is that there is no data at
            # all.  This can also be handled in a generic manner.
            return super(BinTableHDU,self)._calculate_datasum(blocking)

    def _writedata_internal(self, fileobj):
        size = 0

        if self.data is not None:
            size += self._binary_table_byte_swap(fileobj)
            size += self.data.size * self.data.itemsize

        return size

    def _binary_table_byte_swap(self, fileobj):
        swapped = []
        nbytes = 0
        if sys.byteorder == 'little':
            swap_types = ('<', '=')
        else:
            swap_types = ('<',)
        try:
            if not fileobj.simulateonly:
                for idx in range(self.data._nfields):
                    coldata = self.data.field(idx)
                    if isinstance(coldata, chararray.chararray):
                        continue
                    # only swap unswapped
                    # deal with var length table
                    field = np.rec.recarray.field(self.data, idx)
                    if isinstance(coldata, _VLF):
                        for jdx, c in enumerate(coldata):
                            if (not isinstance(c, chararray.chararray) and
                                c.itemsize > 1 and
                                c.dtype.str[0] in swap_types):
                                swapped.append(c)
                            if (field[jdx:jdx+1].dtype.str[0] in swap_types):
                                swapped.append(field[jdx:jdx+1])
                    else:
                        if (coldata.itemsize > 1 and
                            self.data.dtype.descr[idx][1][0] in swap_types):
                            swapped.append(field)

                for obj in swapped:
                    obj.byteswap(True)

                fileobj.writearray(self.data)

                # write out the heap of variable length array
                # columns this has to be done after the
                # "regular" data is written (above)
                fileobj.write(
                    (self.data._gap * '\0').encode('ascii'))

            nbytes = self.data._gap

            for idx in range(self.data._nfields):
                if isinstance(self.data._coldefs._recformats[idx], _FormatP):
                    field = self.data.field(idx)
                    for jdx in range(len(field)):
                        coldata = field[jdx]
                        if len(coldata) > 0:
                            nbytes = nbytes + coldata.nbytes
                            if not fileobj.simulateonly:
                                fileobj.writearray(coldata)

            self.data._heapsize = nbytes - self.data._gap
        finally:
            for obj in swapped:
                obj.byteswap(True)

        return nbytes

    def _populate_table_keywords(self):
        """Populate the new table definition keywords from the header."""

        cols = self.columns
        append = self._header.ascard.append

        for idx, col in enumerate(cols):
            for attr, keyword in zip(KEYWORD_ATTRIBUTES, KEYWORD_NAMES):
                val = getattr(cols, attr + 's')[idx]
                if val:
                    keyword = keyword + str(idx + 1)
                    if attr == 'format':
                        val = cols._recformats[idx]
                        if isinstance(val, _FormatX):
                            val = repr(val._nx) + 'X'
                        elif isinstance(val, _FormatP):
                            VLdata = self.data.field(idx)
                            VLdata._max = max(map(len, VLdata))
                            if val._dtype == 'a':
                                fmt = 'A'
                            else:
                                fmt = _convert_format(val._dtype, reverse=True)
                            val = 'P%s(%d)' % (fmt, VLdata._max)
                        else:
                            val = _convert_format(val, reverse=True)
                    append(Card(keyword, val))

    tdump_file_format = textwrap.dedent("""

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
                        warnings.warn("Overwriting existing file '%s'." % f)
                    else:
                        exist.append(f)

        if exist:
            raise IOError('  '.join(["File '%s' already exists."
                                     for f in exist]))

        # Process the data
        self._dump_data(datafile)

        # Process the column definitions
        if cdfile:
            self._dump_coldefs(cdfile)

        # Process the header parameters
        if hfile:
            self._header.toTxtFile(hfile)

    dump.__doc__ += tdump_file_format.replace('\n', '\n        ')
    tdump = deprecated(name='tdump')(dump)

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
            header.fromTxtFile(hfile, replace)

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
    load.__doc__ += tdump_file_format.replace('\n', '\n        ')
    load = classmethod(load)
    # Have to create a classmethod from this here instead of as a decorator;
    # otherwise we can't update __doc__
    tcreate = deprecated(name='tcreate')(load)

    def _dump_data(self, fileobj):
        """
        Write the table data in the ASCII format read by BinTableHDU.load()
        to fileobj.
        """

        if not fileobj and self._file:
            root, ext = os.path.splitext(self._file.name)
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
                    repeat, dtype, option = _parse_tformat(column.format)
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

        initialpos = fileobj.tell() # We'll be returning here later
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
                recformat = _FormatP(str(length) + recformats[idx])
                # TODO: I shouldn't have to set this magic attribute, but
                # _convert_fits2record does it and apparently _makep requires
                # this.  This spaghetti code needs fixing
                recformat._dtype = recformats[idx]
                recformats[idx] = recformat

        dtype = np.rec.format_parser(recformats, names, None).dtype

        # TODO: In the future maybe enable loading a bit at a time so that we
        # can convert from this format to an actual FITS file on disk without
        # needing enough physical memory to hold the entire thing at once;
        # new_table() could use a similar feature.
        hdu = new_table(np.recarray(shape=1, dtype=dtype), nrows=nrows,
                        fill=True)
        data = hdu.data
        # Unfortunately, the fact that some columns are VLAs is not preserved
        # through the call to new_table(), so we need to fix that up after the
        # fact.
        # TODO: new_table really needs to be fixed so that this mess is
        # unnecessary
        for idx, recformat in enumerate(recformats):
            if not isinstance(recformat, _FormatP):
                continue
            arr = data.columns._arrays[idx]
            hdu.data._convert[idx] = _makep(arr, arr, recformat._dtype)

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


# TODO: Allow tbtype to be either a string or a class; perhaps eventually
# replace this with separate functions for creating tables (possibly in the
# form of a classmethod)  See ticket #60
def new_table(input, header=None, nrows=0, fill=False, tbtype='BinTableHDU'):
    """
    Create a new table from the input column definitions.

    Warning: Creating a new table using this method creates an in-memory *copy*
    of all the column arrays in the input.  This is because if they are
    separate arrays they must be combined into a single contiguous array.

    If the column data is already in a single contiguous array (such as an
    existing record array) it may be better to create a BinTableHDU instance
    directly.  See the PyFITS documentation for more details.

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

    tbtype : str
        Table type to be created ("BinTableHDU" or "TableHDU").
    """

    # construct a table HDU
    # TODO: Something needs to be done about this as part of #60....
    hdu = eval(tbtype)(header=header)

    if isinstance(input, ColDefs):
        # NOTE: This previously raised an error if the tbtype didn't match the
        # tbtype of the input ColDefs. This should no longer be necessary, but
        # just beware.
        columns = hdu.columns = ColDefs(input)
    elif isinstance(input, FITS_rec): # input is a FITS_rec
        # Create a new ColDefs object from the input FITS_rec's ColDefs
        # object and assign it to the ColDefs attribute of the new hdu.
        columns = hdu.columns = ColDefs(input._coldefs, tbtype)
    else: # input is a list of Columns or possibly a recarray
        # Create a new ColDefs object from the input list of Columns and
        # assign it to the ColDefs attribute of the new hdu.
        columns = hdu.columns = ColDefs(input, tbtype)

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
            if (arr is not None):
                dim = arr.shape[0]
            else:
                dim = 0
            if dim > nrows:
                nrows = dim

    if tbtype == 'TableHDU':
        columns = hdu.columns = _ASCIIColDefs(hdu.columns)
        _itemsize = columns.spans[-1] + columns.starts[-1]-1
        dtype = {}

        for j in range(len(columns)):
           data_type = 'S' + str(columns.spans[j])
           dtype[columns.names[j]] = (data_type, columns.starts[j] - 1)

        hdu.data = np.rec.array((' ' * _itemsize * nrows).encode('ascii'),
                                dtype=dtype, shape=nrows).view(FITS_rec)
        hdu.data.setflags(write=True)
    else:
        formats = ','.join(columns._recformats)
        hdu.data = np.rec.array(None, formats=formats,
                                names=columns.names,
                                shape=nrows).view(FITS_rec)

    hdu.data._coldefs = hdu.columns
    hdu.data.formats = hdu.columns.formats

    # Populate data to the new table from the ndarrays in the input ColDefs
    # object.
    for idx in range(len(columns)):
        # For each column in the ColDef object, determine the number
        # of rows in that column.  This will be either the number of
        # rows in the ndarray associated with the column, or the
        # number of rows given in the call to this function, which
        # ever is smaller.  If the input FILL argument is true, the
        # number of rows is set to zero so that no data is copied from
        # the original input data.
        arr = columns._arrays[idx]
        recformat = columns._recformats[idx]

        if arr is None:
            size = 0
        else:
            size = len(arr)

        n = min(size, nrows)
        if fill:
            n = 0

        # Get any scale factors from the FITS_rec
        scale, zero, bscale, bzero, dim = hdu.data._get_scale_factors(idx)[3:]

        field = np.rec.recarray.field(hdu.data, idx)

        if n > 0:
            # Only copy data if there is input data to copy
            # Copy all of the data from the input ColDefs object for this
            # column to the new FITS_rec data array for this column.
            if isinstance(recformat, _FormatX):
                # Data is a bit array
                if arr[:n].shape[-1] == recformat._nx:
                    _wrapx(arr[:n], field[:n], recformat._nx)
                else: # from a table parent data, just pass it
                    field[:n] = arr[:n]
            elif isinstance(recformat, _FormatP):
                hdu.data._convert[idx] = _makep(arr[:n], field,
                                                recformat._dtype, nrows=nrows)
            elif recformat[-2:] == FITS2NUMPY['L'] and arr.dtype == bool:
                # column is boolean
                field[:n] = np.where(arr == False, ord('F'), ord('T'))
            else:
                if tbtype == 'TableHDU':
                    # string no need to convert,
                    if isinstance(arr, chararray.chararray):
                        field[:n] = arr[:n]
                    else:
                        hdu.data._convert[idx] = \
                                np.zeros(nrows, dtype=arr.dtype)
                        if scale or zero:
                            arr = arr.copy()
                        if scale:
                            arr *= bscale
                        if zero:
                            arr += bzero
                        hdu.data._convert[idx][:n] = arr[:n]
                else:
                    field[:n] = arr[:n]

        if n < nrows:
            # If there are additional rows in the new table that were not
            # copied from the input ColDefs object, initialize the new data
            if tbtype == 'BinTableHDU':
                if isinstance(field, np.ndarray):
                    field[n:] = -bzero/bscale
                else:
                    field[n:] = ''
            else:
                field[n:] = ' ' * hdu.data._coldefs.spans[idx]

    # Update the HDU header to match the data
    hdu.update()

    # Make the ndarrays in the Column objects of the ColDefs object of the HDU
    # reference the same ndarray as the HDU's FITS_rec object.
    for idx in range(len(columns)):
        hdu.columns[idx].array = hdu.data.field(idx)

    # Delete the _arrays attribute so that it is recreated to point to the
    # new data placed in the column objects above
    del hdu.columns._arrays

    return hdu
