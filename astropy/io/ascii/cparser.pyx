# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
cimport numpy as np
from numpy import ma
from ...utils.data import get_readable_fileobj
from ...extern import six
from . import core
from distutils import version
import csv

cdef extern from "src/tokenizer.h":
    ctypedef enum tokenizer_state:
        START_LINE
        START_FIELD
        START_QUOTED_FIELD
        FIELD
        QUOTED_FIELD
        QUOTED_FIELD_NEWLINE
        COMMENT

    ctypedef enum err_code:
        NO_ERROR
        INVALID_LINE
        TOO_MANY_COLS
        NOT_ENOUGH_COLS
        CONVERSION_ERROR
        OVERFLOW_ERROR

    ctypedef struct tokenizer_t:
        char *source           # single string containing all of the input
        int source_len         # length of the input
        int source_pos         # current index in source for tokenization
        char delimiter         # delimiter character
        char comment           # comment character
        char quotechar         # quote character
        char *header_output    # string containing header data
        char **output_cols     # array of output strings for each column
        char **col_ptrs        # array of pointers to current output position for each col
        int *output_len        # length of each output column string
        int header_len         # length of the header output string
        int num_cols           # number of table columns
        int num_rows           # number of table rows
        int fill_extra_cols    # represents whether or not to fill rows with too few values
        tokenizer_state state  # current state of the tokenizer
        err_code code          # represents the latest error that has occurred
        int iter_col           # index of the column being iterated over
        char *curr_pos         # current iteration position
        char *buf              # buffer for misc. data
        int strip_whitespace_lines  # whether to strip whitespace at the beginning and end of lines
        int strip_whitespace_fields # whether to strip whitespace at the beginning and end of fields
        # Example input/output
        # --------------------
        # source: "A,B,C\n10,5.,6\n1,2,3"
        # output_cols: ["A\x0010\x001", "B\x005.\x002", "C\x006\x003"]

    tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar, int fill_extra_cols,
                                  int strip_whitespace_lines, int strip_whitespace_fields)
    void delete_tokenizer(tokenizer_t *tokenizer)
    int tokenize(tokenizer_t *self, int start, int end, int header, int *use_cols, int use_cols_len)
    int int_size()
    long str_to_long(tokenizer_t *self, char *str)
    double str_to_double(tokenizer_t *self, char *str)
    void start_iteration(tokenizer_t *self, int col)
    int finished_iteration(tokenizer_t *self)
    char *next_field(tokenizer_t *self)

class CParserError(Exception):
    """
    An instance of this class is thrown when an error occurs
    during C parsing.
    """

ERR_CODES = dict(enumerate([
    "no error",
    "invalid line supplied",
    lambda line: "too many columns found in line {0} of data".format(line),
    lambda line: "not enough columns found in line {0} of data".format(line),
    "type conversion error"
    ]))

cdef class CParser:
    """
    A fast Cython parser class which uses underlying C code
    for tokenization.
    """

    cdef:
        tokenizer_t *tokenizer
        int data_start
        int data_end
        object data_end_obj
        object include_names
        object exclude_names
        object fill_values
        object fill_include_names
        object fill_exclude_names
        object fill_names
        int fill_extra_cols
        np.ndarray use_cols

    cdef public:
        int width
        object names
        bytes source
        object header_start

    def __cinit__(self, source, strip_line_whitespace, strip_line_fields,
                  delimiter=',',
                  comment=None,
                  quotechar='"',
                  header_start=0,
                  data_start=1,
                  data_end=None,
                  names=None,
                  include_names=None,
                  exclude_names=None,
                  fill_values=('', '0'),
                  fill_include_names=None,
                  fill_exclude_names=None,
                  fill_extra_cols=0):

        self.tokenizer = create_tokenizer(ord(delimiter), ord(comment), ord(quotechar), fill_extra_cols,
                                          strip_line_whitespace, strip_line_fields)
        self.source = None
        self.setup_tokenizer(source)
        self.header_start = header_start
        self.data_start = data_start
        self.data_end = -1 # keep reading data until the end
        if data_end is not None and data_end >= 0:
            self.data_end = data_end
        self.data_end_obj = data_end
        self.names = names
        self.include_names = include_names
        self.exclude_names = exclude_names
        self.fill_values = get_fill_values(fill_values)
        self.fill_include_names = fill_include_names
        self.fill_exclude_names = fill_exclude_names
        self.fill_extra_cols = fill_extra_cols
    
    def __dealloc__(self):
        if self.tokenizer:
            delete_tokenizer(self.tokenizer) # perform C memory cleanup

    cdef raise_error(self, msg):
        err_msg = ERR_CODES.get(self.tokenizer.code, "unknown error")

        if callable(err_msg): # error code is lambda function taking current line as input
            err_msg = err_msg(self.tokenizer.num_rows + 1)

        raise CParserError("{0}: {1}".format(msg, err_msg))

    cpdef setup_tokenizer(self, source):
        cdef char *src

        if isinstance(source, six.string_types) or hasattr(source, 'read'):
            if hasattr(source, 'read') or '\n' not in source: # Either filename or file-like object
                with get_readable_fileobj(source) as file_obj:
                    source = file_obj.read()
            # Otherwise, source is the actual data so we leave it be
        else:
            try:
                source = '\n'.join(source) # iterable sequence of lines
            except TypeError:
                raise TypeError('Input "table" must be a file-like object, a string (filename'
                             'or data), or an iterable')
        # Create a reference to the Python object so its char * pointer remains valid
        source_str = source + '\n' # add newline to simplify handling last line of data
        self.source = source_str.encode('ascii') # encode in UTF-8 for char * handling (fixes Python 3 issue)
        src = self.source
        self.tokenizer.source = src
        self.tokenizer.source_len = len(self.source)

    def read_header(self):
        if self.names:
            self.width = len(self.names)

        # header_start is a valid line number
        elif self.header_start is not None and self.header_start >= 0:
            if tokenize(self.tokenizer, self.header_start, -1, 1, <int *> 0, 0) != 0:
                self.raise_error("an error occurred while tokenizing the header line")
            self.names = []
            name = ''
            for i in range(self.tokenizer.header_len):
                c = self.tokenizer.header_output[i] # next char in header string
                if not c: # zero byte -- field terminator
                    if name:
                        self.names.append(name.replace('\x01', '')) # replace empty placeholder with ''
                        name = ''
                    else:
                        break # end of string
                else:
                    name += chr(c)
            self.width = len(self.names)

        else:
            # Get number of columns from first data row
            if tokenize(self.tokenizer, 0, -1, 1, <int *> 0, 0) != 0:
                self.raise_error("an error occurred while tokenizing the first line of data")
            self.width = 0
            for i in range(self.tokenizer.header_len):
                if not self.tokenizer.header_output[i]: # zero byte -- field terminator
                    if i > 0 and self.tokenizer.header_output[i - 1]: # ends valid field
                        self.width += 1
                    else: # end of line
                        break
            if self.width == 0: # no data
                raise core.InconsistentTableError('No data lines found, C reader cannot autogenerate '
                                                  'column names')
            self.names = ['col{0}'.format(i + 1) for i in range(self.width)] # auto-generate names

        size = int_size()
        dtype = np.int16 #TODO: maybe find a better way to do this?
        if size == 64:
            dtype = np.int64
        elif size == 32:
            dtype = np.int32
        self.use_cols = np.ones(self.width, dtype) # "boolean" array denoting whether or not to use each column
        if self.include_names is not None:
            for i, name in enumerate(self.names):
                if name not in self.include_names:
                    self.use_cols[i] = 0
        if self.exclude_names is not None:
            for name in self.exclude_names:
                try:
                    self.use_cols[self.names.index(name)] = 0
                except ValueError: # supplied name is invalid, ignore
                    continue

        # self.names should only contain columns included in output
        self.names = [self.names[i] for i, should_use in enumerate(self.use_cols) if should_use]
        self.width = len(self.names)
        self.tokenizer.num_cols = self.width
            
    def read(self, try_int, try_float, try_string):
        self._tokenize_data()
        self._set_fill_names()
        return self._convert_data(try_int, try_float, try_string)

    cdef _tokenize_data(self):
        if tokenize(self.tokenizer, self.data_start, self.data_end, 0, <int *> self.use_cols.data,
                    len(self.use_cols)) != 0:
            self.raise_error("an error occurred while tokenizing data")

    cdef _set_fill_names(self):
        self.fill_names = set(self.names)
        if self.fill_include_names is not None:
            self.fill_names.intersection_update(self.fill_include_names)
        if self.fill_exclude_names is not None:
            self.fill_names.difference_update(self.fill_exclude_names)

    cdef _convert_data(self, try_int, try_float, try_string):
        cdef int num_rows = self.tokenizer.num_rows
        if self.data_end_obj is not None and self.data_end_obj < 0:
            num_rows += self.data_end_obj # e.g. if data_end = -1, ignore the last row
        cols = {}

        for i, name in enumerate(self.names):
            # Try int first, then float, then string
            try:
                if try_int and not try_int[name]:
                    raise ValueError()
                cols[name] = self._convert_int(i, num_rows)
            except ValueError:
                try:
                    if try_float and not try_float[name]:
                        raise ValueError()
                    cols[name] = self._convert_float(i, num_rows)
                except ValueError:
                    if try_string and not try_string[name]:
                        raise ValueError('Column {0} failed to convert'.format(name))
                    cols[name] = self._convert_str(i, num_rows)

        return cols

    cdef np.ndarray _convert_int(self, int i, int num_rows):
        cdef np.ndarray col = np.empty(num_rows, dtype=np.int_) # intialize ndarray
        cdef long converted
        cdef int row = 0
        cdef int *data = <int *> col.data # pointer to raw data
        cdef bytes field
        cdef bytes new_val
        mask = set() # set of indices for masked values
        start_iteration(self.tokenizer, i) # begin the iteration process in C

        while not finished_iteration(self.tokenizer):
            if row == num_rows: # end prematurely if we aren't using every row
                break
            field = next_field(self.tokenizer) # retrieve the next field in a bytes value

            if field in self.fill_values:
                new_val = str(self.fill_values[field][0]).encode('utf-8')

                # Either this column applies to the field as specified in the fill_values parameter,
                # or no specific columns are specified and this column should apply fill_values
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) or \
                           (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    converted = str_to_long(self.tokenizer, new_val) # try converting the new value
                else:
                    converted = str_to_long(self.tokenizer, field)

            else:
                converted = str_to_long(self.tokenizer, field) # convert the field to long (widest integer type)

            if self.tokenizer.code in (CONVERSION_ERROR, OVERFLOW_ERROR): # no dice
                self.tokenizer.code = NO_ERROR
                raise ValueError()
            col[row] = converted
            row += 1

        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in range(row)]) # convert to masked_array
        else:
            return col

    cdef np.ndarray _convert_float(self, int i, int num_rows):
        # very similar to _convert_int()
        cdef np.ndarray col = np.empty(num_rows, dtype=np.float_)
        cdef double converted
        cdef int row = 0
        cdef float *data = <float *> col.data
        cdef bytes field
        cdef bytes new_val
        mask = set()

        start_iteration(self.tokenizer, i)
        while not finished_iteration(self.tokenizer):
            if row == num_rows:
                break
            field = next_field(self.tokenizer)
            if field in self.fill_values:
                new_val = str(self.fill_values[field][0]).encode('utf-8')
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) or \
                           (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    converted = str_to_double(self.tokenizer, new_val)
                else:
                    converted = str_to_double(self.tokenizer, field)
            else:
                converted = str_to_double(self.tokenizer, field)

            if self.tokenizer.code == CONVERSION_ERROR:
                self.tokenizer.code = NO_ERROR
                raise ValueError()
            elif self.tokenizer.code == OVERFLOW_ERROR:
                self.tokenizer.code = NO_ERROR
                # In numpy < 1.6, using type inference yields a float for overflow values because
                # the error raised is not specific. This replicates the old reading behavior
                # (see #2234).
                if version.LooseVersion(np.__version__) < version.LooseVersion('1.6'):
                    col[row] = new_val if field in self.fill_values else field
                else:
                    raise ValueError()
            else:
                col[row] = converted
            row += 1

        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in range(row)])
        else:
            return col

    cdef np.ndarray _convert_str(self, int i, int num_rows):
        # similar to _convert_int, but no actual conversion
        cdef np.ndarray col = np.empty(num_rows, dtype=object) # TODO: find a faster method here
        cdef int row = 0
        cdef bytes field
        cdef bytes new_val
        cdef int max_len = 0 # greatest length of any element
        mask = set()

        start_iteration(self.tokenizer, i)
        while not finished_iteration(self.tokenizer):
            if row == num_rows:
                break
            field = next_field(self.tokenizer)
            if field in self.fill_values:
                el = str(self.fill_values[field][0])
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) or \
                           (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                else:
                    el = field.decode('utf-8')
            else:
                el = field.decode('utf-8')
            max_len = max(max_len, len(el)) # update max_len with the length of each field
            col[row] = el
            row += 1

        col = col.astype('|S{0}'.format(max_len)) # convert to string with smallest length possible
        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in range(row)])
        else:
            return col

def write(table, output, header_output, output_types,
          delimiter=',',
          comment='# ',
          quotechar='"',
          formats=None,
          strip_whitespace=True,
          names=None, # ignore, already used in _get_writer
          include_names=None,
          exclude_names=None,
          fill_values=[],
          fill_include_names=None,
          fill_exclude_names=None):

    opened_file = False

    if not hasattr(output, 'write'): # output is a filename
        output = open(output, 'w')
        opened_file = True # remember to close file afterwards

    use_names = set(table.colnames)
    if include_names is not None:
        use_names.intersection_update(include_names)
    if exclude_names is not None:
        use_names.difference_update(exclude_names)
    use_names = [x for x in table.colnames if x in use_names]

    fill_values = get_fill_values(fill_values, False)
    fill_names = set(use_names)
    if fill_include_names is not None:
        fill_names.intersection_update(fill_include_names)
    if fill_exclude_names is not None:
        fill_names.difference_update(fill_exclude_names)
    fill_cols = set([i for i, name in enumerate(use_names) if name in fill_names])

    if formats is not None:
        for name in use_names:
            if name in formats:
                table.columns[name].format = formats[name]

    col_iters = [col.iter_str_vals() for col in six.itervalues(table.columns) if col.name in use_names]
    quotechar = None if quotechar is None else str(quotechar)
    delimiter = ' ' if delimiter is None else str(delimiter)
    writer = csv.writer(output, delimiter=delimiter,
                        quotechar=quotechar,
                        lineterminator='\n') # TODO: add more params
    cdef int N = 100 # row chunk size
    cdef int num_cols = len(use_names)
    rows = [[None] * num_cols for i in range(N)]

    if header_output is not None:
        if header_output == 'comment':
            output.write(comment)
        writer.writerow([x.strip() for x in use_names] if strip_whitespace else use_names) # TODO: test this

    if output_types:
        types = ['S' if table.columns[name].dtype.kind == 'S' else 'N' for name in use_names]
        writer.writerow(types)

    for i in range(len(table)):
        for j in range(num_cols):
            field = next(col_iters[j])
            rows[i % N][j] = field
            if table.columns[use_names[j]][i] is np.ma.masked:
                field = core.masked

            if field in fill_values:
                new_val = fill_values[field][0]
                # Either this column applies to the field as specified in the fill_values parameter,
                # or no specific columns are specified and this column should apply fill_values
                if (len(fill_values[field]) > 1 and use_names[j] in fill_values[field][1:]) or \
                           (len(fill_values[field]) == 1 and j in fill_cols):
                    rows[i % N][j] = new_val
            if strip_whitespace:
                rows[i % N][j] = rows[i % N][j].strip()
        if i >= N - 1 and i % N == N - 1:
            writer.writerows(rows)

    if i % N != N - 1:
        writer.writerows(rows[:i % N + 1])
        
    output.write('\n') # append final newline
    # TODO: maybe use os.linesep instead

    if opened_file:
        output.close()

def get_fill_values(fill_values, read=True):
    if len(fill_values) > 0 and isinstance(fill_values[0], six.string_types): # e.g. fill_values=('999', '0')
        fill_values = [fill_values]
    else:
        fill_values = fill_values
    try:
        # Create a dict with the values to be replaced as keys
        if read:
            fill_values = dict([(l[0].encode('utf-8'), l[1:]) for l in fill_values])
        else:
            fill_values = dict([(l[0], l[1:]) for l in fill_values]) # don't worry about unicode for writing

    except IndexError:
        raise ValueError("Format of fill_values must be "
                         "(<bad>, <fill>, <optional col1>, ...)")
    return fill_values
