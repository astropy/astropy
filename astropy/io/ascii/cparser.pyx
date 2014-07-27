# Licensed under a 3-clause BSD style license - see LICENSE.rst

import numpy as np
cimport numpy as np
from numpy import ma
from ...utils.data import get_readable_fileobj
from ...table import pprint
from ...extern import six
from . import core
from libc cimport stdio
from distutils import version
import csv
import os
import math
import multiprocessing
try:
    import Queue
except ImportError: # on python 3, the module is named queue
    import queue as Queue

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
        stdio.FILE *fhandle     # file handle for header reading (if applicable)
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

    tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar,
                                  int fill_extra_cols, int strip_whitespace_lines,
                                  int strip_whitespace_fields)
    tokenizer_t *copy_tokenizer(tokenizer_t *t)
    void delete_tokenizer(tokenizer_t *tokenizer)
    int skip_lines(tokenizer_t *self, int offset, int header)
    int tokenize(tokenizer_t *self, int end, int header,
                 int *use_cols, int use_cols_len)
    long str_to_long(tokenizer_t *self, char *str)
    double str_to_double(tokenizer_t *self, char *str)
    void start_iteration(tokenizer_t *self, int col)
    int finished_iteration(tokenizer_t *self)
    char *next_field(tokenizer_t *self)
    char *read_file_chunk(stdio.FILE *fhandle, int len)
    long file_len(stdio.FILE *fhandle)
    char *get_line(stdio.FILE *fhandle)
    char *read_file_data(stdio.FILE *fhandle, int len)

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
    "type conversion error",
    "overflow error"
    ]))

cdef class FileString:
    """
    A file wrapper class which provides string-like methods.
    """
    cdef:
        stdio.FILE *fhandle
        int pos
        int length

    def __cinit__(self, fname):
        self.fhandle = stdio.fopen(fname, b'r')
        if not self.fhandle:
            raise IOError('File "{0}" could not be opened'.format(fname))
        self.pos = 0
        self.length = -1

    def __dealloc__(self):
        if self.fhandle:
            stdio.fclose(self.fhandle)

    def __len__(self):
        if self.length == -1:
            self.length = file_len(self.fhandle)
        return self.length

    def __getitem__(self, i):
        if i != self.pos + 1:
            stdio.fseek(self.fhandle, i, stdio.SEEK_SET)
        self.pos = i
        return chr(stdio.getc(self.fhandle))

    def splitlines(self):
        stdio.fseek(self.fhandle, 0, stdio.SEEK_SET)
        cdef bytes line
        cdef char *ptr
        while not stdio.feof(self.fhandle):
            ptr = get_line(self.fhandle)
            if not ptr:
                raise IOError("An error occurred while splitting file into lines")
            line = ptr
            yield line[:-1]

cdef class CParser:
    """
    A fast Cython parser class which uses underlying C code
    for tokenization.
    """

    cdef:
        tokenizer_t *tokenizer
        int data_start
        object data_end
        object include_names
        object exclude_names
        dict fill_values
        object fill_include_names
        object fill_exclude_names
        set fill_names
        int fill_extra_cols
        np.ndarray use_cols
        bytes source_bytes
        object parallel

    cdef public:
        int width
        object names
        object source
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
                  fill_extra_cols=0,
                  parallel=False):

        if comment is None:
            comment = '\x00' # tokenizer ignores all comments if comment='\x00'
        self.tokenizer = create_tokenizer(ord(delimiter), ord(comment), ord(quotechar),
                                          fill_extra_cols,
                                          strip_line_whitespace,
                                          strip_line_fields)
        self.source = None
        self.setup_tokenizer(source)
        self.header_start = header_start
        self.data_start = data_start
        self.data_end = data_end
        self.names = names
        self.include_names = include_names
        self.exclude_names = exclude_names
        self.fill_values = get_fill_values(fill_values)
        self.fill_include_names = fill_include_names
        self.fill_exclude_names = fill_exclude_names
        self.fill_extra_cols = fill_extra_cols

        # parallel=True indicates that we should use the CPU count
        if parallel is True:
            parallel = multiprocessing.cpu_count()
        # If parallel = 1 or 0, don't use multiprocessing
        elif parallel is not False and parallel < 2:
            parallel = False
        self.parallel = parallel
    
    def __dealloc__(self):
        if self.tokenizer:
            delete_tokenizer(self.tokenizer) # perform C memory cleanup

    cdef get_error(self, code, num_rows, msg):
        err_msg = ERR_CODES.get(code, "unknown error")

        # error code is lambda function taking current line as input
        if callable(err_msg):
            err_msg = err_msg(num_rows + 1)

        return CParserError("{0}: {1}".format(msg, err_msg))

    cdef raise_error(self, msg):
        raise self.get_error(self.tokenizer.code, self.tokenizer.num_rows, msg)

    cpdef setup_tokenizer(self, source):
        cdef FileString fstring

        if isinstance(source, six.string_types): # filename or data
            if '\n' not in source: # filename
                fstring = FileString(source)
                self.tokenizer.fhandle = fstring.fhandle
                self.source = fstring
                self.tokenizer.source_len = len(fstring)
                return
            # Otherwise, source is the actual data so we leave it be
        elif hasattr(source, 'read'): # file-like object
            with get_readable_fileobj(source) as file_obj:
                source = file_obj.read()
        elif isinstance(source, FileString):
            self.tokenizer.fhandle = (<FileString>source).fhandle
            self.source = source
            self.tokenizer.source_len = len(source)
            return
        else:
            try:
                source = '\n'.join(source) # iterable sequence of lines
            except TypeError:
                raise TypeError('Input "table" must be a file-like object, a '
                                'string (filename or data), or an iterable')
        # Create a reference to the Python object so its char * pointer remains valid
        self.source = source

        # encode in ASCII for char * handling
        self.source_bytes = self.source.encode('ascii')
        self.tokenizer.source = self.source_bytes
        self.tokenizer.source_len = len(self.source_bytes)
        self.tokenizer.fhandle = <stdio.FILE *>0

    cdef _read_file(self):            
        cdef char *ptr = read_file_data(self.tokenizer.fhandle,
                                        self.tokenizer.source_len)
        if not ptr:
            raise IOError("An error occurred while reading whole file into memory")

        self.source_bytes = ptr
        self.source = self.source_bytes.decode('ascii')
        self.tokenizer.fhandle = <stdio.FILE *>0
        self.tokenizer.source = self.source_bytes

    def read_header(self):
        self.tokenizer.source_pos = 0

        if self.names:
            self.width = len(self.names)

        # header_start is a valid line number
        elif self.header_start is not None and self.header_start >= 0:
            if skip_lines(self.tokenizer, self.header_start, 1) != 0:
                self.raise_error("an error occurred while advancing to the "
                                 "first header line")
            if tokenize(self.tokenizer, -1, 1, <int *> 0, 0) != 0:
                self.raise_error("an error occurred while tokenizing the header line")
            self.names = []
            name = ''

            for i in range(self.tokenizer.header_len):
                c = self.tokenizer.header_output[i] # next char in header string
                if not c: # zero byte -- field terminator
                    if name:
                        # replace empty placeholder with ''
                        self.names.append(name.replace('\x01', ''))
                        name = ''
                    else:
                        break # end of string
                else:
                    name += chr(c)
            self.width = len(self.names)

        else:
            # Get number of columns from first data row
            if tokenize(self.tokenizer, -1, 1, <int *> 0, 0) != 0:
                self.raise_error("an error occurred while tokenizing the first line of data")
            self.width = 0
            for i in range(self.tokenizer.header_len):
                # zero byte -- field terminator
                if not self.tokenizer.header_output[i]:
                    # ends valid field
                    if i > 0 and self.tokenizer.header_output[i - 1]:
                        self.width += 1
                    else: # end of line
                        break
            if self.width == 0: # no data
                raise core.InconsistentTableError('No data lines found, C reader '
                                            'cannot autogenerate column names')
            # auto-generate names
            self.names = ['col{0}'.format(i + 1) for i in range(self.width)]

        # "boolean" array denoting whether or not to use each column
        self.use_cols = np.ones(self.width, np.intc)
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
        if self.parallel:
            return self._read_parallel(try_int, try_float, try_string)

        # Read in a single process
        if self.tokenizer.fhandle:
            self._read_file() # get a single string with the entire data

        self.tokenizer.source_pos = 0
        if skip_lines(self.tokenizer, self.data_start, 0) != 0:
            self.raise_error("an error occurred while advancing to the first "
                             "line of data")

        cdef int data_end = -1 # keep reading data until the end
        if self.data_end is not None and self.data_end >= 0:
            data_end = self.data_end - self.data_start #TODO: handle data_start<data_end
        if tokenize(self.tokenizer, data_end, 0, <int *> self.use_cols.data,
                    len(self.use_cols)) != 0:
            self.raise_error("an error occurred while parsing table data")
        elif self.tokenizer.num_rows == 0: # no data
            return [[]] * self.width
        self._set_fill_names()
        cdef int num_rows = self.tokenizer.num_rows
        if self.data_end is not None and self.data_end < 0: # negative indexing
            num_rows += self.data_end
        return self._convert_data(self.tokenizer, try_int, try_float,
                                  try_string, num_rows)

    def _read_parallel(self, try_int, try_float, try_string):
        cdef int source_len = len(self.source)
        self.tokenizer.source_pos = 0

        if self.tokenizer.fhandle:
            stdio.fseek(self.tokenizer.fhandle, 0, stdio.SEEK_SET)

        if skip_lines(self.tokenizer, self.data_start, 0) != 0:
            self.raise_error("an error occurred while advancing to the first "
                             "line of data")

        cdef int N = self.parallel
        queue = multiprocessing.Queue()
        cdef int offset = self.tokenizer.source_pos

        if offset == source_len: # no data
            return dict((name, []) for name in self.names)

        cdef int chunksize = math.ceil((source_len - offset) / float(N))
        cdef list chunkindices = [offset]
        cdef long read_pos = stdio.ftell(self.tokenizer.fhandle) if \
                             self.tokenizer.fhandle else 0

        # This queue is used to signal processes to reconvert if necessary
        reconvert_queue = multiprocessing.Queue()

        for i in range(1, N):
            index = max(offset + chunksize * i, chunkindices[i - 1])
            while index < source_len and self.source[index] != '\n':
                index += 1
            if index < source_len:
                chunkindices.append(index + 1)
            else:
                N = i
                break

        chunkindices.append(source_len) #figure out correct chunkindices
        cdef list source_chunks = []
        cdef char *file_chunk

        if self.tokenizer.fhandle:
            stdio.fseek(self.tokenizer.fhandle, read_pos, stdio.SEEK_SET)

        for i in range(N):
            if self.tokenizer.fhandle:
                read_len = chunkindices[i + 1] - chunkindices[i]
                file_chunk = read_file_chunk(self.tokenizer.fhandle, read_len)
                if not file_chunk:
                    raise IOError('an error occurred while reading file data')
                source_chunks.append(file_chunk)
            else:
                source_chunks.append((chunkindices[i], chunkindices[i + 1]))

        self._set_fill_names()
        cdef list processes = []

        for i in range(N):
            process = multiprocessing.Process(target=self._read_chunk, args=(
                source_chunks[i], try_int, try_float, try_string,
                queue, reconvert_queue, i))
            processes.append(process)
            process.start()

        cdef list chunks = [None] * N
        cdef dict failed_procs = {}

        for i in range(N):
            data, err, proc = queue.get()
            if isinstance(err, Exception):
                for process in processes:
                    process.terminate()
                raise err
            elif err is not None: # err is (error code, error line)
                failed_procs[proc] = err
            chunks[proc] = data

        if failed_procs:
            # find the line number of the error
            line_no = 0
            for i in range(N):
                # ignore errors after data_end
                if i in failed_procs and self.data_end is None or line_no < self.data_end:
                    for process in processes:
                        process.terminate()
                    raise self.get_error(failed_procs[i][0], failed_procs[i][1] + line_no,
                                         "an error occurred while parsing table data")
                line_no += len(chunks[i][self.names[0]])

        seen_str = {}
        seen_numeric = {}
        for name in self.names:
            seen_str[name] = False
            seen_numeric[name] = False

        for chunk in chunks:
            for name in chunk:
                if chunk[name].dtype.kind == 'S': # string values in column
                    seen_str[name] = True
                elif len(chunk[name]) > 0: # ignore empty chunk columns
                    seen_numeric[name] = True

        reconvert_cols = []

        for i, name in enumerate(self.names):
            if seen_str[name] and seen_numeric[name]:
                # Reconvert to str to avoid conversion issues, e.g.
                # 5 (int) -> 5.0 (float) -> 5.0 (string)
                reconvert_cols.append(i)

        reconvert_queue.put(reconvert_cols)
        for process in processes:            
            process.join() # wait for each process to finish
        try:
            while True:
                reconverted, proc, col = queue.get(False)
                chunks[proc][self.names[col]] = reconverted
        except Queue.Empty:
            pass

        if self.data_end is not None:
            if self.data_end < 0:
                # e.g. if data_end = -1, cut the last row
                num_rows = 0
                for chunk in chunks:
                    num_rows += len(chunk[self.names[0]])
                self.data_end += num_rows
            else:
                # TODO: see if this is a problem with RDB
                self.data_end -= 1 # ignore header
            line_no = 0
            for i, chunk in enumerate(chunks):
                num_rows = len(chunk[self.names[0]])
                if line_no + num_rows > self.data_end:
                    for name in self.names:
                        # truncate columns
                        chunk[name] = chunk[name][:self.data_end - line_no]
                    del chunks[i + 1:]
                    break
                line_no += num_rows

        ret = chunks[0]

        for chunk in chunks[1:]:
            for name in chunk:
                if isinstance(ret[name], ma.masked_array) or isinstance(
                        chunk[name], ma.masked_array):
                    ret[name] = ma.concatenate((ret[name], chunk[name]))
                else:
                    ret[name] = np.concatenate((ret[name], chunk[name]))

        return ret

    def _read_chunk(self, source_chunk, try_int, try_float, try_string, queue,
                    reconvert_queue, i):
        cdef tokenizer_t *chunk_tokenizer = copy_tokenizer(self.tokenizer)
        if isinstance(source_chunk, tuple):
            chunk_tokenizer.source = self.source_bytes
            chunk_tokenizer.source_len = source_chunk[1]
            chunk_tokenizer.source_pos = source_chunk[0]
        else:
            chunk_tokenizer.source = source_chunk
            chunk_tokenizer.source_len = len(source_chunk)
        chunk_tokenizer.num_cols = self.width

        data = None
        err = None

        if tokenize(chunk_tokenizer, -1, 0, <int *> self.use_cols.data,
                    len(self.use_cols)) != 0:
            err = (chunk_tokenizer.code, chunk_tokenizer.num_rows)
        if chunk_tokenizer.num_rows == 0: # no data
            data = dict((name, np.array([], np.int_)) for name in self.names)
        else:
            try:
                data = self._convert_data(chunk_tokenizer,
                                          try_int, try_float, try_string, -1)
            except Exception as e:
                delete_tokenizer(chunk_tokenizer)
                queue.put((None, e, i))
                return
        #TODO: check if queue can raise errors
        queue.put((data, err, i))

        reconvert_cols = reconvert_queue.get()
        for col in reconvert_cols:
            queue.put((self._convert_str(chunk_tokenizer, col, -1), i, col))
        reconvert_queue.put(reconvert_cols) # return to the queue for other processes

    cdef _set_fill_names(self):
        self.fill_names = set(self.names)
        if self.fill_include_names is not None:
            self.fill_names.intersection_update(self.fill_include_names)
        if self.fill_exclude_names is not None:
            self.fill_names.difference_update(self.fill_exclude_names)

    cdef _convert_data(self, tokenizer_t *t, try_int, try_float, try_string, num_rows):
        cols = {}

        for i, name in enumerate(self.names):
            # Try int first, then float, then string
            try:
                if try_int and not try_int[name]:
                    raise ValueError()
                cols[name] = self._convert_int(t, i, num_rows)
            except ValueError:
                try:
                    if try_float and not try_float[name]:
                        raise ValueError()
                    cols[name] = self._convert_float(t, i, num_rows)
                except ValueError:
                    if try_string and not try_string[name]:
                        raise ValueError('Column {0} failed to convert'.format(name))
                    cols[name] = self._convert_str(t, i, num_rows)

        return cols

    cdef np.ndarray _convert_int(self, tokenizer_t *t, int i, int nrows):
        cdef int num_rows = t.num_rows
        if nrows != -1:
            num_rows = nrows
        # intialize ndarray
        cdef np.ndarray col = np.empty(num_rows, dtype=np.int_)
        cdef np.ndarray str_col = np.empty(num_rows, dtype=object) 
        cdef long converted
        cdef int row = 0
        cdef long *data = <long *> col.data # pointer to raw data
        cdef bytes field
        cdef bytes new_val
        mask = set() # set of indices for masked values
        start_iteration(t, i) # begin the iteration process in C

        while not finished_iteration(t):
            if row == num_rows: # end prematurely if we aren't using every row
                break
            # retrieve the next field in a bytes value
            field = next_field(t)

            if field in self.fill_values:
                new_val = str(self.fill_values[field][0]).encode('ascii')

                # Either this column applies to the field as specified in the 
                # fill_values parameter, or no specific columns are specified
                # and this column should apply fill_values.
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) \
                   or (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    # try converting the new value
                    converted = str_to_long(t, new_val)
                else:
                    converted = str_to_long(t, field)

            else:
                # convert the field to long (widest integer type)
                converted = str_to_long(t, field)

            if t.code in (CONVERSION_ERROR, OVERFLOW_ERROR):
                # no dice
                t.code = NO_ERROR
                raise ValueError()
            
            data[row] = converted
            row += 1

        if mask:
            # convert to masked_array
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in
                                              range(row)])
        else:
            return col

    cdef np.ndarray _convert_float(self, tokenizer_t *t, int i, int nrows):
        # very similar to _convert_int()
        cdef int num_rows = t.num_rows
        if nrows != -1:
            num_rows = nrows

        cdef np.ndarray col = np.empty(num_rows, dtype=np.float_)
        cdef double converted
        cdef int row = 0
        cdef double *data = <double *> col.data
        cdef bytes field
        cdef bytes new_val
        mask = set()

        start_iteration(t, i)
        while not finished_iteration(t):
            if row == num_rows:
                break
            field = next_field(t)
            if field in self.fill_values:
                new_val = str(self.fill_values[field][0]).encode('ascii')
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) \
                   or (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    converted = str_to_double(t, new_val)
                else:
                    converted = str_to_double(t, field)
            else:
                converted = str_to_double(t, field)

            if t.code == CONVERSION_ERROR:
                t.code = NO_ERROR
                raise ValueError()
            elif t.code == OVERFLOW_ERROR:
                t.code = NO_ERROR
                # In numpy < 1.6, using type inference yields a float for 
                # overflow values because the error raised is not specific.
                # This replicates the old reading behavior (see #2234).
                if version.LooseVersion(np.__version__) < version.LooseVersion('1.6'):
                    col[row] = new_val if field in self.fill_values else field
                else:
                    raise ValueError()
            else:
                data[row] = converted
            row += 1

        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in
                                              range(row)])
        else:
            return col

    cdef np.ndarray _convert_str(self, tokenizer_t *t, int i, int nrows):
        # similar to _convert_int, but no actual conversion
        cdef int num_rows = t.num_rows
        if nrows != -1:
            num_rows = nrows

        cdef np.ndarray col = np.empty(num_rows, dtype=object)
        cdef int row = 0
        cdef bytes field
        cdef bytes new_val
        cdef int max_len = 0 # greatest length of any element
        mask = set()

        start_iteration(t, i)
        while not finished_iteration(t):
            if row == num_rows:
                break
            field = next_field(t)
            if field in self.fill_values:
                el = str(self.fill_values[field][0])
                if (len(self.fill_values[field]) > 1 and self.names[i] in self.fill_values[field][1:]) \
                   or (len(self.fill_values[field]) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                else:
                    el = field.decode('ascii')
            else:
                el = field.decode('ascii')
            # update max_len with the length of each field
            max_len = max(max_len, len(el))
            col[row] = el
            row += 1

        # convert to string with smallest length possible
        col = col.astype('|S{0}'.format(max_len))
        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in
                                              range(row)])
        else:
            return col

cdef class FastWriter:
    """
    A fast Cython writing class for writing tables
    as ASCII data.
    """

    cdef:
        object table
        list use_names
        dict fill_values
        set fill_cols
        list col_iters
        list formats
        list format_funcs
        list types
        str quotechar
        str delimiter
        int strip_whitespace
        object comment

    def __cinit__(self, table,
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

        self.table = table
        self.comment = comment
        self.strip_whitespace = strip_whitespace
        use_names = set(table.colnames)

        # Apply include_names before exclude_names
        if include_names is not None:
            use_names.intersection_update(include_names)
        if exclude_names is not None:
            use_names.difference_update(exclude_names)
        # preserve column ordering via list
        self.use_names = [x for x in table.colnames if x in use_names]

        fill_values = get_fill_values(fill_values, False)
        self.fill_values = fill_values.copy()

        # Add int/float versions of each fill value (if applicable)
        # to the fill_values dict. This prevents the writer from having
        # to call unicode() on every value, which is a major
        # performance hit.
        for key, val in six.iteritems(fill_values):
            try:
                self.fill_values[int(key)] = val
                self.fill_values[float(key)] = val
            except (ValueError, np.ma.MaskError):
                pass

        fill_names = set(self.use_names)
        # Apply fill_include_names before fill_exclude_names
        if fill_include_names is not None:
            fill_names.intersection_update(fill_include_names)
        if fill_exclude_names is not None:
            fill_names.difference_update(fill_exclude_names)
        # Preserve column ordering
        self.fill_cols = set([i for i, name in enumerate(self.use_names) if
                              name in fill_names])

        # formats in user-specified dict should override
        # existing column formats
        if formats is not None:
            for name in self.use_names:
                if name in formats:
                    self.table[name].format = formats[name]

        self.col_iters = []
        self.formats = []
        self.format_funcs = []

        for col in six.itervalues(table.columns):
            if col.name in self.use_names: # iterate over included columns
                if col.format is None:
                    self.format_funcs.append(None)
                else:
                    self.format_funcs.append(pprint._format_funcs.get(col.format,
                                                        pprint._auto_format_func))
                # col is a numpy.ndarray, so we convert it to
                # an ordinary list because csv.writer will call
                # np.array_str() on each numpy value, which is
                # very inefficient
                self.col_iters.append(iter(col.tolist()))
                self.formats.append(col.format)

        self.quotechar = None if quotechar is None else str(quotechar)
        self.delimiter = ' ' if delimiter is None else str(delimiter)
        # 'S' for string types, 'N' for numeric types
        self.types = ['S' if self.table[name].dtype.kind in ('S', 'U') else 'N'
                      for name in self.use_names]

    def _write_header(self, output, writer, header_output, output_types):
        if header_output is not None:
            if header_output == 'comment':
                output.write(self.comment)
            writer.writerow([x.strip() for x in self.use_names] if 
                            self.strip_whitespace else self.use_names)

        if output_types:
            writer.writerow(self.types)

    def write(self, output, header_output, output_types):
        opened_file = False

        if not hasattr(output, 'write'): # output is a filename
            output = open(output, 'w')
            opened_file = True # remember to close file afterwards
        writer = csv.writer(output,
                            delimiter=self.delimiter,
                            doublequote=True,
                            escapechar=None,
                            quotechar=self.quotechar,
                            quoting=csv.QUOTE_MINIMAL,
                            lineterminator=os.linesep)
        self._write_header(output, writer, header_output, output_types)
                                         
        # Split rows into N-sized chunks, since we don't want to
        # store all the rows in memory at one time (inefficient)
        # or fail to take advantage of the speed boost of writerows()
        # over writerow().
        cdef int N = 100
        cdef int num_cols = len(self.use_names)
        cdef int num_rows = len(self.table)
        # cache string columns beforehand
        cdef set string_rows = set([i for i, type in enumerate(self.types) if 
                                    type == 'S'])
        cdef list rows = [[None] * num_cols for i in range(N)]

        for i in range(num_rows):
            for j in range(num_cols):
                orig_field = next(self.col_iters[j]) # get field
                # str_val monitors whether we should check if the field
                # should be stripped
                str_val = True

                if orig_field is None: # tolist() converts ma.masked to None
                    field = core.masked
                    rows[i % N][j] = '--'

                elif self.format_funcs[j] is not None:
                    # TODO: find a better way to format non-numpy types
                    field = self.format_funcs[j](self.formats[j], np.array(
                                                          [orig_field])[0])
                    rows[i % N][j] = field

                else:
                    field = orig_field
                    rows[i % N][j] = field
                    str_val = j in string_rows

                if field in self.fill_values:
                    new_val = self.fill_values[field][0]
                    # Either this column applies to the field as specified in
                    # the fill_values parameter, or no specific columns are
                    # specified and this column should apply fill_values.
                    if (len(self.fill_values[field]) > 1 and self.use_names[j] in self.fill_values[field][1:]) \
                       or (len(self.fill_values[field]) == 1 and j in self.fill_cols):
                        str_val = True
                        rows[i % N][j] = new_val
                        if self.strip_whitespace: # new_val should be a string
                            rows[i % N][j] = rows[i % N][j].strip()

                if str_val and self.strip_whitespace:
                    rows[i % N][j] = rows[i % N][j].strip()

            if i >= N - 1 and i % N == N - 1: # rows is now full
                writer.writerows(rows)

        # Write leftover rows not included in previous chunks
        if i % N != N - 1:
            writer.writerows(rows[:i % N + 1])

        if opened_file:
            output.close()

def get_fill_values(fill_values, read=True):
    if len(fill_values) > 0 and isinstance(fill_values[0], six.string_types):
        # e.g. fill_values=('999', '0')
        fill_values = [fill_values]
    else:
        fill_values = fill_values
    try:
        # Create a dict with the values to be replaced as keys
        if read:
            fill_values = dict([(l[0].encode('ascii'), l[1:]) for l in fill_values])
        else:
            # don't worry about encoding for writing
            fill_values = dict([(l[0], l[1:]) for l in fill_values])

    except IndexError:
        raise ValueError("Format of fill_values must be "
                         "(<bad>, <fill>, <optional col1>, ...)")
    return fill_values
