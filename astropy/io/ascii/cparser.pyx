# Licensed under a 3-clause BSD style license - see LICENSE.rst

import csv
import os
import math
import multiprocessing
import mmap

import numpy as np
cimport numpy as np
from numpy import ma
from libc cimport stdio
from cpython.buffer cimport PyBUF_SIMPLE
from cpython.buffer cimport Py_buffer
from cpython.buffer cimport PyObject_GetBuffer, PyBuffer_Release

from ...utils.data import get_readable_fileobj
from ...table import pprint
from ...extern import six
from . import core

try:
    import Queue
except ImportError: # in python 3, the module is named queue
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
        CARRIAGE_RETURN

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
        char **output_cols     # array of output strings for each column
        char **col_ptrs        # array of pointers to current output position for each col
        int *output_len        # length of each output column string
        int num_cols           # number of table columns
        int num_rows           # number of table rows
        int fill_extra_cols    # represents whether or not to fill rows with too few values
        tokenizer_state state  # current state of the tokenizer
        err_code code          # represents the latest error that has occurred
        int iter_col           # index of the column being iterated over
        char *curr_pos         # current iteration position
        char *buf              # buffer for empty data
        int strip_whitespace_lines  # whether to strip whitespace at the beginning and end of lines
        int strip_whitespace_fields # whether to strip whitespace at the beginning and end of fields
        int use_fast_converter      # whether to use the fast converter for floats
        char *comment_lines    # single null-delimited string containing comment lines
        int comment_lines_len  # length of comment_lines in memory
        int comment_pos        # current index in comment_lines
        # Example input/output
        # --------------------
        # source: "A,B,C\n10,5.,6\n1,2,3"
        # output_cols: ["A\x0010\x001", "B\x005.\x002", "C\x006\x003"]

    ctypedef struct memory_map:
        char *ptr
        int len
        void *file_ptr
        void *handle

    tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar,
                                  int fill_extra_cols, int strip_whitespace_lines,
                                  int strip_whitespace_fields, int use_fast_converter)
    void delete_tokenizer(tokenizer_t *tokenizer)
    int skip_lines(tokenizer_t *self, int offset, int header)
    int tokenize(tokenizer_t *self, int end, int header, int num_cols)
    long str_to_long(tokenizer_t *self, char *str)
    double fast_str_to_double(tokenizer_t *self, char *str)
    double str_to_double(tokenizer_t *self, char *str)
    void start_iteration(tokenizer_t *self, int col)
    char *next_field(tokenizer_t *self, int *size)
    char *get_line(char *ptr, int *len, int map_len)
    void reset_comments(tokenizer_t *self)

cdef extern from "Python.h":
    int PyObject_AsReadBuffer(object obj, const void **buffer, Py_ssize_t *buffer_len)

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
    A wrapper class for a memory-mapped file pointer.
    """
    cdef:
        object fhandle
        object mmap
        void *mmap_ptr
        Py_buffer buf

    def __cinit__(self, fname):
        self.fhandle = open(fname, 'r')
        if not self.fhandle:
            raise IOError('File "{0}" could not be opened'.format(fname))
        self.mmap = mmap.mmap(self.fhandle.fileno(), 0, prot=mmap.PROT_READ)
        cdef Py_ssize_t buf_len = len(self.mmap)
        if six.PY2:
            PyObject_AsReadBuffer(self.mmap, &self.mmap_ptr, &buf_len)
        else:
            PyObject_GetBuffer(self.mmap, &self.buf, PyBUF_SIMPLE)
            self.mmap_ptr = self.buf.buf

    def __dealloc__(self):
        if self.mmap:
            if not six.PY2: # free buffer memory to prevent a resource leak
                PyBuffer_Release(&self.buf)
            self.mmap.close()
            self.fhandle.close()

    def __len__(self):
        return len(self.mmap)

    def __getitem__(self, i):
        return self.mmap[i]

    def splitlines(self):
        """
        Return a generator yielding lines from the memory map.
        """
        cdef char *ptr = <char *>self.mmap_ptr
        cdef char *tmp
        cdef int line_len
        cdef int map_len = len(self.mmap)

        while ptr:
            tmp = get_line(ptr, &line_len, map_len)
            yield ptr[:line_len].decode('ascii')
            ptr = tmp

cdef class CParser:
    """
    A fast Cython parser class which uses underlying C code
    for tokenization.
    """

    cdef:
        tokenizer_t *tokenizer
        object names
        int data_start
        object data_end
        object include_names
        object exclude_names
        object fill_values
        object fill_empty
        object fill_include_names
        object fill_exclude_names
        object fill_names
        int fill_extra_cols
        bytes source_bytes
        char *source_ptr
        object parallel
        set use_cols

    cdef public:
        int width
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
                  fast_reader=True):

        # Handle fast_reader parameter (True or dict specifying options)
        if fast_reader is True or fast_reader == 'force':
            fast_reader = {}
        elif fast_reader is False: # shouldn't happen
            raise core.ParameterError("fast_reader cannot be False for fast readers")
        # parallel and use_fast_reader are False by default
        use_fast_converter = fast_reader.pop('use_fast_converter', False)
        parallel = fast_reader.pop('parallel', False)
        if fast_reader:
            raise core.FastOptionsError("Invalid parameter in fast_reader dict")
        if parallel and os.name == 'nt':
            raise NotImplementedError("Multiprocessing is not yet supported on Windows")

        if comment is None:
            comment = '\x00' # tokenizer ignores all comments if comment='\x00'
        self.tokenizer = create_tokenizer(ord(delimiter), ord(comment), ord(quotechar),
                                          fill_extra_cols,
                                          strip_line_whitespace,
                                          strip_line_fields,
                                          use_fast_converter)
        self.source = None
        if source is not None:
            self.setup_tokenizer(source)
        self.header_start = header_start
        self.data_start = data_start
        self.data_end = data_end
        self.names = names
        self.include_names = include_names
        self.exclude_names = exclude_names
        self.fill_values = fill_values
        self.fill_include_names = fill_include_names
        self.fill_exclude_names = fill_exclude_names
        self.fill_names = None
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
            if '\n' not in source and '\r' not in source: # filename
                fstring = FileString(source)
                self.tokenizer.source = <char *>fstring.mmap_ptr
                self.source_ptr = <char *>fstring.mmap_ptr
                self.source = fstring
                self.tokenizer.source_len = len(fstring)
                return
            # Otherwise, source is the actual data so we leave it be
        elif hasattr(source, 'read'): # file-like object
            with get_readable_fileobj(source) as file_obj:
                source = file_obj.read()
        elif isinstance(source, FileString):
            self.tokenizer.source = <char *>((<FileString>source).mmap_ptr)
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

    def read_header(self):
        self.tokenizer.source_pos = 0

        if self.names:
            self.width = len(self.names)

        # header_start is a valid line number
        elif self.header_start is not None and self.header_start >= 0:
            if skip_lines(self.tokenizer, self.header_start, 1) != 0:
                self.raise_error("an error occurred while advancing to the "
                                 "first header line")
            if tokenize(self.tokenizer, -1, 1, 0) != 0:
                self.raise_error("an error occurred while tokenizing the header line")
            self.names = []
            name = ''

            for i in range(self.tokenizer.output_len[0]): # header is in first col string
                c = self.tokenizer.output_cols[0][i] # next char in header string
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
            if tokenize(self.tokenizer, -1, 1, 0) != 0:
                self.raise_error("an error occurred while tokenizing the first line of data")
            self.width = 0
            for i in range(self.tokenizer.output_len[0]): # header is in first col string
                # zero byte -- field terminator
                if not self.tokenizer.output_cols[0][i]:
                    # ends valid field
                    if i > 0 and self.tokenizer.output_cols[0][i - 1]:
                        self.width += 1
                    else: # end of line
                        break
            if self.width == 0: # no data
                raise core.InconsistentTableError('No data lines found, C reader '
                                            'cannot autogenerate column names')
            # auto-generate names
            self.names = ['col{0}'.format(i + 1) for i in range(self.width)]

        # self.use_cols should only contain columns included in output
        self.use_cols = set(self.names)
        if self.include_names is not None:
            self.use_cols.intersection_update(self.include_names)
        if self.exclude_names is not None:
            self.use_cols.difference_update(self.exclude_names)

        self.width = len(self.names)

    def read(self, try_int, try_float, try_string):
        if self.parallel:
            return self._read_parallel(try_int, try_float, try_string)

        # Read in a single process
        self.tokenizer.source_pos = 0
        if skip_lines(self.tokenizer, self.data_start, 0) != 0:
            self.raise_error("an error occurred while advancing to the first "
                             "line of data")

        cdef int data_end = -1 # keep reading data until the end
        if self.data_end is not None and self.data_end >= 0:
            data_end = max(self.data_end - self.data_start, 0) # read nothing if data_end < 0

        if tokenize(self.tokenizer, data_end, 0, len(self.names)) != 0:
            self.raise_error("an error occurred while parsing table data")
        elif self.tokenizer.num_rows == 0: # no data
            return ([np.array([], dtype=np.int_)] * self.width, [])
        self._set_fill_values()
        cdef int num_rows = self.tokenizer.num_rows
        if self.data_end is not None and self.data_end < 0: # negative indexing
            num_rows += self.data_end
        return self._convert_data(self.tokenizer, try_int, try_float,
                                  try_string, num_rows)

    def _read_parallel(self, try_int, try_float, try_string):
        cdef int source_len = len(self.source)
        self.tokenizer.source_pos = 0

        if skip_lines(self.tokenizer, self.data_start, 0) != 0:
            self.raise_error("an error occurred while advancing to the first "
                             "line of data")

        cdef list line_comments = self._get_comments(self.tokenizer)
        cdef int N = self.parallel
        queue = multiprocessing.Queue()
        cdef int offset = self.tokenizer.source_pos

        if offset == source_len: # no data
            return (dict((name, np.array([], dtype=np.int_)) for name in
                         self.names), [])

        cdef int chunksize = math.ceil((source_len - offset) / float(N))
        cdef list chunkindices = [offset]

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

        self._set_fill_values()
        chunkindices.append(source_len)
        cdef list processes = []

        for i in range(N):
            process = multiprocessing.Process(target=_read_chunk, args=(self,
                chunkindices[i], chunkindices[i + 1],
                try_int, try_float, try_string, queue, reconvert_queue, i))
            processes.append(process)
            process.start()

        cdef list chunks = [None] * N
        cdef dict failed_procs = {}

        for i in range(N):
            queue_ret, err, proc = queue.get()
            if isinstance(err, Exception):
                for process in processes:
                    process.terminate()
                raise err
            elif err is not None: # err is (error code, error line)
                failed_procs[proc] = err
            comments, data = queue_ret
            line_comments.extend(comments)
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
                if chunk[name].dtype.kind in ('S', 'U'):
                    # string values in column
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
                self.data_end -= self.data_start # ignore header

            if self.data_end < 0: # no data
                chunks = [dict((name, []) for name in self.names)]
            else:
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

        ret = {}
        for name in self.get_names():
            col_chunks = [chunk.pop(name) for chunk in chunks]
            if any(isinstance(col_chunk, ma.masked_array) for col_chunk in col_chunks):
                ret[name] = ma.concatenate(col_chunks)
            else:
                ret[name] = np.concatenate(col_chunks)

        for process in processes:
            process.terminate()

        return ret, line_comments

    cdef _set_fill_values(self):
        if self.fill_names is None:
            self.fill_names = set(self.names)
            if self.fill_include_names is not None:
                self.fill_names.intersection_update(self.fill_include_names)
            if self.fill_exclude_names is not None:
                self.fill_names.difference_update(self.fill_exclude_names)
        self.fill_values, self.fill_empty = get_fill_values(self.fill_values)

    cdef _get_comments(self, tokenizer_t *t):
        line_comments = []
        comment = ''
        for i in range(t.comment_pos):
            c = t.comment_lines[i] # next char in comment string
            if not c: # zero byte -- line terminator
                # replace empty placeholder with ''
                line_comments.append(comment.replace('\x01', '').strip())
                comment = ''
            else:
                comment += chr(c)
        return line_comments

    cdef _convert_data(self, tokenizer_t *t, try_int, try_float, try_string, num_rows):
        cols = {}

        for i, name in enumerate(self.names):
            if name not in self.use_cols:
                continue
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

        return cols, self._get_comments(t)

    cdef np.ndarray _convert_int(self, tokenizer_t *t, int i, int nrows):
        cdef int num_rows = t.num_rows
        if nrows != -1:
            num_rows = nrows
        # intialize ndarray
        cdef np.ndarray col = np.empty(num_rows, dtype=np.int_)
        cdef long converted
        cdef int row = 0
        cdef long *data = <long *> col.data # pointer to raw data
        cdef char *field
        cdef char *empty_field = t.buf # memory address of designated empty buffer
        cdef bytes new_value
        mask = set() # set of indices for masked values
        start_iteration(t, i) # begin the iteration process in C

        for row in range(num_rows):
            # retrieve the next field as a C pointer
            field = next_field(t, <int *>0)
            replace_info = None

            if field == empty_field and self.fill_empty:
                replace_info = self.fill_empty
            # hopefully this implicit char * -> byte conversion for fill values
            # checking can be avoided in most cases, since self.fill_values will
            # be empty in the default case (self.fill_empty will do the work
            # instead)
            elif field != empty_field and self.fill_values and field in self.fill_values:
                replace_info = self.fill_values[field]

            if replace_info is not None:
                # Either this column applies to the field as specified in the
                # fill_values parameter, or no specific columns are specified
                # and this column should apply fill_values.
                if (len(replace_info) > 1 and self.names[i] in replace_info[1:]) \
                   or (len(replace_info) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    new_value = str(replace_info[0]).encode('ascii')
                    # try converting the new value
                    converted = str_to_long(t, new_value)
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
        cdef char *field
        cdef char *empty_field = t.buf
        cdef bytes new_value
        cdef int replacing
        mask = set()

        start_iteration(t, i)
        for row in range(num_rows):
            field = next_field(t, <int *>0)
            replace_info = None
            replacing = False

            if field == empty_field and self.fill_empty:
                replace_info = self.fill_empty

            elif field != empty_field and self.fill_values and field in self.fill_values:
                replace_info = self.fill_values[field]

            if replace_info is not None:
                if (len(replace_info) > 1 and self.names[i] in replace_info[1:]) \
                   or (len(replace_info) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    new_value = str(replace_info[0]).encode('ascii')
                    replacing = True
                    converted = str_to_double(t, new_value)
                else:
                    converted = str_to_double(t, field)
            else:
                converted = str_to_double(t, field)

            if t.code == CONVERSION_ERROR:
                t.code = NO_ERROR
                raise ValueError()
            elif t.code == OVERFLOW_ERROR:
                t.code = NO_ERROR
                raise ValueError()
            else:
                data[row] = converted
            row += 1

        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in
                                              range(row)])
        else:
            return col

    cdef _convert_str(self, tokenizer_t *t, int i, int nrows):
        # similar to _convert_int, but no actual conversion
        cdef int num_rows = t.num_rows
        if nrows != -1:
            num_rows = nrows

        cdef int row = 0
        cdef bytes field
        cdef int field_len
        cdef int max_len = 0
        cdef list fields_list = []
        mask = set()

        start_iteration(t, i)
        for row in range(num_rows):
            field = next_field(t, &field_len)
            replace_info = None

            if field_len == 0 and self.fill_empty:
                replace_info = self.fill_empty

            elif field_len > 0 and self.fill_values and field in self.fill_values:
                replace_info = self.fill_values[field]

            if replace_info is not None:
                el = replace_info[0].encode('ascii')
                if (len(replace_info) > 1 and self.names[i] in replace_info[1:]) \
                   or (len(replace_info) == 1 and self.names[i] in self.fill_names):
                    mask.add(row)
                    field = el

            fields_list.append(field)
            if field_len > max_len:
                max_len = field_len
            row += 1

        cdef np.ndarray col = np.array(fields_list, dtype=(np.str, max_len))

        if mask:
            return ma.masked_array(col, mask=[1 if i in mask else 0 for i in
                                              range(row)])
        else:
            return col

    def get_names(self):
        # ignore excluded columns
        return [name for name in self.names if name in self.use_cols]

    def set_names(self, names):
        self.names = names

    def __reduce__(self):
        cdef bytes source_ptr = self.source_ptr if self.source_ptr else b''
        return (_copy_cparser, (source_ptr, self.source_bytes, self.use_cols, self.fill_names,
                                self.fill_values, self.tokenizer.strip_whitespace_lines,
                                self.tokenizer.strip_whitespace_fields,
                                dict(delimiter=chr(self.tokenizer.delimiter),
                                comment=chr(self.tokenizer.comment),
                                quotechar=chr(self.tokenizer.quotechar),
                                header_start=self.header_start,
                                data_start=self.data_start,
                                data_end=self.data_end,
                                names=self.names,
                                include_names=self.include_names,
                                exclude_names=self.exclude_names,
                                fill_values=None,
                                fill_include_names=self.fill_include_names,
                                fill_exclude_names=self.fill_exclude_names,
                                fill_extra_cols=self.tokenizer.fill_extra_cols,
                                use_fast_converter=self.tokenizer.use_fast_converter,
                                parallel=False)))

def _copy_cparser(bytes src_ptr, bytes source_bytes, use_cols, fill_names, fill_values,
                  strip_whitespace_lines, strip_whitespace_fields, kwargs):
    parser = CParser(None, strip_whitespace_lines, strip_whitespace_fields, **kwargs)
    parser.use_cols = use_cols
    parser.fill_names = fill_names
    parser.fill_values = fill_values

    if src_ptr:
        parser.tokenizer.source = src_ptr
    else:
        parser.tokenizer.source = source_bytes
    return parser

def _read_chunk(CParser self, start, end, try_int,
                try_float, try_string, queue, reconvert_queue, i):
    cdef tokenizer_t *chunk_tokenizer = self.tokenizer
    chunk_tokenizer.source_len = end
    chunk_tokenizer.source_pos = start
    reset_comments(chunk_tokenizer)

    data = None
    err = None

    if tokenize(chunk_tokenizer, -1, 0, len(self.names)) != 0:
        err = (chunk_tokenizer.code, chunk_tokenizer.num_rows)
    if chunk_tokenizer.num_rows == 0: # no data
        data = dict((name, np.array([], np.int_)) for name in self.get_names())
        line_comments = self._get_comments(chunk_tokenizer)
    else:
        try:
            data, line_comments = self._convert_data(chunk_tokenizer,
                                      try_int, try_float, try_string, -1)
        except Exception as e:
            delete_tokenizer(chunk_tokenizer)
            queue.put((None, e, i))
            return

    try:
        queue.put(((line_comments, data), err, i))
    except Queue.Full as e:
        # hopefully this shouldn't happen
        delete_tokenizer(chunk_tokenizer)
        queue.pop()
        queue.put((None, e, i))
    reconvert_cols = reconvert_queue.get()
    for col in reconvert_cols:
        queue.put((self._convert_str(chunk_tokenizer, col, -1), i, col))
    delete_tokenizer(chunk_tokenizer)
    reconvert_queue.put(reconvert_cols) # return to the queue for other processes

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
        list line_comments
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
                  fill_exclude_names=None,
                  fast_writer=True):

        if fast_writer is True:
            fast_writer = {}
        # fast_writer might contain custom writing options

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
        self.line_comments = table.meta.get('comments', [])

        for col in six.itervalues(table.columns):
            if col.name in self.use_names: # iterate over included columns
                if col.format is None:
                    self.format_funcs.append(None)
                else:
                    self.format_funcs.append(pprint._format_funcs.get(col.format,
                                                            auto_format_func))
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

    cdef _write_comments(self, output):
        if self.comment is not False:
            for comment_line in self.line_comments:
                output.write(self.comment + comment_line + '\n')

    def _write_header(self, output, writer, header_output, output_types):
        if header_output is not None and header_output == 'comment':
            output.write(self.comment)
            writer.writerow([x.strip() for x in self.use_names] if
                            self.strip_whitespace else self.use_names)
            self._write_comments(output)
        else:
            self._write_comments(output)
            if header_output is not None:
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
                    field = self.format_funcs[j](self.formats[j], orig_field)
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

    # look for an empty replacement to cache for speedy conversion
    fill_empty = None
    for el in fill_values:
        if el[0] == '':
            fill_empty = el[1:]
            break

    try:
        # Create a dict with the values to be replaced as keys
        if read:
            fill_values = dict([(l[0].encode('ascii'), l[1:]) for
                                l in fill_values if l[0] != ''])
        else:
            # don't worry about encoding for writing
            fill_values = dict([(l[0], l[1:]) for l in fill_values])

    except IndexError:
        raise ValueError("Format of fill_values must be "
                         "(<bad>, <fill>, <optional col1>, ...)")
    if read:
        return (fill_values, fill_empty)
    else:
        return fill_values # cache for empty values doesn't matter for writing

def auto_format_func(format_, val):
    """
    Mimics pprint._auto_format_func for non-numpy values.
    """
    if six.callable(format_):
        format_func = lambda format_, val: format_(val)
        try:
            out = format_func(format_, val)
            if not isinstance(out, six.string_types):
                raise ValueError('Format function for value {0} returned {1} instead of string type'
                                 .format(val, type(val)))
        except Exception as err:
            raise ValueError('Format function for value {0} failed: {1}'
                             .format(val, err))
    else:
        try:
            # Convert val to Python object with tolist().  See
            # https://github.com/astropy/astropy/issues/148#issuecomment-3930809
            out = format_.format(val)
            # Require that the format statement actually did something
            if out == format_:
                raise ValueError
            format_func = lambda format_, val: format_.format(val)
        except:  # Not sure what exceptions might be raised
            try:
                out = format_ % val
                if out == format_:
                    raise ValueError
                format_func = lambda format_, val: format_ % val
            except:
                raise ValueError('Unable to parse format string {0}'
                                 .format(format_))
    pprint._format_funcs[format_] = format_func
    return out
