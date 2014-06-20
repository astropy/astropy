# Licensed under a 3-clause BSD style license - see LICENSE.rst

import six
from ...utils.data import get_readable_fileobj

cdef extern from "src/tokenizer.h":
    ctypedef struct tokenizer_t:
        char *source		# single string containing all of the input
	int source_len 		# length of the input
        char delimiter		# delimiter character
        char comment		# comment character
        char quotechar		# quote character
	char **output_cols	# array of output strings for each column
	int *row_positions	# array of indices specifying where each row begins
	"""
	Example input/output
	--------------------

	source: "A,B,C\n10,5.,6\n1,2,3"
	output_cols: ["A101", "B5.2", "C6 3"]
	row_positions: [0, 1, 3]
	"""

    tokenizer_t *create_tokenizer(char delimiter, char comment, char quotechar)
    int tokenize(tokenizer_t *self, int rows)

class CParserError(Exception):
    """
    An instance of this class is thrown when an error occurs
    during C parsing.
    """

cdef class CParser:
    """
    A fast Cython parser class which uses underlying C code
    for tokenization.
    """

    cdef:
        tokenizer_t *tokenizer

    cdef public:
        int width
        object header

    def __cinit__(self, source,
                  delimiter=',',
                  comment=None,
                  quotechar='"',
                  header_start=0,
                  data_start=1,
                  names=None,
                  include_names=None,
                  exclude_names=None,
                  fill_values=('', '0'),
                  fill_include_names=None,
                  fill_exclude_names=None):

        self.tokenizer = create_tokenizer(ord(delimiter), ord(comment), ord(quotechar))
        self.setup_tokenizer(source)
        self.header_start = header_start
        self.data_start = data_start
        self.names = names
        self.include_names = include_names
        self.exclude_names = exclude_names
        self.fill_values = fill_values
        self.fill_include_names = fill_include_names
        self.fill_exclude_names = fill_exclude_names

    cdef setup_tokenizer(self, source):
        cdef char *src

        if isinstance(source, six.string_types):
            if '\n' not in source and '\r' not in source + '': #todo: check else case
                with get_readable_fileobj(source) as file_obj:
                    source = file_obj.read()
        elif hasattr(source, 'read'): # file-like object
            source = source.read()
        else:
            try:
                source = '\n'.join(source) # iterable sequence of lines
            except TypeError:
                raise TypeError('Input "table" must be a file-like object, a string (filename'
		      		       'or data), or an iterable')
	# Create a reference to the Python object so its char * pointer remains valid
        self.source = source
        src = source
        self.parser.source = src
	self.parser.source_len = len(self.source)

    def read_header(self):
        header = []

        # header_start is a valid line number
        if self.header_start is not None and self.header_start >= 0:
            if tokenize(self.parser, 1) != 0:
	       raise CParserError("An error occurred while tokenizing the first line of data")
	    