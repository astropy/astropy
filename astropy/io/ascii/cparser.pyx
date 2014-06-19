# Licensed under a 3-clause BSD style license - see LICENSE.rst

import six
from ...utils.data import get_readable_fileobj

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
            if '\n' not in source and '\r' not in source + '':
                with get_readable_fileobj(source) as file_obj:
                    source = file_obj.read()
        elif hasattr(source, 'read'):
            source = source.read()
        else:
            try:
                source[0]
                source[0:1]
                iter(source)
                source = '\n'.join(source)
            except TypeError:
                raise TypeError('Input "table" must be a string (filename or data) '
                                'or an iterable')
        self.source = source
        src = source
        self.parser.source = src

    def read_header(self):
        header = []

        # header_start is a valid line number
        if self.header_start is not None and self.header_start >= 0:
            pass
