# Licensed under a 3-clause BSD style license - see LICENSE.rst

from . import core
from ...extern import six
from ...table import Table
from . import cparser

@six.add_metaclass(core.MetaBaseReader)
class FastBasic(object):
    """
    This class is intended to handle the same format addressed by the
    ordinary :class:`Basic` writer, but it acts as a wrapper for underlying C
    code and is therefore much faster. Unlike the other ASCII readers and
    writers, this class is not very extensible and is restricted
    by optimization requirements.
    """
    _format_name = 'fast_basic'
    _description = 'Basic table with custom delimiter using the fast C engine'
    _fast = True
    fill_extra_cols = False

    def __init__(self, default_kwargs={}, **user_kwargs):
        kwargs = default_kwargs.copy()
        kwargs.update(user_kwargs) # user kwargs take precedence over defaults
        self.delimiter = str(kwargs.pop('delimiter', ' '))
        self.comment = kwargs.pop('comment', '#')
        if self.comment is not None:
            self.comment = str(self.comment)
        self.quotechar = str(kwargs.pop('quotechar', '"'))
        self.header_start = kwargs.pop('header_start', 0)
        # If data_start is not specified, start reading data right after the header line
        data_start_default = user_kwargs.get('data_start', self.header_start + 1
                                             if self.header_start is not None else 1)
        self.data_start = kwargs.pop('data_start', data_start_default)
        self.kwargs = kwargs
        self.strip_whitespace_lines = True
        self.strip_whitespace_fields = True

    def _read_header(self):
        # Use the tokenizer by default -- this method can be overrided for specialized headers
        self.engine.read_header()

    def read(self, table): # TODO: actually take the parameters from _get_reader()
        if len(self.comment) != 1:
            raise core.ParameterError("The C reader does not support a comment regex")
        elif self.data_start is None:
            raise core.ParameterError("The C reader does not allow data_start to be None")
        elif len(self.delimiter) != 1:
            raise core.ParameterError("The C reader only supports 1-char delimiters")
        elif len(self.quotechar) != 1:
            raise core.ParameterError("The C reader only supports a length-1 quote character")
        elif 'converters' in self.kwargs:
            raise core.ParameterError("The C reader does not support passing "
                                      "specialized converters")
        elif 'Outputter' in self.kwargs:
            raise core.ParameterError("The C reader does not use the Outputter parameter")
        elif 'Inputter' in self.kwargs:
            raise core.ParameterError("The C reader does not use the Inputter parameter")
        elif 'data_Splitter' in self.kwargs or 'header_Splitter' in self.kwargs:
            raise core.ParameterError("The C reader does not use a Splitter class")

        self.engine = cparser.CParser(table, self.strip_whitespace_lines, self.strip_whitespace_fields,
                                      delimiter=self.delimiter,
                                      header_start=self.header_start,
                                      comment=self.comment,
                                      quotechar=self.quotechar,
                                      data_start=self.data_start,
                                      fill_extra_cols=self.fill_extra_cols,
                                      **self.kwargs)
        self._read_header()
        data = self.engine.read()
        return Table(data, names=list(self.engine.names)) # TODO: add masking, units, etc.

class FastCsv(FastBasic):
    """
    A faster version of the ordinary :class:`Csv` writer that uses the optimized
    C parsing engine. Note that this reader will append empty field values to
    the end of any row with not enough columns, while :class:`FastBasic` simply
    raises an error.
    """
    _format_name = 'fast_csv'
    _io_registry_suffix = '.csv'
    _description = 'Comma-separated values table using the fast C engine'
    _fast = True
    fill_extra_cols = True

    def __init__(self, **kwargs):
        FastBasic.__init__(self, {'delimiter': ','}, **kwargs)

class FastTab(FastBasic):
    """
    A faster version of the ordinary :class:`Tab` reader that uses the optimized
    C parsing engine.
    """
    _format_name = 'fast_tab'
    _description = 'Tab-separated values table using the fast C engine'
    _fast = True

    def __init__(self, **kwargs):
        delimiter = kwargs.pop('delimiter', '\t')
        FastBasic.__init__(self, {'delimiter': '\t'}, **kwargs)
        self.strip_whitespace_lines = False
        self.strip_whitespace_fields = False

class FastNoHeader(FastBasic):
    """
    This class uses the fast C engine to read tables with no header line. If the
    names parameter is unspecified, the columns will be autonamed with "col%d".
    """
    _format_name = 'fast_no_header'
    _description = 'Basic table with no headers using the fast C engine'
    _fast = True

    def __init__(self, **kwargs):
        FastBasic.__init__(self, {'header_start': None, 'data_start': 0}, **kwargs)

class FastCommentedHeader(FastBasic):
    """
    A faster version of the :class:`CommentedHeader` reader, which looks for
    column names in a commented line. ``header_start`` denotes the index of
    the header line among all commented lines and is 0 by default.
    """
    _format_name = 'fast_commented_header'
    _description = 'Columns name in a commented line using the fast C engine'
    _fast = True

    def __init__(self, **kwargs):
        FastBasic.__init__(self, {'data_start': 0}, **kwargs)

    def _read_header(self):
        tmp = self.engine.source
        commented_lines = [line.lstrip()[1:] for line in tmp.split('\n') if line and line.lstrip()[0] == '#']
        self.engine.setup_tokenizer([commented_lines[self.header_start]])
        self.engine.header_start = 0
        self.engine.read_header()
        self.engine.setup_tokenizer(tmp)

# TODO: write FastRdb, FastCommentedHeader...will require some changes to tokenizer


