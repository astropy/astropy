# Licensed under a 3-clause BSD style license - see LICENSE.rst

from . import core
from ...extern import six

@six.add_metaclass(core.MetaBaseReader)
class FastBasic(object):
    """This class is intended to handle the same format addressed by the
    ordinary `Basic` writer, but it acts as a wrapper for underlying C
    code and is therefore much faster. Unlike the other readers and writers
    in `io.ascii`, this class is not very extensible and is restricted
    by optimization requirements.
    """
    _format_name = 'fast_basic'
    _description = 'Fast C reader for basic tables with custom delimiters'

    def __init__(self, **kwargs):
        self.delimiter = str(kwargs.get('delimiter', ' '))
        self.comment = kwargs.get('comment', '#')
        if self.comment is not None:
            self.comment = str(self.comment)
        self.kwargs = kwargs

    def read_header(self):
        self.engine.read_header()
        if self.engine.header is not None:
            self.names = list(self.engine.header)
        else:
            self.names = ['col{}'.format(i + 1) for i in range(self.engine.width)]

    def read(self, table):
        if len(self.comment) != 1:
            raise ParameterError("The C reader only supports length-1 comments")
        elif 'converters' in self.kwargs:
            raise ParameterError("The C reader does not support passing "
                                 "specialized converters")
        elif 'Outputter' in self.kwargs:
            raise ParameterError("The C reader does not use the Outputter parameter")
        elif 'Inputter' in self.kwargs:
            raise ParameterError("The C reader does not use the Inputter parameter")
        elif 'data_Splitter' in self.kwargs or 'header_Splitter' in self.kwargs:
            raise ParameterError("The C reader does not use a Splitter class")

        self.engine = parser.CParser(**self.kwargs)
        self.read_header()
        data = self.engine.read()
        return Table(data, names=self.names) # TODO: add masking, units, etc.
