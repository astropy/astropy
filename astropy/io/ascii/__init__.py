# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" An extensible ASCII table reader and writer.

"""

from __future__ import absolute_import, division, print_function

from .core import (InconsistentTableError,
                   NoType, StrType, NumType, FloatType, IntType, AllType,
                   Column,
                   BaseInputter, ContinuationLinesInputter,
                   BaseHeader,
                   BaseData,
                   BaseOutputter, TableOutputter,
                   BaseReader,
                   BaseSplitter, DefaultSplitter, WhitespaceSplitter,
                   convert_numpy,
                   masked
                   )
from .basic import (Basic, BasicHeader, BasicData,
                    Rdb,
                    Csv,
                    Tab,
                    NoHeader,
                    CommentedHeader)
from .cds import Cds
from .latex import Latex, AASTex, latexdicts
from .html import HTML
from .ipac import Ipac
from .daophot import Daophot
from .sextractor import SExtractor
from .fixedwidth import (FixedWidth, FixedWidthNoHeader,
                         FixedWidthTwoLine, FixedWidthSplitter,
                         FixedWidthHeader, FixedWidthData)
from .ui import (set_guess, get_reader, read, get_writer, write)

from . import connect
