# Licensed under a 3-clause BSD style license - see LICENSE.rst
""" An extensible ASCII table reader and writer.

"""
# flake8: noqa

from .core import (InconsistentTableError,
                   ParameterError,
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
from .fastbasic import (FastBasic,
                        FastCsv,
                        FastTab,
                        FastNoHeader,
                        FastCommentedHeader,
                        FastRdb)
from .cds import Cds
from .mrt import Mrt
from .ecsv import Ecsv, EcsvFixedWidth
from .latex import Latex, AASTex, latexdicts
from .html import HTML
from .ipac import Ipac
from .daophot import Daophot
from .qdp import QDP
from .sextractor import SExtractor
from .fixedwidth import (FixedWidth, FixedWidthNoHeader,
                         FixedWidthTwoLine, FixedWidthSplitter,
                         FixedWidthHeader, FixedWidthData)
from .rst import RST
from .ui import (set_guess, get_reader, read, get_writer, write, get_read_trace)

from . import connect
