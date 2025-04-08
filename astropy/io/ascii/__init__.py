# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""An extensible ASCII table reader and writer."""

from astropy import config as _config

from . import connect
from .basic import (
    Basic,
    BasicData,
    BasicHeader,
    CommentedHeader,
    Csv,
    NoHeader,
    Rdb,
    Tab,
)
from .cds import Cds
from .core import (
    AllType,
    BaseData,
    BaseHeader,
    BaseInputter,
    BaseOutputter,
    BaseReader,
    BaseSplitter,
    Column,
    ContinuationLinesInputter,
    DefaultSplitter,
    FloatType,
    InconsistentTableError,
    IntType,
    NoType,
    NumType,
    ParameterError,
    StrType,
    TableOutputter,
    WhitespaceSplitter,
    convert_numpy,
    masked,
)
from .daophot import Daophot
from .ecsv import Ecsv
from .fastbasic import (
    FastBasic,
    FastCommentedHeader,
    FastCsv,
    FastNoHeader,
    FastRdb,
    FastTab,
)
from .fixedwidth import (
    FixedWidth,
    FixedWidthData,
    FixedWidthHeader,
    FixedWidthNoHeader,
    FixedWidthSplitter,
    FixedWidthTwoLine,
)
from .html import HTML
from .ipac import Ipac
from .latex import AASTex, Latex, latexdicts
from .mrt import Mrt
from .qdp import QDP
from .rst import RST
from .sextractor import SExtractor
from .tdat import Tdat
from .ui import get_read_trace, get_reader, get_writer, read, set_guess, write


class Conf(_config.ConfigNamespace):
    """
    Configuration parameters for `astropy.io.ascii`.
    """

    guess_limit_lines = _config.ConfigItem(
        10000,
        "When guessing the format of a table, this is the number of lines that "
        "will be used for initial guessing. If the reading succeeds based on this "
        "number of lines, then reading the full table will be attempted. If the reading "
        "based on the subset of lines fails, the format will no longer be considered. "
        "This can be set to `None` to disable the limit",
    )


conf = Conf()
