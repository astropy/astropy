# This file connects the readers/writers to the astropy.table.Table class

from .ui import read, write
from .basic import Rdb
from .cds import Cds
from .latex import Latex
from .daophot import Daophot
from .ipac import Ipac
from ...table import io_registry

__all__ = []


# CDS
# ===

def read_cds(filename, **kwargs):
    return read(filename, Reader=Cds, **kwargs)

io_registry.register_reader('cds', read_cds)


# DAOPhot
# =======

def read_daophot(filename, **kwargs):
    return read(filename, Reader=Daophot, **kwargs)

io_registry.register_reader('daophot', read_daophot)


# LaTeX
# =====

def read_latex(filename, **kwargs):
    return read(filename, Reader=Latex, **kwargs)


def write_latex(table, filename, **kwargs):
    write(table, filename, Writer=Latex, **kwargs)

io_registry.register_reader('latex', read_latex)
io_registry.register_writer('latex', write_latex)


def is_latex(origin, args, kwargs):
    return args[0].endswith('.tex')

io_registry.register_identifier('latex', is_latex)


# RDB
# ===

def read_rdb(filename, **kwargs):
    return read(filename, Reader=Rdb, **kwargs)


def write_rdb(table, filename, **kwargs):
    write(table, filename, Writer=Rdb, **kwargs)

io_registry.register_reader('rdb', read_rdb)
io_registry.register_writer('rdb', write_rdb)


def is_rdb(origin, args, kwargs):
    return args[0].endswith('.rdb')

io_registry.register_identifier('rdb', is_rdb)
