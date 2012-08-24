# This file connects the readers/writers to the astropy.table.Table class

from .ui import read, write
from .basic import Rdb
from .cds import Cds
from .latex import Latex
from .daophot import Daophot
from .ipac import Ipac
from ...table import io_registry

__all__ = []


# Generic
# =======


def read_asciitable(filename, **kwargs):
    return read(filename, **kwargs)

io_registry.register_reader('ascii', read_asciitable)


def write_asciitable(table, filename, **kwargs):
    return write(table, filename, **kwargs)

io_registry.register_writer('ascii', write_asciitable)


# IPAC
# ====

def read_ipac(filename, **kwargs):
    return read(filename, Reader=Ipac, **kwargs)

io_registry.register_reader('ipac', read_ipac)


def is_ipac(origin, args, kwargs):
    return args[0].endswith('.tbl')

io_registry.register_identifier('ipac', is_ipac)

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

io_registry.register_reader('latex', read_latex)


def write_latex(table, filename, **kwargs):
    write(table, filename, Writer=Latex, **kwargs)

io_registry.register_writer('latex', write_latex)


def is_latex(origin, args, kwargs):
    return args[0].endswith('.tex')

io_registry.register_identifier('latex', is_latex)


# RDB
# ===

def read_rdb(filename, **kwargs):
    return read(filename, Reader=Rdb, **kwargs)

io_registry.register_reader('rdb', read_rdb)


def write_rdb(table, filename, **kwargs):
    write(table, filename, Writer=Rdb, **kwargs)

io_registry.register_writer('rdb', write_rdb)


def is_rdb(origin, args, kwargs):
    return args[0].endswith('.rdb')

io_registry.register_identifier('rdb', is_rdb)
