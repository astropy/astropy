# Licensed under a 3-clause BSD style license - see LICENSE.rst
# This file connects the readers/writers to the astropy.table.Table class



from .. import registry as io_registry
from ...table import Table

__all__ = []


# Generic
# =======


def read_asciitable(filename, **kwargs):
    from .ui import read
    return read(filename, **kwargs)

io_registry.register_reader('ascii', Table, read_asciitable)


def write_asciitable(table, filename, **kwargs):
    from .ui import write
    return write(table, filename, **kwargs)

io_registry.register_writer('ascii', Table, write_asciitable)


# IPAC
# ====

def read_ipac(filename, **kwargs):
    from .ipac import Ipac
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=Ipac, **kwargs)

io_registry.register_reader('ipac', Table, read_ipac)


# CDS
# ===


def read_cds(filename, **kwargs):
    from .cds import Cds
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=Cds, **kwargs)

io_registry.register_reader('cds', Table, read_cds)


# DAOPhot
# =======

def read_daophot(filename, **kwargs):
    from .daophot import Daophot
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=Daophot, **kwargs)

io_registry.register_reader('daophot', Table, read_daophot)

# SExtractor
# =======

def read_sextractor(filename, **kwargs):
    from .sextractor import SExtractor
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=SExtractor, **kwargs)

io_registry.register_reader('sextractor', Table, read_sextractor)

# LaTeX
# =====

def read_latex(filename, **kwargs):
    from .latex import Latex
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=Latex, **kwargs)

io_registry.register_reader('latex', Table, read_latex)


def write_latex(table, filename, **kwargs):
    from .latex import Latex
    from .ui import write
    write(table, filename, Writer=Latex, **kwargs)

io_registry.register_writer('latex', Table, write_latex)


def is_latex(origin, path, fileobj, *args, **kwargs):
    return path is not None and path.endswith('.tex')

io_registry.register_identifier('latex', Table, is_latex)


# RDB
# ===

def read_rdb(filename, **kwargs):
    from .basic import Rdb
    from .ui import read
    if 'guess' not in kwargs:
        kwargs['guess'] = False
    return read(filename, Reader=Rdb, **kwargs)

io_registry.register_reader('rdb', Table, read_rdb)


def write_rdb(table, filename, **kwargs):
    from .basic import Rdb
    from .ui import write
    write(table, filename, Writer=Rdb, **kwargs)

io_registry.register_writer('rdb', Table, write_rdb)


def is_rdb(origin, path, fileobj, *args, **kwargs):
    return path is not None and path.endswith('.rdb')

io_registry.register_identifier('rdb', Table, is_rdb)
