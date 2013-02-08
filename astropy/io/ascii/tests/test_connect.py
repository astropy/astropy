# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from ....tests.helper import pytest

from ....table import Table, Column

from .common import numpy_lt_1p5
ROOT = os.path.abspath(os.path.dirname(__file__))


@pytest.mark.parametrize('filename', ['t/cds.dat', 't/ipac.dat',
                                      't/daophot.dat', 't/latex1.tex'])
def test_read_generic(filename):
    if numpy_lt_1p5 and filename in ['t/cds.dat', 't/ipac.dat']:
        return
    Table.read(os.path.join(ROOT, filename), format='ascii')


def test_write_generic(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1,2,3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    t.write(str(tmpdir.join("test")), format='ascii')


@pytest.mark.xfail('numpy_lt_1p5')
def test_read_ipac():
    Table.read(os.path.join(ROOT, 't/ipac.dat'), format='ipac')


@pytest.mark.xfail('numpy_lt_1p5')
def test_read_cds():
    Table.read(os.path.join(ROOT, 't/cds.dat'), format='cds')


def test_read_dapphot():
    Table.read(os.path.join(ROOT, 't/daophot.dat'), format='daophot')


def test_read_latex():
    Table.read(os.path.join(ROOT, 't/latex1.tex'), format='latex')


def test_read_latex_noformat():
    Table.read(os.path.join(ROOT, 't/latex1.tex'))


def test_write_latex(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1,2,3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path, format='latex')


def test_write_latex_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1,2,3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path)


def test_read_rdb():
    Table.read(os.path.join(ROOT, 't/short.rdb'), format='rdb')


def test_read_rdb_noformat():
    Table.read(os.path.join(ROOT, 't/short.rdb'))


def test_write_rdb(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1,2,3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path, format='rdb')


def test_write_rdb_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1,2,3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path)
