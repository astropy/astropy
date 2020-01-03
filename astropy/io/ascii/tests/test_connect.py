# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os

import pytest

from astropy.table import Table, Column

from astropy.table.table_helpers import simple_table

import numpy as np

ROOT = os.path.abspath(os.path.dirname(__file__))

files = ['data/cds.dat', 'data/ipac.dat', 'data/daophot.dat', 'data/latex1.tex',
         'data/simple_csv.csv']

# Check to see if the BeautifulSoup dependency is present.

try:
    from bs4 import BeautifulSoup  # noqa
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False

try:
    import yaml  # noqa
    HAS_YAML = True
except ImportError:
    HAS_YAML = False

if HAS_BEAUTIFUL_SOUP:
    files.append('data/html.html')


@pytest.mark.parametrize('filename', files)
def test_read_generic(filename):
    Table.read(os.path.join(ROOT, filename), format='ascii')


def test_write_generic(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    t.write(str(tmpdir.join("test")), format='ascii')


def test_read_ipac():
    Table.read(os.path.join(ROOT, 'data/ipac.dat'), format='ipac')


def test_read_cds():
    Table.read(os.path.join(ROOT, 'data/cds.dat'), format='cds')


def test_read_dapphot():
    Table.read(os.path.join(ROOT, 'data/daophot.dat'), format='daophot')


def test_read_latex():
    Table.read(os.path.join(ROOT, 'data/latex1.tex'), format='latex')


def test_read_latex_noformat():
    Table.read(os.path.join(ROOT, 'data/latex1.tex'))


def test_write_latex(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path, format='latex')


def test_write_latex_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path)


@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_read_html():
    Table.read(os.path.join(ROOT, 'data/html.html'), format='html')


@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_read_html_noformat():
    Table.read(os.path.join(ROOT, 'data/html.html'))


def test_write_html(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.html"))
    t.write(path, format='html')


def test_write_html_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.html"))
    t.write(path)


def test_read_rdb():
    Table.read(os.path.join(ROOT, 'data/short.rdb'), format='rdb')


def test_read_rdb_noformat():
    Table.read(os.path.join(ROOT, 'data/short.rdb'))


def test_write_rdb(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path, format='rdb')


def test_write_rdb_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path)


def test_read_csv():
    '''If properly registered, filename should be sufficient to specify format

    #3189
    '''
    Table.read(os.path.join(ROOT, 'data/simple_csv.csv'))


def test_write_csv(tmpdir):
    '''If properly registered, filename should be sufficient to specify format

    #3189
    '''
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.csv"))
    t.write(path)


@pytest.mark.skipif('not HAS_YAML')
def test_auto_identify_ecsv(tmpdir):
    tbl = simple_table()
    tmpfile = str(tmpdir.join('/tmpFile.ecsv'))
    tbl.write(tmpfile)
    tbl2 = Table.read(tmpfile)
    assert np.all(tbl == tbl2)
