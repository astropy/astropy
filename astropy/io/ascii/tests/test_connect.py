# Licensed under a 3-clause BSD style license - see LICENSE.rst
import os
from ....tests.helper import pytest

from ....table import Table, Column

ROOT = os.path.abspath(os.path.dirname(__file__))

files = ['t/cds.dat', 't/ipac.dat', 't/daophot.dat', 't/latex1.tex',
         't/simple_csv.csv']

# Check to see if the BeautifulSoup dependency is present.

try:
    from bs4 import BeautifulSoup
    HAS_BEAUTIFUL_SOUP = True
except ImportError:
    HAS_BEAUTIFUL_SOUP = False

if HAS_BEAUTIFUL_SOUP:
    files.append('t/html.html')

@pytest.mark.parametrize('filename', files)

def test_read_generic(filename):
    Table.read(os.path.join(ROOT, filename), format='ascii')


def test_write_generic(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    t.write(str(tmpdir.join("test")), format='ascii')


def test_read_ipac():
    Table.read(os.path.join(ROOT, 't/ipac.dat'), format='ascii.ipac')


def test_read_cds():
    Table.read(os.path.join(ROOT, 't/cds.dat'), format='ascii.cds')


def test_read_dapphot():
    Table.read(os.path.join(ROOT, 't/daophot.dat'), format='ascii.daophot')


def test_read_latex():
    Table.read(os.path.join(ROOT, 't/latex1.tex'), format='ascii.latex')


def test_read_latex_noformat():
    Table.read(os.path.join(ROOT, 't/latex1.tex'))


def test_write_latex(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path, format='ascii.latex')


def test_write_latex_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.tex"))
    t.write(path)


@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_read_html():
    Table.read(os.path.join(ROOT, 't/html.html'), format='ascii.html')


@pytest.mark.skipif('not HAS_BEAUTIFUL_SOUP')
def test_read_html_noformat():
    Table.read(os.path.join(ROOT, 't/html.html'))


def test_write_html(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.html"))
    t.write(path, format='ascii.html')


def test_write_html_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.html"))
    t.write(path)


def test_read_rdb():
    Table.read(os.path.join(ROOT, 't/short.rdb'), format='ascii.rdb')


def test_read_rdb_noformat():
    Table.read(os.path.join(ROOT, 't/short.rdb'))


def test_write_rdb(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path, format='ascii.rdb')


def test_write_rdb_noformat(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.rdb"))
    t.write(path)


def test_read_csv():
    Table.read(os.path.join(ROOT, 't/simple_csv.csv'), format='ascii.csv')


def test_write_csv(tmpdir):
    t = Table()
    t.add_column(Column(name='a', data=[1, 2, 3]))
    t.add_column(Column(name='b', data=['a', 'b', 'c']))
    path = str(tmpdir.join("data.csv"))
    t.write(path, format='ascii.csv')

