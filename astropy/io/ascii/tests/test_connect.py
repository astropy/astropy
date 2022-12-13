# Licensed under a 3-clause BSD style license - see LICENSE.rst
import numpy as np
import pytest

from astropy.table import Column, Table
from astropy.table.table_helpers import simple_table
from astropy.utils.compat.optional_deps import HAS_BS4
from astropy.utils.data import get_pkg_data_filename

files = [
    "data/cds.dat",
    "data/ipac.dat",
    "data/daophot.dat",
    "data/latex1.tex",
    "data/simple_csv.csv",
]


if HAS_BS4:
    files.append("data/html.html")


@pytest.mark.parametrize("filename", files)
def test_read_generic(filename):
    Table.read(get_pkg_data_filename(filename), format="ascii")


def test_write_generic(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    t.write(tmp_path / "test", format="ascii")


def test_read_ipac():
    Table.read(get_pkg_data_filename("data/ipac.dat"), format="ipac")


def test_read_cds():
    Table.read(get_pkg_data_filename("data/cds.dat"), format="cds")


def test_read_dapphot():
    Table.read(get_pkg_data_filename("data/daophot.dat"), format="daophot")


def test_read_latex():
    Table.read(get_pkg_data_filename("data/latex1.tex"), format="latex")


def test_read_latex_noformat():
    Table.read(get_pkg_data_filename("data/latex1.tex"))


def test_write_latex(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.tex"
    t.write(path, format="latex")


def test_write_latex_noformat(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.tex"
    t.write(path)


@pytest.mark.skipif(not HAS_BS4, reason="requires BeautifulSoup4")
def test_read_html():
    Table.read(get_pkg_data_filename("data/html.html"), format="html")


@pytest.mark.skipif(not HAS_BS4, reason="requires BeautifulSoup4")
def test_read_html_noformat():
    Table.read(get_pkg_data_filename("data/html.html"))


def test_write_html(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.html"
    t.write(path, format="html")


def test_write_html_noformat(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.html"
    t.write(path)


def test_read_rdb():
    Table.read(get_pkg_data_filename("data/short.rdb"), format="rdb")


def test_read_rdb_noformat():
    Table.read(get_pkg_data_filename("data/short.rdb"))


def test_write_rdb(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.rdb"
    t.write(path, format="rdb")


def test_write_rdb_noformat(tmp_path):
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.rdb"
    t.write(path)


def test_read_csv():
    """If properly registered, filename should be sufficient to specify format

    #3189
    """
    Table.read(get_pkg_data_filename("data/simple_csv.csv"))


def test_write_csv(tmp_path):
    """If properly registered, filename should be sufficient to specify format

    #3189
    """
    t = Table()
    t.add_column(Column(name="a", data=[1, 2, 3]))
    t.add_column(Column(name="b", data=["a", "b", "c"]))
    path = tmp_path / "data.csv"
    t.write(path)


def test_auto_identify_ecsv(tmp_path):
    tbl = simple_table()
    tmpfile = tmp_path / "tmpFile.ecsv"
    tbl.write(tmpfile)
    tbl2 = Table.read(tmpfile)
    assert np.all(tbl == tbl2)
