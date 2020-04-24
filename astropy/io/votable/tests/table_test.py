# Licensed under a 3-clause BSD style license - see LICENSE.rst
"""
Test the conversion to/from astropy.table
"""
import io
import os

import pathlib
import pytest
import numpy as np

from astropy.config import set_temp_config, reload_config
from astropy.utils.data import get_pkg_data_filename, get_pkg_data_fileobj
from astropy.io.votable.table import parse, writeto
from astropy.io.votable import tree, conf
from astropy.io.votable.exceptions import VOWarning, W39, E25
from astropy.tests.helper import catch_warnings, raises
from astropy.utils.exceptions import AstropyDeprecationWarning


def test_table(tmpdir):
    # Read the VOTABLE
    votable = parse(get_pkg_data_filename('data/regression.xml'))
    table = votable.get_first_table()
    astropy_table = table.to_table()

    for name in table.array.dtype.names:
        assert np.all(astropy_table.mask[name] == table.array.mask[name])

    votable2 = tree.VOTableFile.from_table(astropy_table)
    t = votable2.get_first_table()

    field_types = [
        ('string_test', {'datatype': 'char', 'arraysize': '*'}),
        ('string_test_2', {'datatype': 'char', 'arraysize': '10'}),
        ('unicode_test', {'datatype': 'unicodeChar', 'arraysize': '*'}),
        ('fixed_unicode_test', {'datatype': 'unicodeChar', 'arraysize': '10'}),
        ('string_array_test', {'datatype': 'char', 'arraysize': '4'}),
        ('unsignedByte', {'datatype': 'unsignedByte'}),
        ('short', {'datatype': 'short'}),
        ('int', {'datatype': 'int'}),
        ('long', {'datatype': 'long'}),
        ('double', {'datatype': 'double'}),
        ('float', {'datatype': 'float'}),
        ('array', {'datatype': 'long', 'arraysize': '2*'}),
        ('bit', {'datatype': 'bit'}),
        ('bitarray', {'datatype': 'bit', 'arraysize': '3x2'}),
        ('bitvararray', {'datatype': 'bit', 'arraysize': '*'}),
        ('bitvararray2', {'datatype': 'bit', 'arraysize': '3x2*'}),
        ('floatComplex', {'datatype': 'floatComplex'}),
        ('doubleComplex', {'datatype': 'doubleComplex'}),
        ('doubleComplexArray', {'datatype': 'doubleComplex', 'arraysize': '*'}),
        ('doubleComplexArrayFixed', {'datatype': 'doubleComplex', 'arraysize': '2'}),
        ('boolean', {'datatype': 'bit'}),
        ('booleanArray', {'datatype': 'bit', 'arraysize': '4'}),
        ('nulls', {'datatype': 'int'}),
        ('nulls_array', {'datatype': 'int', 'arraysize': '2x2'}),
        ('precision1', {'datatype': 'double'}),
        ('precision2', {'datatype': 'double'}),
        ('doublearray', {'datatype': 'double', 'arraysize': '*'}),
        ('bitarray2', {'datatype': 'bit', 'arraysize': '16'})]

    for field, type in zip(t.fields, field_types):
        name, d = type
        assert field.ID == name
        assert field.datatype == d['datatype'], f'{name} expected {d["datatype"]} but get {field.datatype}'  # noqa
        if 'arraysize' in d:
            assert field.arraysize == d['arraysize']

    # W39: Bit values can not be masked
    with pytest.warns(W39):
        writeto(votable2, os.path.join(str(tmpdir), "through_table.xml"))


def test_read_through_table_interface(tmpdir):
    from astropy.table import Table

    with get_pkg_data_fileobj('data/regression.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable', table_id='main_table')

    assert len(t) == 5

    # Issue 8354
    assert t['float'].format is None

    fn = os.path.join(str(tmpdir), "table_interface.xml")

    # W39: Bit values can not be masked
    with pytest.warns(W39):
        t.write(fn, table_id='FOO', format='votable')

    with open(fn, 'rb') as fd:
        t2 = Table.read(fd, format='votable', table_id='FOO')

    assert len(t2) == 5


def test_read_through_table_interface2():
    from astropy.table import Table

    with get_pkg_data_fileobj('data/regression.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable', table_id='last_table')

    assert len(t) == 0


def test_names_over_ids():
    with get_pkg_data_fileobj('data/names.xml', encoding='binary') as fd:
        votable = parse(fd)

    table = votable.get_first_table().to_table(use_names_over_ids=True)

    assert table.colnames == [
        'Name', 'GLON', 'GLAT', 'RAdeg', 'DEdeg', 'Jmag', 'Hmag', 'Kmag',
        'G3.6mag', 'G4.5mag', 'G5.8mag', 'G8.0mag', '4.5mag', '8.0mag',
        'Emag', '24mag', 'f_Name']


def test_explicit_ids():
    with get_pkg_data_fileobj('data/names.xml', encoding='binary') as fd:
        votable = parse(fd)

    table = votable.get_first_table().to_table(use_names_over_ids=False)

    assert table.colnames == [
        'col1', 'col2', 'col3', 'col4', 'col5', 'col6', 'col7', 'col8', 'col9',
        'col10', 'col11', 'col12', 'col13', 'col14', 'col15', 'col16', 'col17']


def test_table_read_with_unnamed_tables():
    """
    Issue #927
    """
    from astropy.table import Table

    with get_pkg_data_fileobj('data/names.xml', encoding='binary') as fd:
        t = Table.read(fd, format='votable')

    assert len(t) == 1


def test_votable_path_object():
    """
    Testing when votable is passed as pathlib.Path object #4412.
    """
    fpath = pathlib.Path(get_pkg_data_filename('data/names.xml'))
    table = parse(fpath).get_first_table().to_table()

    assert len(table) == 1
    assert int(table[0][3]) == 266


def test_from_table_without_mask():
    from astropy.table import Table, Column
    t = Table()
    c = Column(data=[1, 2, 3], name='a')
    t.add_column(c)
    output = io.BytesIO()
    t.write(output, format='votable')


def test_write_with_format():
    from astropy.table import Table, Column
    t = Table()
    c = Column(data=[1, 2, 3], name='a')
    t.add_column(c)

    output = io.BytesIO()
    t.write(output, format='votable', tabledata_format="binary")
    obuff = output.getvalue()
    assert b'VOTABLE version="1.4"' in obuff
    assert b'BINARY' in obuff
    assert b'TABLEDATA' not in obuff

    output = io.BytesIO()
    t.write(output, format='votable', tabledata_format="binary2")
    obuff = output.getvalue()
    assert b'VOTABLE version="1.4"' in obuff
    assert b'BINARY2' in obuff
    assert b'TABLEDATA' not in obuff


def test_empty_table():
    votable = parse(get_pkg_data_filename('data/empty_table.xml'))
    table = votable.get_first_table()
    astropy_table = table.to_table()  # noqa


def test_no_field_not_empty_table():
    votable = parse(get_pkg_data_filename('data/no_field_not_empty_table.xml'))
    table = votable.get_first_table()
    assert len(table.fields) == 0
    assert len(table.infos) == 1


def test_no_field_not_empty_table_exception():
    with pytest.raises(E25):
        votable = parse(get_pkg_data_filename('data/no_field_not_empty_table.xml'), verify='exception')


def test_binary2_masked_strings():
    """
    Issue #8995
    """
    # Read a VOTable which sets the null mask bit for each empty string value.
    votable = parse(get_pkg_data_filename('data/binary2_masked_strings.xml'))
    table = votable.get_first_table()
    astropy_table = table.to_table()

    # Ensure string columns have no masked values and can be written out
    assert not np.any(table.array.mask['epoch_photometry_url'])
    output = io.BytesIO()
    astropy_table.write(output, format='votable')


class TestVerifyOptions:

    # Start off by checking the default (ignore)

    def test_default(self):
        with catch_warnings(VOWarning) as w:
            parse(get_pkg_data_filename('data/gemini.xml'))
        assert len(w) == 0

    # Then try the various explicit options

    def test_verify_ignore(self):
        with catch_warnings(VOWarning) as w:
            parse(get_pkg_data_filename('data/gemini.xml'), verify='ignore')
        assert len(w) == 0

    def test_verify_warn(self):
        with catch_warnings(VOWarning) as w:
            parse(get_pkg_data_filename('data/gemini.xml'), verify='warn')
        assert len(w) == 24

    def test_verify_exception(self):
        with pytest.raises(VOWarning):
            parse(get_pkg_data_filename('data/gemini.xml'), verify='exception')

    # Make sure the pedantic option still works for now (pending deprecation)

    def test_pedantic_false(self):
        with catch_warnings(VOWarning, AstropyDeprecationWarning) as w:
            parse(get_pkg_data_filename('data/gemini.xml'), pedantic=False)
        assert len(w) == 24
        # Make sure we don't yet emit a deprecation warning
        assert not any(isinstance(x.category, AstropyDeprecationWarning) for x in w)

    def test_pedantic_true(self):
        with pytest.raises(VOWarning):
            parse(get_pkg_data_filename('data/gemini.xml'), pedantic=True)

    # Make sure that the default behavior can be set via configuration items

    def test_conf_verify_ignore(self):
        with conf.set_temp('verify', 'ignore'):
            with catch_warnings(VOWarning) as w:
                parse(get_pkg_data_filename('data/gemini.xml'))
            assert len(w) == 0

    def test_conf_verify_warn(self):
        with conf.set_temp('verify', 'warn'):
            with catch_warnings(VOWarning) as w:
                parse(get_pkg_data_filename('data/gemini.xml'))
            assert len(w) == 24

    def test_conf_verify_exception(self):
        with conf.set_temp('verify', 'exception'):
            with pytest.raises(VOWarning):
                parse(get_pkg_data_filename('data/gemini.xml'))

    # And make sure the old configuration item will keep working

    def test_conf_pedantic_false(self, tmpdir):

        with set_temp_config(tmpdir.strpath):

            with open(tmpdir.join('astropy').join('astropy.cfg').strpath, 'w') as f:
                f.write('[io.votable]\npedantic = False')

            reload_config('astropy.io.votable')

            with catch_warnings(VOWarning, AstropyDeprecationWarning) as w:
                parse(get_pkg_data_filename('data/gemini.xml'))
            assert len(w) == 24
            # Make sure we don't yet emit a deprecation warning
            assert not any(isinstance(x.category, AstropyDeprecationWarning) for x in w)

    def test_conf_pedantic_true(self, tmpdir):

        with set_temp_config(tmpdir.strpath):

            with open(tmpdir.join('astropy').join('astropy.cfg').strpath, 'w') as f:
                f.write('[io.votable]\npedantic = True')

            reload_config('astropy.io.votable')

            with pytest.raises(VOWarning):
                parse(get_pkg_data_filename('data/gemini.xml'))
