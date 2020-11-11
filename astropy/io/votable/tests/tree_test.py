# Licensed under a 3-clause BSD style license - see LICENSE.rst

from astropy.io.votable import exceptions
from astropy.io.votable import tree
from astropy.tests.helper import raises


@raises(exceptions.W07)
def test_check_astroyear_fail():
    config = {'verify': 'exception'}
    field = tree.Field(None, name='astroyear', arraysize='1')
    tree.check_astroyear('X2100', field, config)


@raises(exceptions.W08)
def test_string_fail():
    config = {'verify': 'exception'}
    tree.check_string(42, 'foo', config)


def test_make_Fields():
    votable = tree.VOTableFile()
    # ...with one resource...
    resource = tree.Resource()
    votable.resources.append(resource)

    # ... with one table
    table = tree.Table(votable)
    resource.tables.append(table)

    table.fields.extend([tree.Field(
        votable, name='Test', datatype="float", unit="mag")])


def test_unit_format():
    from astropy.io.votable.table import parse
    from astropy.utils.data import get_pkg_data_filename

    data = parse(get_pkg_data_filename('data/irsa-nph-error.xml'))
    assert data._config['version'] == '1.0'
    assert tree._get_default_unit_format(data._config) == 'cds'
    data = parse(get_pkg_data_filename('data/names.xml'))
    assert data._config['version'] == '1.1'
    assert tree._get_default_unit_format(data._config) == 'cds'
    data = parse(get_pkg_data_filename('data/gemini.xml'))
    assert data._config['version'] == '1.2'
    assert tree._get_default_unit_format(data._config) == 'cds'
    data = parse(get_pkg_data_filename('data/binary2_masked_strings.xml'))
    assert data._config['version'] == '1.3'
    assert tree._get_default_unit_format(data._config) == 'cds'
    data = parse(get_pkg_data_filename('data/timesys.xml'))
    assert data._config['version'] == '1.4'
    assert tree._get_default_unit_format(data._config) == 'vounit'
