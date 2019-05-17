import os
from numpy.testing import assert_allclose

from astropy.coordinates import FK5, ICRS
from astropy.table import Table
from astropy.utils.data import get_pkg_data_filename

from ..connect import _votable_meta_to_coo_frames

DATA = 'data'

votables_to_test = ['aj285677t3_votable.xml', 'aj285677t3_votable_2019.vot',
                    'gaia.vot', 'simbad.xml', 'table_from_paper.vot',
                    'vizier-gaia-main-table.vot', 'vizier-gaia-transits.vot',
                    'vizier-panstarrs.vot']


def test_aj285677t3_votable():
    tab = Table.read(get_pkg_data_filename('data/aj285677t3_votable.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_aj285677t3_votable_2019():

    tab = Table.read(get_pkg_data_filename('data/aj285677t3_votable_2019.vot.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_gaia():

    tab = Table.read(get_pkg_data_filename('data/gaia.vot'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['GAIADR2']

    assert isinstance(coosys['GAIADR2'], ICRS)


def test_simbad():

    tab = Table.read(get_pkg_data_filename('data/simbad.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['COOSYS']

    assert isinstance(coosys['COOSYS'], ICRS)


def test_table_from_ivoa_paper():

    tab = Table.read(get_pkg_data_filename('data/table_from_paper.vot'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['system']

    assert isinstance(coosys['system'], ICRS)


def test_vizier_gaia_main_table():

    tab = Table.read(get_pkg_data_filename('data/vizier-gaia-main-table.vot.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 4  # each of the two has a name and an index

    assert coosys[0] is coosys['H_2015.500']

    assert isinstance(coosys['H_2015.500'], ICRS)
    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_vizier_gaia_transits():

    tab = Table.read(get_pkg_data_filename('data/vizier-gaia-transits.vot.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 0


def test_vizier_panstarrs():

    tab = Table.read(get_pkg_data_filename('data/vizier-panstarrs.vot.xml'))

    coosys = _votable_meta_to_coo_frames(tab.meta['votable'])

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)
