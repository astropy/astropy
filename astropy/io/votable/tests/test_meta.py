import os
from numpy.testing import assert_allclose

from astropy.coordinates import FK5, ICRS
from astropy.table import Table

DATA = os.path.dirname(__file__)


def test_aj285677t3_votable():

    tab = Table.read(os.path.join(DATA, 'aj285677t3_votable.xml'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_aj285677t3_votable_2019():

    tab = Table.read(os.path.join(DATA, 'aj285677t3_votable_2019.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_gaia():

    tab = Table.read(os.path.join(DATA, 'gaia.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['GAIADR2']

    assert isinstance(coosys['GAIADR2'], ICRS)


def test_simbad():

    tab = Table.read(os.path.join(DATA, 'simbad.xml'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['COOSYS']

    assert isinstance(coosys['COOSYS'], ICRS)


def test_table_from_ivoa_paper():

    tab = Table.read(os.path.join(DATA, 'table_from_paper.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['system']

    assert isinstance(coosys['system'], ICRS)


def test_vizier_gaia_main_table():

    tab = Table.read(os.path.join(DATA, 'vizier-gaia-main-table.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 3

    assert coosys[0] is coosys['H_2015.500']

    assert isinstance(coosys['H_2015.500'], ICRS)
    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)


def test_vizier_gaia_transits():

    tab = Table.read(os.path.join(DATA, 'vizier-gaia-transits.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 0


def test_vizier_panstarrs():

    tab = Table.read(os.path.join(DATA, 'vizier-panstarrs.vot'))

    coosys = tab.meta['votable']['coosys']

    assert len(coosys) == 2

    assert coosys[0] is coosys['J2000']

    assert isinstance(coosys['J2000'], FK5)
    assert_allclose(coosys['J2000'].equinox.jyear, 2000)
