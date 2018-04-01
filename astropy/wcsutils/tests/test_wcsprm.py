# -*- coding: utf-8 -*-

# Licensed under a 3-clause BSD style license - see LICENSE.rst

import gc
import locale
import re

import pytest
from numpy.testing import assert_array_equal, assert_array_almost_equal
import numpy as np

from ...tests.helper import raises, catch_warnings
from ...io import fits
from .. import wcs
from .. import _wcs
from ...utils.data import get_pkg_data_contents, get_pkg_data_fileobj, get_pkg_data_filename
from ... import units as u


######################################################################


def test_alt():
    w = _wcs.Wcsprm()
    assert w.alt == " "
    w.alt = "X"
    assert w.alt == "X"
    del w.alt
    assert w.alt == " "


@raises(ValueError)
def test_alt_invalid1():
    w = _wcs.Wcsprm()
    w.alt = "$"


@raises(ValueError)
def test_alt_invalid2():
    w = _wcs.Wcsprm()
    w.alt = "  "


def test_axis_types():
    w = _wcs.Wcsprm()
    assert_array_equal(w.axis_types, [0, 0])


def test_cd():
    w = _wcs.Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert w.cd.dtype == float
    assert w.has_cd() is True
    assert_array_equal(w.cd, [[1, 0], [0, 1]])
    del w.cd
    assert w.has_cd() is False


@raises(AttributeError)
def test_cd_missing():
    w = _wcs.Wcsprm()
    assert w.has_cd() is False
    w.cd


@raises(AttributeError)
def test_cd_missing2():
    w = _wcs.Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert w.has_cd() is True
    del w.cd
    assert w.has_cd() is False
    w.cd


@raises(ValueError)
def test_cd_invalid():
    w = _wcs.Wcsprm()
    w.cd = [1, 0, 0, 1]


def test_cdfix():
    w = _wcs.Wcsprm()
    w.cdfix()


def test_cdelt():
    w = _wcs.Wcsprm()
    assert_array_equal(w.cdelt, [1, 1])
    w.cdelt = [42, 54]
    assert_array_equal(w.cdelt, [42, 54])


@raises(TypeError)
def test_cdelt_delete():
    w = _wcs.Wcsprm()
    del w.cdelt


def test_cel_offset():
    w = _wcs.Wcsprm()
    assert w.cel_offset is False
    w.cel_offset = 'foo'
    assert w.cel_offset is True
    w.cel_offset = 0
    assert w.cel_offset is False


def test_celfix():
    # TODO: We need some data with -NCP or -GLS projections to test
    # with.  For now, this is just a smoke test
    w = _wcs.Wcsprm()
    assert w.celfix() == -1


def test_cname():
    w = _wcs.Wcsprm()
    # Test that this works as an iterator
    for x in w.cname:
        assert x == ''
    assert list(w.cname) == ['', '']
    w.cname = [b'foo', 'bar']
    assert list(w.cname) == ['foo', 'bar']


@raises(TypeError)
def test_cname_invalid():
    w = _wcs.Wcsprm()
    w.cname = [42, 54]


def test_colax():
    w = _wcs.Wcsprm()
    assert w.colax.dtype == np.intc
    assert_array_equal(w.colax, [0, 0])
    w.colax = [42, 54]
    assert_array_equal(w.colax, [42, 54])
    w.colax[0] = 0
    assert_array_equal(w.colax, [0, 54])

    with pytest.raises(ValueError):
        w.colax = [1, 2, 3]


def test_colnum():
    w = _wcs.Wcsprm()
    assert w.colnum == 0
    w.colnum = 42
    assert w.colnum == 42

    with pytest.raises(OverflowError):
        w.colnum = 0xffffffffffffffffffff

    with pytest.raises(OverflowError):
        w.colnum = 0xffffffff

    with pytest.raises(TypeError):
        del w.colnum


@raises(TypeError)
def test_colnum_invalid():
    w = _wcs.Wcsprm()
    w.colnum = 'foo'


def test_crder():
    w = _wcs.Wcsprm()
    assert w.crder.dtype == float
    assert np.all(np.isnan(w.crder))
    w.crder[0] = 0
    assert np.isnan(w.crder[1])
    assert w.crder[0] == 0
    w.crder = w.crder


def test_crota():
    w = _wcs.Wcsprm()
    w.crota = [1, 0]
    assert w.crota.dtype == float
    assert w.has_crota() is True
    assert_array_equal(w.crota, [1, 0])
    del w.crota
    assert w.has_crota() is False


@raises(AttributeError)
def test_crota_missing():
    w = _wcs.Wcsprm()
    assert w.has_crota() is False
    w.crota


@raises(AttributeError)
def test_crota_missing2():
    w = _wcs.Wcsprm()
    w.crota = [1, 0]
    assert w.has_crota() is True
    del w.crota
    assert w.has_crota() is False
    w.crota


def test_crpix():
    w = _wcs.Wcsprm()
    assert w.crpix.dtype == float
    assert_array_equal(w.crpix, [0, 0])
    w.crpix = [42, 54]
    assert_array_equal(w.crpix, [42, 54])
    w.crpix[0] = 0
    assert_array_equal(w.crpix, [0, 54])

    with pytest.raises(ValueError):
        w.crpix = [1, 2, 3]


def test_crval():
    w = _wcs.Wcsprm()
    assert w.crval.dtype == float
    assert_array_equal(w.crval, [0, 0])
    w.crval = [42, 54]
    assert_array_equal(w.crval, [42, 54])
    w.crval[0] = 0
    assert_array_equal(w.crval, [0, 54])


def test_csyer():
    w = _wcs.Wcsprm()
    assert w.csyer.dtype == float
    assert np.all(np.isnan(w.csyer))
    w.csyer[0] = 0
    assert np.isnan(w.csyer[1])
    assert w.csyer[0] == 0
    w.csyer = w.csyer


def test_ctype():
    w = _wcs.Wcsprm()
    assert list(w.ctype) == ['', '']
    w.ctype = [b'RA---TAN', 'DEC--TAN']
    assert_array_equal(w.axis_types, [2200, 2201])
    assert w.lat == 1
    assert w.lng == 0
    assert w.lattyp == 'DEC'
    assert w.lngtyp == 'RA'
    assert list(w.ctype) == ['RA---TAN', 'DEC--TAN']
    w.ctype = ['foo', 'bar']
    assert_array_equal(w.axis_types, [0, 0])
    assert list(w.ctype) == ['foo', 'bar']
    assert w.lat == -1
    assert w.lng == -1
    assert w.lattyp == 'DEC'
    assert w.lngtyp == 'RA'


def test_ctype_repr():
    w = _wcs.Wcsprm()
    assert list(w.ctype) == ['', '']
    w.ctype = [b'RA-\t--TAN', 'DEC-\n-TAN']
    assert repr(w.ctype == '["RA-\t--TAN", "DEC-\n-TAN"]')


def test_ctype_index_error():
    w = _wcs.Wcsprm()
    assert list(w.ctype) == ['', '']
    with pytest.raises(IndexError):
        w.ctype[2] = 'FOO'


def test_ctype_invalid_error():
    w = _wcs.Wcsprm()
    assert list(w.ctype) == ['', '']
    with pytest.raises(ValueError):
        w.ctype[0] = 'X' * 100
    with pytest.raises(TypeError):
        w.ctype[0] = True
    with pytest.raises(TypeError):
        w.ctype = ['a', 0]
    with pytest.raises(TypeError):
        w.ctype = None
    with pytest.raises(ValueError):
        w.ctype = ['a', 'b', 'c']
    with pytest.raises(ValueError):
        w.ctype = ['FOO', 'A' * 100]


def test_cubeface():
    w = _wcs.Wcsprm()
    assert w.cubeface == -1
    w.cubeface = 0
    with pytest.raises(OverflowError):
        w.cubeface = -1


def test_cunit():
    w = _wcs.Wcsprm()
    assert list(w.cunit) == [u.Unit(''), u.Unit('')]
    w.cunit = [u.m, 'km']
    assert w.cunit[0] == u.m
    assert w.cunit[1] == u.km


def test_cunit_invalid():
    w = _wcs.Wcsprm()
    with catch_warnings() as warns:
        w.cunit[0] = 'foo'
    assert len(warns) == 1
    assert 'foo' in str(warns[0].message)


def test_cunit_invalid2():
    w = _wcs.Wcsprm()
    with catch_warnings() as warns:
        w.cunit = ['foo', 'bar']
    assert len(warns) == 2
    assert 'foo' in str(warns[0].message)
    assert 'bar' in str(warns[1].message)


def test_unit():
    w = wcs.WCS()
    w.wcs.cunit[0] = u.erg
    assert w.wcs.cunit[0] == u.erg

    assert repr(w.wcs.cunit) == "['erg', '']"


def test_unit2():
    w = wcs.WCS()
    myunit = u.Unit("FOOBAR", parse_strict="warn")
    w.wcs.cunit[0] = myunit


def test_unit3():
    w = wcs.WCS()
    with pytest.raises(IndexError):
        w.wcs.cunit[2] = u.m
    with pytest.raises(ValueError):
        w.wcs.cunit = [u.m, u.m, u.m]


def test_unitfix():
    w = _wcs.Wcsprm()
    w.unitfix()


def test_cylfix():
    # TODO: We need some data with broken cylindrical projections to
    # test with.  For now, this is just a smoke test.
    w = _wcs.Wcsprm()
    assert w.cylfix() == -1

    assert w.cylfix([0, 1]) == -1

    with pytest.raises(ValueError):
        w.cylfix([0, 1, 2])


def test_dateavg():
    w = _wcs.Wcsprm()
    assert w.dateavg == ''
    # TODO: When dateavg is verified, check that it works


def test_dateobs():
    w = _wcs.Wcsprm()
    assert w.dateobs == ''
    # TODO: When dateavg is verified, check that it works


def test_datfix():
    w = _wcs.Wcsprm()
    w.dateobs = '31/12/99'
    assert w.datfix() == 0
    assert w.dateobs == '1999-12-31'
    assert w.mjdobs == 51543.0


def test_equinox():
    w = _wcs.Wcsprm()
    assert np.isnan(w.equinox)
    w.equinox = 0
    assert w.equinox == 0
    del w.equinox
    assert np.isnan(w.equinox)

    with pytest.raises(TypeError):
        w.equinox = None


def test_fix():
    w = _wcs.Wcsprm()
    assert w.fix() == {
        'cdfix': 'No change',
        'cylfix': 'No change',
        'datfix': 'No change',
        'spcfix': 'No change',
        'unitfix': 'No change',
        'celfix': 'No change'}


def test_fix2():
    w = _wcs.Wcsprm()
    w.dateobs = '31/12/99'
    assert w.fix() == {
        'cdfix': 'No change',
        'cylfix': 'No change',
        'datfix': "Changed '31/12/99' to '1999-12-31'",
        'spcfix': 'No change',
        'unitfix': 'No change',
        'celfix': 'No change'}
    assert w.dateobs == '1999-12-31'
    assert w.mjdobs == 51543.0


def test_fix3():
    w = _wcs.Wcsprm()
    w.dateobs = '31/12/F9'
    assert w.fix() == {
        'cdfix': 'No change',
        'cylfix': 'No change',
        'datfix': "Invalid parameter value: invalid date '31/12/F9'",
        'spcfix': 'No change',
        'unitfix': 'No change',
        'celfix': 'No change'}
    assert w.dateobs == '31/12/F9'
    assert np.isnan(w.mjdobs)


def test_fix4():
    w = _wcs.Wcsprm()
    with pytest.raises(ValueError):
        w.fix('X')


def test_fix5():
    w = _wcs.Wcsprm()
    with pytest.raises(ValueError):
        w.fix(naxis=[0, 1, 2])


def test_get_ps():
    # TODO: We need some data with PSi_ma keywords
    w = _wcs.Wcsprm()
    assert len(w.get_ps()) == 0


def test_get_pv():
    # TODO: We need some data with PVi_ma keywords
    w = _wcs.Wcsprm()
    assert len(w.get_pv()) == 0


@raises(AssertionError)
def test_imgpix_matrix():
    w = _wcs.Wcsprm()
    w.imgpix_matrix


@raises(AttributeError)
def test_imgpix_matrix2():
    w = _wcs.Wcsprm()
    w.imgpix_matrix = None


def test_isunity():
    w = _wcs.Wcsprm()
    assert(w.is_unity())


def test_lat():
    w = _wcs.Wcsprm()
    assert w.lat == -1


@raises(AttributeError)
def test_lat_set():
    w = _wcs.Wcsprm()
    w.lat = 0


def test_latpole():
    w = _wcs.Wcsprm()
    assert w.latpole == 90.0
    w.latpole = 45.0
    assert w.latpole == 45.0
    del w.latpole
    assert w.latpole == 90.0


def test_lattyp():
    w = _wcs.Wcsprm()
    print(repr(w.lattyp))
    assert w.lattyp == "    "


@raises(AttributeError)
def test_lattyp_set():
    w = _wcs.Wcsprm()
    w.lattyp = 0


def test_lng():
    w = _wcs.Wcsprm()
    assert w.lng == -1


@raises(AttributeError)
def test_lng_set():
    w = _wcs.Wcsprm()
    w.lng = 0


def test_lngtyp():
    w = _wcs.Wcsprm()
    assert w.lngtyp == "    "


@raises(AttributeError)
def test_lngtyp_set():
    w = _wcs.Wcsprm()
    w.lngtyp = 0


def test_lonpole():
    w = _wcs.Wcsprm()
    assert np.isnan(w.lonpole)
    w.lonpole = 45.0
    assert w.lonpole == 45.0
    del w.lonpole
    assert np.isnan(w.lonpole)


def test_mix():
    w = _wcs.Wcsprm()
    w.ctype = [b'RA---TAN', 'DEC--TAN']
    with pytest.raises(_wcs.InvalidCoordinateError):
        w.mix(1, 1, [240, 480], 1, 5, [0, 2], [54, 32], 1)


def test_mjdavg():
    w = _wcs.Wcsprm()
    assert np.isnan(w.mjdavg)
    w.mjdavg = 45.0
    assert w.mjdavg == 45.0
    del w.mjdavg
    assert np.isnan(w.mjdavg)


def test_mjdobs():
    w = _wcs.Wcsprm()
    assert np.isnan(w.mjdobs)
    w.mjdobs = 45.0
    assert w.mjdobs == 45.0
    del w.mjdobs
    assert np.isnan(w.mjdobs)


def test_name():
    w = _wcs.Wcsprm()
    assert w.name == ''
    w.name = 'foo'
    assert w.name == 'foo'


def test_naxis():
    w = _wcs.Wcsprm()
    assert w.naxis == 2


@raises(AttributeError)
def test_naxis_set():
    w = _wcs.Wcsprm()
    w.naxis = 4


def test_obsgeo():
    w = _wcs.Wcsprm()
    assert np.all(np.isnan(w.obsgeo))
    w.obsgeo = [1, 2, 3]
    assert_array_equal(w.obsgeo, [1, 2, 3])
    del w.obsgeo
    assert np.all(np.isnan(w.obsgeo))


def test_pc():
    w = _wcs.Wcsprm()
    assert w.has_pc()
    assert_array_equal(w.pc, [[1, 0], [0, 1]])
    w.cd = [[1, 0], [0, 1]]
    assert not w.has_pc()
    del w.cd
    assert w.has_pc()
    assert_array_equal(w.pc, [[1, 0], [0, 1]])
    w.pc = w.pc


@raises(AttributeError)
def test_pc_missing():
    w = _wcs.Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert not w.has_pc()
    w.pc


def test_phi0():
    w = _wcs.Wcsprm()
    assert np.isnan(w.phi0)
    w.phi0 = 42.0
    assert w.phi0 == 42.0
    del w.phi0
    assert np.isnan(w.phi0)


@raises(AssertionError)
def test_piximg_matrix():
    w = _wcs.Wcsprm()
    w.piximg_matrix


@raises(AttributeError)
def test_piximg_matrix2():
    w = _wcs.Wcsprm()
    w.piximg_matrix = None


def test_print_contents():
    # In general, this is human-consumable, so we don't care if the
    # content changes, just check the type
    w = _wcs.Wcsprm()
    assert isinstance(str(w), str)


def test_radesys():
    w = _wcs.Wcsprm()
    assert w.radesys == ''
    w.radesys = 'foo'
    assert w.radesys == 'foo'


def test_restfrq():
    w = _wcs.Wcsprm()
    assert w.restfrq == 0.0
    w.restfrq = np.nan
    assert np.isnan(w.restfrq)
    del w.restfrq


def test_restwav():
    w = _wcs.Wcsprm()
    assert w.restwav == 0.0
    w.restwav = np.nan
    assert np.isnan(w.restwav)
    del w.restwav


def test_set_ps():
    w = _wcs.Wcsprm()
    data = [(0, 0, "param1"), (1, 1, "param2")]
    w.set_ps(data)
    assert w.get_ps() == data


def test_set_ps_realloc():
    w = _wcs.Wcsprm()
    w.set_ps([(0, 0, "param1")] * 16)


def test_set_pv():
    w = _wcs.Wcsprm()
    data = [(0, 0, 42.), (1, 1, 54.)]
    w.set_pv(data)
    assert w.get_pv() == data


def test_set_pv_realloc():
    w = _wcs.Wcsprm()
    w.set_pv([(0, 0, 42.)] * 16)


def test_spcfix():
    # TODO: We need some data with broken spectral headers here to
    # really test
    header = get_pkg_data_contents(
        'spectra/orion-velo-1.hdr', encoding='binary')
    w = _wcs.Wcsprm(header)
    assert w.spcfix() == -1


def test_spec():
    w = _wcs.Wcsprm()
    assert w.spec == -1


@raises(AttributeError)
def test_spec_set():
    w = _wcs.Wcsprm()
    w.spec = 0


def test_specsys():
    w = _wcs.Wcsprm()
    assert w.specsys == ''
    w.specsys = 'foo'
    assert w.specsys == 'foo'


def test_sptr():
    # TODO: Write me
    pass


def test_ssysobs():
    w = _wcs.Wcsprm()
    assert w.ssysobs == ''
    w.ssysobs = 'foo'
    assert w.ssysobs == 'foo'


def test_ssyssrc():
    w = _wcs.Wcsprm()
    assert w.ssyssrc == ''
    w.ssyssrc = 'foo'
    assert w.ssyssrc == 'foo'


def test_tab():
    w = _wcs.Wcsprm()
    assert len(w.tab) == 0
    # TODO: Inject some headers that have tables and test


def test_theta0():
    w = _wcs.Wcsprm()
    assert np.isnan(w.theta0)
    w.theta0 = 42.0
    assert w.theta0 == 42.0
    del w.theta0
    assert np.isnan(w.theta0)


def test_toheader():
    w = _wcs.Wcsprm()
    assert isinstance(w.to_header(), str)


def test_velangl():
    w = _wcs.Wcsprm()
    assert np.isnan(w.velangl)
    w.velangl = 42.0
    assert w.velangl == 42.0
    del w.velangl
    assert np.isnan(w.velangl)


def test_velosys():
    w = _wcs.Wcsprm()
    assert np.isnan(w.velosys)
    w.velosys = 42.0
    assert w.velosys == 42.0
    del w.velosys
    assert np.isnan(w.velosys)


def test_velref():
    w = _wcs.Wcsprm()
    assert w.velref == 0.0
    w.velref = 42.0
    assert w.velref == 42.0
    del w.velref
    assert w.velref == 0.0


def test_zsource():
    w = _wcs.Wcsprm()
    assert np.isnan(w.zsource)
    w.zsource = 42.0
    assert w.zsource == 42.0
    del w.zsource
    assert np.isnan(w.zsource)


def test_cd_3d():
    header = get_pkg_data_contents('data/3d_cd.hdr', encoding='binary')
    w = _wcs.Wcsprm(header)
    assert w.cd.shape == (3, 3)
    assert w.get_pc().shape == (3, 3)
    assert w.get_cdelt().shape == (3,)


def test_get_pc():
    header = get_pkg_data_contents('data/3d_cd.hdr', encoding='binary')
    w = _wcs.Wcsprm(header)
    pc = w.get_pc()
    try:
        pc[0, 0] = 42
    except (RuntimeError, ValueError):
        pass
    else:
        raise AssertionError()


@raises(_wcs.SingularMatrixError)
def test_detailed_err():
    w = _wcs.Wcsprm()
    w.pc = [[0, 0], [0, 0]]
    w.set()


def test_header_parse():
    from ...io import fits
    with get_pkg_data_fileobj(
            'data/header_newlines.fits', encoding='binary') as test_file:
        hdulist = fits.open(test_file)
        w = wcs.WCS(hdulist[0].header)
    assert w.wcs.ctype[0] == 'RA---TAN-SIP'


def test_locale():
    orig_locale = locale.getlocale(locale.LC_NUMERIC)[0]

    try:
        locale.setlocale(locale.LC_NUMERIC, 'fr_FR')
    except locale.Error:
        pytest.xfail(
            "Can't set to 'fr_FR' locale, perhaps because it is not installed "
            "on this system")
    try:
        header = get_pkg_data_contents('data/locale.hdr', encoding='binary')
        w = _wcs.Wcsprm(header)
        assert re.search("[0-9]+,[0-9]*", w.to_header()) is None
    finally:
        if orig_locale is None:
            # reset to the default setting
            locale.resetlocale(locale.LC_NUMERIC)
        else:
            # restore to whatever the previous value had been set to for
            # whatever reason
            locale.setlocale(locale.LC_NUMERIC, orig_locale)


@raises(UnicodeEncodeError)
def test_unicode():
    w = _wcs.Wcsprm()
    w.alt = "‰"


def test_sub_segfault():
    # Issue #1960
    header = fits.Header.fromtextfile(
        get_pkg_data_filename('data/sub-segfault.hdr'))
    w = wcs.WCS(header)
    sub = w.sub([wcs.WCSSUB_CELESTIAL])
    gc.collect()


def test_bounds_check():
    w = _wcs.Wcsprm()
    w.bounds_check(False)


def test_wcs_sub_error_message():
    # Issue #1587
    w = _wcs.Wcsprm()
    with pytest.raises(TypeError) as e:
        w.sub('latitude')
    assert str(e).endswith("axes must None, a sequence or an integer")


def test_wcs_sub():
    # Issue #3356
    w = _wcs.Wcsprm()
    w.sub(['latitude'])

    w = _wcs.Wcsprm()
    w.sub([b'latitude'])


def test_compare():
    header = get_pkg_data_contents('data/3d_cd.hdr', encoding='binary')
    w = _wcs.Wcsprm(header)
    w2 = _wcs.Wcsprm(header)

    assert w == w2

    w.equinox = 42
    assert w == w2

    assert not w.compare(w2)
    assert w.compare(w2, _wcs.WCSCOMPARE_ANCILLARY)

    w = _wcs.Wcsprm(header)
    w2 = _wcs.Wcsprm(header)

    w.cdelt[0] = np.float32(0.00416666666666666666666666)
    w2.cdelt[0] = np.float64(0.00416666666666666666666666)

    assert not w.compare(w2)
    assert w.compare(w2, tolerance=1e-6)


def test_radesys_defaults():
    w = _wcs.Wcsprm()
    w.ctype = ['RA---TAN', 'DEC--TAN']
    w.set()
    assert w.radesys == "ICRS"


def test_radesys_defaults_full():

    # As described in Section 3.1 of the FITS standard "Equatorial and ecliptic
    # coordinates", for those systems the RADESYS keyword can be used to
    # indicate the equatorial/ecliptic frame to use. From the standard:

    # "For RADESYSa values of FK4 and FK4-NO-E, any stated equinox is Besselian
    # and, if neither EQUINOXa nor EPOCH are given, a default of 1950.0 is to
    # be taken. For FK5, any stated equinox is Julian and, if neither keyword
    # is given, it defaults to 2000.0.

    # "If the EQUINOXa keyword is given it should always be accompanied by
    # RADESYS a. However, if it should happen to ap- pear by itself then
    # RADESYSa defaults to FK4 if EQUINOXa < 1984.0, or to FK5 if EQUINOXa
    # 1984.0. Note that these defaults, while probably true of older files
    # using the EPOCH keyword, are not required of them.

    # By default RADESYS is empty
    w = _wcs.Wcsprm(naxis=2)
    assert w.radesys == ''
    assert np.isnan(w.equinox)

    # For non-ecliptic or equatorial systems it is still empty
    w = _wcs.Wcsprm(naxis=2)
    for ctype in [('GLON-CAR', 'GLAT-CAR'),
                  ('SLON-SIN', 'SLAT-SIN')]:
        w.ctype = ctype
        w.set()
        assert w.radesys == ''
        assert np.isnan(w.equinox)

    for ctype in [('RA---TAN', 'DEC--TAN'),
                  ('ELON-TAN', 'ELAT-TAN'),
                  ('DEC--TAN', 'RA---TAN'),
                  ('ELAT-TAN', 'ELON-TAN')]:

        # Check defaults for RADESYS
        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.set()
        assert w.radesys == 'ICRS'

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.equinox = 1980
        w.set()
        assert w.radesys == 'FK4'

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.equinox = 1984
        w.set()
        assert w.radesys == 'FK5'

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.radesys = 'foo'
        w.set()
        assert w.radesys == 'foo'

        # Check defaults for EQUINOX
        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.set()
        assert np.isnan(w.equinox)  # frame is ICRS, no equinox

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.radesys = 'ICRS'
        w.set()
        assert np.isnan(w.equinox)

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.radesys = 'FK5'
        w.set()
        assert w.equinox == 2000.

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.radesys = 'FK4'
        w.set()
        assert w.equinox == 1950

        w = _wcs.Wcsprm(naxis=2)
        w.ctype = ctype
        w.radesys = 'FK4-NO-E'
        w.set()
        assert w.equinox == 1950


def test_iteration():
    world = np.array(
        [[-0.58995335, -0.5],
         [0.00664326, -0.5],
         [-0.58995335, -0.25],
         [0.00664326, -0.25],
         [-0.58995335, 0.],
         [0.00664326, 0.],
         [-0.58995335, 0.25],
         [0.00664326, 0.25],
         [-0.58995335, 0.5],
         [0.00664326, 0.5]],
        float
    )

    w = wcs.WCS()
    w.wcs.ctype = ['GLON-CAR', 'GLAT-CAR']
    w.wcs.cdelt = [-0.006666666828, 0.006666666828]
    w.wcs.crpix = [75.907, 74.8485]
    x = w.wcs_world2pix(world, 1)

    expected = np.array(
        [[1.64400000e+02, -1.51498185e-01],
         [7.49105110e+01, -1.51498185e-01],
         [1.64400000e+02, 3.73485009e+01],
         [7.49105110e+01, 3.73485009e+01],
         [1.64400000e+02, 7.48485000e+01],
         [7.49105110e+01, 7.48485000e+01],
         [1.64400000e+02, 1.12348499e+02],
         [7.49105110e+01, 1.12348499e+02],
         [1.64400000e+02, 1.49848498e+02],
         [7.49105110e+01, 1.49848498e+02]],
        float)

    assert_array_almost_equal(x, expected)

    w2 = w.wcs_pix2world(x, 1)

    world[:, 0] %= 360.

    assert_array_almost_equal(w2, world)


def test_invalid_args():
    with pytest.raises(TypeError):
        w = _wcs.Wcsprm(keysel='A')

    with pytest.raises(ValueError):
        w = _wcs.Wcsprm(keysel=2)

    with pytest.raises(ValueError):
        w = _wcs.Wcsprm(colsel=2)

    with pytest.raises(ValueError):
        w = _wcs.Wcsprm(naxis=64)

    header = get_pkg_data_contents(
        'spectra/orion-velo-1.hdr', encoding='binary')
    with pytest.raises(ValueError):
        w = _wcs.Wcsprm(header, relax='FOO')

    with pytest.raises(ValueError):
        w = _wcs.Wcsprm(header, naxis=3)

    with pytest.raises(KeyError):
        w = _wcs.Wcsprm(header, key='A')
