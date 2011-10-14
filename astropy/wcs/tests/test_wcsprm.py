import os
import sys

from nose.tools import raises
from numpy.testing import assert_array_equal
import numpy as np

from astropy import wcs
from astropy.wcs import _wcs


def b(s):
    return s.encode('ascii')

######################################################################

ROOT_DIR = None


def setup():
    global ROOT_DIR
    ROOT_DIR = os.path.dirname(__file__)


def test_alt():
    w = _wcs._Wcsprm()
    assert w.alt == b(" ")
    w.alt = b("X")
    assert w.alt == b("X")
    del w.alt
    assert w.alt == b(" ")


@raises(ValueError)
def test_alt_invalid1():
    w = _wcs._Wcsprm()
    w.alt = b("$")


@raises(ValueError)
def test_alt_invalid2():
    w = _wcs._Wcsprm()
    w.alt = b("  ")


def test_axis_types():
    w = _wcs._Wcsprm()
    assert_array_equal(w.axis_types, [0, 0])


def test_cd():
    w = _wcs._Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert w.cd.dtype == np.float
    assert w.has_cd() == True
    assert_array_equal(w.cd, [[1, 0], [0, 1]])
    del w.cd
    assert w.has_cd() == False


@raises(AttributeError)
def test_cd_missing():
    w = _wcs._Wcsprm()
    assert w.has_cd() == False
    w.cd


@raises(AttributeError)
def test_cd_missing2():
    w = _wcs._Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert w.has_cd() == True
    del w.cd
    assert w.has_cd() == False
    w.cd


@raises(ValueError)
def test_cd_invalid():
    w = _wcs._Wcsprm()
    w.cd = [1, 0, 0, 1]


def test_cdelt():
    w = _wcs._Wcsprm()
    assert_array_equal(w.cdelt, [1, 1])
    w.cdelt = [42, 54]
    assert_array_equal(w.cdelt, [42, 54])


@raises(TypeError)
def test_cdelt_delete():
    w = _wcs._Wcsprm()
    del w.cdelt


def test_cel_offset():
    w = _wcs._Wcsprm()
    assert w.cel_offset is False
    w.cel_offset = 'foo'
    assert w.cel_offset is True
    w.cel_offset = 0
    assert w.cel_offset is False


def test_celfix():
    # TODO: We need some data with -NCP or -GLS projections to test
    # with.  For now, this is just a smoke test
    w = _wcs._Wcsprm()
    assert w.celfix() == -1


def test_cname():
    w = _wcs._Wcsprm()
    # Test that this works as an iterator
    for x in w.cname:
        assert x == b('')
    assert list(w.cname) == [b(''), b('')]
    w.cname = [b('foo'), b('bar')]
    assert list(w.cname) == [b('foo'), b('bar')]


@raises(TypeError)
def test_cname_invalid():
    w = _wcs._Wcsprm()
    w.cname = [42, 54]


def test_colax():
    w = _wcs._Wcsprm()
    assert w.colax.dtype == np.intc
    assert_array_equal(w.colax, [0, 0])
    w.colax = [42, 54]
    assert_array_equal(w.colax, [42, 54])
    w.colax[0] = 0
    assert_array_equal(w.colax, [0, 54])


def test_colnum():
    w = _wcs._Wcsprm()
    assert w.colnum == 0
    w.colnum = 42
    assert w.colnum == 42


@raises(TypeError)
def test_colnum_invalid():
    w = _wcs._Wcsprm()
    w.colnum = 'foo'


def test_crder():
    w = _wcs._Wcsprm()
    assert w.crder.dtype == np.float
    assert np.all(np.isnan(w.crder))
    w.crder[0] = 0
    assert np.isnan(w.crder[1])
    assert w.crder[0] == 0


def test_crota():
    w = _wcs._Wcsprm()
    w.crota = [1, 0]
    assert w.crota.dtype == np.float
    assert w.has_crota() == True
    assert_array_equal(w.crota, [1, 0])
    del w.crota
    assert w.has_crota() == False


@raises(AttributeError)
def test_crota_missing():
    w = _wcs._Wcsprm()
    assert w.has_crota() == False
    w.crota


@raises(AttributeError)
def test_crota_missing2():
    w = _wcs._Wcsprm()
    w.crota = [1, 0]
    assert w.has_crota() == True
    del w.crota
    assert w.has_crota() == False
    w.crota


def test_crpix():
    w = _wcs._Wcsprm()
    assert w.crpix.dtype == np.float
    assert_array_equal(w.crpix, [0, 0])
    w.crpix = [42, 54]
    assert_array_equal(w.crpix, [42, 54])
    w.crpix[0] = 0
    assert_array_equal(w.crpix, [0, 54])


def test_crval():
    w = _wcs._Wcsprm()
    assert w.crval.dtype == np.float
    assert_array_equal(w.crval, [0, 0])
    w.crval = [42, 54]
    assert_array_equal(w.crval, [42, 54])
    w.crval[0] = 0
    assert_array_equal(w.crval, [0, 54])


def test_csyer():
    w = _wcs._Wcsprm()
    assert w.crder.dtype == np.float
    assert np.all(np.isnan(w.crder))
    w.crder[0] = 0
    assert np.isnan(w.crder[1])
    assert w.crder[0] == 0


def test_ctype():
    w = _wcs._Wcsprm()
    assert list(w.ctype) == [b(''), b('')]
    w.ctype = [b('RA---TAN'), b('DEC--TAN')]
    assert_array_equal(w.axis_types, [2200, 2201])
    assert w.lat == 1
    assert w.lng == 0
    assert w.lattyp == b('DEC')
    assert w.lngtyp == b('RA')
    assert list(w.ctype) == [b('RA---TAN'), b('DEC--TAN')]
    w.ctype = [b('foo'), b('bar')]
    assert_array_equal(w.axis_types, [0, 0])
    assert list(w.ctype) == [b('foo'), b('bar')]
    assert w.lat == -1
    assert w.lng == -1
    assert w.lattyp == b('DEC')
    assert w.lngtyp == b('RA')


def test_cubeface():
    w = _wcs._Wcsprm()
    assert w.cubeface == -1


def test_cunit():
    w = _wcs._Wcsprm()
    assert list(w.cunit) == [b(''), b('')]
    w.cunit = [b('m'), b('km')]


@raises(ValueError)
def test_cunit_invalid():
    w = _wcs._Wcsprm()
    w.cunit[0] = b('foo')


@raises(ValueError)
def test_cunit_invalid2():
    w = _wcs._Wcsprm()
    w.cunit = [b('foo'), b('bar')]


def test_cylfix():
    # TODO: We need some data with broken cylindrical projections to
    # test with.  For now, this is just a smoke test.
    w = _wcs._Wcsprm()
    assert w.cylfix() == -1


def test_dateavg():
    w = _wcs._Wcsprm()
    assert w.dateavg == b('')
    # TODO: When dateavg is verified, check that it works


def test_dateobs():
    w = _wcs._Wcsprm()
    assert w.dateobs == b('')
    # TODO: When dateavg is verified, check that it works


def test_datfix():
    w = _wcs._Wcsprm()
    w.dateobs = b('31/12/99')
    assert w.datfix() == 0
    assert w.dateobs == b('1999-12-31')
    assert w.mjdobs == 51543.0


def test_equinox():
    w = _wcs._Wcsprm()
    assert np.isnan(w.equinox)
    w.equinox = 0
    assert w.equinox == 0
    del w.equinox
    assert np.isnan(w.equinox)


def test_fix():
    w = _wcs._Wcsprm()
    assert w.fix() == {
        'cylfix': b('No change'),
        'datfix': b('No change'),
        'spcfix': b('No change'),
        'unitfix': b('No change'),
        'celfix': b('No change')}


def test_fix2():
    w = _wcs._Wcsprm()
    w.dateobs = b('31/12/99')
    assert w.fix() == {
        'cylfix': b('No change'),
        'datfix': b('Success'),
        'spcfix': b('No change'),
        'unitfix': b('No change'),
        'celfix': b('No change')}
    assert w.dateobs == b('1999-12-31')
    assert w.mjdobs == 51543.0


def test_fix3():
    w = _wcs._Wcsprm()
    w.dateobs = b('31/12/F9')
    assert w.fix() == {
        'cylfix': b('No change'),
        'datfix': b("Invalid parameter value: invalid date '31/12/F9'"),
        'spcfix': b('No change'),
        'unitfix': b('No change'),
        'celfix': b('No change')}
    assert w.dateobs == b('31/12/F9')
    assert np.isnan(w.mjdobs)


def test_get_ps():
    # TODO: We need some data with PSi_ma keywords
    w = _wcs._Wcsprm()
    assert len(w.get_ps()) == 0


def test_get_pv():
    # TODO: We need some data with PVi_ma keywords
    w = _wcs._Wcsprm()
    assert len(w.get_pv()) == 0


@raises(AssertionError)
def test_imgpix_matrix():
    w = _wcs._Wcsprm()
    w.imgpix_matrix


@raises(AttributeError)
def test_imgpix_matrix():
    w = _wcs._Wcsprm()
    w.imgpix_matrix = None


def test_isunity():
    w = _wcs._Wcsprm()
    assert(w.is_unity())


def test_lat():
    w = _wcs._Wcsprm()
    assert w.lat == -1


@raises(AttributeError)
def test_lat_set():
    w = _wcs._Wcsprm()
    w.lat = 0


def test_latpole():
    w = _wcs._Wcsprm()
    assert w.latpole == 90.0
    w.latpole = 45.0
    assert w.latpole == 45.0
    del w.latpole
    assert w.latpole == 90.0


def test_lattyp():
    w = _wcs._Wcsprm()
    print(repr(w.lattyp))
    assert w.lattyp == b("    ")


@raises(AttributeError)
def test_lattyp_set():
    w = _wcs._Wcsprm()
    w.lattyp = 0


def test_lng():
    w = _wcs._Wcsprm()
    assert w.lng == -1


@raises(AttributeError)
def test_lng_set():
    w = _wcs._Wcsprm()
    w.lng = 0


def test_lngtyp():
    w = _wcs._Wcsprm()
    assert w.lngtyp == b("    ")


@raises(AttributeError)
def test_lngtyp_set():
    w = _wcs._Wcsprm()
    w.lngtyp = 0


def test_lonpole():
    w = _wcs._Wcsprm()
    assert np.isnan(w.lonpole)
    w.lonpole = 45.0
    assert w.lonpole == 45.0
    del w.lonpole
    assert np.isnan(w.lonpole)


def test_mjdavg():
    w = _wcs._Wcsprm()
    assert np.isnan(w.mjdavg)
    w.mjdavg = 45.0
    assert w.mjdavg == 45.0
    del w.mjdavg
    assert np.isnan(w.mjdavg)


def test_mjdobs():
    w = _wcs._Wcsprm()
    assert np.isnan(w.mjdobs)
    w.mjdobs = 45.0
    assert w.mjdobs == 45.0
    del w.mjdobs
    assert np.isnan(w.mjdobs)


def test_name():
    w = _wcs._Wcsprm()
    assert w.name == b('')
    w.name = b('foo')
    assert w.name == b('foo')


def test_naxis():
    w = _wcs._Wcsprm()
    assert w.naxis == 2


@raises(AttributeError)
def test_naxis_set():
    w = _wcs._Wcsprm()
    w.naxis = 4


def test_obsgeo():
    w = _wcs._Wcsprm()
    assert np.all(np.isnan(w.obsgeo))
    w.obsgeo = [1, 2, 3]
    assert_array_equal(w.obsgeo, [1, 2, 3])
    del w.obsgeo
    assert np.all(np.isnan(w.obsgeo))


def test_pc():
    w = _wcs._Wcsprm()
    assert w.has_pc()
    assert_array_equal(w.pc, [[1, 0], [0, 1]])
    w.cd = [[1, 0], [0, 1]]
    assert not w.has_pc()
    del w.cd
    assert w.has_pc()
    assert_array_equal(w.pc, [[1, 0], [0, 1]])


@raises(AttributeError)
def test_pc_missing():
    w = _wcs._Wcsprm()
    w.cd = [[1, 0], [0, 1]]
    assert not w.has_pc()
    w.pc


def test_phi0():
    w = _wcs._Wcsprm()
    assert np.isnan(w.phi0)
    w.phi0 = 42.0
    assert w.phi0 == 42.0
    del w.phi0
    assert np.isnan(w.phi0)


@raises(AssertionError)
def test_piximg_matrix():
    w = _wcs._Wcsprm()
    w.piximg_matrix


@raises(AttributeError)
def test_piximg_matrix():
    w = _wcs._Wcsprm()
    w.piximg_matrix = None


def test_print_contents():
    # In general, this is human-consumable, so we don't care if the
    # content changes, just check the type
    w = _wcs._Wcsprm()
    assert isinstance(str(w), str)


def test_radesys():
    w = _wcs._Wcsprm()
    assert w.radesys == b('')
    w.radesys = b('foo')
    assert w.radesys == b('foo')


def test_restfrq():
    w = _wcs._Wcsprm()
    assert w.restfrq == 0.0
    w.restfrq = np.nan
    assert np.isnan(w.restfrq)


def test_restwav():
    w = _wcs._Wcsprm()
    assert w.restwav == 0.0
    w.restwav = np.nan
    assert np.isnan(w.restwav)


def test_set_ps():
    w = _wcs._Wcsprm()
    data = [(0, 0, "param1"), (1, 1, "param2")]
    w.set_ps(data)
    assert w.get_ps() == data


def test_set_ps_realloc():
    w = _wcs._Wcsprm()
    w.set_ps([(0, 0, "param1")] * 16)


def test_set_pv():
    w = _wcs._Wcsprm()
    data = [(0, 0, 42.), (1, 1, 54.)]
    w.set_pv(data)
    assert w.get_pv() == data


def test_set_pv_realloc():
    w = _wcs._Wcsprm()
    w.set_pv([(0, 0, 42.)] * 16)


def test_spcfix():
    # TODO: We need some data with broken spectral headers here to
    # really test
    header = open(os.path.join(ROOT_DIR, 'spectra', 'orion-velo-1.hdr'),
                  'rb').read()
    w = _wcs._Wcsprm(header)
    print w.spcfix()
    assert w.spcfix() == 0


def test_spec():
    w = _wcs._Wcsprm()
    assert w.spec == -1


@raises(AttributeError)
def test_spec_set():
    w = _wcs._Wcsprm()
    w.spec = 0


def test_specsys():
    w = _wcs._Wcsprm()
    assert w.specsys == b('')
    w.specsys = b('foo')
    assert w.specsys == b('foo')


def test_sptr():
    #TODO: Write me
    pass


def test_ssysobs():
    w = _wcs._Wcsprm()
    assert w.ssysobs == b('')
    w.ssysobs = b('foo')
    assert w.ssysobs == b('foo')


def test_ssyssrc():
    w = _wcs._Wcsprm()
    assert w.ssyssrc == b('')
    w.ssyssrc = b('foo')
    assert w.ssyssrc == b('foo')


def test_tab():
    w = _wcs._Wcsprm()
    assert len(w.tab) == 0
    # TODO: Inject some headers that have tables and test


def test_theta0():
    w = _wcs._Wcsprm()
    assert np.isnan(w.theta0)
    w.theta0 = 42.0
    assert w.theta0 == 42.0
    del w.theta0
    assert np.isnan(w.theta0)


def test_toheader():
    w = _wcs._Wcsprm()
    if sys.version_info[0] >= 3:
        assert isinstance(w.to_header(), bytes)
    else:
        assert isinstance(w.to_header(), str)


def test_velangl():
    w = _wcs._Wcsprm()
    assert w.velangl == 0.0
    w.velangl = 42.0
    assert w.velangl == 42.0
    del w.velangl
    assert np.isnan(w.velangl)


def test_velosys():
    w = _wcs._Wcsprm()
    assert np.isnan(w.velosys)
    w.velosys = 42.0
    assert w.velosys == 42.0
    del w.velosys
    assert np.isnan(w.velosys)


def test_zsource():
    w = _wcs._Wcsprm()
    assert np.isnan(w.zsource)
    w.zsource = 42.0
    assert w.zsource == 42.0
    del w.zsource
    assert np.isnan(w.zsource)


def test_cd_3d():
    header = open(os.path.join(ROOT_DIR, 'data', '3d_cd.hdr'), 'rb').read()
    w = _wcs._Wcsprm(header)
    assert w.cd.shape == (3, 3)
    assert w.get_pc().shape == (3, 3)
    assert w.get_cdelt().shape == (3,)


@raises(RuntimeError)
def test_get_pc():
    header = open(os.path.join(ROOT_DIR, 'data', '3d_cd.hdr'), 'rb').read()
    w = _wcs._Wcsprm(header)
    w.get_pc()[0, 0] = 42


@raises(_wcs.SingularMatrixError)
def test_detailed_err():
    w = _wcs._Wcsprm()
    w.pc = [[0, 0], [0, 0]]
    w.set()
