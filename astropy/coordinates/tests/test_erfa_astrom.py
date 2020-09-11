import numpy as np
import pytest

import astropy.units as u
from astropy.time import Time
from astropy.utils.exceptions import AstropyWarning
from astropy.coordinates import EarthLocation, AltAz, GCRS, SkyCoord

from astropy.coordinates.erfa_astrom import (
    erfa_astrom, ErfaAstrom, ErfaAstromInterpolator
)


def test_science_state():

    assert erfa_astrom.get().__class__ is ErfaAstrom

    res = 300 * u.s
    with erfa_astrom.set(ErfaAstromInterpolator(res)):
        assert isinstance(erfa_astrom.get(), ErfaAstromInterpolator)
        erfa_astrom.get().mjd_resolution == res.to_value(u.day)

    # context manager should have switched it back
    assert erfa_astrom.get().__class__ is ErfaAstrom

    # must be a subclass of BaseErfaAstrom
    with pytest.raises(TypeError):
        erfa_astrom.set('foo')


def test_warnings():

    with pytest.warns(AstropyWarning):
        with erfa_astrom.set(ErfaAstromInterpolator(9 * u.us)):
            pass


def test_erfa_astrom():
    # I was having a pretty hard time in coming
    # up with a unit test only testing the astrom provider
    # that would not just test its implementation with its implementation
    # so we test a coordinate conversion using it

    location = EarthLocation(
        lon=-17.891105 * u.deg,
        lat=28.761584 * u.deg,
        height=2200 * u.m,
    )
    obstime = Time('2020-01-01T18:00') + np.linspace(0, 1, 100) * u.hour

    altaz = AltAz(location=location, obstime=obstime)
    coord = SkyCoord(ra=83.63308333, dec=22.0145, unit=u.deg)

    # do the reference transformation, no interpolation
    ref = coord.transform_to(altaz)

    with erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)):
        interp_300s = coord.transform_to(altaz)

    # make sure they are actually different
    assert np.any(ref.separation(interp_300s) > 0.01 * u.microarcsecond)

    # make sure the resolution is as good as we expect
    assert np.all(ref.separation(interp_300s) < 1 * u.microarcsecond)


def test_interpolation_nd():
    '''
    Test that the interpolation also works for nd-arrays
    '''

    fact = EarthLocation(
        lon=-17.891105 * u.deg,
        lat=28.761584 * u.deg,
        height=2200 * u.m,
    )

    interp_provider = ErfaAstromInterpolator(300 * u.s)
    provider = ErfaAstrom()

    for shape in [tuple(), (1, ), (10, ), (3, 2), (2, 10, 5), (4, 5, 3, 2)]:
        # create obstimes of the desired shapes
        delta_t = np.linspace(0, 12, np.prod(shape, dtype=int)) * u.hour
        obstime = (Time('2020-01-01T18:00') + delta_t).reshape(shape)

        altaz = AltAz(location=fact, obstime=obstime)
        gcrs = GCRS(obstime=obstime)

        for frame, tcode in zip([altaz, altaz, gcrs], ['apio13', 'apci', 'apcs']):
            without_interp = getattr(provider, tcode)(frame)
            assert without_interp.shape == shape

            with_interp = getattr(interp_provider, tcode)(frame)
            assert with_interp.shape == shape


def test_interpolation_broadcasting():
    from astropy.coordinates.tests.utils import randomly_sample_sphere
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import Time
    import astropy.units as u

    from astropy.coordinates.erfa_astrom import erfa_astrom, ErfaAstromInterpolator

    # 1000 random locations on the sky
    ra, dec, _ = randomly_sample_sphere(100)
    coord = SkyCoord(ra, dec)

    # 30 times over the space of 1 hours
    times = Time('2020-01-01T20:00') + np.linspace(-0.5, 0.5, 30) * u.hour

    lst1 = EarthLocation(
        lon=-17.891498 * u.deg,
        lat=28.761443 * u.deg,
        height=2200 * u.m,
    )

    # note the use of broadcasting so that 300 times are broadcast against 1000 positions
    aa_frame = AltAz(obstime=times[:, np.newaxis], location=lst1)
    aa_coord = coord.transform_to(aa_frame)

    with erfa_astrom.set(ErfaAstromInterpolator(300 * u.s)):
        aa_coord_interp = coord.transform_to(aa_frame)

    assert aa_coord.shape == aa_coord_interp.shape
    assert np.all(aa_coord.separation(aa_coord_interp) < 1 * u.microarcsecond)
