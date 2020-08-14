import numpy as np


def test_astrom_manager():
    # I was having a pretty hard time in coming
    # up with a unit test only testing the get_astrom function
    # that would not just test its implementation with its implementation
    # so we test a coordinate conversion using it

    from astropy.coordinates.astrom_manager import astrom_interpolation_resolution
    from astropy.coordinates import SkyCoord, EarthLocation, AltAz
    from astropy.time import Time
    import astropy.units as u

    location = EarthLocation(
        lon=-17.891105 * u.deg,
        lat=28.761584 * u.deg,
        height=2200 * u.m,
    )
    obstime = Time('2020-01-01T18:00') + np.linspace(0, 12, 1000) * u.hour

    altaz = AltAz(location=location, obstime=obstime)
    coord = SkyCoord(ra=83.63308333, dec=22.0145, unit=u.deg)

    # do the reference transformation, no interpolation
    ref = coord.transform_to(altaz)

    with astrom_interpolation_resolution.set(300 * u.s):
        interp_300s = coord.transform_to(altaz)

    # make sure they are actually different
    assert np.any(ref.separation(interp_300s) > u.Quantity(0.01, u.microarcsecond))

    # make sure the resolution is as good as we expect
    assert np.all(ref.separation(interp_300s) < u.Quantity(1, u.microarcsecond))


def test_interpolation_2d():
    # I was having a pretty hard time in coming
    # up with a unit test only testing the get_astrom function
    # that would not just test its implementation with its implementation
    # so we test a coordinate conversion using it

    from astropy.coordinates.astrom_manager import get_astrom

    from astropy.coordinates import EarthLocation, AltAz, GCRS
    from astropy.time import Time
    import astropy.units as u

    location = EarthLocation(
        lon=-17.891105 * u.deg,
        lat=28.761584 * u.deg,
        height=2200 * u.m,
    )
    obstime = Time('2020-01-01T18:00') + np.linspace(0, 12, 1000) * u.hour
    obstime = obstime.reshape((25, 40))

    altaz = AltAz(location=location, obstime=obstime)
    gcrs = GCRS(obstime=obstime)

    for frame, tcode in zip([altaz, altaz, gcrs], ['apio13', 'apci', 'apcs']):
        without_interp = get_astrom(frame, tcode, interpolation_resolution=0 * u.s)

        with_interp = get_astrom(frame, tcode)
        assert without_interp.shape == with_interp.shape
