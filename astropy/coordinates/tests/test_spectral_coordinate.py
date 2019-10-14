from astropy.coordinates import SkyCoord, EarthLocation, ICRS, GCRS, Galactic, CartesianDifferential, SphericalRepresentation, RadialDifferential, get_body_barycentric_posvel
from astropy import time
import numpy as np
import astropy.units as u

from ..spectra.spectral_coordinate import SpectralCoord


def test_create_spectral_coord():
    keck = EarthLocation.of_site('keck')
    obstime = time.Time('2018-12-13 9:00')

    observer_gcrs = keck.get_gcrs(obstime)

    spectral_axis = np.linspace(500, 2500, 1000) * u.AA
    spec_coord = SpectralCoord(spectral_axis)

    assert isinstance(spec_coord, u.Quantity)
    assert len(spec_coord) == 1000

    spec_coord.rest = 1215 * u.AA


