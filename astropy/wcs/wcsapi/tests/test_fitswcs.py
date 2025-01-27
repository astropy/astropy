# Note that we test the main astropy.wcs.WCS class directly rather than testing
# the mix-in class on its own (since it's not functional without being used as
# a mix-in)

import warnings
from itertools import product

import numpy as np
import pytest
from numpy.testing import assert_allclose, assert_array_equal, assert_equal
from packaging.version import Version

from astropy import units as u
from astropy.coordinates import (
    FK5,
    ICRS,
    ITRS,
    EarthLocation,
    Galactic,
    SkyCoord,
    SpectralCoord,
    StokesCoord,
)
from astropy.io import fits
from astropy.io.fits import Header
from astropy.io.fits.verify import VerifyWarning
from astropy.tests.helper import assert_quantity_allclose
from astropy.time import Time
from astropy.units import Quantity, UnitsWarning
from astropy.utils import iers
from astropy.utils.data import get_pkg_data_filename
from astropy.utils.exceptions import AstropyUserWarning
from astropy.wcs._wcs import __version__ as wcsver
from astropy.wcs.wcs import WCS, FITSFixedWarning, NoConvergence, Sip
from astropy.wcs.wcsapi.fitswcs import VELOCITY_FRAMES, custom_ctype_to_ucd_mapping

###############################################################################
# The following example is the simplest WCS with default values
###############################################################################


WCS_EMPTY = WCS(naxis=1)
WCS_EMPTY.wcs.crpix = [1]


def test_empty():
    wcs = WCS_EMPTY

    # Low-level API

    assert wcs.pixel_n_dim == 1
    assert wcs.world_n_dim == 1
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == [None]
    assert wcs.world_axis_units == [""]
    assert wcs.pixel_axis_names == [""]
    assert wcs.world_axis_names == [""]

    assert_equal(wcs.axis_correlation_matrix, True)

    assert wcs.world_axis_object_components == [("world", 0, "value")]

    assert wcs.world_axis_object_classes["world"][0] is Quantity
    assert wcs.world_axis_object_classes["world"][1] == ()
    assert wcs.world_axis_object_classes["world"][2]["unit"] is u.one

    assert_allclose(wcs.pixel_to_world_values(29), 29)
    assert_allclose(wcs.array_index_to_world_values(29), 29)

    assert np.ndim(wcs.pixel_to_world_values(29)) == 0
    assert np.ndim(wcs.array_index_to_world_values(29)) == 0

    assert_allclose(wcs.world_to_pixel_values(29), 29)
    assert_equal(wcs.world_to_array_index_values(29), (29,))

    assert np.ndim(wcs.world_to_pixel_values(29)) == 0
    assert np.ndim(wcs.world_to_array_index_values(29)) == 0

    # High-level API

    coord = wcs.pixel_to_world(29)
    assert_quantity_allclose(coord, 29 * u.one)
    assert np.ndim(coord) == 0

    coord = wcs.array_index_to_world(29)
    assert_quantity_allclose(coord, 29 * u.one)
    assert np.ndim(coord) == 0

    coord = 15 * u.one

    x = wcs.world_to_pixel(coord)
    assert_allclose(x, 15.0)
    assert np.ndim(x) == 0

    i = wcs.world_to_array_index(coord)
    assert_equal(i, 15)
    assert np.ndim(i) == 0


###############################################################################
# The following example is a simple 2D image with celestial coordinates
###############################################################################

HEADER_SIMPLE_CELESTIAL = """
WCSAXES = 2
CTYPE1  = RA---TAN
CTYPE2  = DEC--TAN
CRVAL1  = 10
CRVAL2  = 20
CRPIX1  = 30
CRPIX2  = 40
CDELT1  = -0.1
CDELT2  =  0.1
CROTA2  = 0.
CUNIT1  = deg
CUNIT2  = deg
"""

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_SIMPLE_CELESTIAL = WCS(Header.fromstring(HEADER_SIMPLE_CELESTIAL, sep="\n"))


def test_simple_celestial():
    wcs = WCS_SIMPLE_CELESTIAL

    # Low-level API

    assert wcs.pixel_n_dim == 2
    assert wcs.world_n_dim == 2
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == ["pos.eq.ra", "pos.eq.dec"]
    assert wcs.world_axis_units == ["deg", "deg"]
    assert wcs.pixel_axis_names == ["", ""]
    assert wcs.world_axis_names == ["", ""]

    assert_equal(wcs.axis_correlation_matrix, True)

    assert wcs.world_axis_object_components == [
        ("celestial", 0, "spherical.lon.degree"),
        ("celestial", 1, "spherical.lat.degree"),
    ]

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], ICRS)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] == (u.deg, u.deg)

    assert_allclose(wcs.pixel_to_world_values(29, 39), (10, 20))
    assert_allclose(wcs.array_index_to_world_values(39, 29), (10, 20))

    assert_allclose(wcs.world_to_pixel_values(10, 20), (29.0, 39.0))
    assert_equal(wcs.world_to_array_index_values(10, 20), (39, 29))

    # High-level API

    coord = wcs.pixel_to_world(29, 39)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert_allclose(coord.ra.deg, 10)
    assert_allclose(coord.dec.deg, 20)

    coord = wcs.array_index_to_world(39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert_allclose(coord.ra.deg, 10)
    assert_allclose(coord.dec.deg, 20)

    coord = SkyCoord(10, 20, unit="deg", frame="icrs")

    x, y = wcs.world_to_pixel(coord)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)

    i, j = wcs.world_to_array_index(coord)
    assert_equal(i, 39)
    assert_equal(j, 29)

    # Check that if the coordinates are passed in a different frame things still
    # work properly

    coord_galactic = coord.galactic

    x, y = wcs.world_to_pixel(coord_galactic)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)

    i, j = wcs.world_to_array_index(coord_galactic)
    assert_equal(i, 39)
    assert_equal(j, 29)

    # Check that we can actually index the array

    data = np.arange(3600).reshape((60, 60))

    coord = SkyCoord(10, 20, unit="deg", frame="icrs")
    index = wcs.world_to_array_index(coord)
    assert_equal(data[index], 2369)

    coord = SkyCoord([10, 12], [20, 22], unit="deg", frame="icrs")
    index = wcs.world_to_array_index(coord)
    assert_equal(data[index], [2369, 3550])


###############################################################################
# The following example is a spectral cube with axes in an unusual order
###############################################################################

HEADER_SPECTRAL_CUBE = """
WCSAXES = 3
CTYPE1  = GLAT-CAR
CTYPE2  = FREQ
CTYPE3  = GLON-CAR
CNAME1  = Latitude
CNAME2  = Frequency
CNAME3  = Longitude
CRVAL1  = 10
CRVAL2  = 20
CRVAL3  = 25
CRPIX1  = 30
CRPIX2  = 40
CRPIX3  = 45
CDELT1  = -0.1
CDELT2  =  0.5
CDELT3  =  0.1
CUNIT1  = deg
CUNIT2  = Hz
CUNIT3  = deg
"""

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_SPECTRAL_CUBE = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep="\n"))


def test_spectral_cube():
    # Spectral cube with a weird axis ordering

    wcs = WCS_SPECTRAL_CUBE

    # Low-level API

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, False, True], [False, True, False], [True, False, True]],
    )

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] == (u.deg, u.deg)

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29.0, 39.0, 44.0))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    # High-level API

    coord, spec = wcs.pixel_to_world(29, 39, 44)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert_allclose(coord.l.deg, 25)
    assert_allclose(coord.b.deg, 10)
    assert isinstance(spec, SpectralCoord)
    assert_allclose(spec.to_value(u.Hz), 20)

    coord, spec = wcs.array_index_to_world(44, 39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert_allclose(coord.l.deg, 25)
    assert_allclose(coord.b.deg, 10)
    assert isinstance(spec, SpectralCoord)
    assert_allclose(spec.to_value(u.Hz), 20)

    coord = SkyCoord(25, 10, unit="deg", frame="galactic")
    spec = 20 * u.Hz

    with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
        x, y, z = wcs.world_to_pixel(coord, spec)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)
    assert_allclose(z, 44.0)

    # Order of world coordinates shouldn't matter
    with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
        x, y, z = wcs.world_to_pixel(spec, coord)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)
    assert_allclose(z, 44.0)

    with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
        i, j, k = wcs.world_to_array_index(coord, spec)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)

    # Order of world coordinates shouldn't matter
    with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
        i, j, k = wcs.world_to_array_index(spec, coord)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)


HEADER_SPECTRAL_CUBE_NONALIGNED = (
    HEADER_SPECTRAL_CUBE.strip()
    + "\n"
    + """
PC2_3 = -0.5
PC3_2 = +0.5
"""
)

with warnings.catch_warnings():
    warnings.simplefilter("ignore", VerifyWarning)
    WCS_SPECTRAL_CUBE_NONALIGNED = WCS(
        Header.fromstring(HEADER_SPECTRAL_CUBE_NONALIGNED, sep="\n")
    )


def test_spectral_cube_nonaligned():
    # Make sure that correlation matrix gets adjusted if there are non-identity
    # CD matrix terms.

    wcs = WCS_SPECTRAL_CUBE_NONALIGNED

    assert wcs.world_axis_physical_types == [
        "pos.galactic.lat",
        "em.freq",
        "pos.galactic.lon",
    ]
    assert wcs.world_axis_units == ["deg", "Hz", "deg"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["Latitude", "Frequency", "Longitude"]

    assert_equal(
        wcs.axis_correlation_matrix,
        [
            [True, True, True],
            [False, True, True],
            [True, True, True],
        ],
    )

    # NOTE: we check world_axis_object_components and world_axis_object_classes
    # again here because in the past this failed when non-aligned axes were
    # present, so this serves as a regression test.

    assert len(wcs.world_axis_object_components) == 3
    assert wcs.world_axis_object_components[0] == (
        "celestial",
        1,
        "spherical.lat.degree",
    )
    assert wcs.world_axis_object_components[1][:2] == ("spectral", 0)
    assert wcs.world_axis_object_components[2] == (
        "celestial",
        0,
        "spherical.lon.degree",
    )

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], Galactic)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] == (u.deg, u.deg)

    assert wcs.world_axis_object_classes["spectral"][0] is Quantity
    assert wcs.world_axis_object_classes["spectral"][1] == ()
    assert wcs.world_axis_object_classes["spectral"][2] == {}


###############################################################################
# The following example is from Rots et al (2015), Table 5. It represents a
# cube with two spatial dimensions and one time dimension
###############################################################################


HEADER_TIME_CUBE = """
SIMPLE  = T / Fits standard
BITPIX  = -32 / Bits per pixel
NAXIS   = 3 / Number of axes
NAXIS1  = 2048 / Axis length
NAXIS2  = 2048 / Axis length
NAXIS3  = 11 / Axis length
DATE    = '2008-10-28T14:39:06' / Date FITS file was generated
OBJECT  = '2008 TC3' / Name of the object observed
EXPTIME = 1.0011 / Integration time
MJD-OBS = 54746.02749237 / Obs start
DATE-OBS= '2008-10-07T00:39:35.3342' / Observing date
TELESCOP= 'VISTA' / ESO Telescope Name
INSTRUME= 'VIRCAM' / Instrument used.
TIMESYS = 'UTC' / From Observatory Time System
TREFPOS = 'TOPOCENT' / Topocentric
MJDREF  = 54746.0 / Time reference point in MJD
RADESYS = 'ICRS' / Not equinoctal
CTYPE2  = 'RA---ZPN' / Zenithal Polynomial Projection
CRVAL2  = 2.01824372640628 / RA at ref pixel
CUNIT2  = 'deg' / Angles are degrees always
CRPIX2  = 2956.6 / Pixel coordinate at ref point
CTYPE1  = 'DEC--ZPN' / Zenithal Polynomial Projection
CRVAL1  = 14.8289418840003 / Dec at ref pixel
CUNIT1  = 'deg' / Angles are degrees always
CRPIX1  = -448.2 / Pixel coordinate at ref point
CTYPE3  = 'UTC' / linear time (UTC)
CRVAL3  = 2375.341 / Relative time of first frame
CUNIT3  = 's' / Time unit
CRPIX3  = 1.0 / Pixel coordinate at ref point
CTYPE3A = 'TT' / alternative linear time (TT)
CRVAL3A = 2440.525 / Relative time of first frame
CUNIT3A = 's' / Time unit
CRPIX3A = 1.0 / Pixel coordinate at ref point
OBSGEO-B= -24.6157 / [deg] Tel geodetic latitude (=North)+
OBSGEO-L= -70.3976 / [deg] Tel geodetic longitude (=East)+
OBSGEO-H= 2530.0000 / [m] Tel height above reference ellipsoid
CRDER3  = 0.0819 / random error in timings from fit
CSYER3  = 0.0100 / absolute time error
PC1_1   = 0.999999971570892 / WCS transform matrix element
PC1_2   = 0.000238449608932 / WCS transform matrix element
PC2_1   = -0.000621542859395 / WCS transform matrix element
PC2_2   = 0.999999806842218 / WCS transform matrix element
CDELT1  = -9.48575432499806E-5 / Axis scale at reference point
CDELT2  = 9.48683176211164E-5 / Axis scale at reference point
CDELT3  = 13.3629 / Axis scale at reference point
PV1_1   = 1. / ZPN linear term
PV1_3   = 42. / ZPN cubic term
"""

with warnings.catch_warnings():
    warnings.simplefilter("ignore", (VerifyWarning, FITSFixedWarning))
    WCS_TIME_CUBE = WCS(Header.fromstring(HEADER_TIME_CUBE, sep="\n"))


def test_time_cube():
    # Spectral cube with a weird axis ordering

    wcs = WCS_TIME_CUBE

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (11, 2048, 2048)
    assert wcs.pixel_shape == (2048, 2048, 11)
    assert wcs.world_axis_physical_types == ["pos.eq.dec", "pos.eq.ra", "time"]
    assert wcs.world_axis_units == ["deg", "deg", "s"]
    assert wcs.pixel_axis_names == ["", "", ""]
    assert wcs.world_axis_names == ["", "", ""]

    assert_equal(
        wcs.axis_correlation_matrix,
        [[True, True, False], [True, True, False], [False, False, True]],
    )

    components = wcs.world_axis_object_components
    assert components[0] == ("celestial", 1, "spherical.lat.degree")
    assert components[1] == ("celestial", 0, "spherical.lon.degree")
    assert components[2][:2] == ("time", 0)
    assert callable(components[2][2])

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], ICRS)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] == (u.deg, u.deg)

    assert wcs.world_axis_object_classes["time"][0] is Time
    assert wcs.world_axis_object_classes["time"][1] == ()
    assert wcs.world_axis_object_classes["time"][2] == {}
    assert callable(wcs.world_axis_object_classes["time"][3])

    assert_allclose(
        wcs.pixel_to_world_values(-449.2, 2955.6, 0),
        (14.8289418840003, 2.01824372640628, 2375.341),
    )

    assert_allclose(
        wcs.array_index_to_world_values(0, 2955.6, -449.2),
        (14.8289418840003, 2.01824372640628, 2375.341),
    )

    assert_allclose(
        wcs.world_to_pixel_values(14.8289418840003, 2.01824372640628, 2375.341),
        (-449.2, 2955.6, 0),
    )
    assert_equal(
        wcs.world_to_array_index_values(14.8289418840003, 2.01824372640628, 2375.341),
        (0, 2956, -449),
    )

    # High-level API

    coord, time = wcs.pixel_to_world(29, 39, 44)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert_allclose(coord.ra.deg, 1.7323356692202325)
    assert_allclose(coord.dec.deg, 14.783516054817797)
    assert isinstance(time, Time)
    assert_allclose(time.mjd, 54746.03429755324)

    coord, time = wcs.array_index_to_world(44, 39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert_allclose(coord.ra.deg, 1.7323356692202325)
    assert_allclose(coord.dec.deg, 14.783516054817797)
    assert isinstance(time, Time)
    assert_allclose(time.mjd, 54746.03429755324)

    x, y, z = wcs.world_to_pixel(coord, time)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)
    assert_allclose(z, 44.0)

    # Order of world coordinates shouldn't matter
    x, y, z = wcs.world_to_pixel(time, coord)
    assert_allclose(x, 29.0)
    assert_allclose(y, 39.0)
    assert_allclose(z, 44.0)

    i, j, k = wcs.world_to_array_index(coord, time)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)

    # Order of world coordinates shouldn't matter
    i, j, k = wcs.world_to_array_index(time, coord)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)


###############################################################################
# The following tests are to make sure that Time objects are constructed
# correctly for a variety of combinations of WCS keywords
###############################################################################


HEADER_TIME_1D = """
SIMPLE  = T
BITPIX  = -32
NAXIS   = 1
NAXIS1  = 2048
TIMESYS = 'UTC'
TREFPOS = 'TOPOCENT'
MJDREF  = 50002.6
CTYPE1  = 'UTC'
CRVAL1  = 5
CUNIT1  = 's'
CRPIX1  = 1.0
CDELT1  = 2
OBSGEO-L= -20
OBSGEO-B= -70
OBSGEO-H= 2530
"""

if Version(wcsver) >= Version("7.1"):
    HEADER_TIME_1D += "DATEREF = '1995-10-12T14:24:00'\n"


@pytest.fixture
def header_time_1d():
    return Header.fromstring(HEADER_TIME_1D, sep="\n")


def assert_time_at(header, position, jd1, jd2, scale, format):
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)
    time = wcs.pixel_to_world(position)
    assert_allclose(time.jd1, jd1, rtol=1e-10)
    assert_allclose(time.jd2, jd2, rtol=1e-10)
    assert time.format == format
    assert time.scale == scale


@pytest.mark.parametrize(
    "scale", ("tai", "tcb", "tcg", "tdb", "tt", "ut1", "utc", "local")
)
def test_time_1d_values(header_time_1d, scale):
    # Check that Time objects are instantiated with the correct values,
    # scales, and formats.

    header_time_1d["CTYPE1"] = scale.upper()
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, scale, "mjd")


def test_time_1d_values_gps(header_time_1d):
    # Special treatment for GPS scale
    header_time_1d["CTYPE1"] = "GPS"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + (7 + 19) / 3600 / 24, "tai", "mjd")


def test_time_1d_values_deprecated(header_time_1d):
    # Deprecated (in FITS) scales
    header_time_1d["CTYPE1"] = "TDT"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "tt", "mjd")
    header_time_1d["CTYPE1"] = "IAT"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "tai", "mjd")
    header_time_1d["CTYPE1"] = "GMT"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "utc", "mjd")
    header_time_1d["CTYPE1"] = "ET"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "tt", "mjd")


def test_time_1d_values_time(header_time_1d):
    header_time_1d["CTYPE1"] = "TIME"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "utc", "mjd")
    header_time_1d["TIMESYS"] = "TAI"
    assert_time_at(header_time_1d, 1, 2450003, 0.1 + 7 / 3600 / 24, "tai", "mjd")


@pytest.mark.remote_data
@pytest.mark.parametrize("scale", ("tai", "tcb", "tcg", "tdb", "tt", "ut1", "utc"))
def test_time_1d_roundtrip(header_time_1d, scale):
    # Check that coordinates round-trip

    pixel_in = np.arange(3, 10)

    header_time_1d["CTYPE1"] = scale.upper()

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header_time_1d)

    # Simple test
    time = wcs.pixel_to_world(pixel_in)
    pixel_out = wcs.world_to_pixel(time)
    assert_allclose(pixel_in, pixel_out)

    # Test with an intermediate change to a different scale/format
    time = wcs.pixel_to_world(pixel_in).tdb
    time.format = "isot"
    pixel_out = wcs.world_to_pixel(time)
    assert_allclose(pixel_in, pixel_out)


def test_time_1d_high_precision(header_time_1d):
    # Case where the MJDREF is split into two for high precision
    del header_time_1d["MJDREF"]
    header_time_1d["MJDREFI"] = 52000.0
    header_time_1d["MJDREFF"] = 1e-11

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header_time_1d)

    time = wcs.pixel_to_world(10)

    # Here we have to use a very small rtol to really test that MJDREFF is
    # taken into account
    assert_allclose(time.jd1, 2452001.0, rtol=1e-12)
    assert_allclose(time.jd2, -0.5 + 25 / 3600 / 24 + 1e-11, rtol=1e-13)


def test_time_1d_location_geodetic(header_time_1d):
    # Make sure that the location is correctly returned (geodetic case)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header_time_1d)

    time = wcs.pixel_to_world(10)

    lon, lat, alt = time.location.to_geodetic()

    # FIXME: alt won't work for now because ERFA doesn't implement the IAU 1976
    # ellipsoid (https://github.com/astropy/astropy/issues/9420)
    assert_allclose(lon.degree, -20)
    assert_allclose(lat.degree, -70)
    # assert_allclose(alt.to_value(u.m), 2530.)


@pytest.fixture
def header_time_1d_no_obs():
    header = Header.fromstring(HEADER_TIME_1D, sep="\n")
    del header["OBSGEO-L"]
    del header["OBSGEO-B"]
    del header["OBSGEO-H"]
    return header


def test_time_1d_location_geocentric(header_time_1d_no_obs):
    # Make sure that the location is correctly returned (geocentric case)

    header = header_time_1d_no_obs

    header["OBSGEO-X"] = 10
    header["OBSGEO-Y"] = -20
    header["OBSGEO-Z"] = 30

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)

    time = wcs.pixel_to_world(10)

    x, y, z = time.location.to_geocentric()

    assert_allclose(x.to_value(u.m), 10)
    assert_allclose(y.to_value(u.m), -20)
    assert_allclose(z.to_value(u.m), 30)


def test_time_1d_location_geocenter(header_time_1d_no_obs):
    header_time_1d_no_obs["TREFPOS"] = "GEOCENTER"

    wcs = WCS(header_time_1d_no_obs)
    time = wcs.pixel_to_world(10)

    x, y, z = time.location.to_geocentric()

    assert_allclose(x.to_value(u.m), 0)
    assert_allclose(y.to_value(u.m), 0)
    assert_allclose(z.to_value(u.m), 0)


def test_time_1d_location_missing(header_time_1d_no_obs):
    # Check what happens when no location is present

    wcs = WCS(header_time_1d_no_obs)
    with pytest.warns(
        UserWarning,
        match=(
            "Missing or incomplete observer location "
            "information, setting location in Time to None"
        ),
    ):
        time = wcs.pixel_to_world(10)

    assert time.location is None


def test_time_1d_location_incomplete(header_time_1d_no_obs):
    # Check what happens when location information is incomplete

    header_time_1d_no_obs["OBSGEO-L"] = 10.0

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header_time_1d_no_obs)

    with pytest.warns(
        UserWarning,
        match=(
            "Missing or incomplete observer location "
            "information, setting location in Time to None"
        ),
    ):
        time = wcs.pixel_to_world(10)

    assert time.location is None


def test_time_1d_location_unsupported(header_time_1d_no_obs):
    # Check what happens when TREFPOS is unsupported

    header_time_1d_no_obs["TREFPOS"] = "BARYCENTER"

    wcs = WCS(header_time_1d_no_obs)
    with pytest.warns(
        UserWarning,
        match=(
            "Observation location 'barycenter' is not "
            "supported, setting location in Time to None"
        ),
    ):
        time = wcs.pixel_to_world(10)

    assert time.location is None


def test_time_1d_unsupported_ctype(header_time_1d_no_obs):
    # For cases that we don't support yet, e.g. UT(...), use Time and drop sub-scale

    # Case where the MJDREF is split into two for high precision
    header_time_1d_no_obs["CTYPE1"] = "UT(WWV)"

    wcs = WCS(header_time_1d_no_obs)
    with (
        pytest.warns(
            UserWarning,
            match="Dropping unsupported sub-scale WWV from scale UT",
        ),
        pytest.warns(
            UserWarning,
            match="Missing or incomplete observer location information",
        ),
    ):
        time = wcs.pixel_to_world(10)

    assert isinstance(time, Time)


###############################################################################
# Extra corner cases
###############################################################################


def test_unrecognized_unit():
    # TODO: Determine whether the following behavior is desirable
    wcs = WCS(naxis=1)
    with pytest.warns(UnitsWarning):
        wcs.wcs.cunit = ["bananas // sekonds"]
        assert wcs.world_axis_units == ["bananas // sekonds"]


def test_distortion_correlations():
    filename = get_pkg_data_filename("../../tests/data/sip.fits")
    with pytest.warns(FITSFixedWarning):
        w = WCS(filename)
    assert_equal(w.axis_correlation_matrix, True)

    # Changing PC to an identity matrix doesn't change anything since
    # distortions are still present.
    w.wcs.pc = [[1, 0], [0, 1]]
    assert_equal(w.axis_correlation_matrix, True)

    # Nor does changing the name of the axes to make them non-celestial
    w.wcs.ctype = ["X", "Y"]
    assert_equal(w.axis_correlation_matrix, True)

    # However once we turn off the distortions the matrix changes
    w.sip = None
    assert_equal(w.axis_correlation_matrix, [[True, False], [False, True]])

    # If we go back to celestial coordinates then the matrix is all True again
    w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
    assert_equal(w.axis_correlation_matrix, True)

    # Or if we change to X/Y but have a non-identity PC
    w.wcs.pc = [[0.9, -0.1], [0.1, 0.9]]
    w.wcs.ctype = ["X", "Y"]
    assert_equal(w.axis_correlation_matrix, True)


def test_custom_ctype_to_ucd_mappings():
    wcs = WCS(naxis=1)
    wcs.wcs.ctype = ["SPAM"]

    assert wcs.world_axis_physical_types == [None]

    # Check simple behavior

    with custom_ctype_to_ucd_mapping({"APPLE": "food.fruit"}):
        assert wcs.world_axis_physical_types == [None]

    with custom_ctype_to_ucd_mapping({"APPLE": "food.fruit", "SPAM": "food.spam"}):
        assert wcs.world_axis_physical_types == ["food.spam"]

    # Check nesting

    with custom_ctype_to_ucd_mapping({"SPAM": "food.spam"}):
        with custom_ctype_to_ucd_mapping({"APPLE": "food.fruit"}):
            assert wcs.world_axis_physical_types == ["food.spam"]

    with custom_ctype_to_ucd_mapping({"APPLE": "food.fruit"}):
        with custom_ctype_to_ucd_mapping({"SPAM": "food.spam"}):
            assert wcs.world_axis_physical_types == ["food.spam"]

    # Check priority in nesting

    with custom_ctype_to_ucd_mapping({"SPAM": "notfood"}):
        with custom_ctype_to_ucd_mapping({"SPAM": "food.spam"}):
            assert wcs.world_axis_physical_types == ["food.spam"]

    with custom_ctype_to_ucd_mapping({"SPAM": "food.spam"}):
        with custom_ctype_to_ucd_mapping({"SPAM": "notfood"}):
            assert wcs.world_axis_physical_types == ["notfood"]


def test_caching_components_and_classes():
    # Make sure that when we change the WCS object, the classes and components
    # are updated (we use a cache internally, so we need to make sure the cache
    # is invalidated if needed)

    wcs = WCS_SIMPLE_CELESTIAL.deepcopy()

    assert wcs.world_axis_object_components == [
        ("celestial", 0, "spherical.lon.degree"),
        ("celestial", 1, "spherical.lat.degree"),
    ]

    assert wcs.world_axis_object_classes["celestial"][0] is SkyCoord
    assert wcs.world_axis_object_classes["celestial"][1] == ()
    assert isinstance(wcs.world_axis_object_classes["celestial"][2]["frame"], ICRS)
    assert wcs.world_axis_object_classes["celestial"][2]["unit"] == (u.deg, u.deg)

    wcs.wcs.radesys = "FK5"

    frame = wcs.world_axis_object_classes["celestial"][2]["frame"]
    assert isinstance(frame, FK5)
    assert frame.equinox.jyear == 2000.0

    wcs.wcs.equinox = 2010

    frame = wcs.world_axis_object_classes["celestial"][2]["frame"]
    assert isinstance(frame, FK5)
    assert frame.equinox.jyear == 2010.0


def test_sub_wcsapi_attributes():
    # Regression test for a bug that caused some of the WCS attributes to be
    # incorrect when using WCS.sub or WCS.celestial (which is an alias for sub
    # with lon/lat types).

    wcs = WCS_SPECTRAL_CUBE.deepcopy()
    wcs.pixel_shape = (30, 40, 50)
    wcs.pixel_bounds = [(-1, 11), (-2, 18), (5, 15)]

    # Use celestial shortcut

    wcs_sub1 = wcs.celestial

    assert wcs_sub1.pixel_n_dim == 2
    assert wcs_sub1.world_n_dim == 2
    assert wcs_sub1.array_shape == (50, 30)
    assert wcs_sub1.pixel_shape == (30, 50)
    assert wcs_sub1.pixel_bounds == [(-1, 11), (5, 15)]
    assert wcs_sub1.world_axis_physical_types == [
        "pos.galactic.lat",
        "pos.galactic.lon",
    ]
    assert wcs_sub1.world_axis_units == ["deg", "deg"]
    assert wcs_sub1.world_axis_names == ["Latitude", "Longitude"]

    # Try adding axes

    wcs_sub2 = wcs.sub([0, 2, 0])

    assert wcs_sub2.pixel_n_dim == 3
    assert wcs_sub2.world_n_dim == 3
    assert wcs_sub2.array_shape == (None, 40, None)
    assert wcs_sub2.pixel_shape == (None, 40, None)
    assert wcs_sub2.pixel_bounds == [None, (-2, 18), None]
    assert wcs_sub2.world_axis_physical_types == [None, "em.freq", None]
    assert wcs_sub2.world_axis_units == ["", "Hz", ""]
    assert wcs_sub2.world_axis_names == ["", "Frequency", ""]

    # Use strings

    wcs_sub3 = wcs.sub(["longitude", "latitude"])

    assert wcs_sub3.pixel_n_dim == 2
    assert wcs_sub3.world_n_dim == 2
    assert wcs_sub3.array_shape == (30, 50)
    assert wcs_sub3.pixel_shape == (50, 30)
    assert wcs_sub3.pixel_bounds == [(5, 15), (-1, 11)]
    assert wcs_sub3.world_axis_physical_types == [
        "pos.galactic.lon",
        "pos.galactic.lat",
    ]
    assert wcs_sub3.world_axis_units == ["deg", "deg"]
    assert wcs_sub3.world_axis_names == ["Longitude", "Latitude"]

    # Now try without CNAME set

    wcs.wcs.cname = [""] * wcs.wcs.naxis
    wcs_sub4 = wcs.sub(["longitude", "latitude"])

    assert wcs_sub4.pixel_n_dim == 2
    assert wcs_sub4.world_n_dim == 2
    assert wcs_sub4.array_shape == (30, 50)
    assert wcs_sub4.pixel_shape == (50, 30)
    assert wcs_sub4.pixel_bounds == [(5, 15), (-1, 11)]
    assert wcs_sub4.world_axis_physical_types == [
        "pos.galactic.lon",
        "pos.galactic.lat",
    ]
    assert wcs_sub4.world_axis_units == ["deg", "deg"]
    assert wcs_sub4.world_axis_names == ["", ""]


###############################################################################
# Spectral transformations
###############################################################################

HEADER_SPECTRAL_FRAMES = """
BUNIT   = 'Jy/beam'
EQUINOX =      2.000000000E+03
CTYPE1  = 'RA---SIN'
CRVAL1  =    2.60108333333E+02
CDELT1  =     -2.777777845E-04
CRPIX1  =                  1.0
CUNIT1  = 'deg'
CTYPE2  = 'DEC--SIN'
CRVAL2  =   -9.75000000000E-01
CDELT2  =      2.777777845E-04
CRPIX2  =                  1.0
CUNIT2  = 'deg'
CTYPE3  = 'FREQ'
CRVAL3  =    1.37835117405E+09
CDELT3  =      9.765625000E+04
CRPIX3  =                 32.0
CUNIT3  = 'Hz'
SPECSYS = 'TOPOCENT'
RESTFRQ =      1.420405752E+09 / [Hz]
RADESYS = 'FK5'
"""


@pytest.fixture
def header_spectral_frames():
    return Header.fromstring(HEADER_SPECTRAL_FRAMES, sep="\n")


def test_spectralcoord_frame(header_spectral_frames):
    # This is a test to check the numerical results of transformations between
    # different velocity frames. We simply make sure that the returned
    # SpectralCoords are in the right frame but don't check the transformations
    # since this is already done in test_spectralcoord_accuracy
    # in astropy.coordinates.

    with iers.conf.set_temp("auto_download", False):
        obstime = Time("2009-05-04T04:44:23", scale="utc")

        header = header_spectral_frames.copy()
        header["MJD-OBS"] = obstime.mjd
        header["CRVAL1"] = 16.33211
        header["CRVAL2"] = -34.2221
        header["OBSGEO-L"] = 144.2
        header["OBSGEO-B"] = -20.2
        header["OBSGEO-H"] = 0.0

        # We start off with a WCS defined in topocentric frequency
        with pytest.warns(FITSFixedWarning):
            wcs_topo = WCS(header)

        # We convert a single pixel coordinate to world coordinates and keep only
        # the second high level object - a SpectralCoord:
        sc_topo = wcs_topo.pixel_to_world(0, 0, 31)[1]

        # We check that this is in topocentric frame with zero velocities
        assert isinstance(sc_topo, SpectralCoord)
        assert isinstance(sc_topo.observer, ITRS)
        assert sc_topo.observer.obstime.isot == obstime.isot
        assert_equal(sc_topo.observer.data.differentials["s"].d_xyz.value, 0)

        observatory = (
            EarthLocation.from_geodetic(144.2, -20.2)
            .get_itrs(obstime=obstime)
            .transform_to(ICRS())
        )
        assert (
            observatory.separation_3d(sc_topo.observer.transform_to(ICRS())) < 1 * u.km
        )

        for specsys, expected_frame in VELOCITY_FRAMES.items():
            header["SPECSYS"] = specsys
            with pytest.warns(FITSFixedWarning):
                wcs = WCS(header)
            sc = wcs.pixel_to_world(0, 0, 31)[1]

            # Now transform to the expected velocity frame, which should leave
            # the spectral coordinate unchanged
            sc_check = sc.with_observer_stationary_relative_to(expected_frame)
            assert_quantity_allclose(sc.quantity, sc_check.quantity)


@pytest.mark.parametrize(
    ("ctype3", "observer"),
    product(["ZOPT", "BETA", "VELO", "VRAD", "VOPT"], [False, True]),
)
def test_different_ctypes(header_spectral_frames, ctype3, observer):
    header = header_spectral_frames.copy()
    header["CTYPE3"] = ctype3
    header["CRVAL3"] = 0.1
    header["CDELT3"] = 0.001

    if ctype3[0] == "V":
        header["CUNIT3"] = "m s-1"
    else:
        header["CUNIT3"] = ""

    header["RESTWAV"] = 1.420405752e09
    header["MJD-OBS"] = 55197

    if observer:
        header["OBSGEO-L"] = 144.2
        header["OBSGEO-B"] = -20.2
        header["OBSGEO-H"] = 0.0
        header["SPECSYS"] = "BARYCENT"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)

    skycoord, spectralcoord = wcs.pixel_to_world(0, 0, 31)

    assert isinstance(spectralcoord, SpectralCoord)

    if observer:
        pix = wcs.world_to_pixel(skycoord, spectralcoord)
    else:
        with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
            pix = wcs.world_to_pixel(skycoord, spectralcoord)

    assert_allclose(pix, [0, 0, 31], rtol=1e-6, atol=1e-9)


def test_non_convergence_warning():
    """Test case for issue #11446
    Since we can't define a target accuracy when plotting a WCS `all_world2pix`
    should not error but only warn when the default accuracy can't be reached.
    """
    # define a minimal WCS where convergence fails for certain image positions
    wcs = WCS(naxis=2)
    crpix = [0, 0]
    a = b = ap = bp = np.zeros((4, 4))
    a[3, 0] = -1.20116753e-07

    test_pos_x = [1000, 1]
    test_pos_y = [0, 2]

    wcs.sip = Sip(a, b, ap, bp, crpix)
    # first make sure the WCS works when using a low accuracy
    expected = wcs.all_world2pix(test_pos_x, test_pos_y, 0, tolerance=1e-3)

    # then check that it fails when using the default accuracy
    with pytest.raises(NoConvergence):
        wcs.all_world2pix(test_pos_x, test_pos_y, 0)

    # at last check that world_to_pixel_values raises a warning but returns
    # the same 'low accuray' result
    with pytest.warns(UserWarning):
        assert_allclose(wcs.world_to_pixel_values(test_pos_x, test_pos_y), expected)


HEADER_SPECTRAL_1D = """
CTYPE1  = 'FREQ'
CRVAL1  =    1.37835117405E+09
CDELT1  =      9.765625000E+04
CRPIX1  =                 32.0
CUNIT1  = 'Hz'
SPECSYS = 'TOPOCENT'
RESTFRQ =      1.420405752E+09 / [Hz]
RADESYS = 'FK5'
"""


@pytest.fixture
def header_spectral_1d():
    return Header.fromstring(HEADER_SPECTRAL_1D, sep="\n")


@pytest.mark.parametrize(
    ("ctype1", "observer"),
    product(["ZOPT", "BETA", "VELO", "VRAD", "VOPT"], [False, True]),
)
def test_spectral_1d(header_spectral_1d, ctype1, observer):
    # This is a regression test for issues that happened with 1-d WCS
    # where the target is not defined but observer is.

    header = header_spectral_1d.copy()
    header["CTYPE1"] = ctype1
    header["CRVAL1"] = 0.1
    header["CDELT1"] = 0.001

    if ctype1[0] == "V":
        header["CUNIT1"] = "m s-1"
    else:
        header["CUNIT1"] = ""

    header["RESTWAV"] = 1.420405752e09
    header["MJD-OBS"] = 55197

    if observer:
        header["OBSGEO-L"] = 144.2
        header["OBSGEO-B"] = -20.2
        header["OBSGEO-H"] = 0.0
        header["SPECSYS"] = "BARYCENT"

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", FITSFixedWarning)
        wcs = WCS(header)

    # First ensure that transformations round-trip

    spectralcoord = wcs.pixel_to_world(31)

    assert isinstance(spectralcoord, SpectralCoord)
    assert spectralcoord.target is None
    assert (spectralcoord.observer is not None) is observer

    if observer:
        expected_message = "No target defined on SpectralCoord"
    else:
        expected_message = "No observer defined on WCS"

    with pytest.warns(AstropyUserWarning, match=expected_message):
        pix = wcs.world_to_pixel(spectralcoord)

    assert_allclose(pix, [31], rtol=1e-6)

    # Also make sure that we can convert a SpectralCoord on which the observer
    # is not defined but the target is.

    with pytest.warns(AstropyUserWarning, match="No velocity defined on frame"):
        spectralcoord_no_obs = SpectralCoord(
            spectralcoord.quantity,
            doppler_rest=spectralcoord.doppler_rest,
            doppler_convention=spectralcoord.doppler_convention,
            target=ICRS(10 * u.deg, 20 * u.deg, distance=1 * u.kpc),
        )

    if observer:
        expected_message = "No observer defined on SpectralCoord"
    else:
        expected_message = "No observer defined on WCS"

    with pytest.warns(AstropyUserWarning, match=expected_message):
        pix2 = wcs.world_to_pixel(spectralcoord_no_obs)
    assert_allclose(pix2, [31], rtol=1e-6)

    # And finally check case when both observer and target are defined on the
    # SpectralCoord

    with pytest.warns(AstropyUserWarning, match="No velocity defined on frame"):
        spectralcoord_no_obs = SpectralCoord(
            spectralcoord.quantity,
            doppler_rest=spectralcoord.doppler_rest,
            doppler_convention=spectralcoord.doppler_convention,
            observer=ICRS(10 * u.deg, 20 * u.deg, distance=0 * u.kpc),
            target=ICRS(10 * u.deg, 20 * u.deg, distance=1 * u.kpc),
        )

    if observer:
        pix3 = wcs.world_to_pixel(spectralcoord_no_obs)
    else:
        with pytest.warns(AstropyUserWarning, match="No observer defined on WCS"):
            pix3 = wcs.world_to_pixel(spectralcoord_no_obs)

    assert_allclose(pix3, [31], rtol=1e-6)


HEADER_SPECTRAL_WITH_TIME = """
WCSAXES = 3
CTYPE1  = 'RA---TAN'
CTYPE2  = 'DEC--TAN'
CTYPE3  = 'WAVE'
CRVAL1  = 98.83153
CRVAL2  = -66.818
CRVAL3  = 6.4205
CRPIX1  = 21.
CRPIX2  = 22.
CRPIX3  = 1.
CDELT1  = 3.6111E-05
CDELT2  = 3.6111E-05
CDELT3  = 0.001
CUNIT1  = 'deg'
CUNIT2  = 'deg'
CUNIT3  = 'um'
MJD-AVG = 59045.41466
RADESYS = 'ICRS'
SPECSYS = 'BARYCENT'
TIMESYS = 'UTC'
"""


@pytest.fixture
def header_spectral_with_time():
    return Header.fromstring(HEADER_SPECTRAL_WITH_TIME, sep="\n")


def test_spectral_with_time_kw(header_spectral_with_time):
    def check_wcs(header):
        assert_allclose(w.all_pix2world(*w.wcs.crpix, 1), w.wcs.crval)
        sky, spec = w.pixel_to_world(*w.wcs.crpix)
        assert_allclose(
            (sky.spherical.lon.degree, sky.spherical.lat.degree, spec.value),
            w.wcs.crval,
            rtol=1e-3,
        )

    # Check with MJD-AVG and TIMESYS
    hdr = header_spectral_with_time.copy()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (VerifyWarning, FITSFixedWarning))
        w = WCS(hdr)
        # Make sure the correct keyword is used in a test
        assert ~np.isnan(w.wcs.mjdavg)
        assert np.isnan(w.wcs.mjdobs)

    check_wcs(w)

    # Check fall back to MJD-OBS
    hdr["MJD-OBS"] = hdr["MJD-AVG"]
    del hdr["MJD-AVG"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (VerifyWarning, FITSFixedWarning))
        w = WCS(hdr)
        # Make sure the correct keyword is used in a test
        assert ~np.isnan(w.wcs.mjdobs)
        assert np.isnan(w.wcs.mjdavg)
    check_wcs(w)

    # Check fall back to DATE--OBS
    hdr["DATE-OBS"] = "2020-07-15"
    del hdr["MJD-OBS"]
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", (VerifyWarning, FITSFixedWarning))
        w = WCS(hdr)
        w.wcs.mjdobs = np.nan
        # Make sure the correct keyword is used in a test
        assert np.isnan(w.wcs.mjdobs)
        assert np.isnan(w.wcs.mjdavg)
        assert w.wcs.dateobs != ""
    check_wcs(hdr)

    # Check fall back to scale='utc'
    del hdr["TIMESYS"]
    check_wcs(hdr)


def test_fits_tab_time_and_units():
    """
    This test is a regression test for https://github.com/astropy/astropy/issues/12095

    It checks the following:
      - If a spatial WCS isn't converted to units of deg by wcslib it still works.
      - If TIMESYS is upper case we parse it correctly
      - If a TIME CTYPE key has an algorithm code (in this case -TAB) it still works.

    The file used here was generated by gWCS and then edited to add the TIMESYS key.
    """
    with (
        fits.open(get_pkg_data_filename("data/example_4d_tab.fits")) as hdul,
        pytest.warns(FITSFixedWarning),
    ):
        w = WCS(header=hdul[0].header, fobj=hdul)

    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", message=r".*dubious year \(Note \d\)")
        world = w.pixel_to_world(0, 0, 0, 0)

    assert isinstance(world[0], SkyCoord)
    assert world[0].data.lat.unit is u.arcsec
    assert world[0].data.lon.unit is u.arcsec
    assert u.allclose(world[0].l, 0.06475506 * u.deg)
    assert u.allclose(world[0].b, -0.02430561 * u.deg)
    assert isinstance(world[1], SpectralCoord)
    assert u.allclose(world[1], 24.96 * u.Hz)
    assert isinstance(world[2], Time)
    assert world[2].scale == "utc"
    assert u.allclose(world[2].mjd, 0.00032986111111110716)


################################################################################
# Tests with Stokes
################################################################################


HEADER_POLARIZED = """
CTYPE1  = 'HPLT-TAN'
CTYPE2  = 'HPLN-TAN'
CTYPE3  = 'STOKES'
"""


@pytest.fixture
def header_polarized():
    return Header.fromstring(HEADER_POLARIZED, sep="\n")


@pytest.fixture
def wcs_polarized(header_polarized):
    return WCS(header_polarized)


def test_phys_type_polarization(wcs_polarized):
    w = wcs_polarized
    assert w.world_axis_physical_types[2] == "phys.polarization.stokes"


def test_pixel_to_world_stokes(wcs_polarized):
    w = wcs_polarized
    world = w.pixel_to_world(0, 0, 0)
    assert world[2] == 1
    assert isinstance(world[2], StokesCoord)
    assert_equal(world[2].symbol, "I")

    world = w.pixel_to_world(0, 0, [0, 1, 2, 3])
    assert isinstance(world[2], StokesCoord)
    assert_array_equal(world[2], [1, 2, 3, 4])
    assert_array_equal(world[2].symbol, ["I", "Q", "U", "V"])


@pytest.mark.parametrize("direction", ("world_to_pixel", "pixel_to_world"))
def test_out_of_bounds(direction):
    # Make sure that we correctly deal with any out-of-bound values in the
    # low-level API.

    wcs = WCS(naxis=2)
    wcs.wcs.crpix = (1, 1)
    wcs.wcs.set()

    func = (
        wcs.world_to_pixel_values
        if direction == "world_to_pixel"
        else wcs.pixel_to_world_values
    )

    xp = np.arange(5) + 1
    yp = np.arange(5) + 1

    # Before setting bounds

    # Python Scalars
    xw, yw = func(1, 1)
    assert_array_equal(xw, 1)
    assert_array_equal(yw, 1)

    # Numpy Scalars
    xw, yw = func(xp[0], yp[0])
    assert_array_equal(xw, 1)
    assert_array_equal(yw, 1)

    # Arrays
    xw, yw = func(xp, yp)
    assert_array_equal(xw, [1, 2, 3, 4, 5])
    assert_array_equal(yw, [1, 2, 3, 4, 5])

    # Mixed
    xw, yw = func(xp[0], yp)
    assert_array_equal(xw, 1)
    assert_array_equal(yw, [1, 2, 3, 4, 5])

    # Setting bounds on one dimension

    wcs.pixel_bounds = [(-0.5, 3.5), None]

    # Python Scalars

    xw, yw = func(1, 1)
    assert_array_equal(xw, 1)
    assert_array_equal(yw, 1)

    xw, yw = func(5, 5)
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, 5)

    # Numpy Scalars

    xw, yw = func(xp[0], yp[0])
    assert_array_equal(xw, 1)
    assert_array_equal(yw, 1)

    xw, yw = func(xp[-1], yp[-1])
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, 5)

    # Arrays
    xw, yw = func(xp, yp)
    assert_array_equal(xw, [1, 2, 3, np.nan, np.nan])
    assert_array_equal(yw, [1, 2, 3, 4, 5])

    # Mixed
    xw, yw = func(xp[-1], yp)
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, [1, 2, 3, 4, 5])

    # Setting bounds on both dimensions

    wcs.pixel_bounds = [(-0.5, 3.5), (2.5, 5.5)]

    # Python Scalars

    xw, yw = func(1, 1)
    assert_array_equal(xw, 1)
    assert_array_equal(yw, np.nan)

    xw, yw = func(5, 5)
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, 5)

    # Numpy Scalars

    xw, yw = func(xp[0], yp[0])
    assert_array_equal(xw, 1)
    assert_array_equal(yw, np.nan)

    xw, yw = func(xp[-1], yp[-1])
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, 5)

    # Arrays
    xw, yw = func(xp, yp)
    assert_array_equal(xw, [1, 2, 3, np.nan, np.nan])
    assert_array_equal(yw, [np.nan, np.nan, 3, 4, 5])

    # Mixed
    xw, yw = func(xp[-1], yp)
    assert_array_equal(xw, np.nan)
    assert_array_equal(yw, [np.nan, np.nan, 3, 4, 5])
