# Note that we test the main astropy.wcs.WCS class directly rather than testing
# the mix-in class on its own (since it's not functional without being used as
# a mix-in)

import warnings

import numpy as np
from numpy.testing import assert_equal, assert_allclose

from .... import units as u
from ....tests.helper import assert_quantity_allclose
from ....units import Quantity
from ....coordinates import ICRS, FK5, Galactic, SkyCoord
from ....io.fits import Header
from ....utils.data import get_pkg_data_filename
from ...wcs import WCS
from ..fitswcs import custom_ctype_to_ucd_mapping

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
    assert wcs.world_axis_units == ['']

    assert_equal(wcs.axis_correlation_matrix, True)

    assert wcs.world_axis_object_components == [('world', 0, 'value')]

    assert wcs.world_axis_object_classes['world'][0] is Quantity
    assert wcs.world_axis_object_classes['world'][1] == ()
    assert wcs.world_axis_object_classes['world'][2]['unit'] is u.one

    assert_allclose(wcs.pixel_to_world_values(29), 29)
    assert_allclose(wcs.array_index_to_world_values(29), 29)

    assert_allclose(wcs.world_to_pixel_values(29), 29)
    assert_equal(wcs.world_to_array_index_values(29), (29,))

    # High-level API

    coord = wcs.pixel_to_world(29)
    assert_quantity_allclose(coord, 29 * u.one)

    coord = 15 * u.one

    x = wcs.world_to_pixel(coord)
    assert_allclose(x, 15.)

    i = wcs.world_to_array_index(coord)
    assert_equal(i, (15,))


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

WCS_SIMPLE_CELESTIAL = WCS(Header.fromstring(HEADER_SIMPLE_CELESTIAL, sep='\n'))


def test_simple_celestial():

    wcs = WCS_SIMPLE_CELESTIAL

    # Low-level API

    assert wcs.pixel_n_dim == 2
    assert wcs.world_n_dim == 2
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == ['pos.eq.ra', 'pos.eq.dec']
    assert wcs.world_axis_units == ['deg', 'deg']

    assert_equal(wcs.axis_correlation_matrix, True)

    assert wcs.world_axis_object_components == [('celestial', 0, 'spherical.lon.degree'),
                                                  ('celestial', 1, 'spherical.lat.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert_allclose(wcs.pixel_to_world_values(29, 39), (10, 20))
    assert_allclose(wcs.array_index_to_world_values(39, 29), (10, 20))

    assert_allclose(wcs.world_to_pixel_values(10, 20), (29., 39.))
    assert_equal(wcs.world_to_array_index_values(10, 20), (39, 29))

    # High-level API

    coord = wcs.pixel_to_world(29, 39)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert coord.ra.deg == 10
    assert coord.dec.deg == 20

    coord = wcs.array_index_to_world(39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, ICRS)
    assert coord.ra.deg == 10
    assert coord.dec.deg == 20

    coord = SkyCoord(10, 20, unit='deg', frame='icrs')

    x, y = wcs.world_to_pixel(coord)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)

    i, j = wcs.world_to_array_index(coord)
    assert_equal(i, 39)
    assert_equal(j, 29)

    # Check that if the coordinates are passed in a different frame things still
    # work properly

    coord_galactic = coord.galactic

    x, y = wcs.world_to_pixel(coord_galactic)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)

    i, j = wcs.world_to_array_index(coord_galactic)
    assert_equal(i, 39)
    assert_equal(j, 29)

    # Check that we can actually index the array

    data = np.arange(3600).reshape((60, 60))

    coord = SkyCoord(10, 20, unit='deg', frame='icrs')
    index = wcs.world_to_array_index(coord)
    assert_equal(data[index], 2369)

    coord = SkyCoord([10, 12], [20, 22], unit='deg', frame='icrs')
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

WCS_SPECTRAL_CUBE = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE, sep='\n'))


def test_spectral_cube():

    # Spectral cube with a weird axis ordering

    wcs = WCS_SPECTRAL_CUBE

    # Low-level API

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape is None
    assert wcs.pixel_shape is None
    assert wcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert wcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(wcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('freq', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['freq'][0] is Quantity
    assert wcs.world_axis_object_classes['freq'][1] == ()
    assert wcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(wcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(wcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(wcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(wcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))

    # High-level API

    coord, spec = wcs.pixel_to_world(29, 39, 44)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert coord.l.deg == 25
    assert coord.b.deg == 10
    assert isinstance(spec, Quantity)
    assert spec.to_value(u.Hz) == 20

    coord, spec = wcs.array_index_to_world(44, 39, 29)
    assert isinstance(coord, SkyCoord)
    assert isinstance(coord.frame, Galactic)
    assert coord.l.deg == 25
    assert coord.b.deg == 10
    assert isinstance(spec, Quantity)
    assert spec.to_value(u.Hz) == 20

    coord = SkyCoord(25, 10, unit='deg', frame='galactic')
    spec = 20 * u.Hz

    x, y, z = wcs.world_to_pixel(coord, spec)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)
    assert_allclose(z, 44.)

    # Order of world coordinates shouldn't matter
    x, y, z = wcs.world_to_pixel(spec, coord)
    assert_allclose(x, 29.)
    assert_allclose(y, 39.)
    assert_allclose(z, 44.)

    i, j, k = wcs.world_to_array_index(coord, spec)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)

    # Order of world coordinates shouldn't matter
    i, j, k = wcs.world_to_array_index(spec, coord)
    assert_equal(i, 44)
    assert_equal(j, 39)
    assert_equal(k, 29)


HEADER_SPECTRAL_CUBE_NONALIGNED = HEADER_SPECTRAL_CUBE.strip() + '\n' + """
PC2_3 = -0.5
PC3_2 = +0.5
"""

WCS_SPECTRAL_CUBE_NONALIGNED = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE_NONALIGNED, sep='\n'))


def test_spectral_cube_nonaligned():

    # Make sure that correlation matrix gets adjusted if there are non-identity
    # CD matrix terms.

    wcs = WCS_SPECTRAL_CUBE_NONALIGNED

    assert_equal(wcs.axis_correlation_matrix, [[True, True, True], [False, True, True], [True, True, True]])


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
OBSGEO-B= -24.6157 / [deg] Tel geodetic latitute (=North)+
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

WCS_TIME_CUBE = WCS(Header.fromstring(HEADER_TIME_CUBE, sep='\n'))


def test_time_cube():

    # Spectral cube with a weird axis ordering

    wcs = WCS_TIME_CUBE

    assert wcs.pixel_n_dim == 3
    assert wcs.world_n_dim == 3
    assert wcs.array_shape == (11, 2048, 2048)
    assert wcs.pixel_shape == (2048, 2048, 11)
    assert wcs.world_axis_physical_types == ['pos.eq.dec', 'pos.eq.ra', 'time']
    assert wcs.world_axis_units == ['deg', 'deg', 's']

    assert_equal(wcs.axis_correlation_matrix, [[True, True, False], [True, True, False], [False, False, True]])

    assert wcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('celestial', 0, 'spherical.lon.degree'),
                                                  ('utc', 0, 'value')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert wcs.world_axis_object_classes['utc'][0] is Quantity
    assert wcs.world_axis_object_classes['utc'][1] == ()
    assert wcs.world_axis_object_classes['utc'][2] == {'unit': 's'}

    assert_allclose(wcs.pixel_to_world_values(-449.2, 2955.6, 0),
                    (14.8289418840003, 2.01824372640628, 2375.341))

    assert_allclose(wcs.array_index_to_world_values(0, 2955.6, -449.2),
                    (14.8289418840003, 2.01824372640628, 2375.341))

    assert_allclose(wcs.world_to_pixel_values(14.8289418840003, 2.01824372640628, 2375.341),
                    (-449.2, 2955.6, 0))
    assert_equal(wcs.world_to_array_index_values(14.8289418840003, 2.01824372640628, 2375.341),
                 (0, 2956, -449))

    # High-level API

    # Make sure that we get a FutureWarning about the time column
    with warnings.catch_warnings(record=True) as warning_entries:
        warnings.resetwarnings()
        coord, time = wcs.pixel_to_world(29, 39, 44)

    # Check that there's at least one warning of the right category/message
    assert len(warning_entries) > 0
    found_warning = False
    for w in warning_entries:
        msg = 'In future, times will be represented by the Time class'
        if w.category is FutureWarning and str(w.message).startswith(msg):
            found_warning = True
    assert found_warning

    assert isinstance(coord, SkyCoord)
    assert isinstance(time, Quantity)

###############################################################################
# Extra corner cases
###############################################################################


def test_unrecognized_unit():
    # TODO: Determine whether the following behavior is desirable
    wcs = WCS(naxis=1)
    wcs.wcs.cunit = ['bananas // sekonds']
    assert wcs.world_axis_units == ['bananas // sekonds']


def test_distortion_correlations():

    filename = get_pkg_data_filename('../../tests/data/sip.fits')
    w = WCS(filename)
    assert_equal(w.axis_correlation_matrix, True)

    # Changing PC to an identity matrix doesn't change anything since
    # distortions are still present.
    w.wcs.pc = [[1, 0], [0, 1]]
    assert_equal(w.axis_correlation_matrix, True)

    # Nor does changing the name of the axes to make them non-celestial
    w.wcs.ctype = ['X', 'Y']
    assert_equal(w.axis_correlation_matrix, True)

    # However once we turn off the distortions the matrix changes
    w.sip = None
    assert_equal(w.axis_correlation_matrix, [[True, False], [False, True]])

    # If we go back to celestial coordinates then the matrix is all True again
    w.wcs.ctype = ['RA---TAN', 'DEC--TAN']
    assert_equal(w.axis_correlation_matrix, True)

    # Or if we change to X/Y but have a non-identity PC
    w.wcs.pc = [[0.9, -0.1], [0.1, 0.9]]
    w.wcs.ctype = ['X', 'Y']
    assert_equal(w.axis_correlation_matrix, True)


def test_custom_ctype_to_ucd_mappings():

    wcs = WCS(naxis=1)
    wcs.wcs.ctype = ['SPAM']

    assert wcs.world_axis_physical_types == [None]

    # Check simple behavior

    with custom_ctype_to_ucd_mapping({'APPLE': 'food.fruit'}):
        assert wcs.world_axis_physical_types == [None]

    with custom_ctype_to_ucd_mapping({'APPLE': 'food.fruit', 'SPAM': 'food.spam'}):
        assert wcs.world_axis_physical_types == ['food.spam']

    # Check nesting

    with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
        with custom_ctype_to_ucd_mapping({'APPLE': 'food.fruit'}):
            assert wcs.world_axis_physical_types == ['food.spam']

    with custom_ctype_to_ucd_mapping({'APPLE': 'food.fruit'}):
        with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
            assert wcs.world_axis_physical_types == ['food.spam']

    # Check priority in nesting

    with custom_ctype_to_ucd_mapping({'SPAM': 'notfood'}):
        with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
            assert wcs.world_axis_physical_types == ['food.spam']

    with custom_ctype_to_ucd_mapping({'SPAM': 'food.spam'}):
        with custom_ctype_to_ucd_mapping({'SPAM': 'notfood'}):
            assert wcs.world_axis_physical_types == ['notfood']


def test_caching_components_and_classes():

    # Make sure that when we change the WCS object, the classes and components
    # are updated (we use a cache internally, so we need to make sure the cache
    # is invalidated if needed)

    wcs = WCS_SIMPLE_CELESTIAL

    assert wcs.world_axis_object_components == [('celestial', 0, 'spherical.lon.degree'),
                                                ('celestial', 1, 'spherical.lat.degree')]

    assert wcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert wcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(wcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert wcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    wcs.wcs.radesys = 'FK5'

    frame = wcs.world_axis_object_classes['celestial'][2]['frame']
    assert isinstance(frame, FK5)
    assert frame.equinox.jyear == 2000.

    wcs.wcs.equinox = 2010

    frame = wcs.world_axis_object_classes['celestial'][2]['frame']
    assert isinstance(frame, FK5)
    assert frame.equinox.jyear == 2010.
