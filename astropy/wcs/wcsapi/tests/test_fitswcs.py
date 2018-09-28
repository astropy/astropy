# Note that we test the main astropy.wcs.WCS class directly rather than testing
# the mix-in class on its own (since it's not functional without being used as
# a mix-in)

from numpy.testing import assert_equal, assert_allclose

from .... import units as u
from ....units import Quantity
from ....coordinates import ICRS, Galactic, SkyCoord
from ....io.fits import Header
from ...wcs import WCS

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

    llwcs = WCS_SIMPLE_CELESTIAL

    assert llwcs.pixel_n_dim == 2
    assert llwcs.world_n_dim == 2
    assert llwcs.world_axis_physical_types == ['pos.eq.ra', 'pos.eq.dec']
    assert llwcs.world_axis_units == ['deg', 'deg']

    assert_equal(llwcs.axis_correlation_matrix, True)

    assert llwcs.world_axis_object_components == [('celestial', 0, 'spherical.lon.degree'),
                                                  ('celestial', 1, 'spherical.lat.degree')]

    assert llwcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert llwcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(llwcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert llwcs.world_axis_object_classes['celestial'][2]['unit'], 'deg'

    assert_allclose(llwcs.pixel_to_world_values(29, 39), (10, 20))
    assert_allclose(llwcs.array_index_to_world_values(39, 29), (10, 20))

    assert_allclose(llwcs.world_to_pixel_values(10, 20), (29., 39.))
    assert_equal(llwcs.world_to_array_index_values(10, 20), (39, 29))


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

    llwcs = WCS_SPECTRAL_CUBE

    assert llwcs.pixel_n_dim == 3
    assert llwcs.world_n_dim == 3
    assert llwcs.world_axis_physical_types == ['pos.galactic.lat', 'em.freq', 'pos.galactic.lon']
    assert llwcs.world_axis_units == ['deg', 'Hz', 'deg']

    assert_equal(llwcs.axis_correlation_matrix, [[True, False, True], [False, True, False], [True, False, True]])

    assert llwcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('freq', 0, 'value'),
                                                  ('celestial', 0, 'spherical.lon.degree')]

    assert llwcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert llwcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(llwcs.world_axis_object_classes['celestial'][2]['frame'], Galactic)
    assert llwcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert llwcs.world_axis_object_classes['freq'][0] is Quantity
    assert llwcs.world_axis_object_classes['freq'][1] == ()
    assert llwcs.world_axis_object_classes['freq'][2] == {'unit': 'Hz'}

    assert_allclose(llwcs.pixel_to_world_values(29, 39, 44), (10, 20, 25))
    assert_allclose(llwcs.array_index_to_world_values(44, 39, 29), (10, 20, 25))

    assert_allclose(llwcs.world_to_pixel_values(10, 20, 25), (29., 39., 44.))
    assert_equal(llwcs.world_to_array_index_values(10, 20, 25), (44, 39, 29))


HEADER_SPECTRAL_CUBE_NONALIGNED = HEADER_SPECTRAL_CUBE.strip() + '\n' + """
PC2_3 = -0.5
PC3_2 = +0.5
"""

WCS_SPECTRAL_CUBE_NONALIGNED = WCS(Header.fromstring(HEADER_SPECTRAL_CUBE_NONALIGNED, sep='\n'))


def test_spectral_cube_nonaligned():

    # Make sure that correlation matrix gets adjusted if there are non-identity
    # CD matrix terms.

    llwcs = WCS_SPECTRAL_CUBE_NONALIGNED

    assert_equal(llwcs.axis_correlation_matrix, [[True, True, True], [False, True, True], [True, True, True]])


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

    llwcs = WCS_TIME_CUBE

    assert llwcs.pixel_n_dim == 3
    assert llwcs.world_n_dim == 3
    assert llwcs.world_axis_physical_types == ['pos.eq.dec', 'pos.eq.ra', 'time']
    assert llwcs.world_axis_units == ['deg', 'deg', 's']

    assert_equal(llwcs.axis_correlation_matrix, [[True, True, False], [True, True, False], [False, False, True]])

    assert llwcs.world_axis_object_components == [('celestial', 1, 'spherical.lat.degree'),
                                                  ('celestial', 0, 'spherical.lon.degree'),
                                                  ('utc', 0, 'value')]

    assert llwcs.world_axis_object_classes['celestial'][0] is SkyCoord
    assert llwcs.world_axis_object_classes['celestial'][1] == ()
    assert isinstance(llwcs.world_axis_object_classes['celestial'][2]['frame'], ICRS)
    assert llwcs.world_axis_object_classes['celestial'][2]['unit'] is u.deg

    assert llwcs.world_axis_object_classes['utc'][0] is Quantity
    assert llwcs.world_axis_object_classes['utc'][1] == ()
    assert llwcs.world_axis_object_classes['utc'][2] == {'unit': 's'}

    assert_allclose(llwcs.pixel_to_world_values(-449.2, 2955.6, 0),
                    (14.8289418840003, 2.01824372640628, 2375.341))

    assert_allclose(llwcs.array_index_to_world_values(0, 2955.6, -449.2),
                    (14.8289418840003, 2.01824372640628, 2375.341))

    assert_allclose(llwcs.world_to_pixel_values(14.8289418840003, 2.01824372640628, 2375.341),
                    (-449.2, 2955.6, 0))
    assert_equal(llwcs.world_to_array_index_values(14.8289418840003, 2.01824372640628, 2375.341),
                 (0, 2956, -449))
