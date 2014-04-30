# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

# import numpy as np
# from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import pytest
from ...coordinates import (ICRS, FK4, FK5, Galactic, SkyCoord, SphericalRepresentation)
from ...time import Time

RA = 1.0 * u.deg
DEC = 2.0 * u.deg
C_ICRS = ICRS(RA, DEC)
C_FK5 = C_ICRS.transform_to(FK5)
J2001 = Time('J2001', scale='utc')


def tst_transform_to():
    for frame in (FK5, FK5(equinox=Time('J1975.0')),
                  FK4, FK4(equinox=Time('J1975.0')),
                  Galactic):
        c_frame = C_ICRS.transform_to(frame)
        s_icrs = SkyCoord(RA, DEC, frame='icrs')
        s_frame = s_icrs.transform_to(frame)
        assert c_frame.ra == s_frame.ra
        assert c_frame.dec == s_frame.dec


def test_frame_init():
    """
    Different ways of providing the frame.
    """
    sc = SkyCoord(RA, DEC, frame='icrs')
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, frame=ICRS)
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, 'icrs')
    assert sc.frame == 'icrs'

    sc = SkyCoord(RA, DEC, ICRS)
    assert sc.frame == 'icrs'

    sc = SkyCoord('icrs', RA, DEC)
    assert sc.frame == 'icrs'

    sc = SkyCoord(ICRS, RA, DEC)
    assert sc.frame == 'icrs'

    sc = SkyCoord(sc)
    assert sc.frame == 'icrs'

    sc = SkyCoord(C_ICRS)
    assert sc.frame == 'icrs'

    SkyCoord(C_ICRS, frame='icrs')
    assert sc.frame == 'icrs'

    with pytest.raises(ValueError) as err:
        SkyCoord(C_ICRS, frame='galactic')
    assert 'Cannot override frame=' in str(err)


def test_attr_inheritance():
    """
    When initializing from an existing coord the preferred attrs like
    equinox should be inherited to the SkyCoord.  If there is a conflict
    then raise an exception.
    """
    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc)
    assert sc2.equinox == sc.equinox
    assert sc2.obstime == sc.obstime

    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    sc2 = SkyCoord(sc._coord)  # Doesn't have equinox there
    assert sc2.equinox != sc.equinox
    assert sc2.obstime != sc.obstime


def test_frame_attr_getattr():
    """
    When accessing frame attributes like equinox, the value should come
    from self._coord when that object has the relevant attribute, otherwise
    from self.
    """
    sc = SkyCoord('icrs', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == 'J1999'  # Just the raw value
    assert sc.obstime == 'J2001'

    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999', obstime='J2001')
    assert sc.equinox == Time('J1999')  # Coming from the self._coord object
    assert sc.obstime == Time('J2001')

    sc = SkyCoord('fk4', 1, 2, unit='deg', equinox='J1999')
    assert sc.equinox == Time('J1999')
    assert sc.obstime == Time('J1999')


def test_api():
    """
    Verify that the API take-2 examples run
    """

    # NOT YET
    if False:
        sc = SkyCoord(SphericalRepresentation(lon=8 * u.hour, lat=5 * u.deg, distance=1 * u.kpc),
                      frame='icrs')

    sc = SkyCoord(ra=8 * u.hour, dec=5 * u.deg, frame='icrs')
    sc = SkyCoord(l=120 * u.deg, b=5 * u.deg, frame='galactic')

    # High-level classes can also be initialized directly from low-level objects
    sc = SkyCoord(ICRS(ra=8 * u.hour, dec=5 * u.deg))

    # The next example raises an error because the high-level class must always
    # have position data.
    if False:
        with pytest.raises(ValueError):
            sc = SkyCoord(FK5(equinox=J2001))  # raises ValueError

    # similarly, the low-level object can always be accessed

    # NOT YET.  NEVER?
    # assert str(sframe) == '<ICRS RA=120.000 deg, Dec=5.00000 deg>'

    # Should (eventually) support a variety of possible complex string formats
    sc = SkyCoord('8h00m00s +5d00m00.0s', frame='icrs')

    # In the next example, the unit is only needed b/c units are ambiguous.  In
    # general, we *never* accept ambiguity
    sc = SkyCoord('8:00:00 +5:00:00.0', unit=(u.hour, u.deg), frame='icrs')

    # The next one would yield length-2 array coordinates, because of the comma

    # NOT YET
    # sc = SkyCoord(['8h 5d', '2°5\'12.3" 0.3rad'], frame='icrs')

    # It should also interpret common designation styles as a coordinate
    # NOT YET
    # sc = SkyCoord('SDSS J123456.89-012345.6', frame='icrs')

    # the string representation can be inherited from the low-level class.

    # NOT YET
    # assert str(sc) == '<SkyCoord (ICRS) RA=120.000 deg, Dec=5.00000 deg>'

    # but it should also be possible to provide formats for outputting to strings,
    # similar to `Time`.  This can be added right away or at a later date.

    # transformation is done the same as for low-level classes, which it delegates to

    # NOT YET
    # scfk5_j2001 = stransform_to(FK5(equinox=J2001))

    # The key difference is that the high-level class remembers frame information
    # necessary for round-tripping, unlike the low-level classes:
    sc1 = SkyCoord(ra=8 * u.hour, dec=5 * u.deg, equinox=J2001, frame='fk5')
    sc2 = sc1.transform_to('icrs')
    # The next assertion succeeds, but it doesn't mean anything for ICRS, as ICRS
    # isn't defined in terms of an equinox
    assert sc2.equinox == J2001
    # But it *is* necessary once we transform to FK5
    sc3 = sc2.transform_to('fk5')
    assert sc3.equinox == J2001
    assert sc1.ra == sc3.ra
    # Note that this did *not* work in the low-level class example shown above,
    # because the ICRS low-level class does not store `equinox`.

    # `SkyCoord` will also include the attribute-style access that is in the
    # v0.2/0.3 coordinate objects.  This will *not* be in the low-level classes
    sc = SkyCoord(ra=8 * u.hour, dec=5 * u.deg, frame='icrs')
    scgal = sc.galactic
    assert str(scgal).startswith('<Galactic SkyCoord: l=216.317')

    # the existing `from_name` and `match_to_catalog_*` methods will be moved to the
    # high-level class as convenience functionality.

    if False:
        m31icrs = SkyCoord.from_name('M31', frame='icrs')
        assert str(m31icrs) == '<SkyCoord (ICRS) RA=10.68471 deg, Dec=41.26875 deg>'

        cat1 = SkyCoord(ra=1 * u.hr, dec=2 * u.deg, distance=3 * u.kpc)
        cat2 = SkyCoord(ra=1 * u.hr, dec=2 * u.deg, distance=3 * u.kpc)
        idx2, sep2d, dist3d = cat1.match_to_catalog_sky(cat2)
        idx2, sep2d, dist3d = cat1.match_to_catalog_3d(cat2)

        # additional convenience functionality for the future should be added as methods
        # on `SkyCoord`, *not* the low-level classes.
