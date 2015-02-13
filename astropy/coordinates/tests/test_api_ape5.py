# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

"""
This is the APE5 coordinates API document re-written to work as a series of test
functions.

Note that new tests for coordinates functionality should generally *not* be
added to this file - instead, add them to other appropriate test modules  in
this package, like ``test_sky_coord.py``, ``test_frames.py``, or
``test_representation.py``.  This file is instead meant mainly to keep track of
deviations from the original APE5 plan.
"""

import numpy as np
from numpy.random import randn
from numpy import testing as npt
from ...tests.helper import pytest
raises = pytest.raises

from ...extern import six

from ... import units as u
from ... import time
from ... import coordinates as coords
from ..errors import *

try:
    import scipy
except ImportError:
    HAS_SCIPY = False
else:
    HAS_SCIPY = True


def test_representations_api():
    from ..representation import SphericalRepresentation, \
        UnitSphericalRepresentation, PhysicsSphericalRepresentation, \
        CartesianRepresentation
    from ... coordinates import Angle, Longitude, Latitude, Distance

    #<-----------------Classes for representation of coordinate data--------------->
    # These classes inherit from a common base class and internally contain Quantity
    # objects, which are arrays (although they may act as scalars, like numpy's
    # length-0  "arrays")

    # They can be initialized with a variety of ways that make intuitive sense.
    # Distance is optional.
    UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg)
    UnitSphericalRepresentation(lon=8*u.hourangle, lat=5*u.deg)
    SphericalRepresentation(lon=8*u.hourangle, lat=5*u.deg, distance=10*u.kpc)

    # In the initial implementation, the lat/lon/distance arguments to the
    # initializer must be in order. A *possible* future change will be to allow
    # smarter guessing of the order.  E.g. `Latitude` and `Longitude` objects can be
    # given in any order.
    UnitSphericalRepresentation(Longitude(8, u.hour), Latitude(5, u.deg))
    SphericalRepresentation(Longitude(8, u.hour), Latitude(5, u.deg), Distance(10, u.kpc))

    # Arrays of any of the inputs are fine
    UnitSphericalRepresentation(lon=[8, 9]*u.hourangle, lat=[5, 6]*u.deg)

    # Default is to copy arrays, but optionally, it can be a reference
    UnitSphericalRepresentation(lon=[8, 9]*u.hourangle, lat=[5, 6]*u.deg, copy=False)

    # strings are parsed by `Latitude` and `Longitude` constructors, so no need to
    # implement parsing in the Representation classes
    UnitSphericalRepresentation(lon=Angle('2h6m3.3s'), lat=Angle('0.1rad'))

    # Or, you can give `Quantity`s with keywords, and they will be internally
    # converted to Angle/Distance
    c1 = SphericalRepresentation(lon=8*u.hourangle, lat=5*u.deg, distance=10*u.kpc)

    # Can also give another representation object with the `reprobj` keyword.
    c2 = SphericalRepresentation.from_representation(c1)

    #  distance, lat, and lon typically will just match in shape
    SphericalRepresentation(lon=[8, 9]*u.hourangle, lat=[5, 6]*u.deg, distance=[10, 11]*u.kpc)
    # if the inputs are not the same, if possible they will be broadcast following
    # numpy's standard broadcasting rules.
    c2 = SphericalRepresentation(lon=[8, 9]*u.hourangle, lat=[5, 6]*u.deg, distance=10*u.kpc)
    assert len(c2.distance) == 2
    #when they can't be broadcast, it is a ValueError (same as Numpy)
    with raises(ValueError):
        c2 = UnitSphericalRepresentation(lon=[8, 9, 10]*u.hourangle, lat=[5, 6]*u.deg)

    # It's also possible to pass in scalar quantity lists with mixed units. These
    # are converted to array quantities following the same rule as `Quantity`: all
    # elements are converted to match the first element's units.
    c2 = UnitSphericalRepresentation(lon=Angle([8*u.hourangle, 135*u.deg]),
                                     lat=Angle([5*u.deg, (6*np.pi/180)*u.rad]))
    assert c2.lat.unit == u.deg and c2.lon.unit == u.hourangle
    npt.assert_almost_equal(c2.lon[1].value, 9)

    # The Quantity initializer itself can also be used to force the unit even if the
    # first element doesn't have the right unit
    lon = u.Quantity([120*u.deg, 135*u.deg], u.hourangle)
    lat = u.Quantity([(5*np.pi/180)*u.rad, 0.4*u.hourangle], u.deg)
    c2 = UnitSphericalRepresentation(lon, lat)

    # regardless of how input, the `lat` and `lon` come out as angle/distance
    assert isinstance(c1.lat, Angle)
    assert isinstance(c1.lat, Latitude)  # `Latitude` is an `Angle` subclass
    assert isinstance(c1.distance, Distance)

    # but they are read-only, as representations are immutable once created
    with raises(AttributeError):
        c1.lat = Latitude(5, u.deg)
    # Note that it is still possible to modify the array in-place, but this is not
    # sanctioned by the API, as this would prevent things like caching.
    c2.lat[:] = [0] * u.deg  # possible, but NOT SUPPORTED

    # To address the fact that there are various other conventions for how spherical
    # coordinates are defined, other conventions can be included as new classes.
    # Later there may be other conventions that we implement - for now just the
    # physics convention, as it is one of the most common cases.
    c3 = PhysicsSphericalRepresentation(phi=120*u.deg, theta=85*u.deg, r=3*u.kpc)

    # first dimension must be length-3 if a lone `Quantity` is passed in.
    c1 = CartesianRepresentation(randn(3, 100) * u.kpc)
    assert c1.xyz.shape[0] == 3
    assert c1.xyz.unit == u.kpc
    assert c1.x.shape[0] == 100
    assert c1.y.shape[0] == 100
    assert c1.z.shape[0] == 100
    # can also give each as separate keywords
    CartesianRepresentation(x=randn(100)*u.kpc, y=randn(100)*u.kpc, z=randn(100)*u.kpc)
    # if the units don't match but are all distances, they will automatically be
    # converted to match `x`
    xarr, yarr, zarr = randn(3, 100)
    c1 = CartesianRepresentation(x=xarr*u.kpc, y=yarr*u.kpc, z=zarr*u.kpc)
    c2 = CartesianRepresentation(x=xarr*u.kpc, y=yarr*u.kpc, z=zarr*u.pc)
    assert c1.xyz.unit ==  c2.xyz.unit == u.kpc
    npt.assert_allclose((c1.z / 1000) - c2.z, 0, atol=1e-10)

    # representations convert into other representations via  `represent_as`
    srep = SphericalRepresentation(lon=90*u.deg, lat=0*u.deg, distance=1*u.pc)
    crep = srep.represent_as(CartesianRepresentation)
    npt.assert_allclose(crep.x.value, 0, atol=1e-10)
    npt.assert_allclose(crep.y.value, 1, atol=1e-10)
    npt.assert_allclose(crep.z.value, 0, atol=1e-10)
    # The functions that actually do the conversion are defined via methods on the
    # representation classes. This may later be expanded into a full registerable
    # transform graph like the coordinate frames, but initially it will be a simpler
    # method system


def test_frame_api():
    from ..representation import SphericalRepresentation, \
                                 UnitSphericalRepresentation
    from ..builtin_frames import ICRS, FK5
    #<---------------------Reference Frame/"Low-level" classes--------------------->
    # The low-level classes have a dual role: they act as specifiers of coordinate
    # frames and they *may* also contain data as one of the representation objects,
    # in which case they are the actual coordinate objects themselves.

    # They can always accept a representation as a first argument
    icrs = ICRS(UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg))

    # which is stored as the `data` attribute
    assert icrs.data.lat == 5*u.deg
    assert icrs.data.lon == 8*u.hourangle

    # Frames that require additional information like equinoxs or obstimes get them
    # as keyword parameters to the frame constructor.  Where sensible, defaults are
    # used. E.g., FK5 is almost always J2000 equinox
    fk5 = FK5(UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg))
    J2000 = time.Time('J2000', scale='utc')
    fk5_2000 = FK5(UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg), equinox=J2000)
    assert fk5.equinox == fk5_2000.equinox

    # the information required to specify the frame is immutable
    J2001 = time.Time('J2001', scale='utc')
    with raises(AttributeError):
        fk5.equinox = J2001

    # Similar for the representation data.
    with raises(AttributeError):
        fk5.data = UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg)

    # There is also a class-level attribute that lists the attributes needed to
    # identify the frame.  These include attributes like `equinox` shown above.
    assert all([nm in ('equinox', 'obstime') for nm in fk5.get_frame_attr_names()])

    # the result of `get_frame_attr_names` is called for particularly in  the
    # high-level class (discussed below) to allow round-tripping between various
    # frames.  It is also part of the public API for other similar developer /
    # advanced users' use.

    # The actual position information is accessed via the representation objects
    npt.assert_allclose(icrs.represent_as(SphericalRepresentation).lat.to(u.deg), 5*u.deg)
    npt.assert_allclose(icrs.spherical.lat.to(u.deg), 5*u.deg)  # shorthand for the above
    assert icrs.cartesian.z.value > 0

    # Many frames have a "default" representation, the one in which they are
    # conventionally described, often with a special name for some of the
    # coordinates. E.g., most equatorial coordinate systems are spherical with RA and
    # Dec. This works simply as a shorthand for the longer form above

    npt.assert_allclose(icrs.dec.to(u.deg), 5*u.deg)
    npt.assert_allclose(fk5.ra.to(u.hourangle), 8*u.hourangle)

    assert icrs.representation == SphericalRepresentation

    # low-level classes can also be initialized with names valid for that representation
    # and frame:
    icrs_2 = ICRS(ra=8*u.hour, dec=5*u.deg, distance=1*u.kpc)
    npt.assert_allclose(icrs.ra.to(u.deg), icrs_2.ra.to(u.deg))

    # and these are taken as the default if keywords are not given:
    #icrs_nokwarg = ICRS(8*u.hour, 5*u.deg, distance=1*u.kpc)
    #assert icrs_nokwarg.ra == icrs_2.ra and icrs_nokwarg.dec == icrs_2.dec

    # they also are capable of computing on-sky or 3d separations from each other,
    # which will be a direct port of the existing methods:
    coo1 = ICRS(ra=0*u.hour, dec=0*u.deg)
    coo2 = ICRS(ra=0*u.hour, dec=1*u.deg)
    # `separation` is the on-sky separation
    assert coo1.separation(coo2).degree == 1.0

    # while `separation_3d` includes the 3D distance information
    coo3 = ICRS(ra=0*u.hour, dec=0*u.deg, distance=1*u.kpc)
    coo4 = ICRS(ra=0*u.hour, dec=0*u.deg, distance=2*u.kpc)
    assert coo3.separation_3d(coo4).kpc == 1.0

    # The next example fails because `coo1` and `coo2` don't have distances
    with raises(ValueError):
        assert coo1.separation_3d(coo2).kpc == 1.0

    # repr/str also shows info, with frame and data
    #assert repr(fk5) == ''


def test_transform_api():
    from ..representation import UnitSphericalRepresentation
    from ..builtin_frames import ICRS, FK5
    from ..baseframe import frame_transform_graph, BaseCoordinateFrame
    from ..transformations import DynamicMatrixTransform
    #<-------------------------Transformations------------------------------------->
    # Transformation functionality is the key to the whole scheme: they transform
    # low-level classes from one frame to another.

    #(used below but defined above in the API)
    fk5 = FK5(ra=8*u.hour, dec=5*u.deg)

    # If no data (or `None`) is given, the class acts as a specifier of a frame, but
    # without any stored data.
    J2001 = time.Time('J2001', scale='utc')
    fk5_J2001_frame = FK5(equinox=J2001)

    # if they do not have data, the string instead is the frame specification
    assert repr(fk5_J2001_frame) == "<FK5 Frame (equinox=J2001.000)>"

    #  Note that, although a frame object is immutable and can't have data added, it
    #  can be used to create a new object that does have data by giving the
    # `realize_frame` method a representation:
    srep = UnitSphericalRepresentation(lon=8*u.hour, lat=5*u.deg)
    fk5_j2001_with_data = fk5_J2001_frame.realize_frame(srep)
    assert fk5_j2001_with_data.data is not None
    # Now `fk5_j2001_with_data` is in the same frame as `fk5_J2001_frame`, but it
    # is an actual low-level coordinate, rather than a frame without data.

    # These frames are primarily useful for specifying what a coordinate should be
    # transformed *into*, as they are used by the `transform_to` method
    # E.g., this snippet precesses the point to the new equinox
    newfk5 = fk5.transform_to(fk5_J2001_frame)
    assert newfk5.equinox == J2001

    # classes can also be given to `transform_to`, which then uses the defaults for
    # the frame information:
    samefk5 = fk5.transform_to(FK5)
    # `fk5` was initialized using default `obstime` and `equinox`, so:
    npt.assert_allclose(samefk5.ra - fk5.ra, 0, atol=1e-10)
    npt.assert_allclose(samefk5.dec - fk5.dec, 0, atol=1e-10)

    # transforming to a new frame necessarily loses framespec information if that
    # information is not applicable to the new frame.  This means transforms are not
    # always round-trippable:
    fk5_2 = FK5(ra=8*u.hour, dec=5*u.deg, equinox=J2001)
    ic_trans = fk5_2.transform_to(ICRS)

    # `ic_trans` does not have an `equinox`, so now when we transform back to FK5,
    # it's a *different* RA and Dec
    fk5_trans = ic_trans.transform_to(FK5)
    assert not np.allclose(fk5_2.ra.to(u.deg), fk5_trans.ra.to(u.deg), rtol=0, atol=1e-10)

    # But if you explicitly give the right equinox, all is fine
    fk5_trans_2 = fk5_2.transform_to(FK5(equinox=J2001))
    npt.assert_allclose(fk5_2.ra.to(u.deg), fk5_trans_2.ra.to(u.deg), rtol=0, atol=1e-10)

    # Trying to tansforming a frame with no data is of course an error:
    with raises(ValueError):
        FK5(equinox=J2001).transform_to(ICRS)


    # To actually define a new transformation, the same scheme as in the
    # 0.2/0.3 coordinates framework can be re-used - a graph of transform functions
    # connecting various coordinate classes together.  The main changes are:
    # 1) The transform functions now get the frame object they are transforming the
    #    current data into.
    # 2) Frames with additional information need to have a way to transform between
    #    objects of the same class, but with different framespecinfo values

    # An example transform function:
    class SomeNewSystem(BaseCoordinateFrame):
        pass
    @frame_transform_graph.transform(DynamicMatrixTransform, SomeNewSystem, FK5)
    def new_to_fk5(newobj, fk5frame):
        ot = newobj.obstime
        eq = fk5frame.equinox
        # ... build a *cartesian* transform matrix using `eq` that transforms from
        # the `newobj` frame as observed at `ot` to FK5 an equinox `eq`
        matrix = np.eye(3)
        return matrix

    # Other options for transform functions include one that simply returns the new
    # coordinate object, and one that returns a cartesian matrix but does *not*
    # require `newobj` or `fk5frame` - this allows optimization of the transform.


def test_highlevel_api():
    J2001 = time.Time('J2001', scale='utc')

    #<---------------------------"High-level" class-------------------------------->
    # The "high-level" class is intended to wrap the lower-level classes in such a
    # way that they can be round-tripped, as well as providing a variety of
    # convenience functionality.  This document is not intended to show *all* of the
    # possible high-level functionality, rather how the high-level classes are
    # initialized and interact with the low-level classes

    # this creates an object that contains an `ICRS` low-level class, initialized
    # identically to the first ICRS example further up.

    sc = coords.SkyCoord(coords.SphericalRepresentation(lon=8 * u.hour,
                         lat=5 * u.deg, distance=1 * u.kpc), frame='icrs')

    # Other representations and `system` keywords delegate to the appropriate
    # low-level class. The already-existing registry for user-defined coordinates
    # will be used by `SkyCoordinate` to figure out what various the `system`
    # keyword actually means.

    sc = coords.SkyCoord(ra=8 * u.hour, dec=5 * u.deg, frame='icrs')
    sc = coords.SkyCoord(l=120 * u.deg, b=5 * u.deg, frame='galactic')

    # High-level classes can also be initialized directly from low-level objects
    sc = coords.SkyCoord(coords.ICRS(ra=8 * u.hour, dec=5 * u.deg))

    # The next example raises an error because the high-level class must always
    # have position data.
    with pytest.raises(ValueError):
        sc = coords.SkyCoord(coords.FK5(equinox=J2001))  # raises ValueError

    # similarly, the low-level object can always be accessed

    #this is how it's supposed to look, but sometimes the numbers get rounded in
    #funny ways
    #assert repr(sc.frame) == '<ICRS Coordinate: ra=120.0 deg, dec=5.0 deg>'
    rscf = repr(sc.frame)
    assert rscf.startswith('<ICRS Coordinate: (ra, dec) in deg')

    # and  the string representation will be inherited from the low-level class.

    # same deal, should loook like this, but different archituectures/ python
    # versions may round the numbers differently
    #assert repr(sc) == '<SkyCoord (ICRS): ra=120.0 deg, dec=5.0 deg>'
    rsc = repr(sc)
    assert rsc.startswith('<SkyCoord (ICRS): (ra, dec) in deg')

    # Supports a variety of possible complex string formats
    sc = coords.SkyCoord('8h00m00s +5d00m00.0s', frame='icrs')

    # In the next example, the unit is only needed b/c units are ambiguous.  In
    # general, we *never* accept ambiguity
    sc = coords.SkyCoord('8:00:00 +5:00:00.0', unit=(u.hour, u.deg), frame='icrs')

    # The next one would yield length-2 array coordinates, because of the comma

    sc = coords.SkyCoord(['8h 5d', '2°2′3″ 0.3rad'], frame='icrs')

    # It should also interpret common designation styles as a coordinate
    # NOT YET
    # sc = coords.SkyCoord('SDSS J123456.89-012345.6', frame='icrs')

    # but it should also be possible to provide formats for outputting to strings,
    # similar to `Time`.  This can be added right away or at a later date.

    # transformation is done the same as for low-level classes, which it delegates to

    sc_fk5_j2001 = sc.transform_to(coords.FK5(equinox=J2001))
    assert sc_fk5_j2001.equinox == J2001

    # The key difference is that the high-level class remembers frame information
    # necessary for round-tripping, unlike the low-level classes:
    sc1 = coords.SkyCoord(ra=8 * u.hour, dec=5 * u.deg, equinox=J2001, frame='fk5')
    sc2 = sc1.transform_to('icrs')

    # The next assertion succeeds, but it doesn't mean anything for ICRS, as ICRS
    # isn't defined in terms of an equinox
    assert sc2.equinox == J2001

    # But it *is* necessary once we transform to FK5
    sc3 = sc2.transform_to('fk5')
    assert sc3.equinox == J2001
    npt.assert_allclose(sc1.ra, sc3.ra)

    # `SkyCoord` will also include the attribute-style access that is in the
    # v0.2/0.3 coordinate objects.  This will *not* be in the low-level classes
    sc = coords.SkyCoord(ra=8 * u.hour, dec=5 * u.deg, frame='icrs')
    scgal = sc.galactic
    assert str(scgal).startswith('<SkyCoord (Galactic): (l, b)')

    # the existing `from_name` and `match_to_catalog_*` methods will be moved to the
    # high-level class as convenience functionality.

    #in remote-data test below!
    #m31icrs = coords.SkyCoord.from_name('M31', frame='icrs')
    #assert str(m31icrs) == '<SkyCoord (ICRS) RA=10.68471 deg, Dec=41.26875 deg>'

    if HAS_SCIPY:
        cat1 = coords.SkyCoord(ra=[1, 2]*u.hr, dec=[3, 4.01]*u.deg, distance=[5, 6]*u.kpc, frame='icrs')
        cat2 = coords.SkyCoord(ra=[1, 2, 2.01]*u.hr, dec=[3, 4, 5]*u.deg, distance=[5, 200, 6]*u.kpc, frame='icrs')
        idx1, sep2d1, dist3d1 = cat1.match_to_catalog_sky(cat2)
        idx2, sep2d2, dist3d2 = cat1.match_to_catalog_3d(cat2)

        assert np.any(idx1 != idx2)

    # additional convenience functionality for the future should be added as methods
    # on `SkyCoord`, *not* the low-level classes.


@pytest.mark.remote_data
def test_highlevel_api_remote():
    m31icrs = coords.SkyCoord.from_name('M31', frame='icrs')

    assert str(m31icrs) == '<SkyCoord (ICRS): (ra, dec) in deg\n    (10.6847083, 41.26875)>'

    m31fk4 = coords.SkyCoord.from_name('M31', frame='fk4')

    assert m31icrs.frame != m31fk4.frame
    assert np.abs(m31icrs.ra - m31fk4.ra) > .5*u.deg
