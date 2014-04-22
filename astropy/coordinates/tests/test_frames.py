# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst

from __future__ import (absolute_import, division, print_function,
                        unicode_literals)

import numpy as np
from numpy.testing import assert_allclose

from ... import units as u
from ...tests.helper import pytest
from .. import  representation

def test_create_data_frames():
    from ..builtin_frames import ICRS

    #from repr
    i1 = ICRS(representation.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc))
    i2 = ICRS(representation.UnitSphericalRepresentation(lon=1*u.deg, lat=2*u.deg))

    #from preferred name
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)
    i4 = ICRS(ra=1*u.deg, dec=2*u.deg)

    assert i1.data.lat == i3.data.lat
    assert i1.data.lon == i3.data.lon
    assert i1.data.distance == i3.data.distance

    assert i2.data.lat == i4.data.lat
    assert i2.data.lon == i4.data.lon

    #now make sure the preferred names work as properties
    assert_allclose(i1.ra, i3.ra)
    assert_allclose(i2.ra, i4.ra)
    assert_allclose(i1.distance, i3.distance)

    with pytest.raises(AttributeError):
        i1.ra = [11.]*u.deg


def test_create_nodata_frames():
    from ..builtin_frames import ICRS, FK4, FK5

    i = ICRS()
    assert len(i.frame_attr_names) == 0

    f5 = FK5()
    assert f5.equinox == FK5.frame_attr_names['equinox']

    f4 = FK4()
    assert f4.equinox == FK4.frame_attr_names['equinox']
    assert f4.obstime == FK4.frame_attr_names['obstime']


def test_frame_repr():
    from ..builtin_frames import ICRS, FK5

    i = ICRS()
    assert repr(i) == '<ICRS frame>'

    f5 = FK5()
    assert repr(f5).startswith('<FK5 frame: equinox=')

    i2 = ICRS(ra=1*u.deg, dec=2*u.deg)
    i3 = ICRS(ra=1*u.deg, dec=2*u.deg, distance=3*u.kpc)

    assert repr(i2).startswith('<ICRS coordinate: ')
    assert 'ra=' in repr(i2)
    assert 'distance=' in repr(i3)

def test_realizing():
    from ..builtin_frames import ICRS, FK5
    from ...time import Time

    rep = representation.SphericalRepresentation(1*u.deg, 2*u.deg, 3*u.kpc)

    i = ICRS()
    i2 = i.realize_frame(rep)

    assert not i.has_data
    assert i2.has_data

    f = FK5(equinox=Time('J2001', scale='utc'))
    f2 = f.realize_frame(rep)

    assert not f.has_data
    assert f2.has_data

    assert f2.equinox == f.equinox
    assert f2.equinox != FK5.frame_attr_names['equinox']


def test_getitem():
    from ..builtin_frames import ICRS

    rep = representation.SphericalRepresentation(
        [1, 2, 3]*u.deg, [4, 5, 6]*u.deg, [7, 8, 9]*u.kpc)

    i = ICRS(rep)
    assert len(i.ra) == 3

    iidx = i[1:]
    assert len(iidx.ra) == 2

    iidx2 = i[0]
    assert iidx2.ra.isscalar

def test_transform():
    """
    This test just makes sure the transform architecture works, but does *not*
    actually test all the builtin transforms themselves are accurate
    """
    from ..builtin_frames import ICRS, FK4, FK5
    from ...time import Time

    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ == representation.UnitSphericalRepresentation

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)


    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg, distance=[5, 6]*u.kpc)
    f = i.transform_to(FK5)
    i2 = f.transform_to(ICRS)

    assert i2.data.__class__ != representation.UnitSphericalRepresentation



    #make sure self-transforms also work
    i = ICRS(ra=[1, 2]*u.deg, dec=[3, 4]*u.deg)
    i2 = i.transform_to(ICRS)

    assert_allclose(i.ra, i2.ra)
    assert_allclose(i.dec, i2.dec)

