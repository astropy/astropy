"""
This file tests the behaviour of subclasses of Representation and Frames
"""

from copy import deepcopy
from collections import OrderedDict

from astropy.coordinates import Longitude, Latitude
from astropy.coordinates.representation import (REPRESENTATION_CLASSES,
                                                SphericalRepresentation,
                                                UnitSphericalRepresentation)
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.coordinates.transformations import FunctionTransform
from astropy.coordinates import ICRS
from astropy.coordinates.baseframe import RepresentationMapping

import astropy.units as u

import astropy.coordinates

# Classes setup, borrowed from SunPy.

# Here we define the classes *inside* the tests to make sure that we can wipe
# the slate clean when the tests have finished running.


def setup_function(func):
    func.REPRESENTATION_CLASSES_ORIG = deepcopy(REPRESENTATION_CLASSES)


def teardown_function(func):
    REPRESENTATION_CLASSES.clear()
    REPRESENTATION_CLASSES.update(func.REPRESENTATION_CLASSES_ORIG)


def test_unit_representation_subclass():

    class Longitude180(Longitude):
        def __new__(cls, angle, unit=None, wrap_angle=180*u.deg, **kwargs):
            self = super(Longitude180, cls).__new__(cls, angle, unit=unit,
                                                    wrap_angle=wrap_angle, **kwargs)
            return self

    class UnitSphericalWrap180Representation(UnitSphericalRepresentation):
        attr_classes = OrderedDict([('lon', Longitude180),
                                    ('lat', Latitude)])
        recommended_units = {'lon': u.deg, 'lat': u.deg}

    class SphericalWrap180Representation(SphericalRepresentation):
        attr_classes = OrderedDict([('lon', Longitude180),
                                    ('lat', Latitude),
                                    ('distance', u.Quantity)])
        recommended_units = {'lon': u.deg, 'lat': u.deg}

        _unit_representation = UnitSphericalWrap180Representation

    class myframe(ICRS):
        default_representation = SphericalWrap180Representation
        frame_specific_representation_info = {
            'spherical': [RepresentationMapping('lon', 'ra'),
                          RepresentationMapping('lat', 'dec')]
        }
        frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['unitsphericalwrap180'] = \
        frame_specific_representation_info['sphericalwrap180'] = \
            frame_specific_representation_info['spherical']

    @frame_transform_graph.transform(FunctionTransform,
                                     myframe, astropy.coordinates.ICRS)
    def myframe_to_icrs(myframe_coo, icrs):
        return icrs.realize_frame(myframe_coo._data)

    f = myframe(10*u.deg, 10*u.deg)
    assert isinstance(f._data, UnitSphericalWrap180Representation)
    assert isinstance(f.ra, Longitude180)

    g = f.transform_to(astropy.coordinates.ICRS)
    assert isinstance(g, astropy.coordinates.ICRS)
    assert isinstance(g._data, UnitSphericalWrap180Representation)

    frame_transform_graph.remove_transform(myframe,
                                           astropy.coordinates.ICRS,
                                           None)
