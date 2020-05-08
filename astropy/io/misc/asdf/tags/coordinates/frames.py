# Licensed under a 3-clause BSD style license - see LICENSE.rst
# -*- coding: utf-8 -*-
import os
import glob

from asdf import tagged

import astropy.units as u
import astropy.coordinates
from astropy.coordinates.baseframe import frame_transform_graph
from astropy.units import Quantity
from astropy.coordinates import ICRS, Longitude, Latitude, Angle

from astropy.io.misc.asdf.types import AstropyType


__all__ = ['CoordType']

SCHEMA_PATH = os.path.abspath(
    os.path.join(os.path.dirname(__file__), '..', '..', 'data', 'schemas', 'astropy.org', 'astropy'))


def _get_frames():
    """
    By reading the schema files, get the list of all the frames we can
    save/load.
    """
    search = os.path.join(SCHEMA_PATH, 'coordinates', 'frames', '*.yaml')
    files = glob.glob(search)

    names = []
    for fpath in files:
        path, fname = os.path.split(fpath)
        frame, _ = fname.split('-')
        # Skip baseframe because we cannot directly save / load it.
        # Skip icrs because we have an explicit tag for it because there are
        # two versions.
        if frame not in ['baseframe', 'icrs']:
            names.append(frame)

    return names


class BaseCoordType:
    """
    This defines the base methods for coordinates, without defining anything
    related to asdf types. This allows subclasses with different types and
    schemas to use this without confusing the metaclass machinery.
    """
    @staticmethod
    def _tag_to_frame(tag):
        """
        Extract the frame name from the tag.
        """
        tag = tag[tag.rfind('/')+1:]
        tag = tag[:tag.rfind('-')]
        return frame_transform_graph.lookup_name(tag)

    @classmethod
    def _frame_name_to_tag(cls, frame_name):
        return cls.make_yaml_tag(cls._tag_prefix + frame_name)

    @classmethod
    def from_tree_tagged(cls, node, ctx):

        frame = cls._tag_to_frame(node._tag)

        data = node.get('data', None)
        if data is not None:
            return frame(node['data'], **node['frame_attributes'])

        return frame(**node['frame_attributes'])

    @classmethod
    def to_tree_tagged(cls, frame, ctx):
        if type(frame) not in frame_transform_graph.frame_set:
            raise ValueError("Can only save frames that are registered with the "
                             "transformation graph.")

        node = {}
        if frame.has_data:
            node['data'] = frame.data
        frame_attributes = {}
        for attr in frame.frame_attributes.keys():
            value = getattr(frame, attr, None)
            if value is not None:
                frame_attributes[attr] = value
        node['frame_attributes'] = frame_attributes

        return tagged.tag_object(cls._frame_name_to_tag(frame.name), node, ctx=ctx)

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(new, type(old))
        if new.has_data:
            assert u.allclose(new.data.lon, old.data.lon)
            assert u.allclose(new.data.lat, old.data.lat)


class CoordType(BaseCoordType, AstropyType):
    _tag_prefix = "coordinates/frames/"
    name = ["coordinates/frames/" + f for f in _get_frames()]
    types = [astropy.coordinates.BaseCoordinateFrame]
    handle_dynamic_subclasses = True
    requires = ['astropy']
    version = "1.0.0"


class ICRSType(CoordType):
    """
    Define a special tag for ICRS so we can make it version 1.1.0.
    """
    name = "coordinates/frames/icrs"
    types = ['astropy.coordinates.ICRS']
    version = "1.1.0"


class ICRSType10(AstropyType):
    name = "coordinates/frames/icrs"
    types = [astropy.coordinates.ICRS]
    requires = ['astropy']
    version = "1.0.0"

    @classmethod
    def from_tree(cls, node, ctx):
        wrap_angle = Angle(node['ra']['wrap_angle'])
        ra = Longitude(
            node['ra']['value'],
            unit=node['ra']['unit'],
            wrap_angle=wrap_angle)
        dec = Latitude(node['dec']['value'], unit=node['dec']['unit'])

        return ICRS(ra=ra, dec=dec)

    @classmethod
    def to_tree(cls, frame, ctx):
        node = {}

        wrap_angle = Quantity(frame.ra.wrap_angle)
        node['ra'] = {
            'value': frame.ra.value,
            'unit': frame.ra.unit.to_string(),
            'wrap_angle': wrap_angle
        }
        node['dec'] = {
            'value': frame.dec.value,
            'unit': frame.dec.unit.to_string()
        }

        return node

    @classmethod
    def assert_equal(cls, old, new):
        assert isinstance(old, ICRS)
        assert isinstance(new, ICRS)
        assert u.allclose(new.ra, old.ra)
        assert u.allclose(new.dec, old.dec)
