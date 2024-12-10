# Licensed under a 3-clause BSD style license - see LICENSE.rst


import sys
import warnings

import numpy as np
from matplotlib.patches import Polygon

from astropy import units as u
from astropy.coordinates import SkyCoord
from astropy.coordinates.matrix_utilities import rotation_matrix
from astropy.coordinates.representation import (
    SphericalRepresentation,
    UnitSphericalRepresentation,
)
from astropy.utils.exceptions import AstropyUserWarning

__all__ = ["Quadrangle", "SphericalCircle"]

# Monkey-patch the docs to fix CapStyle and JoinStyle subs.
# TODO! delete when upstream fix matplotlib/matplotlib#19839
if sys.flags.optimize < 2:
    Polygon.__init__.__doc__ = Polygon.__init__.__doc__.replace(
        "`.CapStyle`", "``matplotlib._enums.CapStyle``"
    )
    Polygon.__init__.__doc__ = Polygon.__init__.__doc__.replace(
        "`.JoinStyle`", "``matplotlib._enums.JoinStyle``"
    )
    Polygon.set_capstyle.__doc__ = Polygon.set_capstyle.__doc__.replace(
        "`.CapStyle`", "``matplotlib._enums.CapStyle``"
    )
    Polygon.set_joinstyle.__doc__ = Polygon.set_joinstyle.__doc__.replace(
        "`.JoinStyle`", "``matplotlib._enums.JoinStyle``"
    )


def _rotate_polygon(lon, lat, lon0, lat0):
    """
    Given a polygon with vertices defined by (lon, lat), rotate the polygon
    such that the North pole of the spherical coordinates is now at (lon0,
    lat0). Therefore, to end up with a polygon centered on (lon0, lat0), the
    polygon should initially be drawn around the North pole.
    """
    # Create a representation object
    polygon = UnitSphericalRepresentation(lon=lon, lat=lat)

    # Determine rotation matrix to make it so that the circle is centered
    # on the correct longitude/latitude.
    transform_matrix = rotation_matrix(-lon0, axis="z") @ rotation_matrix(
        -(0.5 * np.pi * u.radian - lat0), axis="y"
    )

    # Apply 3D rotation
    polygon = polygon.to_cartesian()
    polygon = polygon.transform(transform_matrix)
    polygon = UnitSphericalRepresentation.from_cartesian(polygon)

    return polygon.lon, polygon.lat


class SphericalCircle(Polygon):
    """
    Create a patch representing a spherical circle - that is, a circle that is
    formed of all the points that are within a certain angle of the central
    coordinates on a sphere. Here we assume that latitude goes from -90 to +90.

    This class is needed in cases where the user wants to add a circular patch
    to a celestial image, since otherwise the circle will be distorted, because
    a fixed interval in longitude corresponds to a different angle on the sky
    depending on the latitude.

    Parameters
    ----------
    center : tuple or `~astropy.units.Quantity` ['angle']
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements
        or a `~astropy.coordinates.SkyCoord` object.
    radius : `~astropy.units.Quantity` ['angle']
        The radius of the circle
    resolution : int, optional
        The number of points that make up the circle - increase this to get a
        smoother circle.
    vertex_unit : `~astropy.units.Unit`
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """

    def __init__(self, center, radius, resolution=100, vertex_unit=u.degree, **kwargs):
        # Extract longitude/latitude, either from a SkyCoord object, or
        # from a tuple of two quantities or a single 2-element Quantity.
        # The SkyCoord is converted to SphericalRepresentation, if not already.
        if isinstance(center, SkyCoord):
            rep_type = center.representation_type
            if not issubclass(
                rep_type, (SphericalRepresentation, UnitSphericalRepresentation)
            ):
                warnings.warn(
                    f"Received `center` of representation type {rep_type} "
                    "will be converted to SphericalRepresentation ",
                    AstropyUserWarning,
                )
            longitude, latitude = center.spherical.lon, center.spherical.lat
        else:
            longitude, latitude = center

        # Start off by generating the circle around the North pole
        lon = np.linspace(0.0, 2 * np.pi, resolution + 1)[:-1] * u.radian
        lat = np.repeat(0.5 * np.pi - radius.to_value(u.radian), resolution) * u.radian

        lon, lat = _rotate_polygon(lon, lat, longitude, latitude)

        # Extract new longitude/latitude in the requested units
        lon = lon.to_value(vertex_unit)
        lat = lat.to_value(vertex_unit)

        # Create polygon vertices
        vertices = np.array([lon, lat]).transpose()

        super().__init__(vertices, **kwargs)


class Quadrangle(Polygon):
    """
    Create a patch representing a latitude-longitude quadrangle.

    The edges of the quadrangle lie on two lines of constant longitude and two
    lines of constant latitude (or the equivalent component names in the
    coordinate frame of interest, such as right ascension and declination).
    Note that lines of constant latitude are not great circles.

    Unlike `matplotlib.patches.Rectangle`, the edges of this patch will render
    as curved lines if appropriate for the WCS transformation.

    Parameters
    ----------
    anchor : tuple or `~astropy.units.Quantity` ['angle']
        This can be either a tuple of two `~astropy.units.Quantity` objects, or
        a single `~astropy.units.Quantity` array with two elements.
    width : `~astropy.units.Quantity` ['angle']
        The width of the quadrangle in longitude (or, e.g., right ascension)
    height : `~astropy.units.Quantity` ['angle']
        The height of the quadrangle in latitude (or, e.g., declination)
    resolution : int, optional
        The number of points that make up each side of the quadrangle -
        increase this to get a smoother quadrangle.
    vertex_unit : `~astropy.units.Unit` ['angle']
        The units in which the resulting polygon should be defined - this
        should match the unit that the transformation (e.g. the WCS
        transformation) expects as input.

    Notes
    -----
    Additional keyword arguments are passed to `~matplotlib.patches.Polygon`
    """

    def __init__(
        self, anchor, width, height, resolution=100, vertex_unit=u.degree, **kwargs
    ):
        # Extract longitude/latitude, either from a tuple of two quantities, or
        # a single 2-element Quantity.
        longitude, latitude = u.Quantity(anchor).to_value(vertex_unit)

        # Convert the quadrangle dimensions to the appropriate units
        width = width.to_value(vertex_unit)
        height = height.to_value(vertex_unit)

        # Create progressions in longitude and latitude
        lon_seq = longitude + np.linspace(0, width, resolution + 1)
        lat_seq = latitude + np.linspace(0, height, resolution + 1)

        # Trace the path of the quadrangle
        lon = np.concatenate(
            [
                lon_seq[:-1],
                np.repeat(lon_seq[-1], resolution),
                np.flip(lon_seq[1:]),
                np.repeat(lon_seq[0], resolution),
            ]
        )
        lat = np.concatenate(
            [
                np.repeat(lat_seq[0], resolution),
                lat_seq[:-1],
                np.repeat(lat_seq[-1], resolution),
                np.flip(lat_seq[1:]),
            ]
        )

        # Create polygon vertices
        vertices = np.array([lon, lat]).transpose()

        super().__init__(vertices, **kwargs)
