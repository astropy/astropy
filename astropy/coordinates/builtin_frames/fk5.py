# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ...time import Time
from ..representation import SphericalRepresentation
from ..baseframe import (BaseCoordinateFrame, frame_transform_graph,
                         TimeFrameAttribute, RepresentationMapping)
from ..transformations import DynamicMatrixTransform

from .. import earth_orientation as earth
from .icrs import ICRS


# The UTC time scale is not properly defined prior to 1960, so Time('B1950',
# scale='utc') will emit a warning. Instead, we use Time('B1950', scale='tai')
# which is equivalent, but does not emit a warning.
_EQUINOX_J2000 = Time('J2000', scale='utc')
_EQUINOX_B1950 = Time('B1950', scale='tai')


class FK5(BaseCoordinateFrame):
    """
    A coordinate or frame in the FK5 system.

    Parameters
    ----------
    representation : `BaseRepresentation` or None
        A representation object or None to have no data (or use the other keywords)
    ra : `Angle`, optional, must be keyword
        The RA for this object (``dec`` must also be given and ``representation``
        must be None).
    dec : `Angle`, optional, must be keyword
        The Declination for this object (``ra`` must also be given and
        ``representation`` must be None).
    distance : `~astropy.units.Quantity`, optional, must be keyword
        The Distance for this object along the line-of-sight.
        (``representation`` must be None).
    equinox : `~astropy.time.Time`, optional, must be keyword
        The equinox of this frame.
    """
    frame_specific_representation_info = {
        'spherical': [RepresentationMapping('lon', 'ra'),
                      RepresentationMapping('lat', 'dec')]
    }
    frame_specific_representation_info['unitspherical'] = \
        frame_specific_representation_info['spherical']

    default_representation = SphericalRepresentation
    equinox = TimeFrameAttribute(default=_EQUINOX_J2000)

    @staticmethod
    def _precession_matrix(oldequinox, newequinox):
        """
        Compute and return the precession matrix for FK5 based on Capitaine et
        al. 2003/IAU2006.  Used inside some of the transformation functions.

        Parameters
        ----------
        oldequinox : `~astropy.time.Time`
            The equinox to precess from.
        newequinox : `~astropy.time.Time`
            The equinox to precess to.

        Returns
        -------
        newcoord : array
            The precession matrix to transform to the new equinox
        """
        return earth.precession_matrix_Capitaine(oldequinox, newequinox)


# Has to be defined at module level, because `transform` needs an FK5 reference
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK5)
def fk5_to_fk5(fk5coord1, fk5frame2):
    return fk5coord1._precession_matrix(fk5coord1.equinox, fk5frame2.equinox)


# ICRS to/from FK5 -------------------------------->
@frame_transform_graph.transform(DynamicMatrixTransform, ICRS, FK5)
def icrs_to_fk5(icrscoord, fk5frame):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5frame._precession_matrix(_EQUINOX_J2000, fk5frame.equinox)
    return pmat * icrscoord._ICRS_TO_FK5_J2000_MAT


# can't be static because the equinox is needed
@frame_transform_graph.transform(DynamicMatrixTransform, FK5, ICRS)
def fk5_to_icrs(fk5coord, icrsframe):
    # ICRS is by design very close to J2000 equinox
    pmat = fk5coord._precession_matrix(fk5coord.equinox, _EQUINOX_J2000)
    return icrsframe._ICRS_TO_FK5_J2000_MAT.T * pmat
