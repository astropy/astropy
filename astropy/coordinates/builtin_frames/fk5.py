# -*- coding: utf-8 -*-
# Licensed under a 3-clause BSD style license - see LICENSE.rst
from __future__ import (absolute_import, unicode_literals, division,
                        print_function)

from ..baseframe import frame_transform_graph
from ..attributes import TimeAttribute
from ..transformations import DynamicMatrixTransform
from .. import earth_orientation as earth

from .baseradec import _base_radec_docstring, BaseRADecFrame
from .utils import EQUINOX_J2000


class FK5(BaseRADecFrame):
    """
    A coordinate or frame in the FK5 system.

    Note that this is a barycentric version of FK5 - that is, the origin for
    this frame is the Solar System Barycenter, *not* the Earth geocenter.

    The frame attributes are listed under **Other Parameters**.

    {params}

    Other parameters
    ----------------
    equinox : `~astropy.time.Time`
        The equinox of this frame.
    """

    equinox = TimeAttribute(default=EQUINOX_J2000)

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


FK5.__doc__ = FK5.__doc__.format(params=_base_radec_docstring)

# This is the "self-transform".  Defined at module level because the decorator
#  needs a reference to the FK5 class


@frame_transform_graph.transform(DynamicMatrixTransform, FK5, FK5)
def fk5_to_fk5(fk5coord1, fk5frame2):
    return fk5coord1._precession_matrix(fk5coord1.equinox, fk5frame2.equinox)
